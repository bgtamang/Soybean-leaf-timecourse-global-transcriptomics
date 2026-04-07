# JAG1 Targets module - Target gene browser

jag1_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Filter Targets"),
          card_body(
            checkboxGroupInput(ns("tiers"), "Confidence Tier:",
              choices = c("Gold", "Silver", "Bronze"),
              selected = c("Gold", "Silver")
            ),
            hr(),
            sliderInput(ns("logfc_min"), "Min |Mean logFC|:",
              min = 0, max = 10, value = 1, step = 0.5
            ),
            sliderInput(ns("score_min"), "Min Confidence Score:",
              min = 0, max = 1, value = 0, step = 0.1
            ),
            hr(),
            selectInput(ns("pattern"), "Expression Pattern:",
              choices = c("All", "Persistent", "TP1_Specific", "Late_Onset", "Variable"),
              selected = "All"
            ),
            hr(),
            checkboxInput(ns("jag1_correlated"), "JAG1-correlated only", FALSE)
          )
        ),
        card(
          card_header("Statistics"),
          card_body(
            uiOutput(ns("stats_panel"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Target Table",
            card(
              card_body(
                DT::DTOutput(ns("target_table")),
                hr(),
                downloadButton(ns("download_targets"), "Download Filtered Targets", class = "btn-primary")
              )
            )
          ),
          tabPanel("Tier Overview",
            card(
              card_body(
                fluidRow(
                  column(6, plotlyOutput(ns("tier_scatter"), height = "400px")),
                  column(6, plotlyOutput(ns("tier_bar"), height = "400px"))
                )
              )
            )
          ),
          tabPanel("Expression Patterns",
            card(
              card_body(
                plotOutput(ns("pattern_heatmap"), height = "500px")
              )
            )
          ),
          tabPanel("Selected Gene",
            card(
              card_body(
                uiOutput(ns("gene_detail"))
              )
            )
          )
        )
      )
    )
  )
}

jag1_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Filtered targets
    filtered_targets <- reactive({
      result <- jag1_targets %>%
        filter(Confidence_Tier %in% input$tiers)

      if (input$logfc_min > 0) {
        result <- result %>% filter(abs(Mean_logFC_Pairwise) >= input$logfc_min)
      }

      if (input$score_min > 0) {
        result <- result %>% filter(Confidence_Score >= input$score_min)
      }

      if (input$pattern != "All") {
        result <- result %>% filter(Pattern == input$pattern)
      }

      if (input$jag1_correlated) {
        result <- result %>% filter(Correlated_with_JAG1 == TRUE)
      }

      result
    })

    # Statistics panel
    output$stats_panel <- renderUI({
      targets <- filtered_targets()

      tagList(
        p(strong("Filtered Targets:"), nrow(targets)),
        tags$hr(),
        tags$table(class = "table table-sm",
          tags$tr(
            tags$td(span(class = "badge", style = "background-color: #FFD700", "Gold")),
            tags$td(sum(targets$Confidence_Tier == "Gold"))
          ),
          tags$tr(
            tags$td(span(class = "badge", style = "background-color: #C0C0C0", "Silver")),
            tags$td(sum(targets$Confidence_Tier == "Silver"))
          ),
          tags$tr(
            tags$td(span(class = "badge", style = "background-color: #CD7F32", "Bronze")),
            tags$td(sum(targets$Confidence_Tier == "Bronze"))
          )
        )
      )
    })

    # Target table
    output$target_table <- DT::renderDT({
      targets <- filtered_targets()

      display_cols <- c("Rank", "GeneID", "Confidence_Tier", "Mean_logFC_Pairwise",
                        "Confidence_Score", "Pattern", "Best_hit_arabi_name")
      display_cols <- display_cols[display_cols %in% names(targets)]

      targets %>%
        select(all_of(display_cols)) %>%
        mutate(
          Mean_logFC_Pairwise = round(Mean_logFC_Pairwise, 3),
          Confidence_Score = round(Confidence_Score, 3)
        ) %>%
        datatable(
          filter = "top",
          selection = "single",
          options = list(
            pageLength = 20,
            scrollX = TRUE,
            order = list(list(0, "asc"))
          ),
          rownames = FALSE
        ) %>%
        formatStyle(
          "Confidence_Tier",
          backgroundColor = styleEqual(
            c("Gold", "Silver", "Bronze"),
            c("#FFD70050", "#C0C0C050", "#CD7F3250")
          )
        )
    })

    # Handle table selection
    observeEvent(input$target_table_rows_selected, {
      targets <- filtered_targets()
      if (!is.null(input$target_table_rows_selected)) {
        gene <- targets$GeneID[input$target_table_rows_selected]
        selected_gene(gene)
      }
    })

    # Tier scatter plot (publication quality)
    output$tier_scatter <- renderPlotly({
      targets <- filtered_targets()

      p <- ggplot(targets, aes(x = Mean_logFC_Pairwise, y = Confidence_Score,
                               color = Confidence_Tier,
                               text = paste("Gene:", GeneID,
                                           "<br>Tier:", Confidence_Tier,
                                           "<br>logFC:", round(Mean_logFC_Pairwise, 2)))) +
        geom_point(alpha = 0.8, size = 3, stroke = 0.3) +
        geom_hline(yintercept = c(0.33, 0.66), linetype = "dashed",
                   color = "grey60", linewidth = 0.4) +
        scale_color_manual(values = tier_colors, name = "Confidence\nTier") +
        theme_publication(base_size = 12) +
        labs(
          title = "JAG1 Target Distribution",
          subtitle = paste(nrow(targets), "targets shown"),
          x = "Mean log\u2082 Fold Change",
          y = "Confidence Score"
        )

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    })

    # Tier bar chart (publication quality)
    output$tier_bar <- renderPlotly({
      tier_counts <- filtered_targets() %>%
        group_by(Confidence_Tier) %>%
        summarise(Count = n(), .groups = "drop") %>%
        mutate(Confidence_Tier = factor(Confidence_Tier, levels = c("Gold", "Silver", "Bronze")))

      p <- ggplot(tier_counts, aes(x = Confidence_Tier, y = Count, fill = Confidence_Tier)) +
        geom_col(width = 0.65, color = "black", linewidth = 0.4) +
        geom_text(aes(label = Count), vjust = -0.5, size = 4.5, fontface = "bold") +
        scale_fill_manual(values = tier_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        theme_publication(base_size = 12) +
        theme(
          legend.position = "none",
          panel.grid.major.x = element_blank()
        ) +
        labs(
          title = "Target Count by Confidence Tier",
          x = NULL,
          y = "Number of Genes"
        )

      ggplotly(p, tooltip = c("x", "y"))
    })

    # Pattern heatmap
    output$pattern_heatmap <- renderPlot({
      targets <- filtered_targets()
      if (nrow(targets) == 0) {
        plot.new()
        text(0.5, 0.5, "No targets match filters", cex = 1.5)
        return()
      }

      # Get FC columns for timepoints
      fc_cols <- c("FC_TP1", "FC_TP2", "FC_TP3", "FC_TP4", "FC_TP5")
      fc_cols <- fc_cols[fc_cols %in% names(targets)]

      if (length(fc_cols) == 0) {
        plot.new()
        text(0.5, 0.5, "Temporal FC data not available", cex = 1.5)
        return()
      }

      # Limit to top 50 for visibility
      top_targets <- targets %>%
        arrange(Rank) %>%
        head(50)

      mat <- as.matrix(top_targets[, fc_cols])
      rownames(mat) <- top_targets$GeneID

      # Annotation
      row_anno <- data.frame(
        Tier = top_targets$Confidence_Tier,
        row.names = top_targets$GeneID
      )

      # Publication-quality heatmap
      gg_heatmap(
        mat,
        title = "JAG1 Target Expression Across Timepoints (Top 50)",
        row_anno = row_anno,
        cluster_rows = TRUE,
        fontsize_row = 7,
        legend_title = expression(log[2]*FC)
      )
    })

    # Gene detail panel
    output$gene_detail <- renderUI({
      gene <- selected_gene()
      if (is.null(gene)) {
        return(p(class = "text-muted", "Select a gene from the table to see details"))
      }

      target <- jag1_targets %>% filter(GeneID == gene)
      if (nrow(target) == 0) {
        return(p(class = "text-warning", "Gene not found in JAG1 targets"))
      }

      target <- target[1, ]

      tagList(
        h4(target$GeneID),
        if (!is.na(target$Best_hit_arabi_name)) {
          p(tags$em(target$Best_hit_arabi_name))
        },
        hr(),
        fluidRow(
          column(6,
            tags$dl(
              tags$dt("Confidence Tier"),
              tags$dd(span(class = paste0("badge bg-", tolower(target$Confidence_Tier)), target$Confidence_Tier)),
              tags$dt("Rank"),
              tags$dd(target$Rank),
              tags$dt("Mean logFC"),
              tags$dd(round(target$Mean_logFC_Pairwise, 3)),
              tags$dt("Confidence Score"),
              tags$dd(round(target$Confidence_Score, 3))
            )
          ),
          column(6,
            tags$dl(
              tags$dt("Pattern"),
              tags$dd(target$Pattern),
              tags$dt("JAG1 Correlated"),
              tags$dd(ifelse(target$Correlated_with_JAG1, "Yes", "No")),
              if (!is.na(target$JAG1_Correlation_TP1)) tagList(
                tags$dt("JAG1 Correlation (TP1)"),
                tags$dd(round(target$JAG1_Correlation_TP1, 3))
              )
            )
          )
        ),
        if (!is.na(target$Best_hit_arabi_defline)) {
          tagList(
            hr(),
            h5("Annotation"),
            p(target$Best_hit_arabi_defline)
          )
        }
      )
    })

    # Download handler
    output$download_targets <- downloadHandler(
      filename = function() {
        paste0("JAG1_targets_filtered_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(filtered_targets(), file, row.names = FALSE)
      }
    )
  })
}
