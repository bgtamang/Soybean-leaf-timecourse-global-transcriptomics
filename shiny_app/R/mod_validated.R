# Validated Genes module - 79 high-confidence JAG1 targets (Figure 6E)

validated_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Filter Validated Genes"),
          card_body(
            checkboxGroupInput(ns("de_tier"), "DE Tier:",
              choices = c("Gold", "Silver", "Bronze"),
              selected = c("Gold", "Silver", "Bronze")
            ),
            hr(),
            selectInput(ns("binding_class"), "Binding Class:",
              choices = c("All", "Both", "ChIP-only", "DAP-only"),
              selected = "All"
            ),
            hr(),
            sliderInput(ns("min_pheno_r"), "Min Phenotype |r|:",
              min = 0, max = 1, value = 0, step = 0.05
            ),
            sliderInput(ns("min_logfc"), "Min |logFC|:",
              min = 0, max = 10, value = 0, step = 0.5
            )
          )
        ),
        card(
          card_header("Evidence Summary"),
          card_body(
            uiOutput(ns("evidence_panel"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Scatter Plot",
            card(
              card_body(
                plotlyOutput(ns("scatter_plot"), height = "500px")
              )
            )
          ),
          tabPanel("Binding Classes",
            card(
              card_body(
                plotlyOutput(ns("binding_bar"), height = "400px")
              )
            )
          ),
          tabPanel("Gene Table",
            card(
              card_body(
                DT::DTOutput(ns("gene_table")),
                hr(),
                downloadButton(ns("download_validated"), "Download Validated Genes",
                               class = "btn-primary")
              )
            )
          )
        )
      )
    )
  )
}

validated_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Filtered validated genes
    filtered_data <- reactive({
      if (is.null(validated_genes)) return(NULL)

      data <- validated_genes %>%
        filter(DE_Tier %in% input$de_tier)

      if (input$binding_class != "All") {
        data <- data %>% filter(Binding_Class == input$binding_class)
      }

      if (input$min_pheno_r > 0 && "Pheno_r" %in% names(data)) {
        data <- data %>% filter(abs(Pheno_r) >= input$min_pheno_r)
      }

      if (input$min_logfc > 0 && "logFC" %in% names(data)) {
        data <- data %>% filter(abs(logFC) >= input$min_logfc)
      }

      data
    })

    # Evidence summary sidebar
    output$evidence_panel <- renderUI({
      if (!is.null(evidence_summary)) {
        tags$table(class = "table table-sm",
          lapply(seq_len(nrow(evidence_summary)), function(i) {
            tags$tr(
              tags$td(class = "small", evidence_summary$Metric[i]),
              tags$td(tags$b(evidence_summary$Count[i]))
            )
          })
        )
      } else if (!is.null(validated_genes)) {
        tagList(
          p(tags$b("Total validated:"), nrow(validated_genes)),
          p(tags$b("With binding:"), sum(validated_genes$Has_Binding, na.rm = TRUE)),
          p(tags$b("With phenotype:"), sum(validated_genes$Has_Pheno, na.rm = TRUE))
        )
      } else {
        p("Evidence data not available")
      }
    })

    # Scatter plot: logFC vs phenotype correlation
    output$scatter_plot <- renderPlotly({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No data matches filters"))
      }

      # Check required columns
      has_pheno <- "Pheno_r" %in% names(data) || "Cor_LW_Ratio" %in% names(data)
      has_logfc <- "logFC" %in% names(data) || "Mean_logFC" %in% names(data)

      if (!has_pheno || !has_logfc) {
        return(plotly_empty() %>% layout(title = "Required columns not available"))
      }

      # Use available column names
      pheno_col <- if ("Pheno_r" %in% names(data)) "Pheno_r" else "Cor_LW_Ratio"
      fc_col <- if ("logFC" %in% names(data)) "logFC" else "Mean_logFC"

      data <- data %>%
        mutate(
          pheno_r = .data[[pheno_col]],
          log_fc = .data[[fc_col]],
          evidence_n = if ("Evidence_Count" %in% names(.)) Evidence_Count else 1
        ) %>%
        filter(!is.na(pheno_r), !is.na(log_fc))

      if (nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No genes with both logFC and phenotype data"))
      }

      p <- ggplot(data, aes(x = log_fc, y = pheno_r,
                            color = DE_Tier,
                            size = evidence_n,
                            text = paste("Gene:", GeneID,
                                        "<br>DE Tier:", DE_Tier,
                                        "<br>logFC:", round(log_fc, 2),
                                        "<br>Pheno r:", round(pheno_r, 3),
                                        "<br>Binding:", Binding_Class,
                                        "<br>Evidence:", evidence_n))) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_point(alpha = 0.7, stroke = 0.3) +
        scale_color_manual(values = tier_colors, name = "DE Tier") +
        scale_size_continuous(range = c(3, 8), name = "Evidence\nCount") +
        theme_publication(base_size = 12) +
        labs(
          x = "log2 Fold Change (Narrow vs Broad)",
          y = "Phenotype Correlation (r)",
          title = paste("79 Validated JAG1 Targets:", nrow(data), "shown"),
          subtitle = "Sized by evidence count, colored by DE tier"
        )

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    })

    # Binding class bar chart
    output$binding_bar <- renderPlotly({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No data"))
      }

      if (!"Binding_Class" %in% names(data)) {
        return(plotly_empty() %>% layout(title = "Binding class data not available"))
      }

      binding_counts <- data %>%
        group_by(Binding_Class) %>%
        summarise(Count = n(), .groups = "drop")

      bind_colors <- c("Both" = "#27AE60", "ChIP-only" = "#3498DB",
                        "DAP-only" = "#E67E22", "None" = "#95A5A6")

      p <- ggplot(binding_counts, aes(x = reorder(Binding_Class, Count),
                                      y = Count,
                                      fill = Binding_Class)) +
        geom_col(width = 0.65, color = "black", linewidth = 0.4) +
        geom_text(aes(label = Count), hjust = -0.3, size = 4.5, fontface = "bold") +
        coord_flip() +
        scale_fill_manual(values = bind_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        theme_publication(base_size = 12) +
        theme(legend.position = "none", panel.grid.major.y = element_blank()) +
        labs(x = NULL, y = "Number of Genes",
             title = "Validated Genes by Binding Evidence")

      ggplotly(p, tooltip = c("x", "y"))
    })

    # Gene table
    output$gene_table <- DT::renderDT({
      data <- filtered_data()
      if (is.null(data)) return(NULL)

      display_cols <- c("GeneID", "DE_Tier", "Mean_logFC", "Binding_Class",
                        "Evidence_Count", "Priority_Score", "Final_Tier",
                        "Cor_LW_Ratio", "Pval_LW_Ratio", "Has_ChIPseq",
                        "Has_DAPseq", "Has_Pheno", "Arabi_Defline")
      display_cols <- display_cols[display_cols %in% names(data)]

      data %>%
        select(all_of(display_cols)) %>%
        datatable(
          filter = "top",
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        ) %>%
        formatStyle(
          "DE_Tier",
          backgroundColor = styleEqual(
            c("Gold", "Silver", "Bronze"),
            c("#FFD70050", "#C0C0C050", "#CD7F3250")
          )
        )
    })

    # Handle table selection → update selected_gene
    observeEvent(input$gene_table_rows_selected, {
      data <- filtered_data()
      if (!is.null(data) && !is.null(input$gene_table_rows_selected)) {
        gene <- data$GeneID[input$gene_table_rows_selected]
        selected_gene(gene)
      }
    })

    # Download handler
    output$download_validated <- downloadHandler(
      filename = function() {
        paste0("validated_79_genes_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    )
  })
}
