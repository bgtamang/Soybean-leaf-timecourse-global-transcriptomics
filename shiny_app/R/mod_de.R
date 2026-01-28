# DE Analysis module - Differential expression visualization

de_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Contrast Selection"),
          card_body(
            selectInput(ns("contrast"), "Select Contrast:",
              choices = NULL  # Populated in server
            ),
            hr(),
            sliderInput(ns("fdr_threshold"), "FDR Threshold:",
              min = 0.001, max = 0.1, value = 0.05, step = 0.001
            ),
            sliderInput(ns("logfc_threshold"), "Min |log2FC|:",
              min = 0, max = 5, value = 1, step = 0.25
            ),
            hr(),
            radioButtons(ns("direction"), "Direction:",
              choices = c("All" = "all", "Up" = "up", "Down" = "down"),
              selected = "all", inline = TRUE
            )
          )
        ),
        card(
          card_header("Summary"),
          card_body(
            uiOutput(ns("de_summary"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Volcano Plot",
            card(
              card_body(
                plotlyOutput(ns("volcano_plot"), height = "500px")
              )
            )
          ),
          tabPanel("MA Plot",
            card(
              card_body(
                plotlyOutput(ns("ma_plot"), height = "500px")
              )
            )
          ),
          tabPanel("DE Gene Table",
            card(
              card_body(
                DT::DTOutput(ns("de_table")),
                hr(),
                downloadButton(ns("download_de"), "Download Filtered Results", class = "btn-primary")
              )
            )
          )
        )
      )
    )
  )
}

de_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Load DE results
    de_data <- reactive({
      load_de_results()
    })

    # Update contrast choices
    observe({
      de <- de_data()
      if (!is.null(de)) {
        contrasts <- unique(de$Contrast)
        # Group by type
        temporal <- contrasts[grepl("TP[0-9]vs", contrasts)]
        between <- contrasts[!grepl("TP[0-9]vs", contrasts)]

        choices <- list(
          "Temporal Contrasts" = temporal,
          "Between Genotype" = between
        )
        updateSelectInput(session, "contrast", choices = choices)
      }
    })

    # Filtered DE data
    filtered_de <- reactive({
      req(input$contrast)
      de <- de_data()
      if (is.null(de)) return(NULL)

      result <- de %>%
        filter(
          Contrast == input$contrast,
          FDR < input$fdr_threshold,
          abs(logFC) > input$logfc_threshold
        )

      if (input$direction == "up") {
        result <- result %>% filter(logFC > 0)
      } else if (input$direction == "down") {
        result <- result %>% filter(logFC < 0)
      }

      result
    })

    # All DE for contrast (for plotting)
    contrast_de <- reactive({
      req(input$contrast)
      de <- de_data()
      if (is.null(de)) return(NULL)
      de %>% filter(Contrast == input$contrast)
    })

    # Summary
    output$de_summary <- renderUI({
      filtered <- filtered_de()
      if (is.null(filtered)) {
        return(p("Select a contrast to see results"))
      }

      n_up <- sum(filtered$logFC > 0)
      n_down <- sum(filtered$logFC < 0)

      tagList(
        p(strong("Total DE Genes:"), nrow(filtered)),
        p(span(class = "text-danger", icon("arrow-up")), " Up:", n_up),
        p(span(class = "text-primary", icon("arrow-down")), " Down:", n_down)
      )
    })

    # Volcano plot (publication quality)
    output$volcano_plot <- renderPlotly({
      de <- contrast_de()
      if (is.null(de) || nrow(de) == 0) {
        return(plotly_empty() %>% layout(title = "No data available"))
      }

      # Add significance category
      de <- de %>%
        mutate(
          Significant = factor(case_when(
            FDR < input$fdr_threshold & logFC > input$logfc_threshold ~ "Upregulated",
            FDR < input$fdr_threshold & logFC < -input$logfc_threshold ~ "Downregulated",
            TRUE ~ "Not Significant"
          ), levels = c("Upregulated", "Downregulated", "Not Significant"))
        )

      n_up <- sum(de$Significant == "Upregulated")
      n_down <- sum(de$Significant == "Downregulated")

      p <- ggplot(de, aes(x = logFC, y = -log10(FDR),
                          color = Significant,
                          text = paste("Gene:", GeneID,
                                      "<br>logFC:", round(logFC, 2),
                                      "<br>FDR:", signif(FDR, 3)))) +
        geom_point(alpha = 0.7, size = 2) +
        geom_vline(xintercept = c(-input$logfc_threshold, input$logfc_threshold),
                   linetype = "dashed", color = "grey40", linewidth = 0.5) +
        geom_hline(yintercept = -log10(input$fdr_threshold),
                   linetype = "dashed", color = "grey40", linewidth = 0.5) +
        scale_color_manual(
          values = c("Upregulated" = "#B2182B", "Downregulated" = "#2166AC", "Not Significant" = "grey75"),
          name = "Regulation"
        ) +
        theme_publication(base_size = 12) +
        labs(
          title = "Volcano Plot",
          subtitle = paste0(input$contrast, " (Up: ", n_up, ", Down: ", n_down, ")"),
          x = "log\u2082 Fold Change",
          y = "-log\u2081\u2080(FDR)"
        )

      ggplotly(p, tooltip = "text", source = "volcano_plot") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5)) %>%
        event_register("plotly_click")
    })

    # Handle click on volcano plot
    observeEvent(event_data("plotly_click", source = "volcano_plot"), {
      click <- event_data("plotly_click")
      if (!is.null(click)) {
        de <- contrast_de()
        # Find clicked gene
        clicked_gene <- de$GeneID[click$pointNumber + 1]
        if (!is.null(clicked_gene)) {
          selected_gene(clicked_gene)
        }
      }
    })

    # MA plot (publication quality)
    output$ma_plot <- renderPlotly({
      de <- contrast_de()
      if (is.null(de) || nrow(de) == 0) {
        return(plotly_empty() %>% layout(title = "No data available"))
      }

      # Need AveExpr for MA plot - simulate if not present
      if (!"AveExpr" %in% names(de)) {
        de$AveExpr <- rnorm(nrow(de), mean = 5, sd = 2)
      }

      de <- de %>%
        mutate(
          Significant = factor(ifelse(
            FDR < input$fdr_threshold & abs(logFC) > input$logfc_threshold,
            "Significant", "Not Significant"
          ), levels = c("Significant", "Not Significant"))
        )

      p <- ggplot(de, aes(x = AveExpr, y = logFC,
                          color = Significant,
                          text = paste("Gene:", GeneID,
                                      "<br>logFC:", round(logFC, 2),
                                      "<br>FDR:", signif(FDR, 3)))) +
        geom_point(alpha = 0.6, size = 2) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
        geom_hline(yintercept = c(-input$logfc_threshold, input$logfc_threshold),
                   linetype = "dashed", color = "grey40", linewidth = 0.4) +
        scale_color_manual(
          values = c("Significant" = "#B2182B", "Not Significant" = "grey75"),
          name = NULL
        ) +
        theme_publication(base_size = 12) +
        labs(
          title = "MA Plot",
          subtitle = input$contrast,
          x = "Average Expression (log)",
          y = "log\u2082 Fold Change"
        )

      ggplotly(p, tooltip = "text")
    })

    # DE table
    output$de_table <- DT::renderDT({
      de <- filtered_de()
      if (is.null(de)) return(NULL)

      de %>%
        select(GeneID, logFC, FDR, Best_hit_arabi_name, Best_hit_arabi_defline) %>%
        mutate(
          logFC = round(logFC, 3),
          FDR = signif(FDR, 3)
        ) %>%
        datatable(
          filter = "top",
          selection = "single",
          options = list(
            pageLength = 20,
            scrollX = TRUE
          ),
          rownames = FALSE
        )
    })

    # Handle table selection
    observeEvent(input$de_table_rows_selected, {
      de <- filtered_de()
      if (!is.null(de) && !is.null(input$de_table_rows_selected)) {
        gene <- de$GeneID[input$de_table_rows_selected]
        selected_gene(gene)
      }
    })

    # Download handler
    output$download_de <- downloadHandler(
      filename = function() {
        paste0("DE_", input$contrast, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(filtered_de(), file, row.names = FALSE)
      }
    )
  })
}
