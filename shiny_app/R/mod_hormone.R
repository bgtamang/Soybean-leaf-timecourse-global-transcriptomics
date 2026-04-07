# Hormone Pathway module - Hormone enrichment among JAG1 targets (Figure 4E)

hormone_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Pathway Filters"),
          card_body(
            checkboxGroupInput(ns("pathways"), "Hormone Pathways:",
              choices = c("SA", "Auxin", "ABA", "JA", "Ethylene",
                          "Brassinosteroid", "Cytokinin", "Gibberellin"),
              selected = c("SA", "Auxin", "ABA", "JA", "Ethylene",
                           "Brassinosteroid", "Cytokinin", "Gibberellin")
            ),
            hr(),
            checkboxInput(ns("sig_only"), "Significant only (P < 0.05)", FALSE)
          )
        ),
        card(
          card_header("Key Finding"),
          card_body(
            p(tags$b("SA and Auxin pathways"), " are significantly enriched",
              " among JAG1 targets:"),
            tags$ul(class = "small",
              tags$li(tags$b("SA:"), " 3.99x enrichment (P = 0.016)"),
              tags$li(tags$b("Auxin:"), " 1.78x enrichment (P = 0.022)")
            ),
            hr(),
            p(class = "small text-muted",
              "34 hormone-related genes among 1,567 JAG1 targets, ",
              "from 459 total hormone-annotated genes in soybean.")
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Enrichment Plot",
            card(
              card_body(
                plotlyOutput(ns("enrichment_plot"), height = "450px")
              )
            )
          ),
          tabPanel("Enrichment Table",
            card(
              card_body(
                DT::DTOutput(ns("enrichment_table"))
              )
            )
          ),
          tabPanel("Hormone JAG1 Targets",
            card(
              card_body(
                p(class = "text-muted", "34 hormone-related genes that are JAG1 targets"),
                DT::DTOutput(ns("gene_table"))
              )
            )
          )
        )
      )
    )
  )
}

hormone_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Filtered enrichment data
    filtered_enrichment <- reactive({
      if (is.null(hormone_enrichment)) return(NULL)

      data <- hormone_enrichment %>%
        filter(Hormone %in% input$pathways)

      if (input$sig_only) {
        data <- data %>% filter(P_value < 0.05)
      }

      data
    })

    # Filtered gene data
    filtered_genes <- reactive({
      if (is.null(hormone_genes_jag1)) return(NULL)

      # Filter by selected pathways if Matched_KO column exists
      hormone_genes_jag1
    })

    # Enrichment bar plot
    output$enrichment_plot <- renderPlotly({
      data <- filtered_enrichment()
      if (is.null(data) || nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No enrichment data available"))
      }

      data <- data %>%
        mutate(
          neg_log10_p = -log10(P_value),
          sig_label = ifelse(P_value < 0.05, "*", "")
        )

      p <- ggplot(data, aes(x = reorder(Hormone, Fold_enrichment),
                            y = Fold_enrichment,
                            fill = neg_log10_p,
                            text = paste("Pathway:", Hormone,
                                        "<br>Fold enrichment:", round(Fold_enrichment, 2),
                                        "<br>P-value:", signif(P_value, 3),
                                        "<br>FDR:", signif(FDR, 3),
                                        "<br>Genes in test:", Genes_in_test_set,
                                        "of", Genes_in_background))) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
        geom_text(aes(label = sig_label), hjust = -0.5, size = 6, color = "red") +
        coord_flip() +
        scale_fill_gradient(low = "#FEE0D2", high = "#CB181D",
                            name = "-log10(P)") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme_publication(base_size = 12) +
        theme(panel.grid.major.y = element_blank()) +
        labs(x = NULL, y = "Fold Enrichment",
             title = "Hormone Pathway Enrichment Among JAG1 Targets",
             subtitle = "Dashed line = no enrichment (fold = 1), * = P < 0.05")

      ggplotly(p, tooltip = "text")
    })

    # Enrichment summary table
    output$enrichment_table <- DT::renderDT({
      data <- filtered_enrichment()
      if (is.null(data)) return(NULL)

      data %>%
        mutate(
          Fold_enrichment = round(Fold_enrichment, 2),
          P_value = signif(P_value, 3),
          Odds_ratio = round(Odds_ratio, 2),
          FDR = signif(FDR, 3)
        ) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        ) %>%
        formatStyle(
          "P_value",
          backgroundColor = styleInterval(0.05, c("#FADBD8", "white"))
        )
    })

    # Hormone gene table
    output$gene_table <- DT::renderDT({
      data <- filtered_genes()
      if (is.null(data)) return(NULL)

      display_cols <- c("GeneID", "Confidence_Tier", "Matched_KO",
                        "Mean_logFC_Pairwise", "NarrowvsBroad_TP1_logFC",
                        "NarrowvsBroad_TP1_FDR", "Best_hit_arabi_defline",
                        "KEGG_Arabi_defline")
      display_cols <- display_cols[display_cols %in% names(data)]

      data %>%
        select(all_of(display_cols)) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE),
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
  })
}
