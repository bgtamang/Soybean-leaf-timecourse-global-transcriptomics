# GO Enrichment module

go_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Enrichment Settings"),
          card_body(
            selectInput(ns("target_set"), "Target Set:",
              choices = c("All JAG1 Targets" = "all",
                          "Gold Tier" = "gold",
                          "Silver Tier" = "silver",
                          "Bronze Tier" = "bronze")
            ),
            hr(),
            sliderInput(ns("fdr_threshold"), "FDR Threshold:",
              min = 0.001, max = 0.1, value = 0.05, step = 0.005
            ),
            sliderInput(ns("top_n"), "Top N Terms:",
              min = 10, max = 50, value = 20, step = 5
            )
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Dot Plot",
            card(
              card_body(
                plotOutput(ns("dot_plot"), height = "600px")
              )
            )
          ),
          tabPanel("GO Table",
            card(
              card_body(
                DT::DTOutput(ns("go_table")),
                hr(),
                downloadButton(ns("download_go"), "Download GO Results", class = "btn-primary")
              )
            )
          ),
          tabPanel("Category Summary",
            card(
              card_body(
                plotlyOutput(ns("category_plot"), height = "400px")
              )
            )
          )
        )
      )
    )
  )
}

go_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Get GO data based on selection (uses pre-loaded global data)
    go_data <- reactive({
      get_go_data(input$target_set)
    })

    # Filtered GO results
    filtered_go <- reactive({
      go <- go_data()
      if (is.null(go)) return(NULL)

      # Filter by FDR
      fdr_col <- intersect(c("FDR", "p.adjust", "qvalue", "adj.P.Val"), names(go))
      if (length(fdr_col) > 0) {
        go <- go %>% filter(.data[[fdr_col[1]]] < input$fdr_threshold)
      }

      # Sort and limit
      if ("Fold_Enrichment" %in% names(go)) {
        go <- go %>% arrange(desc(Fold_Enrichment))
      } else if ("Count" %in% names(go)) {
        go <- go %>% arrange(desc(Count))
      }

      head(go, input$top_n)
    })

    # Dot plot
    output$dot_plot <- renderPlot({
      go <- filtered_go()
      if (is.null(go) || nrow(go) == 0) {
        plot.new()
        text(0.5, 0.5, "No significant GO terms found", cex = 1.5)
        return()
      }

      # Determine columns
      term_col <- intersect(c("Description", "Term", "GO_term", "ID"), names(go))[1]
      count_col <- intersect(c("Count", "GeneCount", "geneCount"), names(go))[1]
      fdr_col <- intersect(c("FDR", "p.adjust", "qvalue"), names(go))[1]
      enrich_col <- intersect(c("Fold_Enrichment", "GeneRatio", "enrichmentScore"), names(go))[1]

      if (is.na(term_col) || is.na(count_col)) {
        plot.new()
        text(0.5, 0.5, "Required columns not found in GO data", cex = 1.5)
        return()
      }

      # Truncate long terms for display
      go$Term_short <- substr(go[[term_col]], 1, 45)
      go$Term_short <- ifelse(nchar(go[[term_col]]) > 45,
                              paste0(go$Term_short, "..."),
                              go$Term_short)

      # Publication-quality dot plot
      p <- ggplot(go, aes(x = .data[[count_col]],
                          y = reorder(Term_short, .data[[count_col]])))

      if (!is.na(fdr_col) && !is.na(enrich_col)) {
        p <- p + geom_point(aes(size = .data[[count_col]], color = -log10(.data[[fdr_col]])),
                            alpha = 0.85)
      } else {
        p <- p + geom_point(aes(size = .data[[count_col]]), color = "#2171B5", alpha = 0.85)
      }

      p <- p +
        scale_color_gradientn(
          colors = c("#4575B4", "#74ADD1", "#FEE090", "#F46D43", "#D73027"),
          name = expression(-log[10](FDR)),
          guide = guide_colorbar(barwidth = 1, barheight = 6, frame.colour = "black")
        ) +
        scale_size_continuous(name = "Gene\nCount", range = c(3, 10),
                              guide = guide_legend(override.aes = list(alpha = 1))) +
        theme_publication(base_size = 12) +
        theme(
          axis.text.y = element_text(size = 10),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.4)
        ) +
        labs(
          title = paste("GO Enrichment Analysis"),
          subtitle = switch(input$target_set,
                           "all" = "All JAG1 Targets",
                           "gold" = "Gold Tier Targets",
                           "silver" = "Silver Tier Targets",
                           "bronze" = "Bronze Tier Targets"),
          x = "Gene Count",
          y = NULL
        )

      p
    })

    # GO table
    output$go_table <- DT::renderDT({
      go <- filtered_go()
      if (is.null(go)) return(NULL)

      datatable(go,
        filter = "top",
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    # Category summary - use Module if available, otherwise Ontology
    output$category_plot <- renderPlotly({
      go <- go_data()
      if (is.null(go) || nrow(go) == 0) {
        return(plotly_empty() %>% layout(title = "No data available"))
      }

      # Check for Module or Ontology/Category column
      cat_col <- intersect(c("Module", "Ontology", "Category", "ONTOLOGY"), names(go))

      if (length(cat_col) == 0) {
        return(plotly_empty() %>% layout(title = "Category data not available"))
      }

      cat_summary <- go %>%
        group_by(.data[[cat_col[1]]]) %>%
        summarise(Count = n(), .groups = "drop") %>%
        arrange(desc(Count)) %>%
        head(15)  # Top 15 categories

      p <- ggplot(cat_summary, aes(x = reorder(.data[[cat_col[1]]], Count), y = Count, fill = .data[[cat_col[1]]])) +
        geom_col(width = 0.7) +
        geom_text(aes(label = Count), hjust = -0.2) +
        coord_flip() +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = NULL, y = "Number of Terms",
             title = paste("GO Terms by", ifelse(cat_col[1] == "Module", "Module", "Category")))

      ggplotly(p)
    })

    # Download handler
    output$download_go <- downloadHandler(
      filename = function() {
        paste0("GO_enrichment_", input$target_set, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(go_data(), file, row.names = FALSE)
      }
    )
  })
}
