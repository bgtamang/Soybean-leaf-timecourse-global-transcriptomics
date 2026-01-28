# WGCNA module - Co-expression network analysis

wgcna_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Module Selection"),
          card_body(
            selectInput(ns("selected_module"), "Select Module:",
              choices = NULL  # Populated in server
            ),
            hr(),
            checkboxInput(ns("show_jag1_enriched"), "JAG1-enriched modules only", FALSE)
          )
        ),
        card(
          card_header("Module Info"),
          card_body(
            uiOutput(ns("module_info"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Module-Trait Heatmap",
            card(
              card_body(
                plotOutput(ns("trait_heatmap"), height = "500px")
              )
            )
          ),
          tabPanel("Module Summary",
            card(
              card_body(
                DT::DTOutput(ns("module_table"))
              )
            )
          ),
          tabPanel("JAG1 Enrichment",
            card(
              card_body(
                plotlyOutput(ns("enrichment_plot"), height = "400px"),
                hr(),
                DT::DTOutput(ns("enrichment_table"))
              )
            )
          ),
          tabPanel("Hub Genes",
            card(
              card_body(
                DT::DTOutput(ns("hub_table"))
              )
            )
          )
        )
      )
    )
  )
}

wgcna_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Update module choices
    observe({
      modules <- module_sizes$Module
      if (!is.null(modules)) {
        # Format nicely
        choices <- setNames(
          modules,
          paste0(modules, " (n=", module_sizes$Size, ")")
        )
        updateSelectInput(session, "selected_module", choices = choices)
      }
    })

    # Module info panel
    output$module_info <- renderUI({
      req(input$selected_module)

      mod_info <- module_sizes %>% filter(Module == input$selected_module)
      trait_info <- module_traits %>% filter(Module == input$selected_module)

      if (nrow(mod_info) == 0) {
        return(p("Module not found"))
      }

      tagList(
        p(strong("Size:"), mod_info$Size, "genes"),
        if (nrow(trait_info) > 0 && "Cor_Narrow" %in% names(trait_info)) {
          tagList(
            hr(),
            p(strong("Leaf Type Correlation:")),
            p("Narrow: ", round(trait_info$Cor_Narrow, 3)),
            p("Broad: ", round(trait_info$Cor_Broad, 3))
          )
        }
      )
    })

    # Module-trait heatmap
    output$trait_heatmap <- renderPlot({
      if (is.null(module_traits) || nrow(module_traits) == 0) {
        plot.new()
        text(0.5, 0.5, "Module-trait data not available", cex = 1.5)
        return()
      }

      # Select correlation columns
      cor_cols <- grep("^Cor_", names(module_traits), value = TRUE)
      if (length(cor_cols) == 0) {
        plot.new()
        text(0.5, 0.5, "Correlation data not available", cex = 1.5)
        return()
      }

      mat <- as.matrix(module_traits[, cor_cols])
      rownames(mat) <- module_traits$Module
      colnames(mat) <- gsub("Cor_", "", colnames(mat))

      # Publication-quality heatmap
      gg_heatmap(
        mat,
        title = "Module-Trait Correlation Heatmap",
        cluster_rows = TRUE,
        fontsize_row = 9,
        legend_title = "Correlation\n(r)"
      )
    })

    # Module summary table
    output$module_table <- DT::renderDT({
      if (is.null(module_traits)) return(NULL)

      display_data <- module_traits %>%
        left_join(module_sizes, by = "Module") %>%
        select(Module, Size, everything())

      datatable(display_data,
        filter = "top",
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      ) %>%
        formatRound(columns = grep("Cor_", names(display_data), value = TRUE), digits = 3)
    })

    # JAG1 enrichment plot
    output$enrichment_plot <- renderPlotly({
      if (is.null(jag1_module_enrichment)) {
        return(plotly_empty() %>% layout(title = "Enrichment data not available"))
      }

      enrich_data <- jag1_module_enrichment

      if (nrow(enrich_data) == 0) {
        return(plotly_empty() %>% layout(title = "No enrichment data"))
      }

      # Assume columns: Module, Enrichment, Pvalue or similar
      if ("Fold_Enrichment" %in% names(enrich_data)) {
        p <- ggplot(enrich_data, aes(x = reorder(Module, Fold_Enrichment),
                                     y = Fold_Enrichment,
                                     fill = Fold_Enrichment > 1)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#27AE60")) +
          theme_minimal() +
          theme(legend.position = "none") +
          labs(x = NULL, y = "Fold Enrichment", title = "JAG1 Target Enrichment by Module")

        ggplotly(p)
      } else {
        plotly_empty() %>% layout(title = "Enrichment metric not found")
      }
    })

    # Enrichment table
    output$enrichment_table <- DT::renderDT({
      if (is.null(jag1_module_enrichment)) return(NULL)
      datatable(jag1_module_enrichment, options = list(pageLength = 10), rownames = FALSE)
    })

    # Hub genes table
    output$hub_table <- DT::renderDT({
      if (is.null(hub_genes)) return(NULL)
      req(input$selected_module)

      hub_filtered <- hub_genes %>%
        filter(Module == input$selected_module)

      if (nrow(hub_filtered) == 0) {
        return(NULL)
      }

      datatable(hub_filtered,
        options = list(pageLength = 20),
        rownames = FALSE
      )
    })
  })
}
