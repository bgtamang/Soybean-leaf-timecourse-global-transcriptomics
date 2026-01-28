# Phenotype module - Leaf trait correlations

phenotype_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Filter Options"),
          card_body(
            sliderInput(ns("cor_threshold"), "Min |Correlation|:",
              min = 0, max = 0.9, value = 0.3, step = 0.1
            ),
            selectInput(ns("direction"), "Direction:",
              choices = c("All", "Positive", "Negative"),
              selected = "All"
            ),
            hr(),
            checkboxInput(ns("jag1_only"), "JAG1 targets only", FALSE)
          )
        ),
        card(
          card_header("Phenotype Summary"),
          card_body(
            uiOutput(ns("pheno_summary"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("L/W Ratio Comparison",
            card(
              card_body(
                plotlyOutput(ns("lw_plot"), height = "400px")
              )
            )
          ),
          tabPanel("Gene-Phenotype Correlations",
            card(
              card_body(
                plotlyOutput(ns("cor_scatter"), height = "400px"),
                hr(),
                DT::DTOutput(ns("cor_table"))
              )
            )
          ),
          tabPanel("Module-Phenotype",
            card(
              card_body(
                plotOutput(ns("module_pheno_heatmap"), height = "400px")
              )
            )
          )
        )
      )
    )
  )
}

phenotype_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Load phenotype correlation data (uses pre-loaded global data)
    pheno_cors <- reactive({
      gene_phenotype_cors
    })

    # Filtered correlations
    filtered_cors <- reactive({
      cors <- pheno_cors()
      if (is.null(cors)) return(NULL)

      # Find correlation column (Cor_LW_Ratio is the actual column name)
      cor_col <- intersect(c("Cor_LW_Ratio", "Correlation", "Cor", "r", "cor"), names(cors))[1]
      if (is.na(cor_col)) return(cors)

      result <- cors %>%
        filter(abs(.data[[cor_col]]) >= input$cor_threshold)

      if (input$direction == "Positive") {
        result <- result %>% filter(.data[[cor_col]] > 0)
      } else if (input$direction == "Negative") {
        result <- result %>% filter(.data[[cor_col]] < 0)
      }

      if (input$jag1_only) {
        result <- result %>% filter(GeneID %in% jag1_targets$GeneID)
      }

      result
    })

    # Phenotype summary (uses pre-loaded global data)
    output$pheno_summary <- renderUI({
      if (is.null(phenotype_traits)) {
        return(p("Phenotype data not available"))
      }

      pheno <- phenotype_traits

      tagList(
        h5("Leaf L/W Ratio"),
        tags$table(class = "table table-sm",
          tags$tr(
            tags$th(""),
            tags$th("V1"),
            tags$th("V2")
          ),
          tags$tr(
            tags$td("Broad"),
            tags$td("1.26"),
            tags$td("1.24")
          ),
          tags$tr(
            tags$td("Narrow"),
            tags$td("1.77"),
            tags$td("2.53")
          )
        ),
        hr(),
        p(class = "small text-muted", "Values represent mean L/W ratio per leaf type")
      )
    })

    # L/W ratio plot (publication quality)
    output$lw_plot <- renderPlotly({
      lw_data <- data.frame(
        Leaf_Type = factor(rep(c("Broad", "Narrow"), each = 2), levels = c("Broad", "Narrow")),
        Stage = factor(rep(c("V1", "V2"), 2), levels = c("V1", "V2")),
        LW_Ratio = c(1.26, 1.24, 1.77, 2.53)
      )

      p <- ggplot(lw_data, aes(x = Stage, y = LW_Ratio, fill = Leaf_Type, group = Leaf_Type)) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7,
                 color = "black", linewidth = 0.4) +
        geom_text(aes(label = LW_Ratio), position = position_dodge(width = 0.75),
                  vjust = -0.5, size = 4, fontface = "bold") +
        scale_fill_manual(values = leaftype_colors, name = "Leaf Type") +
        scale_y_continuous(limits = c(0, 3.2), expand = expansion(mult = c(0, 0.12))) +
        theme_publication(base_size = 12) +
        theme(
          panel.grid.major.x = element_blank(),
          legend.position = "top"
        ) +
        labs(
          title = "Leaf Length/Width Ratio by Developmental Stage",
          x = "Developmental Stage",
          y = "L/W Ratio"
        )

      ggplotly(p) %>% layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.1))
    })

    # Correlation scatter
    output$cor_scatter <- renderPlotly({
      cors <- filtered_cors()
      if (is.null(cors) || nrow(cors) == 0) {
        return(plotly_empty() %>% layout(title = "No correlations meet threshold"))
      }

      # Find correlation column (Cor_LW_Ratio is the actual column name)
      cor_col <- intersect(c("Cor_LW_Ratio", "Correlation", "Cor", "r", "cor"), names(cors))[1]
      if (is.na(cor_col)) {
        return(plotly_empty() %>% layout(title = "Correlation column not found"))
      }

      # Check for p-value column (Pval_LW_Ratio or FDR_LW_Ratio)
      p_col <- intersect(c("Pval_LW_Ratio", "FDR_LW_Ratio", "P_value", "Pvalue", "p.value", "P"), names(cors))

      if (length(p_col) > 0) {
        p <- ggplot(cors, aes(x = .data[[cor_col]], y = -log10(.data[[p_col[1]]]),
                              text = paste("Gene:", GeneID,
                                          "<br>r:", round(.data[[cor_col]], 3),
                                          "<br>-log10(P):", round(-log10(.data[[p_col[1]]]), 2)))) +
          geom_point(alpha = 0.7, color = "#2171B5", size = 2.5)
      } else {
        p <- ggplot(cors, aes(x = .data[[cor_col]], y = seq_len(nrow(cors)),
                              text = paste("Gene:", GeneID,
                                          "<br>r:", round(.data[[cor_col]], 3)))) +
          geom_point(alpha = 0.7, color = "#2171B5", size = 2.5)
      }

      p <- p +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
        theme_publication(base_size = 12) +
        labs(
          title = "Gene-Phenotype Correlations",
          subtitle = paste(nrow(cors), "genes shown"),
          x = "Correlation with L/W Ratio (r)",
          y = ifelse(length(p_col) > 0, "-log\u2081\u2080(P-value)", "Rank")
        )

      ggplotly(p, tooltip = "text")
    })

    # Correlation table
    output$cor_table <- DT::renderDT({
      cors <- filtered_cors()
      if (is.null(cors)) return(NULL)

      datatable(cors,
        filter = "top",
        selection = "single",
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    # Handle table selection
    observeEvent(input$cor_table_rows_selected, {
      cors <- filtered_cors()
      if (!is.null(cors) && !is.null(input$cor_table_rows_selected)) {
        gene <- cors$GeneID[input$cor_table_rows_selected]
        selected_gene(gene)
      }
    })

    # Module-phenotype heatmap
    output$module_pheno_heatmap <- renderPlot({
      if (is.null(module_traits) || nrow(module_traits) == 0) {
        plot.new()
        text(0.5, 0.5, "Module-trait data not available", cex = 1.5)
        return()
      }

      # Select phenotype-related columns
      pheno_cols <- grep("(Narrow|Broad|LW|Ratio)", names(module_traits), value = TRUE)
      pheno_cols <- intersect(pheno_cols, grep("^Cor_", names(module_traits), value = TRUE))

      if (length(pheno_cols) == 0) {
        pheno_cols <- grep("^Cor_", names(module_traits), value = TRUE)[1:2]
      }

      if (length(pheno_cols) == 0) {
        plot.new()
        text(0.5, 0.5, "Phenotype correlation columns not found", cex = 1.5)
        return()
      }

      mat <- as.matrix(module_traits[, pheno_cols, drop = FALSE])
      rownames(mat) <- module_traits$Module
      colnames(mat) <- gsub("Cor_", "", colnames(mat))

      # Publication-quality heatmap
      gg_heatmap(
        mat,
        title = "Module-Phenotype Correlation Heatmap",
        cluster_rows = TRUE,
        fontsize_row = 9,
        legend_title = "Correlation\n(r)"
      )
    })
  })
}
