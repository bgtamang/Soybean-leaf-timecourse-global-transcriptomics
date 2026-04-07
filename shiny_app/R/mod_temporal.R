# Temporal module - Expression dynamics over time

temporal_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Temporal Filters"),
          card_body(
            selectInput(ns("pattern_type"), "Expression Pattern:",
              choices = c("All", "Persistent", "TP1_Specific", "Late_Onset", "Variable",
                          "Housekeeping (Machado et al.)" = "Housekeeping"),
              selected = "All"
            ),
            hr(),
            checkboxInput(ns("jag1_only"), "JAG1 targets only", TRUE),
            helpText(class = "small text-muted",
                     "Note: Uncheck for Housekeeping genes"),
            hr(),
            sliderInput(ns("n_genes"), "Number of genes to show:",
              min = 10, max = 100, value = 30, step = 10
            )
          )
        ),
        card(
          card_header("Pattern Summary"),
          card_body(
            uiOutput(ns("pattern_summary"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Temporal Heatmap",
            card(
              card_body(
                plotOutput(ns("temporal_heatmap"), height = "550px")
              )
            )
          ),
          tabPanel("Pattern Distribution",
            card(
              card_body(
                plotlyOutput(ns("pattern_bar"), height = "400px")
              )
            )
          ),
          tabPanel("Gene Trajectories",
            card(
              card_body(
                textInput(ns("gene_search"), "Search gene:", placeholder = "e.g., Glyma.20G116200"),
                plotlyOutput(ns("trajectory_plot"), height = "400px")
              )
            )
          ),
          tabPanel("Temporal Table",
            card(
              card_body(
                DT::DTOutput(ns("temporal_table"))
              )
            )
          )
        )
      )
    )
  )
}

temporal_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Load temporal data (uses pre-loaded global data)
    temporal_data <- reactive({
      if (!is.null(temporal_classification)) {
        return(temporal_classification)
      }
      # Fallback to JAG1 targets with pattern info
      if (!is.null(jag1_targets) && "Pattern" %in% names(jag1_targets)) {
        return(jag1_targets)
      }
      NULL
    })

    # Auto-uncheck JAG1 only when Housekeeping is selected
    observeEvent(input$pattern_type, {
      if (input$pattern_type == "Housekeeping") {
        updateCheckboxInput(session, "jag1_only", value = FALSE)
      }
    })

    # Filtered data
    filtered_data <- reactive({
      # Handle Housekeeping genes separately
      if (input$pattern_type == "Housekeeping") {
        if (is.null(housekeeping_genes)) return(NULL)
        # Return housekeeping genes with pattern label
        hk_data <- housekeeping_genes %>%
          mutate(
            GeneID = Gene,
            Pattern = "Housekeeping",
            Confidence_Tier = NA_character_
          ) %>%
          select(GeneID, Pattern, Mean_Expression, CV_percent, everything(), -Gene)
        return(hk_data)
      }

      data <- temporal_data()
      if (is.null(data)) return(NULL)

      if (input$pattern_type != "All" && "Pattern" %in% names(data)) {
        data <- data %>% filter(Pattern == input$pattern_type)
      }

      if (input$jag1_only) {
        data <- data %>% filter(GeneID %in% jag1_targets$GeneID)
      }

      data
    })

    # Pattern summary
    output$pattern_summary <- renderUI({
      data <- temporal_data()

      # Build pattern counts
      pattern_list <- list()

      if (!is.null(data) && "Pattern" %in% names(data)) {
        pattern_counts <- data %>%
          group_by(Pattern) %>%
          summarise(n = n(), .groups = "drop")
        for (i in seq_len(nrow(pattern_counts))) {
          pattern_list[[pattern_counts$Pattern[i]]] <- pattern_counts$n[i]
        }
      }

      # Add housekeeping count
      if (!is.null(housekeeping_genes)) {
        pattern_list[["Housekeeping"]] <- nrow(housekeeping_genes)
      }

      if (length(pattern_list) == 0) {
        return(p("Pattern data not available"))
      }

      tags$table(class = "table table-sm",
        lapply(names(pattern_list), function(pat) {
          tags$tr(
            tags$td(pat),
            tags$td(pattern_list[[pat]])
          )
        })
      )
    })

    # Temporal heatmap
    output$temporal_heatmap <- renderPlot({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) {
        plot.new()
        text(0.5, 0.5, "No data matches filters", cex = 1.5)
        return()
      }

      # Handle Housekeeping genes differently - show expression across samples
      if (input$pattern_type == "Housekeeping") {
        expr <- expression_matrix
        if (is.null(expr)) {
          plot.new()
          text(0.5, 0.5, "Expression matrix not available", cex = 1.5)
          return()
        }

        # Get top N housekeeping genes
        hk_genes <- head(data$GeneID, input$n_genes)

        # Direct match - gene-level IDs should match directly
        valid_genes <- hk_genes[hk_genes %in% rownames(expr)]

        if (length(valid_genes) < 2) {
          plot.new()
          text(0.5, 0.5, "Not enough housekeeping genes found in expression matrix", cex = 1.5)
          return()
        }

        # Get expression matrix subset
        mat <- as.matrix(expr[valid_genes, ])

        # Z-score scale for visualization
        mat <- t(scale(t(mat)))

        # Row annotation with CV%
        cv_vals <- data$CV_percent[match(valid_genes, data$GeneID)]
        row_anno <- data.frame(
          CV_pct = round(cv_vals, 1),
          row.names = valid_genes
        )

        p <- gg_heatmap(
          mat,
          title = paste("Housekeeping Gene Expression (n =", length(valid_genes), "most stable)"),
          row_anno = row_anno,
          cluster_rows = FALSE,  # Keep ordered by stability
          fontsize_row = if(length(valid_genes) <= 30) 8 else 6,
          legend_title = "Z-score"
        )
        print(p)
        return()
      }

      # Get FC columns for other patterns
      fc_cols <- grep("^FC_TP[1-5]$", names(data), value = TRUE)
      if (length(fc_cols) == 0) {
        fc_cols <- c("FC_TP1", "FC_TP2", "FC_TP3", "FC_TP4", "FC_TP5")
        fc_cols <- fc_cols[fc_cols %in% names(data)]
      }

      if (length(fc_cols) == 0) {
        plot.new()
        text(0.5, 0.5, "Temporal FC columns not found", cex = 1.5)
        return()
      }

      # Limit genes
      plot_data <- data %>%
        head(input$n_genes)

      mat <- as.matrix(plot_data[, fc_cols])
      rownames(mat) <- plot_data$GeneID

      # Annotation
      row_anno <- NULL
      if ("Pattern" %in% names(plot_data)) {
        row_anno <- data.frame(Pattern = plot_data$Pattern, row.names = plot_data$GeneID)
      }

      # Publication-quality temporal heatmap
      gg_heatmap(
        mat,
        title = paste("Temporal Expression Dynamics (n =", nrow(plot_data), "genes)"),
        row_anno = row_anno,
        cluster_rows = TRUE,
        fontsize_row = 7,
        legend_title = expression(log[2]*FC)
      )
    })

    # Pattern distribution bar
    output$pattern_bar <- renderPlotly({
      data <- temporal_data()
      if (is.null(data) || !"Pattern" %in% names(data)) {
        return(plotly_empty() %>% layout(title = "Pattern data not available"))
      }

      pattern_counts <- data %>%
        group_by(Pattern) %>%
        summarise(Count = n(), .groups = "drop") %>%
        mutate(Pattern = factor(Pattern))

      p <- ggplot(pattern_counts, aes(x = reorder(Pattern, Count), y = Count, fill = Pattern)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        geom_text(aes(label = Count), hjust = -0.3, size = 4, fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        coord_flip() +
        theme_publication(base_size = 12) +
        theme(legend.position = "none", panel.grid.major.y = element_blank()) +
        labs(x = NULL, y = "Number of Genes", title = "Genes by Temporal Pattern")

      ggplotly(p)
    })

    # Gene trajectory plot
    output$trajectory_plot <- renderPlotly({
      req(input$gene_search)
      search_term <- trimws(input$gene_search)

      if (nchar(search_term) < 5) {
        return(plotly_empty() %>% layout(title = "Enter a gene ID"))
      }

      # Find in JAG1 targets
      gene_data <- jag1_targets %>%
        filter(grepl(search_term, GeneID, ignore.case = TRUE))

      if (nrow(gene_data) == 0) {
        return(plotly_empty() %>% layout(title = paste("Gene", search_term, "not found")))
      }

      gene_data <- gene_data[1, ]

      # Get FC values
      fc_cols <- c("FC_TP1", "FC_TP2", "FC_TP3", "FC_TP4", "FC_TP5")
      fc_cols <- fc_cols[fc_cols %in% names(gene_data)]

      if (length(fc_cols) == 0) {
        return(plotly_empty() %>% layout(title = "Temporal data not available for this gene"))
      }

      fc_values <- as.numeric(gene_data[, fc_cols])
      timepoints <- gsub("FC_", "", fc_cols)

      plot_df <- data.frame(
        Timepoint = factor(timepoints, levels = timepoints),
        logFC = fc_values
      )

      # Publication-quality trajectory plot
      p <- ggplot(plot_df, aes(x = Timepoint, y = logFC, group = 1)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
        geom_line(color = "#C0392B", linewidth = 1.5) +
        geom_point(color = "#C0392B", size = 5, shape = 21, fill = "#E74C3C", stroke = 1.2) +
        geom_area(alpha = 0.15, fill = "#E74C3C") +
        theme_publication(base_size = 13) +
        labs(
          title = paste("Expression Trajectory:", gene_data$GeneID),
          subtitle = if (!is.na(gene_data$Best_hit_arabi_name[1])) gene_data$Best_hit_arabi_name[1] else NULL,
          x = "Developmental Timepoint",
          y = "log\u2082 Fold Change (Narrow vs Broad)"
        )

      ggplotly(p)
    })

    # Temporal table
    output$temporal_table <- DT::renderDT({
      data <- filtered_data()
      if (is.null(data)) return(NULL)

      # Select relevant columns based on pattern type
      if (input$pattern_type == "Housekeeping") {
        display_cols <- c("GeneID", "Pattern", "Mean_Expression", "SD", "CV_percent")
      } else {
        display_cols <- c("GeneID", "Pattern", "Confidence_Tier",
                          "FC_TP1", "FC_TP2", "FC_TP3", "FC_TP4", "FC_TP5",
                          "Best_hit_arabi_name")
      }
      display_cols <- display_cols[display_cols %in% names(data)]

      data %>%
        select(all_of(display_cols)) %>%
        datatable(
          filter = "top",
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        )
    })

    # Handle table selection
    observeEvent(input$temporal_table_rows_selected, {
      data <- filtered_data()
      if (!is.null(data) && !is.null(input$temporal_table_rows_selected)) {
        gene <- data$GeneID[input$temporal_table_rows_selected]
        selected_gene(gene)
        updateTextInput(session, "gene_search", value = gene)
      }
    })
  })
}
