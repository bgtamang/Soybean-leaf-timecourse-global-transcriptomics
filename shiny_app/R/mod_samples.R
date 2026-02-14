# Samples module - PCA and QC visualization

samples_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Controls"),
          card_body(
            selectInput(ns("color_by"), "Color by:",
              choices = c("Timepoint", "Leaf_type", "Line", "Batch"),
              selected = "Timepoint"
            ),
            selectInput(ns("shape_by"), "Shape by:",
              choices = c("None", "Leaf_type", "Timepoint", "Line", "Batch"),
              selected = "Leaf_type"
            ),
            hr(),
            checkboxGroupInput(ns("timepoints"), "Timepoints:",
              choices = c("TP0", "TP1", "TP2", "TP3", "TP4"),
              selected = c("TP0", "TP1", "TP2", "TP3", "TP4")
            ),
            hr(),
            checkboxGroupInput(ns("leaf_types"), "Leaf Types:",
              choices = c("Broad", "Narrow"),
              selected = c("Broad", "Narrow")
            )
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("PCA Plot",
            card(
              card_body(
                plotlyOutput(ns("pca_plot"), height = "500px")
              )
            )
          ),
          tabPanel("Sample Table",
            card(
              card_body(
                DT::DTOutput(ns("sample_table"))
              )
            )
          ),
          tabPanel("QC Summary",
            card(
              card_body(
                fluidRow(
                  column(6, plotOutput(ns("lib_size_plot"), height = "350px")),
                  column(6, plotOutput(ns("detected_genes_plot"), height = "350px"))
                )
              )
            )
          )
        )
      )
    )
  )
}

samples_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Filter samples
    filtered_samples <- reactive({
      experimental_design %>%
        filter(
          Timepoint %in% input$timepoints,
          Leaf_type %in% input$leaf_types
        )
    })

    # PCA plot
    output$pca_plot <- renderPlotly({
      pca_data <- load_pca_data()

      if (is.null(pca_data)) {
        return(plotly_empty() %>% layout(title = "PCA data not available"))
      }

      # Merge with sample info - Sample column matches experimental design Sample
      plot_data <- pca_data$coords %>%
        left_join(experimental_design, by = "Sample")

      plot_data <- plot_data %>%
        filter(
          Timepoint %in% input$timepoints,
          Leaf_type %in% input$leaf_types
        )

      if (nrow(plot_data) == 0) {
        return(plotly_empty() %>% layout(title = "No samples match filters"))
      }

      # Leaf type shapes matching publication (circle=Broad, triangle=Narrow)
      leaf_shapes <- c("Broad" = 21, "Narrow" = 24)

      # Build publication-quality PCA plot (matches Fig 2A)
      p <- ggplot(plot_data, aes(x = PC1, y = PC2,
                                  text = paste("Sample:", Sample,
                                             "<br>Line:", Line,
                                             "<br>Timepoint:", Timepoint)))

      if (input$shape_by != "None") {
        p <- p + geom_point(aes(fill = .data[[input$color_by]],
                                shape = .data[[input$shape_by]]),
                            size = 3.5, stroke = 0.3, color = "black", alpha = 0.9)
      } else {
        p <- p + geom_point(aes(fill = .data[[input$color_by]]),
                            shape = 21, size = 3.5, stroke = 0.3,
                            color = "black", alpha = 0.9)
      }

      p <- p +
        theme_publication(base_size = 13) +
        labs(
          x = paste0("PC1 (", round(pca_data$variance[1], 1), "%)"),
          y = paste0("PC2 (", round(pca_data$variance[2], 1), "%)"),
          fill = gsub("_", " ", input$color_by),
          shape = if (input$shape_by != "None") gsub("_", " ", input$shape_by) else NULL
        ) +
        theme(
          aspect.ratio = 1,
          legend.position = "right",
          legend.box = "vertical"
        )

      # Color scale
      if (input$color_by == "Leaf_type") {
        p <- p + scale_fill_manual(values = leaftype_colors, name = "Leaf Type")
      } else if (input$color_by == "Timepoint") {
        p <- p + scale_fill_manual(values = timepoint_colors, name = "Timepoint")
      }

      # Shape scale
      if (input$shape_by == "Leaf_type") {
        p <- p + scale_shape_manual(values = leaf_shapes, name = "Leaf Type",
                                    labels = c("Broad (Ln)", "Narrow (ln)"))
      }

      # Legend overrides for filled shapes
      p <- p + guides(
        fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 3)),
        shape = guide_legend(order = 2, override.aes = list(size = 3, fill = "gray50"))
      )

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    })

    # Sample table
    output$sample_table <- DT::renderDT({
      datatable(
        filtered_samples(),
        filter = "top",
        options = list(
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })

    # Library size plot
    output$lib_size_plot <- renderPlot({
      if (is.null(qc_metrics)) {
        plot.new()
        text(0.5, 0.5, "QC data not available", cex = 1.5)
        return()
      }

      qc_data <- qc_metrics %>%
        filter(
          Timepoint %in% input$timepoints,
          Leaf_type %in% input$leaf_types
        )

      if (!"Library_size" %in% names(qc_data)) {
        plot.new()
        text(0.5, 0.5, "Library size data not available", cex = 1.5)
        return()
      }

      # Use Sample column if Sample_name doesn't exist
      sample_col <- if ("Sample_name" %in% names(qc_data)) "Sample_name" else "Sample"
      qc_data$display_name <- qc_data[[sample_col]]

      ggplot(qc_data, aes(x = reorder(display_name, Library_size), y = Library_size / 1e6, fill = Leaf_type)) +
        geom_col() +
        scale_fill_manual(values = leaftype_colors) +
        coord_flip() +
        theme_minimal() +
        labs(x = NULL, y = "Library Size (millions)", fill = "Leaf Type") +
        theme(axis.text.y = element_text(size = 7))
    })

    # Detected genes plot
    output$detected_genes_plot <- renderPlot({
      # Load expression data
      expr <- tryCatch(load_expression(), error = function(e) NULL)

      if (is.null(expr)) {
        plot.new()
        text(0.5, 0.5, "Expression data not available", cex = 1.2)
        return()
      }

      # Count detected genes per sample (logCPM > 0, i.e., CPM > 1)
      detected <- colSums(expr > 0, na.rm = TRUE)

      det_df <- data.frame(
        Sample = names(detected),
        Detected = as.numeric(detected)
      )

      # Merge with experimental design
      det_df <- det_df %>%
        left_join(experimental_design, by = "Sample")

      det_df <- det_df %>%
        filter(
          Timepoint %in% input$timepoints,
          Leaf_type %in% input$leaf_types
        )

      if (nrow(det_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No samples match filters", cex = 1.2)
        return()
      }

      ggplot(det_df, aes(x = reorder(Sample, Detected), y = Detected / 1000, fill = Leaf_type)) +
        geom_col() +
        scale_fill_manual(values = leaftype_colors) +
        coord_flip() +
        theme_minimal() +
        labs(x = NULL, y = "Detected Genes (thousands, CPM > 1)", fill = "Leaf Type") +
        theme(axis.text.y = element_text(size = 7))
    })
  })
}
