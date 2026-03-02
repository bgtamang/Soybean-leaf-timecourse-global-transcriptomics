# Single-Cell Integration module - snRNA-seq spatial validation (Figure S6)

single_cell_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        card(
          card_body(class = "bg-light",
            fluidRow(
              column(3,
                value_box(
                  title = "LP1 JAG1+ Cells",
                  value = "35.1%",
                  theme = "danger",
                  showcase = icon("bullseye")
                )
              ),
              column(3,
                value_box(
                  title = "LP1 Enrichment",
                  value = "6.7x",
                  theme = "warning",
                  showcase = icon("chart-line")
                )
              ),
              column(3,
                value_box(
                  title = "LP1 Cluster Cells",
                  value = "515",
                  theme = "info",
                  showcase = icon("braille")
                )
              ),
              column(3,
                value_box(
                  title = "Data Source",
                  value = "Fan et al. 2025",
                  theme = "primary",
                  showcase = icon("database")
                )
              )
            )
          )
        )
      )
    ),
    fluidRow(
      column(3,
        card(
          card_header("Visualization Options"),
          card_body(
            selectInput(ns("view"), "Select View:",
              choices = c(
                "JAG1 by Cluster" = "jag_bar",
                "Component Heatmap" = "heatmap",
                "Co-localization" = "coloc",
                "Target Scores" = "target_score"
              ),
              selected = "jag_bar"
            ),
            hr(),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'coloc'", ns("view")),
              checkboxGroupInput(ns("coloc_tiers"), "Filter Tiers:",
                choices = c("Gold", "Silver", "Bronze"),
                selected = c("Gold", "Silver", "Bronze")
              )
            )
          )
        ),
        card(
          card_header("LP1 Component Summary"),
          card_body(
            uiOutput(ns("lp1_stats"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Plot",
            card(
              card_body(
                plotlyOutput(ns("main_plot"), height = "500px")
              )
            )
          ),
          tabPanel("Cluster Data",
            card(
              card_body(
                DT::DTOutput(ns("cluster_table"))
              )
            )
          ),
          tabPanel("JAG1/JAG2 Per Cluster",
            card(
              card_body(
                DT::DTOutput(ns("jag_table"))
              )
            )
          )
        )
      )
    )
  )
}

single_cell_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # LP1 stats sidebar
    output$lp1_stats <- renderUI({
      if (is.null(sc_cluster_summary)) {
        return(p("Single-cell data not available"))
      }

      lp1 <- sc_cluster_summary %>% filter(grepl("LP1|Leaf [Pp]rimordium\\s*1", cluster, ignore.case = TRUE))

      if (nrow(lp1) == 0) {
        # Try first row or row with highest jag1_pct
        lp1 <- sc_cluster_summary %>% slice_max(jag1_pct, n = 1)
      }

      if (nrow(lp1) == 0) return(p("LP1 data not found"))

      lp1 <- lp1[1, ]

      component_cols <- c("jag1_pct", "tpr_mean_pct", "hda_mean_pct",
                           "krp_mean_pct", "cyclin_target_mean_pct", "cdk_mean_pct")
      component_labels <- c("JAG1", "TPR", "HDAC", "KRP", "Cyclin (targets)", "CDK")

      available <- component_cols[component_cols %in% names(lp1)]
      labels <- component_labels[component_cols %in% names(lp1)]

      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th("Component"), tags$th("% in LP1"))),
        lapply(seq_along(available), function(i) {
          val <- round(as.numeric(lp1[[available[i]]]), 1)
          tags$tr(
            tags$td(labels[i]),
            tags$td(tags$b(paste0(val, "%")))
          )
        })
      )
    })

    # Main plot dispatch
    output$main_plot <- renderPlotly({
      switch(input$view,
        "jag_bar"      = render_jag_bar(),
        "heatmap"      = render_component_heatmap(),
        "coloc"        = render_coloc_scatter(),
        "target_score" = render_target_score()
      )
    })

    # JAG1 % expressing across clusters
    render_jag_bar <- function() {
      if (is.null(sc_cluster_summary)) {
        return(plotly_empty() %>% layout(title = "Single-cell data not available"))
      }

      data <- sc_cluster_summary %>%
        mutate(
          is_lp1 = grepl("LP1|Leaf [Pp]rimordium\\s*1", cluster, ignore.case = TRUE),
          highlight = ifelse(is_lp1, "LP1", "Other")
        )

      # If no LP1 label found, highlight the max jag1_pct cluster
      if (!any(data$is_lp1)) {
        max_idx <- which.max(data$jag1_pct)
        data$highlight[max_idx] <- "LP1"
        data$is_lp1[max_idx] <- TRUE
      }

      p <- ggplot(data, aes(x = reorder(cluster, jag1_pct), y = jag1_pct,
                            fill = highlight,
                            text = paste("Cluster:", cluster,
                                        "<br>JAG1+ cells:", round(jag1_pct, 1), "%",
                                        "<br>N cells:", n_cells))) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        coord_flip() +
        scale_fill_manual(values = c("LP1" = "#E74C3C", "Other" = "#BDC3C7"),
                          name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_publication(base_size = 12) +
        theme(panel.grid.major.y = element_blank(),
              legend.position = "none") +
        labs(x = NULL, y = "% Cells Expressing JAG1",
             title = "JAG1 Expression Across snRNA-seq Clusters",
             subtitle = "LP1 (Leaf Primordium 1) highlighted in red")

      ggplotly(p, tooltip = "text")
    }

    # Component heatmap (z-scaled)
    render_component_heatmap <- function() {
      if (is.null(sc_cluster_summary)) {
        return(plotly_empty() %>% layout(title = "Single-cell data not available"))
      }

      component_cols <- c("jag1_pct", "tpr_mean_pct", "hda_mean_pct",
                           "krp_mean_pct", "cyclin_target_mean_pct", "cdk_mean_pct")
      available_cols <- component_cols[component_cols %in% names(sc_cluster_summary)]

      if (length(available_cols) < 2) {
        return(plotly_empty() %>% layout(title = "Not enough component data"))
      }

      mat <- as.matrix(sc_cluster_summary[, available_cols])
      rownames(mat) <- sc_cluster_summary$cluster

      # Z-scale each column
      mat_z <- apply(mat, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
      mat_z[is.na(mat_z)] <- 0

      # Clean column names
      colnames(mat_z) <- gsub("_mean_pct|_pct", "", colnames(mat_z))
      colnames(mat_z) <- gsub("hda", "HDAC", colnames(mat_z))
      colnames(mat_z) <- toupper(colnames(mat_z))
      colnames(mat_z) <- gsub("CYCLIN_TARGET", "Cyclin\n(targets)", colnames(mat_z))

      # Build heatmap with plotly for interactivity
      plot_ly(
        x = colnames(mat_z),
        y = rownames(mat_z),
        z = mat_z,
        type = "heatmap",
        colors = colorRamp(c("#2166AC", "#F7F7F7", "#B2182B")),
        hovertemplate = "Cluster: %{y}<br>Component: %{x}<br>Z-score: %{z:.2f}<extra></extra>"
      ) %>%
        layout(
          title = "Cell Cycle Component Expression Across Clusters (Z-scaled)",
          xaxis = list(title = ""),
          yaxis = list(title = "", autorange = "reversed"),
          margin = list(l = 120)
        )
    }

    # Co-localization scatter
    render_coloc_scatter <- function() {
      if (is.null(sc_target_coloc)) {
        return(plotly_empty() %>% layout(title = "Co-localization data not available"))
      }

      data <- sc_target_coloc

      if ("Tier" %in% names(data)) {
        data <- data %>% filter(Tier %in% input$coloc_tiers)
      }

      if (nrow(data) == 0) {
        return(plotly_empty() %>% layout(title = "No data matches selected tiers"))
      }

      p <- ggplot(data, aes(x = JAG1_Pct_Expressing, y = Target_Mean_Pct_Expressing,
                            color = if ("Tier" %in% names(.)) Tier else Cluster,
                            size = N_Cells,
                            text = paste("Cluster:", Cluster,
                                        "<br>JAG1%:", round(JAG1_Pct_Expressing, 1),
                                        "<br>Target mean%:", round(Target_Mean_Pct_Expressing, 1),
                                        "<br>N targets:", N_Targets_Detected,
                                        "<br>N cells:", N_Cells))) +
        geom_point(alpha = 0.7, stroke = 0.3) +
        scale_size_continuous(range = c(3, 10), name = "N Cells") +
        theme_publication(base_size = 12) +
        labs(
          x = "% Cells Expressing JAG1",
          y = "Mean % Cells Expressing Targets",
          title = "JAG1-Target Co-localization Across Clusters"
        )

      if ("Tier" %in% names(data)) {
        p <- p + scale_color_manual(values = tier_colors, name = "Tier")
      }

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    }

    # Target detection score by cluster
    render_target_score <- function() {
      if (is.null(sc_target_score)) {
        return(plotly_empty() %>% layout(title = "Target score data not available"))
      }

      data <- sc_target_score %>%
        mutate(
          is_lp1 = grepl("LP1|Leaf [Pp]rimordium\\s*1", cluster, ignore.case = TRUE),
          highlight = ifelse(is_lp1, "LP1", "Other")
        )

      # If no LP1 label found, highlight the max jag1_pct cluster
      if (!any(data$is_lp1) && "jag1_pct" %in% names(data)) {
        max_idx <- which.max(data$jag1_pct)
        data$highlight[max_idx] <- "LP1"
      }

      score_col <- if ("target_mean_pct_expressing" %in% names(data)) {
        "target_mean_pct_expressing"
      } else if ("pct_targets_detected" %in% names(data)) {
        "pct_targets_detected"
      } else {
        names(data)[3]  # fallback
      }

      p <- ggplot(data, aes(x = reorder(cluster, .data[[score_col]]),
                            y = .data[[score_col]],
                            fill = highlight,
                            text = paste("Cluster:", cluster,
                                        "<br>Score:", round(.data[[score_col]], 1),
                                        "<br>N cells:", n_cells))) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        coord_flip() +
        scale_fill_manual(values = c("LP1" = "#E74C3C", "Other" = "#BDC3C7"),
                          name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_publication(base_size = 12) +
        theme(panel.grid.major.y = element_blank(), legend.position = "none") +
        labs(x = NULL, y = "Mean % Expressing Targets",
             title = "JAG1 Target Detection Score by Cluster",
             subtitle = "LP1 highlighted in red")

      ggplotly(p, tooltip = "text")
    }

    # Cluster summary table
    output$cluster_table <- DT::renderDT({
      if (is.null(sc_cluster_summary)) return(NULL)

      sc_cluster_summary %>%
        mutate(across(where(is.numeric), ~ round(., 2))) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        )
    })

    # JAG per cluster table
    output$jag_table <- DT::renderDT({
      if (is.null(sc_jag_per_cluster)) return(NULL)

      sc_jag_per_cluster %>%
        mutate(across(where(is.numeric), ~ round(., 2))) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        )
    })
  })
}
