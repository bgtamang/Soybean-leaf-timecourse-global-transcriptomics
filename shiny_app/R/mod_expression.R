# Expression module - Gene search and expression visualization

expression_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Gene Search"),
          card_body(
            textInput(ns("gene_search"), "Search Gene:",
              placeholder = "e.g., Glyma.20G116200 or JAG1"
            ),
            helpText("Search by gene ID or Arabidopsis homolog name"),
            hr(),
            selectInput(ns("plot_type"), "Plot Type:",
              choices = c("Line Plot" = "line", "Box Plot" = "box", "Heatmap" = "heatmap")
            ),
            hr(),
            selectInput(ns("group_by"), "Group by:",
              choices = c("Timepoint", "Leaf_type", "Line"),
              selected = "Timepoint"
            ),
            selectInput(ns("color_by"), "Color by:",
              choices = c("Leaf_type", "Timepoint", "Line"),
              selected = "Leaf_type"
            )
          )
        ),
        card(
          card_header("Gene Info"),
          card_body(
            uiOutput(ns("gene_info"))
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Expression Plot",
            card(
              card_body(
                plotlyOutput(ns("expr_plot"), height = "450px")
              )
            )
          ),
          tabPanel("Expression Table",
            card(
              card_body(
                DT::DTOutput(ns("expr_table"))
              )
            )
          ),
          tabPanel("Multi-Gene Heatmap",
            card(
              card_body(
                fluidRow(
                  column(6,
                    selectInput(ns("gene_set"), "Select Gene Set:",
                      choices = c(
                        "Custom (enter below)" = "custom",
                        "Housekeeping - Top 20 (most stable)" = "hk_top20",
                        "Housekeeping - Top 50" = "hk_top50",
                        "Housekeeping - All" = "hk_all",
                        "JAG1 Targets - Gold" = "jag1_gold",
                        "JAG1 Targets - Silver" = "jag1_silver",
                        "JAG1 Targets - Bronze" = "jag1_bronze",
                        "Temporal - Persistent" = "temp_persistent",
                        "Temporal - TP1-Specific" = "temp_tp1",
                        "Temporal - Late Onset" = "temp_late"
                      ),
                      selected = "custom"
                    )
                  ),
                  column(6,
                    sliderInput(ns("heatmap_ngenes"), "Max genes to display:",
                      min = 10, max = 100, value = 50, step = 10
                    )
                  )
                ),
                textAreaInput(ns("gene_list"), "Enter Gene IDs (one per line):",
                  rows = 4,
                  placeholder = "Glyma.20G116200\nGlyma.10G273800"
                ),
                actionButton(ns("plot_heatmap"), "Generate Heatmap", class = "btn-primary mb-3"),
                plotOutput(ns("heatmap"), height = "500px")
              )
            )
          )
        )
      )
    )
  )
}

expression_server <- function(id, selected_gene) {
  moduleServer(id, function(input, output, session) {

    # Load expression data reactively
    expr_matrix <- reactive({
      load_expression()
    })

    # Search result
    search_result <- reactive({
      req(input$gene_search)
      search_term <- trimws(input$gene_search)

      if (nchar(search_term) < 3) return(NULL)

      # Check if it's in JAG1 targets
      matches <- jag1_targets %>%
        filter(
          grepl(search_term, GeneID, ignore.case = TRUE) |
          grepl(search_term, Best_hit_arabi_name, ignore.case = TRUE)
        )

      if (nrow(matches) > 0) {
        return(matches[1, ])
      }

      # Check expression matrix
      expr <- expr_matrix()
      if (search_term %in% rownames(expr)) {
        return(data.frame(GeneID = search_term, Best_hit_arabi_name = NA))
      }

      NULL
    }) %>% debounce(500)

    # Update selected gene
    observeEvent(search_result(), {
      result <- search_result()
      if (!is.null(result)) {
        selected_gene(result$GeneID)
      }
    })

    # Respond to external gene selection
    observeEvent(selected_gene(), {
      gene <- selected_gene()
      if (!is.null(gene) && gene != input$gene_search) {
        updateTextInput(session, "gene_search", value = gene)
      }
    }, ignoreInit = TRUE)

    # Gene info panel
    output$gene_info <- renderUI({
      result <- search_result()
      if (is.null(result)) {
        return(p(class = "text-muted", "Enter a gene ID to see details"))
      }

      tagList(
        tags$dl(
          tags$dt("Gene ID"),
          tags$dd(result$GeneID),
          if (!is.na(result$Best_hit_arabi_name)) tagList(
            tags$dt("Arabidopsis Homolog"),
            tags$dd(result$Best_hit_arabi_name)
          ),
          if ("Confidence_Tier" %in% names(result) && !is.na(result$Confidence_Tier)) tagList(
            tags$dt("JAG1 Target Tier"),
            tags$dd(
              span(class = paste0("badge bg-", tolower(result$Confidence_Tier)),
                   result$Confidence_Tier)
            )
          )
        )
      )
    })

    # Expression plot
    output$expr_plot <- renderPlotly({
      result <- search_result()
      if (is.null(result)) {
        return(plotly_empty() %>% layout(title = "Search for a gene to view expression"))
      }

      gene_id <- result$GeneID
      expr <- expr_matrix()

      if (!gene_id %in% rownames(expr)) {
        return(plotly_empty() %>% layout(title = paste("Gene", gene_id, "not found in expression matrix")))
      }

      # Get expression values
      gene_expr <- data.frame(
        Sample = colnames(expr),
        Expression = as.numeric(expr[gene_id, ])
      ) %>%
        left_join(experimental_design, by = "Sample")

      # Build publication-quality plot based on type
      if (input$plot_type == "line") {
        p <- ggplot(gene_expr, aes(x = Timepoint, y = Expression,
                                    color = .data[[input$color_by]],
                                    group = interaction(Line, .data[[input$color_by]]))) +
          stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
          stat_summary(fun = mean, geom = "point", size = 4) +
          stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8)

      } else if (input$plot_type == "box") {
        p <- ggplot(gene_expr, aes(x = .data[[input$group_by]], y = Expression,
                                    fill = .data[[input$color_by]])) +
          geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.6) +
          geom_jitter(width = 0.15, alpha = 0.6, size = 2, shape = 21,
                      fill = "white", stroke = 0.5)
      }

      p <- p +
        theme_publication(base_size = 13) +
        labs(
          title = paste("Expression Profile:", gene_id),
          subtitle = if (!is.na(result$Best_hit_arabi_name)) result$Best_hit_arabi_name else NULL,
          y = "Expression (log2-CPM)",
          x = NULL,
          color = gsub("_", " ", input$color_by),
          fill = gsub("_", " ", input$color_by)
        )

      if (input$color_by == "Leaf_type") {
        p <- p + scale_color_manual(values = leaftype_colors, name = "Leaf Type") +
          scale_fill_manual(values = leaftype_colors, name = "Leaf Type")
      } else if (input$color_by == "Timepoint") {
        p <- p + scale_color_manual(values = timepoint_colors, name = "Timepoint") +
          scale_fill_manual(values = timepoint_colors, name = "Timepoint")
      }

      ggplotly(p) %>% layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    })

    # Expression table
    output$expr_table <- DT::renderDT({
      result <- search_result()
      if (is.null(result)) return(NULL)

      gene_id <- result$GeneID
      expr <- expr_matrix()

      if (!gene_id %in% rownames(expr)) return(NULL)

      expr_df <- data.frame(
        Sample = colnames(expr),
        log2CPM = round(as.numeric(expr[gene_id, ]), 2)
      ) %>%
        left_join(experimental_design, by = "Sample") %>%
        select(Sample, Line, Leaf_type, Timepoint, log2CPM)

      datatable(expr_df, options = list(pageLength = 15), rownames = FALSE)
    })

    # Get genes based on selected gene set
    get_gene_set <- reactive({
      gene_set <- input$gene_set
      max_genes <- input$heatmap_ngenes

      genes <- switch(gene_set,
        "custom" = {
          g <- strsplit(input$gene_list, "\n")[[1]]
          g <- trimws(g)
          g[g != ""]
        },
        "hk_top20" = {
          if (!is.null(housekeeping_genes)) head(housekeeping_genes$Gene, 20) else character(0)
        },
        "hk_top50" = {
          if (!is.null(housekeeping_genes)) head(housekeeping_genes$Gene, 50) else character(0)
        },
        "hk_all" = {
          if (!is.null(housekeeping_genes)) head(housekeeping_genes$Gene, max_genes) else character(0)
        },
        "jag1_gold" = {
          if (!is.null(jag1_targets)) head(jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Gold"], max_genes) else character(0)
        },
        "jag1_silver" = {
          if (!is.null(jag1_targets)) head(jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Silver"], max_genes) else character(0)
        },
        "jag1_bronze" = {
          if (!is.null(jag1_targets)) head(jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Bronze"], max_genes) else character(0)
        },
        "temp_persistent" = {
          if (!is.null(jag1_targets) && "Pattern" %in% names(jag1_targets)) {
            head(jag1_targets$GeneID[jag1_targets$Pattern == "Persistent"], max_genes)
          } else character(0)
        },
        "temp_tp0" = {
          if (!is.null(jag1_targets) && "Pattern" %in% names(jag1_targets)) {
            head(jag1_targets$GeneID[jag1_targets$Pattern == "TP1_Specific"], max_genes)
          } else character(0)
        },
        "temp_late" = {
          if (!is.null(jag1_targets) && "Pattern" %in% names(jag1_targets)) {
            head(jag1_targets$GeneID[jag1_targets$Pattern == "Late_Onset"], max_genes)
          } else character(0)
        },
        character(0)
      )

      genes
    })

    # Update gene list text area when gene set changes
    observeEvent(input$gene_set, {
      if (input$gene_set != "custom") {
        genes <- get_gene_set()
        if (length(genes) > 0) {
          updateTextAreaInput(session, "gene_list", value = paste(genes, collapse = "\n"))
        }
      }
    })

    # Multi-gene heatmap
    output$heatmap <- renderPlot({
      req(input$plot_heatmap)

      isolate({
        genes <- get_gene_set()

        if (length(genes) < 2) {
          plot.new()
          text(0.5, 0.5, "Enter at least 2 gene IDs or select a gene set", cex = 1.5)
          return()
        }

        expr <- expr_matrix()

        # Direct match - gene-level IDs should match directly
        valid_genes <- genes[genes %in% rownames(expr)]

        if (length(valid_genes) < 2) {
          plot.new()
          text(0.5, 0.5, "Not enough valid genes found in expression matrix", cex = 1.5)
          return()
        }

        expr_sub <- expr[valid_genes, ]

        # Scale by row
        expr_scaled <- t(scale(t(expr_sub)))

        # Annotation
        col_anno <- experimental_design %>%
          filter(Sample %in% colnames(expr_scaled)) %>%
          select(Sample, Leaf_type, Timepoint) %>%
          as.data.frame()
        rownames(col_anno) <- col_anno$Sample
        col_anno$Sample <- NULL

        # Get title based on gene set
        title <- switch(input$gene_set,
          "custom" = "Custom Gene Set Heatmap",
          "hk_top20" = "Housekeeping Genes - Top 20 Most Stable",
          "hk_top50" = "Housekeeping Genes - Top 50 Most Stable",
          "hk_all" = paste("Housekeeping Genes (n =", length(valid_genes), ")"),
          "jag1_gold" = paste("JAG1 Gold Tier Targets (n =", length(valid_genes), ")"),
          "jag1_silver" = paste("JAG1 Silver Tier Targets (n =", length(valid_genes), ")"),
          "jag1_bronze" = paste("JAG1 Bronze Tier Targets (n =", length(valid_genes), ")"),
          "temp_persistent" = paste("Persistent Pattern Genes (n =", length(valid_genes), ")"),
          "temp_tp1" = paste("TP1-Specific Pattern Genes (n =", length(valid_genes), ")"),
          "temp_late" = paste("Late Onset Pattern Genes (n =", length(valid_genes), ")"),
          "Multi-Gene Expression Heatmap"
        )

        # Publication-quality multi-gene heatmap
        gg_heatmap(
          expr_scaled,
          title = title,
          cluster_rows = TRUE,
          show_rownames = length(valid_genes) <= 50,
          fontsize_row = if(length(valid_genes) <= 30) 9 else 7,
          colors = c("#2166AC", "#F7F7F7", "#B2182B"),
          legend_title = "Z-score"
        )
      })
    })
  })
}
