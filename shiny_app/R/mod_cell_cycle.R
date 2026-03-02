# Cell Cycle module - TPR/HDAC, KRP, Cyclin, CDK analysis (Figure 4)

cell_cycle_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(3,
        card(
          card_header("Cell Cycle Component"),
          card_body(
            selectInput(ns("component"), "Select Component:",
              choices = c(
                "TPR/HDAC Machinery" = "machinery",
                "KRP Inhibitors" = "krp",
                "Cyclins (JAG1 Targets)" = "cyclin",
                "CDK Kinases" = "cdk"
              ),
              selected = "machinery"
            ),
            hr(),
            uiOutput(ns("component_info"))
          )
        ),
        card(
          card_header("Key Finding"),
          card_body(
            p(tags$b("JAG1 represses D-type cyclins"),
              " — a mechanism distinct from Arabidopsis,",
              " where JAG promotes KRP expression."),
            hr(),
            tags$ul(class = "small",
              tags$li("7 cyclins are JAG1 targets (all D-type)"),
              tags$li("0 of 9 KRPs are JAG1 targets"),
              tags$li("0 of 46 CDKs are JAG1 targets"),
              tags$li("TPR/HDAC: co-repressor complex for JAG1")
            )
          )
        )
      ),
      column(9,
        tabsetPanel(
          tabPanel("Visualization",
            card(
              card_body(
                plotlyOutput(ns("main_plot"), height = "500px")
              )
            )
          ),
          tabPanel("Component Table",
            card(
              card_body(
                DT::DTOutput(ns("component_table"))
              )
            )
          ),
          tabPanel("Full Cyclin List",
            card(
              card_body(
                p(class = "text-muted", "All 101 soybean cyclins (PF00134/PF02984 validated)"),
                DT::DTOutput(ns("cyclin_full_table"))
              )
            )
          ),
          tabPanel("Full CDK List",
            card(
              card_body(
                p(class = "text-muted", "All 46 soybean CDKs (kinase-domain validated)"),
                DT::DTOutput(ns("cdk_full_table"))
              )
            )
          )
        )
      )
    )
  )
}

cell_cycle_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Component info sidebar
    output$component_info <- renderUI({
      switch(input$component,
        "machinery" = tagList(
          p(tags$b("TPR/HDAC Complex")),
          p("12 TPR + 15 HDAC genes form the co-repressor complex ",
            "recruited by JAG1's EAR motif."),
          if (!is.null(cell_cycle_machinery)) {
            stats <- cell_cycle_machinery$stats
            tags$table(class = "table table-sm",
              lapply(seq_len(nrow(stats)), function(i) {
                tags$tr(
                  tags$td(stats$Gene_Family[i]),
                  tags$td(paste0(stats$N_Genes[i], " genes")),
                  tags$td(paste0("mean ", round(stats$Mean_logCPM[i], 1), " logCPM"))
                )
              })
            )
          }
        ),
        "krp" = tagList(
          p(tags$b("KRP Inhibitors")),
          p("9 KRP (Kip-Related Protein) CDK inhibitors. ",
            "None are JAG1 targets — unlike Arabidopsis."),
          p(class = "small text-muted", "5 of 9 show higher expression in narrow-leaf at TP0")
        ),
        "cyclin" = tagList(
          p(tags$b("JAG1-Target Cyclins")),
          p("7 cyclins identified as JAG1 transcriptional targets. ",
            "All are D-type cyclins, repressed in narrow-leaf genotypes."),
          p(class = "small text-muted",
            "101 total cyclins in soybean; 7 are JAG1 targets")
        ),
        "cdk" = tagList(
          p(tags$b("CDK Kinases")),
          p("9 core CDKs (CDKA + CDKB) shown at TP0. ",
            "None are JAG1 targets — constitutive cell cycle machinery."),
          p(class = "small text-muted",
            "46 total CDKs in soybean; 0 are JAG1 targets or DE")
        )
      )
    })

    # Main visualization
    output$main_plot <- renderPlotly({
      switch(input$component,
        "machinery" = render_machinery_plot(),
        "krp"       = render_krp_plot(),
        "cyclin"    = render_cyclin_plot(),
        "cdk"       = render_cdk_plot()
      )
    })

    # TPR/HDAC horizontal bar chart
    render_machinery_plot <- function() {
      if (is.null(cell_cycle_machinery)) {
        return(plotly_empty() %>% layout(title = "Machinery data not available"))
      }
      genes <- cell_cycle_machinery$genes

      p <- ggplot(genes, aes(x = reorder(Name, Mean_Expression),
                             y = Mean_Expression,
                             fill = Gene_Family,
                             text = paste("Gene:", GeneID,
                                         "<br>Name:", Name,
                                         "<br>Mean Expression:", round(Mean_Expression, 2)))) +
        geom_col(width = 0.7, color = "black", linewidth = 0.3) +
        coord_flip() +
        scale_fill_manual(values = c("TPR" = "#3498DB", "HDAC" = "#E74C3C"),
                          name = "Family") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_publication(base_size = 12) +
        theme(panel.grid.major.y = element_blank()) +
        labs(x = NULL, y = "Mean Expression (log2-CPM)",
             title = "TPR/HDAC Co-repressor Complex Expression at TP0")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    }

    # KRP grouped bar chart
    render_krp_plot <- function() {
      if (is.null(cell_cycle_krp)) {
        return(plotly_empty() %>% layout(title = "KRP data not available"))
      }
      # Pivot Broad/Narrow columns to long format
      krp_long <- cell_cycle_krp %>%
        select(Name, Broad, Narrow) %>%
        pivot_longer(cols = c(Broad, Narrow), names_to = "Leaf_Type", values_to = "Expression")

      p <- ggplot(krp_long, aes(x = Name, y = Expression, fill = Leaf_Type,
                                text = paste("Gene:", Name,
                                            "<br>Leaf Type:", Leaf_Type,
                                            "<br>Expression:", round(Expression, 2)))) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.3) +
        scale_fill_manual(values = leaftype_colors, name = "Leaf Type") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_publication(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank()) +
        labs(x = NULL, y = "Expression (log2-CPM)",
             title = "KRP Expression at TP0: Broad vs Narrow")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    }

    # Cyclin lollipop chart
    render_cyclin_plot <- function() {
      if (is.null(cell_cycle_cyclin_targets)) {
        return(plotly_empty() %>% layout(title = "Cyclin data not available"))
      }

      # Use locusName as gene ID column
      gene_col <- if ("locusName" %in% names(cell_cycle_cyclin_targets)) "locusName" else "GeneID"

      df <- cell_cycle_cyclin_targets %>%
        mutate(Gene = .data[[gene_col]],
               Direction = ifelse(Mean_logFC_Pairwise > 0, "Up in Narrow", "Down in Narrow"))

      p <- ggplot(df, aes(x = reorder(Gene, Mean_logFC_Pairwise),
                          y = Mean_logFC_Pairwise,
                          color = Direction,
                          text = paste("Gene:", Gene,
                                      "<br>Cyclin type:", cyclin_type,
                                      "<br>Tier:", Confidence_Tier,
                                      "<br>logFC:", round(Mean_logFC_Pairwise, 2)))) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_segment(aes(xend = Gene, y = 0, yend = Mean_logFC_Pairwise), linewidth = 1) +
        geom_point(size = 4) +
        coord_flip() +
        scale_color_manual(values = c("Up in Narrow" = "#E74C3C", "Down in Narrow" = "#2E86AB")) +
        theme_publication(base_size = 12) +
        theme(panel.grid.major.y = element_blank()) +
        labs(x = NULL, y = "log2 Fold Change (Narrow vs Broad)",
             title = "JAG1-Target Cyclins at TP0")

      ggplotly(p, tooltip = "text")
    }

    # CDK grouped bar chart
    render_cdk_plot <- function() {
      if (is.null(cell_cycle_cdk)) {
        return(plotly_empty() %>% layout(title = "CDK data not available"))
      }

      # Use CDK_Label if available, fall back to GeneID
      label_col <- if ("CDK_Label" %in% names(cell_cycle_cdk)) "CDK_Label" else "GeneID"

      cdk_long <- cell_cycle_cdk %>%
        mutate(Label = .data[[label_col]]) %>%
        select(Label, Broad, Narrow) %>%
        pivot_longer(cols = c(Broad, Narrow), names_to = "Leaf_Type", values_to = "Expression")

      p <- ggplot(cdk_long, aes(x = Label, y = Expression, fill = Leaf_Type,
                                text = paste("Gene:", Label,
                                            "<br>Leaf Type:", Leaf_Type,
                                            "<br>Expression:", round(Expression, 2)))) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 color = "black", linewidth = 0.3) +
        scale_fill_manual(values = leaftype_colors, name = "Leaf Type") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_publication(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank()) +
        labs(x = NULL, y = "Expression (log2-CPM)",
             title = "Core CDK Expression at TP0: Broad vs Narrow")

      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "v", x = 1.02, y = 0.5))
    }

    # Component-specific data table
    output$component_table <- DT::renderDT({
      data <- switch(input$component,
        "machinery" = {
          if (!is.null(cell_cycle_machinery)) cell_cycle_machinery$genes else NULL
        },
        "krp" = cell_cycle_krp,
        "cyclin" = cell_cycle_cyclin_targets,
        "cdk" = cell_cycle_cdk
      )

      if (is.null(data)) return(NULL)

      datatable(data,
        filter = "top",
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      )
    })

    # Full cyclin browsable table
    output$cyclin_full_table <- DT::renderDT({
      if (is.null(cyclin_comprehensive)) return(NULL)

      display_cols <- c("locusName", "cyclin_type", "mean_CPM", "is_expressed",
                        "expression_level", "is_jag1_target", "Confidence_Tier",
                        "Mean_logFC_Pairwise", "Pattern", "has_binding",
                        "Best-hit-arabi-name", "Best-hit-arabi-defline")
      display_cols <- display_cols[display_cols %in% names(cyclin_comprehensive)]

      cyclin_comprehensive %>%
        select(all_of(display_cols)) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        )
    })

    # Full CDK browsable table
    output$cdk_full_table <- DT::renderDT({
      if (is.null(cdk_comprehensive)) return(NULL)

      display_cols <- c("locusName", "cdk_type", "cdk_function", "mean_CPM",
                        "is_expressed", "expression_level", "is_jag1_target",
                        "Confidence_Tier", "Mean_logFC_Pairwise", "Pattern",
                        "has_binding", "Best-hit-arabi-name", "Best-hit-arabi-defline")
      display_cols <- display_cols[display_cols %in% names(cdk_comprehensive)]

      cdk_comprehensive %>%
        select(all_of(display_cols)) %>%
        datatable(
          filter = "top",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE
        )
    })
  })
}
