# Download module - Data export center

download_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        card(
          card_header(icon("download"), " Download Center"),
          card_body(
            p("Download processed data tables and analysis results from this study."),
            p(class = "text-muted small",
              "All data is provided in CSV format for easy import into R, Python, or Excel.")
          )
        )
      )
    ),
    fluidRow(
      column(4,
        card(
          card_header("JAG1 Target Data"),
          card_body(
            p("Download JAG1 target gene lists"),
            uiOutput(ns("tier_selector")),
            downloadButton(ns("download_targets"), "Download Targets", class = "btn-primary btn-block w-100")
          )
        )
      ),
      column(4,
        card(
          card_header("Expression Data"),
          card_body(
            p("Download expression matrices"),
            selectInput(ns("expr_type"), "Select Data:",
              choices = c("Sample Metadata" = "metadata",
                          "logCPM Matrix (Batch-corrected)" = "logcpm")
            ),
            downloadButton(ns("download_expression"), "Download Expression", class = "btn-info btn-block w-100")
          )
        )
      ),
      column(4,
        card(
          card_header("Analysis Results"),
          card_body(
            p("Download analysis summaries"),
            selectInput(ns("analysis_type"), "Select Analysis:",
              choices = c("GO Enrichment" = "go",
                          "WGCNA Modules" = "wgcna",
                          "Module-Trait Correlations" = "module_traits")
            ),
            downloadButton(ns("download_analysis"), "Download Analysis", class = "btn-success btn-block w-100")
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Custom Gene Set Export"),
          card_body(
            fluidRow(
              column(8,
                textAreaInput(ns("gene_list"), "Enter Gene IDs (one per line):",
                  rows = 5,
                  placeholder = "Glyma.20G116200\nGlyma.10G273800\n..."
                )
              ),
              column(4,
                br(),
                p("Enter gene IDs to export their expression values and annotations."),
                downloadButton(ns("download_custom"), "Export Custom Set", class = "btn-secondary")
              )
            )
          )
        )
      )
    ),
    fluidRow(
      column(12,
        card(
          card_header("Data Citation"),
          card_body(
            p("If you use data from this dashboard in your research, please cite:"),
            tags$blockquote(class = "blockquote",
              p(class = "mb-0",
                "Tamang, B.G., Kramer, C., and Ainsworth, E.A. (2026). Natural leaf shape variation reveals ",
                "diverse transcriptional targets of GmJAG1 during soybean leaf development. ",
                tags$em("The Plant Journal"), " 127, e71088. ",
                tags$a(href = "https://doi.org/10.1111/tpj.71088",
                       "https://doi.org/10.1111/tpj.71088", target = "_blank")
              )
            ),
            hr(),
            p(class = "small text-muted",
              "Raw sequencing data: ",
              tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317596",
                     "NCBI GEO GSE317596", target = "_blank")),
            p(class = "small text-muted",
              "Analysis code: ",
              tags$a(href = "https://github.com/bgtamang/Soybean-leaf-timecourse-global-transcriptomics",
                     "GitHub", target = "_blank"))
          )
        )
      )
    )
  )
}

download_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Dynamic tier selector with counts from data
    output$tier_selector <- renderUI({
      ns <- session$ns
      n_all <- project_info$n_targets
      n_gold <- project_info$n_gold
      n_silver <- project_info$n_silver
      n_bronze <- project_info$n_bronze
      choice_vals <- c("all", "gold", "silver", "bronze", "top50")
      choice_names <- c(
        paste0("All Targets (", format(n_all, big.mark = ","), ")"),
        paste0("Gold Tier (", n_gold, ")"),
        paste0("Silver Tier (", n_silver, ")"),
        paste0("Bronze Tier (", format(n_bronze, big.mark = ","), ")"),
        "Top 50"
      )
      selectInput(ns("target_tier"), "Select Tier:",
        choices = setNames(choice_vals, choice_names)
      )
    })

    # Download JAG1 targets
    output$download_targets <- downloadHandler(
      filename = function() {
        paste0("JAG1_targets_", input$target_tier, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        data <- switch(input$target_tier,
          "all" = jag1_targets,
          "gold" = jag1_targets %>% filter(Confidence_Tier == "Gold"),
          "silver" = jag1_targets %>% filter(Confidence_Tier == "Silver"),
          "bronze" = jag1_targets %>% filter(Confidence_Tier == "Bronze"),
          "top50" = jag1_targets %>% arrange(Rank) %>% head(50)
        )
        write.csv(data, file, row.names = FALSE)
      }
    )

    # Download expression data
    output$download_expression <- downloadHandler(
      filename = function() {
        paste0("expression_", input$expr_type, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        if (input$expr_type == "metadata") {
          write.csv(experimental_design, file, row.names = FALSE)
        } else if (input$expr_type == "logcpm") {
          expr <- load_expression()
          expr_df <- as.data.frame(expr)
          expr_df$GeneID <- rownames(expr_df)
          expr_df <- expr_df[, c("GeneID", setdiff(names(expr_df), "GeneID"))]
          write.csv(expr_df, file, row.names = FALSE)
        }
      }
    )

    # Download analysis results
    output$download_analysis <- downloadHandler(
      filename = function() {
        paste0(input$analysis_type, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        data <- switch(input$analysis_type,
          "go" = go_enrichment,
          "wgcna" = module_sizes,
          "module_traits" = module_traits
        )
        write.csv(data, file, row.names = FALSE)
      }
    )

    # Download custom gene set
    output$download_custom <- downloadHandler(
      filename = function() {
        paste0("custom_gene_set_", Sys.Date(), ".csv")
      },
      content = function(file) {
        genes <- strsplit(input$gene_list, "\n")[[1]]
        genes <- trimws(genes)
        genes <- genes[genes != ""]

        if (length(genes) == 0) {
          write.csv(data.frame(Message = "No genes provided"), file, row.names = FALSE)
          return()
        }

        # Get data for these genes
        expr <- load_expression()
        valid_genes <- genes[genes %in% rownames(expr)]

        if (length(valid_genes) == 0) {
          write.csv(data.frame(Message = "No valid genes found"), file, row.names = FALSE)
          return()
        }

        expr_sub <- as.data.frame(t(expr[valid_genes, , drop = FALSE]))
        expr_sub$Sample <- rownames(expr_sub)
        expr_sub <- expr_sub[, c("Sample", valid_genes)]

        # Add sample info (try both possible column names)
        join_col <- if ("Sample_name" %in% names(experimental_design)) "Sample_name" else "Sample"
        expr_out <- expr_sub %>%
          left_join(experimental_design, by = c("Sample" = join_col))

        write.csv(expr_out, file, row.names = FALSE)
      }
    )
  })
}
