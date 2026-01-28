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
            selectInput(ns("target_tier"), "Select Tier:",
              choices = c("All Targets (1,511)" = "all",
                          "Gold Tier (74)" = "gold",
                          "Silver Tier (65)" = "silver",
                          "Bronze Tier (1,372)" = "bronze",
                          "Top 50" = "top50")
            ),
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
                          "TPM Matrix (Top 20K genes)" = "tpm")
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
                "Tamang, B. et al. (2026). Transcriptome-wide identification of GmJAG1 ",
                "target genes in soybean leaf development. [Journal TBD]"
              )
            ),
            hr(),
            p(class = "small text-muted",
              "Raw sequencing data available at NCBI GEO: GSExxxxxx")
          )
        )
      )
    )
  )
}

download_server <- function(id) {
  moduleServer(id, function(input, output, session) {

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
        } else if (input$expr_type == "tpm") {
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

        # Add sample info
        expr_out <- expr_sub %>%
          left_join(experimental_design, by = c("Sample" = "Sample_name"))

        write.csv(expr_out, file, row.names = FALSE)
      }
    )
  })
}
