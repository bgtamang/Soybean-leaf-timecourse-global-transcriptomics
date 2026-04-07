# Home module - Project overview

home_ui <- function(id) {
  ns <- NS(id)

 tagList(
    # Hero section with image
    div(
      class = "bg-dark text-white p-4 rounded-3 mb-4",
      fluidRow(
        column(7,
          div(class = "p-3",
            h1("GmJAG1 RNA-Seq Explorer", class = "display-5 fw-bold"),
            p(class = "lead mb-3",
              "Interactive exploration of soybean leaf development transcriptome data"),
            p("This dashboard provides comprehensive visualization and analysis tools for ",
              "RNA-seq data from the GmJAG1 soybean leaf development study.")
          )
        ),
        column(5,
          div(class = "d-flex align-items-start justify-content-center gap-3",
            div(class = "text-center",
              img(src = "Narrow_Broad_Soybean.png",
                  alt = "Narrow vs Broad Soybean Leaves",
                  class = "img-fluid rounded shadow",
                  style = "max-height: 250px; border: 2px solid #444;"),
              p(class = "text-muted small mt-2 mb-0",
                "Narrow (left) vs Broad (right) soybean trifoliate leaves")
            ),
            div(class = "text-center",
              img(src = "JAG1_structure.png",
                  alt = "JAG1 AlphaFold Structure",
                  class = "img-fluid rounded shadow",
                  style = "max-height: 250px; border: 2px solid #444;")
            )
          )
        )
      )
    ),

    # Stats cards
    fluidRow(
      column(3, value_box(
        title = "Total Samples",
        value = textOutput(ns("n_samples")),
        theme = "primary",
        showcase = icon("vials")
      )),
      column(3, value_box(
        title = "Genes Analyzed",
        value = textOutput(ns("n_genes")),
        theme = "info",
        showcase = icon("dna")
      )),
      column(3, value_box(
        title = "JAG1 Targets",
        value = textOutput(ns("n_targets")),
        theme = "success",
        showcase = icon("bullseye")
      )),
      column(3, value_box(
        title = "WGCNA Modules",
        value = textOutput(ns("n_modules")),
        theme = "warning",
        showcase = icon("project-diagram")
      ))
    ),

    # Experimental design
    fluidRow(
      column(6,
        card(
          card_header("Experimental Design"),
          card_body(
            h5("Genotypes"),
            tags$ul(
              tags$li(tags$b("Broad-leaf:"), " PI 532462A, LD11-2170 (functional GmJAG1)"),
              tags$li(tags$b("Narrow-leaf:"), " PI 547745, PI 612713B (D9H mutation)")
            ),
            h5("Timepoints"),
            tags$ul(
              tags$li(tags$b("TP1:"), " Shoot apex/meristem"),
              tags$li(tags$b("TP2-TP3:"), " Early leaf development"),
              tags$li(tags$b("TP4-TP5:"), " Mid to late leaf expansion")
            ),
            h5("Replicates"),
            p("3 biological replicates per genotype-timepoint combination")
          )
        )
      ),
      column(6,
        card(
          card_header("JAG1 Target Summary"),
          card_body(
            plotOutput(ns("tier_plot"), height = "250px"),
            hr(),
            p(class = "text-muted",
              "Targets classified into Gold (4/4 DE comparisons), ",
              "Silver (3/4), and Bronze (2/4 or pooled) confidence tiers.")
          )
        )
      )
    ),

    # Quick navigation
    fluidRow(
      column(12,
        card(
          card_header("Quick Navigation"),
          card_body(
            fluidRow(
              column(4,
                h5(icon("chart-line"), " Expression"),
                p("Search genes and visualize expression patterns across samples and timepoints."),
                actionButton(ns("go_expression"), "Explore Expression", class = "btn-outline-primary btn-sm")
              ),
              column(4,
                h5(icon("bullseye"), " JAG1 Targets"),
                p("Browse identified JAG1 target genes with tier classification."),
                actionButton(ns("go_jag1"), "View Targets", class = "btn-outline-success btn-sm")
              ),
              column(4,
                h5(icon("download"), " Download"),
                p("Export data tables, figures, and analysis results."),
                actionButton(ns("go_download"), "Download Data", class = "btn-outline-secondary btn-sm")
              )
            )
          )
        )
      )
    ),

    # Citation
    fluidRow(
      column(12,
        card(
          card_header("Citation"),
          card_body(
            p("If you use data from this dashboard, please cite:"),
            p(class = "fst-italic",
              "Tamang, B.G., Kramer, C., and Ainsworth, E.A. (2026). Natural leaf shape variation reveals ",
              "diverse transcriptional targets of GmJAG1 during soybean leaf development."),
            p(class = "text-muted small",
              "Data deposited in NCBI GEO: ",
              tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317596",
                     "GSE317596", target = "_blank"),
              " | Code: ",
              tags$a(href = "https://github.com/bgtamang/Soybean-leaf-timecourse-global-transcriptomics",
                     "GitHub", target = "_blank"))
          )
        )
      )
    )
  )
}

home_server <- function(id, parent_session = NULL) {
  moduleServer(id, function(input, output, session) {

    output$n_samples <- renderText({
      project_info$n_samples
    })

    output$n_genes <- renderText({
      format(project_info$n_genes, big.mark = ",")
    })

    output$n_targets <- renderText({
      format(project_info$n_targets, big.mark = ",")
    })

    output$n_modules <- renderText({
      project_info$n_modules
    })

    output$tier_plot <- renderPlot({
      tier_data <- data.frame(
        Tier = factor(c("Gold", "Silver", "Bronze"), levels = c("Gold", "Silver", "Bronze")),
        Count = c(project_info$n_gold, project_info$n_silver, project_info$n_bronze)
      )

      ggplot(tier_data, aes(x = Tier, y = Count, fill = Tier)) +
        geom_col(width = 0.7, color = "black", linewidth = 0.4) +
        geom_text(aes(label = Count), vjust = -0.7, size = 5, fontface = "bold") +
        scale_fill_manual(values = tier_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        theme_publication(base_size = 14) +
        theme(
          legend.position = "none",
          panel.grid.major.x = element_blank()
        ) +
        labs(x = NULL, y = "Number of Targets")
    })

    # Navigation buttons - use parent session if available
    nav_session <- if (!is.null(parent_session)) parent_session else session

    observeEvent(input$go_expression, {
      updateNavbarPage(nav_session, "main_nav", selected = "Expression")
    })

    observeEvent(input$go_jag1, {
      updateNavbarPage(nav_session, "main_nav", selected = "JAG1 Targets")
    })

    observeEvent(input$go_download, {
      updateNavbarPage(nav_session, "main_nav", selected = "Download")
    })
  })
}
