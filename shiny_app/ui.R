# Main UI structure

ui <- page_navbar(
  title = "GmJAG1 RNA-Seq Explorer",
  id = "main_nav",

  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    secondary = "#27AE60",
    success = "#27AE60",
    info = "#3498DB",
    warning = "#F39C12",
    danger = "#E74C3C"
  ),

  # Custom CSS
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style(HTML("
      .value-box { margin-bottom: 15px; }
      .nav-link { font-weight: 500; }
      .card { margin-bottom: 20px; }
      .dataTables_wrapper { font-size: 0.9em; }
    "))
  ),

  # Tab 1: Home
  nav_panel(
    title = "Home",
    icon = icon("home"),
    home_ui("home")
  ),


  # Tab 2: Samples
  nav_panel(
    title = "Samples",
    icon = icon("vials"),
    samples_ui("samples")
  ),

  # Tab 3: Phenotype
  nav_panel(
    title = "Phenotype",
    icon = icon("leaf"),
    phenotype_ui("phenotype")
  ),

  # Tab 4: Expression
  nav_panel(
    title = "Expression",
    icon = icon("chart-line"),
    expression_ui("expression")
  ),

  # Tab 5: DE Analysis
  nav_panel(
    title = "DE Analysis",
    icon = icon("chart-area"),
    de_ui("de")
  ),

  # Tab 6: JAG1 Targets
  nav_panel(
    title = "JAG1 Targets",
    icon = icon("bullseye"),
    jag1_ui("jag1")
  ),

  # Tab 7: WGCNA
  nav_panel(
    title = "WGCNA",
    icon = icon("project-diagram"),
    wgcna_ui("wgcna")
  ),

  # Tab 8: Cell Cycle
  nav_panel(
    title = "Cell Cycle",
    icon = icon("sync"),
    cell_cycle_ui("cell_cycle")
  ),

  # Tab 9: Hormones
  nav_panel(
    title = "Hormones",
    icon = icon("seedling"),
    hormone_ui("hormone")
  ),

  # Tab 10: GO Enrichment
  nav_panel(
    title = "GO Enrichment",
    icon = icon("sitemap"),
    go_ui("go")
  ),

  # Tab 11: Temporal
  nav_panel(
    title = "Temporal",
    icon = icon("clock"),
    temporal_ui("temporal")
  ),

  # Tab 12: 79 Validated
  nav_panel(
    title = "79 Validated",
    icon = icon("check-double"),
    validated_ui("validated")
  ),

  # Tab 13: Single-Cell
  nav_panel(
    title = "Single-Cell",
    icon = icon("braille"),
    single_cell_ui("single_cell")
  ),

  # Tab 14: Download
  nav_panel(
    title = "Download",
    icon = icon("download"),
    download_ui("download")
  ),

  # Tab 15: Methods
  nav_panel(
    title = "Methods",
    icon = icon("book"),
    methods_ui("methods")
  ),

  # Right-aligned items

  nav_spacer(),
  nav_item(
    tags$a(
      icon("github"),
      "GitHub",
      href = "https://github.com/",
      target = "_blank",
      class = "nav-link"
    )
  )
)
