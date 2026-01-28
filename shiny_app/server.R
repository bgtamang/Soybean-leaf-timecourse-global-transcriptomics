# Main server logic

server <- function(input, output, session) {

  # Reactive value for cross-tab gene selection
  selected_gene <- reactiveVal(NULL)

  # Call all module servers
  home_server("home", parent_session = session)
  samples_server("samples")
  expression_server("expression", selected_gene)
  de_server("de", selected_gene)
  jag1_server("jag1", selected_gene)
  wgcna_server("wgcna")
  go_server("go")
  phenotype_server("phenotype", selected_gene)
  temporal_server("temporal", selected_gene)
  download_server("download")
  methods_server("methods")

  # Cross-tab navigation when gene is selected
  observeEvent(selected_gene(), {
    gene <- selected_gene()
    if (!is.null(gene) && nchar(gene) > 0) {
      # Show notification
      showNotification(
        paste("Selected gene:", gene),
        type = "message",
        duration = 3
      )
    }
  })

  # Session info for debugging
  session$onSessionEnded(function() {
    cat("Session ended\n")
  })
}
