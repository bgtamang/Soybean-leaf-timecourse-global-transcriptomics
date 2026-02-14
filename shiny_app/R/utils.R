# Utility functions for the Shiny dashboard

# Format p-values for display
format_pvalue <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "<0.001",
      ifelse(p < 0.01, sprintf("%.3f", p),
        sprintf("%.2f", p))))
}

# Format large numbers with commas
format_number <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}

# Truncate long strings
truncate_string <- function(s, max_length = 50) {
  ifelse(nchar(s) > max_length,
    paste0(substr(s, 1, max_length - 3), "..."),
    s)
}

# Create empty plotly plot with message
plotly_empty <- function(message = "No data available") {
  plot_ly() %>%
    layout(
      title = list(text = message, y = 0.5),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

# Safe data loading with error handling
safe_load_rds <- function(file_path, default = NULL) {
  tryCatch({
    if (file.exists(file_path)) {
      readRDS(file_path)
    } else {
      default
    }
  }, error = function(e) {
    warning(paste("Error loading", file_path, ":", e$message))
    default
  })
}

# Color functions (must match global.R palettes)
get_tier_color <- function(tier) {
  colors <- c("Gold" = "#D4AF37", "Silver" = "#8E8E8E", "Bronze" = "#CD7F32")
  colors[tier]
}

get_leaftype_color <- function(leaf_type) {
  colors <- c("Broad" = "#2E86AB", "Narrow" = "#E94F37")
  colors[leaf_type]
}

# Gene search helper
search_genes <- function(search_term, gene_table, columns = c("GeneID", "Best_hit_arabi_name")) {
  if (is.null(search_term) || nchar(search_term) < 3) {
    return(NULL)
  }

  search_term <- tolower(search_term)

  for (col in columns) {
    if (col %in% names(gene_table)) {
      matches <- gene_table %>%
        filter(grepl(search_term, tolower(.data[[col]])))
      if (nrow(matches) > 0) {
        return(matches)
      }
    }
  }

  NULL
}

# Create value box (compatible with both shinydashboard and bslib)
create_value_box <- function(value, title, icon_name = "info", color = "primary") {
  if (requireNamespace("bslib", quietly = TRUE)) {
    bslib::value_box(
      title = title,
      value = value,
      theme = color,
      showcase = icon(icon_name)
    )
  } else {
    div(class = paste0("value-box bg-", color, " text-white p-3 rounded"),
      h5(title),
      h2(value),
      icon(icon_name)
    )
  }
}
