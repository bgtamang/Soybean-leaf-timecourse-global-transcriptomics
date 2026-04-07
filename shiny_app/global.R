# Global data loading and packages
# Optimized: All data loaded once at startup for efficiency

# Required packages
library(shiny)
library(bslib)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Optional packages (with graceful fallback)
if (requireNamespace("pheatmap", quietly = TRUE)) library(pheatmap)
if (requireNamespace("waiter", quietly = TRUE)) library(waiter)
if (requireNamespace("shinyWidgets", quietly = TRUE)) library(shinyWidgets)

# Data directory
data_dir <- "data"

# Helper function to safely load RDS files
safe_load <- function(filename) {
  path <- file.path(data_dir, filename)
  if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(e) NULL)
  } else {
    NULL
  }
}

# ============================================================
# LOAD ALL DATA ONCE AT STARTUP
# This is efficient because:
# 1. Data is read from disk only once
# 2. All sessions share the same data in memory
# 3. No repeated file I/O during user interactions
# ============================================================

message("Loading dashboard data...")

# Core data (always needed)
experimental_design <- safe_load("experimental_design.rds")
jag1_targets <- safe_load("jag1_targets.rds")
tier_summary <- safe_load("tier_summary.rds")

# Module/WGCNA data
module_traits <- safe_load("module_traits.rds")
module_sizes <- safe_load("module_sizes.rds")
hub_genes <- safe_load("hub_genes.rds")
jag1_module_enrichment <- safe_load("jag1_module_enrichment.rds")

# GO enrichment data
go_enrichment <- safe_load("go_enrichment.rds")
go_gold <- safe_load("go_gold.rds")
go_silver <- safe_load("go_silver.rds")
go_bronze <- safe_load("go_bronze.rds")

# Expression and DE data (GENE-LEVEL, not transcript-level)
# Run data_prep/prepare_data.R if you need to recreate from checkpoints
expression_matrix <- safe_load("expression_matrix.rds")
de_results <- safe_load("de_results.rds")
contrast_list <- safe_load("contrast_list.rds")

# Sample/QC data
pca_data <- safe_load("pca_data.rds")
qc_metrics <- safe_load("qc_metrics.rds")

# Phenotype data
phenotype_traits <- safe_load("phenotype_traits.rds")
gene_phenotype_cors <- safe_load("gene_phenotype_cors.rds")

# Temporal data
temporal_classification <- safe_load("temporal_classification.rds")

# Other data
binding_integration <- safe_load("binding_integration.rds")
gene_list <- safe_load("gene_list.rds")

# Cell cycle data (Figure 4)
cell_cycle_machinery      <- safe_load("cell_cycle_machinery.rds")
cell_cycle_krp            <- safe_load("cell_cycle_krp.rds")
cell_cycle_cyclin_targets <- safe_load("cell_cycle_cyclin_targets.rds")
cell_cycle_cdk            <- safe_load("cell_cycle_cdk.rds")
cyclin_comprehensive      <- safe_load("cyclin_comprehensive.rds")
cdk_comprehensive         <- safe_load("cdk_comprehensive.rds")

# Hormone pathway data (Figure 4E)
hormone_enrichment  <- safe_load("hormone_enrichment.rds")
hormone_genes_jag1  <- safe_load("hormone_genes_jag1.rds")

# Validated genes - 79 high-confidence targets (Figure 6E)
validated_genes  <- safe_load("validated_genes.rds")
evidence_summary <- safe_load("evidence_summary.rds")

# Single-cell integration data (Fan et al. 2025 snRNA-seq)
sc_cluster_summary <- safe_load("sc_cluster_summary.rds")
sc_jag_per_cluster <- safe_load("sc_jag_per_cluster.rds")
sc_target_coloc    <- safe_load("sc_target_coloc.rds")
sc_target_score    <- safe_load("sc_target_score.rds")

# Housekeeping genes (Machado et al. 2020)
housekeeping_genes <- tryCatch({
  # Try multiple possible paths
  possible_paths <- c(
    file.path(data_dir, "housekeeping_stability.csv"),
    "../03_results/tables/housekeeping_stability.csv",
    file.path(dirname(getwd()), "03_results", "tables", "housekeeping_stability.csv")
  )
  hk_path <- NULL
  for (p in possible_paths) {
    if (file.exists(p)) {
      hk_path <- p
      break
    }
  }
  if (!is.null(hk_path)) {
    hk <- read.csv(hk_path, stringsAsFactors = FALSE)
    message("  Loaded housekeeping genes from: ", hk_path, " (", nrow(hk), " genes)")
    hk
  } else {
    message("  WARNING: Housekeeping genes file not found")
    NULL
  }
}, error = function(e) {
  message("  ERROR loading housekeeping genes: ", e$message)
  NULL
})

message("Data loading complete!")

# ============================================================
# ACCESSOR FUNCTIONS (for backward compatibility)
# These now just return the pre-loaded data
# ============================================================

load_expression <- function() {
  expression_matrix
}

load_de_results <- function() {
  de_results
}

load_pca_data <- function() {
  pca_data
}

load_qc_metrics <- function() {
  qc_metrics
}

load_temporal_data <- function() {
  temporal_classification
}

load_phenotype_cors <- function() {
  gene_phenotype_cors
}

# Get GO data by tier
get_go_data <- function(tier = "all") {
  switch(tier,
    "gold" = go_gold,
    "silver" = go_silver,
    "bronze" = go_bronze,
    go_enrichment  # default: all

)
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

format_pvalue <- function(p) {
  ifelse(p < 0.001, "<0.001",
         ifelse(p < 0.01, sprintf("%.3f", p),
                sprintf("%.2f", p)))
}

# ggplot2-based heatmap function (publication quality)
gg_heatmap <- function(mat, title = "", row_anno = NULL,
                       colors = NULL, cluster_rows = TRUE,
                       show_rownames = TRUE, fontsize_row = 9,
                       legend_title = "Value") {

  # Use professional diverging palette by default
 if (is.null(colors)) {
    colors <- c("#2166AC", "#F7F7F7", "#B2182B")
  }

  # Convert matrix to long format
  df <- as.data.frame(mat)
  df$Row <- rownames(mat)

  # Row ordering (hierarchical clustering or original)
  if (cluster_rows && nrow(mat) > 1) {
    row_order <- tryCatch({
      hc <- hclust(dist(mat))
      hc$labels[hc$order]
    }, error = function(e) rownames(mat))
  } else {
    row_order <- rownames(mat)
  }

  df_long <- df %>%
    pivot_longer(cols = -Row, names_to = "Col", values_to = "Value") %>%
    mutate(
      Row = factor(Row, levels = rev(row_order)),
      Col = factor(Col, levels = colnames(mat))
    )

  # Add annotation column if provided
  if (!is.null(row_anno) && ncol(row_anno) > 0) {
    anno_col <- names(row_anno)[1]
    df_long <- df_long %>%
      left_join(
        data.frame(Row = rownames(row_anno), Anno = row_anno[[anno_col]]),
        by = "Row"
      )
  }

  # Calculate appropriate midpoint
  val_range <- range(df_long$Value, na.rm = TRUE)
  midpoint <- if (val_range[1] < 0 && val_range[2] > 0) 0 else mean(val_range)

  # Build publication-quality plot
  p <- ggplot(df_long, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = colors[1], mid = colors[2], high = colors[3],
      midpoint = midpoint, name = legend_title,
      guide = guide_colorbar(
        barwidth = 1, barheight = 8,
        title.position = "top", title.hjust = 0.5,
        frame.colour = "black", ticks.colour = "black"
      )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_heatmap(base_size = 11) +
    theme(
      axis.text.y = element_text(size = fontsize_row)
    ) +
    labs(title = title)

  # Add facet for annotation if present
  if (!is.null(row_anno) && ncol(row_anno) > 0) {
    p <- p +
      facet_grid(Anno ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 10)
      )
  }

  if (!show_rownames) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  p
}

# ============================================================
# PUBLICATION-QUALITY THEME
# ============================================================

# Professional theme inspired by GraphPad Prism / Nature journals
theme_publication <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),

      # Axis
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(colour = "black", size = rel(0.9)),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.text.y = element_text(margin = margin(r = 4)),
      axis.title = element_text(colour = "black", face = "bold", size = rel(1.0)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10), angle = 90),

      # Legend
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.text = element_text(size = rel(0.85)),
      legend.position = "right",
      legend.margin = margin(5, 5, 5, 5),

      # Title and subtitle
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5,
                                   margin = margin(b = 10), colour = "grey40"),
      plot.caption = element_text(size = rel(0.75), hjust = 1, colour = "grey50"),
      plot.margin = margin(15, 15, 15, 15),

      # Strip (facets)
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = rel(0.9), margin = margin(5, 5, 5, 5))
    )
}

# Clean theme for heatmaps
theme_heatmap <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                  colour = "black", size = rel(0.95)),
      axis.text.y = element_text(colour = "black", size = rel(0.85)),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.5,
                                margin = margin(b = 12)),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(0.4, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================
# PROFESSIONAL COLOR PALETTES
# ============================================================

# Tier colors (metallic, professional)
tier_colors <- c(
  "Gold" = "#D4AF37",
  "Silver" = "#8E8E8E",
  "Bronze" = "#CD7F32"
)

# Leaf type colors (colorblind-friendly, high contrast)
leaftype_colors <- c(
  "Broad" = "#2E86AB",   # Steel blue
  "Narrow" = "#E94F37"   # Vermillion red
)

# Timepoint colors (viridis gradient — matches publication Fig 2A)
timepoint_colors <- c(
  "TP1" = "#440154",  # Dark purple (shoot apex)
  "TP2" = "#3B528B",  # Blue
  "TP3" = "#21908C",  # Teal
  "TP4" = "#5DC863",  # Green
  "TP5" = "#FDE725"   # Yellow (mature)
)

# Scientific color palettes
palette_diverging <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
                       "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
palette_sequential <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
                        "#6BAED6", "#4292C6", "#2171B5", "#084594")
palette_categorical <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# Heatmap colors (diverging, publication standard)
heatmap_colors <- colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7",
                                      "#EF8A62", "#B2182B"))(100)

# ============================================================
# PROJECT INFO
# ============================================================

project_info <- list(
  title = "GmJAG1 RNA-Seq Explorer",
  subtitle = "Soybean Leaf Development Transcriptome",
  n_samples = if (!is.null(experimental_design)) nrow(experimental_design) else 0,
  n_genes = if (!is.null(expression_matrix)) nrow(expression_matrix) else 0,
  n_targets = if (!is.null(jag1_targets)) nrow(jag1_targets) else 0,
  n_gold = if (!is.null(jag1_targets)) sum(jag1_targets$Confidence_Tier == "Gold") else 0,
  n_silver = if (!is.null(jag1_targets)) sum(jag1_targets$Confidence_Tier == "Silver") else 0,
  n_bronze = if (!is.null(jag1_targets)) sum(jag1_targets$Confidence_Tier == "Bronze") else 0,
  n_modules = if (!is.null(module_sizes)) nrow(module_sizes) else 0
)

# ============================================================
# LOAD MODULE FILES
# ============================================================

for (mod_file in list.files("R", pattern = "^mod_.*\\.R$", full.names = TRUE)) {
  source(mod_file)
}

# Load utility functions
if (file.exists("R/utils.R")) source("R/utils.R")
