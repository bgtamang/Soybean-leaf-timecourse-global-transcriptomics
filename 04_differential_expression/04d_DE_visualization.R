
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 13: DE VISUALIZATION\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/13_DE_viz", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "gplots",
  "UpSetR",
  "reshape2"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
load("03_results/checkpoints/12_DE_organized.RData")
cat("Loaded DE results\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded expression data\n\n")

# Get expression matrix
expr_mat <- v_primary$E
targets <- targets_primary

# Ensure Sample column exists
if (!"Sample" %in% colnames(targets)) {
  targets$Sample <- rownames(targets)
}
# Function to create volcano plot
create_volcano <- function(results, title, fdr_threshold = 0.05, logfc_threshold = 1) {
  results$Significance <- "NS"
  results$Significance[results$FDR < fdr_threshold & results$logFC > logfc_threshold] <- "Up"
  results$Significance[results$FDR < fdr_threshold & results$logFC < -logfc_threshold] <- "Down"
  results$Significance <- factor(results$Significance, levels = c("Down", "NS", "Up"))

  # Count significant genes
  n_up <- sum(results$Significance == "Up")
  n_down <- sum(results$Significance == "Down")

  p <- ggplot(results, aes(x = logFC, y = -log10(FDR), color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "#4575B4", "NS" = "gray70", "Up" = "#D73027")) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "gray40") +
    theme_bw(base_size = 12) +
    labs(title = title,
         subtitle = paste0("Up: ", n_up, " | Down: ", n_down),
         x = "log2 Fold Change",
         y = "-log10(FDR)") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")

  return(p)
}

# Create volcano plots for key contrasts
cat("Creating volcano plots for key contrasts...\n")

# Key contrasts (Narrow vs Broad at TP0)
for (contrast_name in key_contrasts) {
  if (contrast_name %in% names(de_results)) {
    results <- de_results[[contrast_name]]
    p <- create_volcano(results, paste0("Volcano: ", contrast_name))
    ggsave(paste0("03_results/figures/13_DE_viz/volcano_", contrast_name, ".png"),
           p, width = 8, height = 6, dpi = 300)
    cat("  Saved: volcano_", contrast_name, ".png\n", sep = "")
  }
}

# Pooled comparison
if ("NarrowvsBroad_TP0" %in% names(de_results)) {
  results <- de_results[["NarrowvsBroad_TP0"]]
  p <- create_volcano(results, "Volcano: Narrow vs Broad (Pooled) at TP0")
  ggsave("03_results/figures/13_DE_viz/volcano_NarrowvsBroad_pooled.png",
         p, width = 8, height = 6, dpi = 300)
  cat("  Saved: volcano_NarrowvsBroad_pooled.png\n")
}


cat("\n========================================\n")

# Function to create MD plot
create_md <- function(results, title, fdr_threshold = 0.05, logfc_threshold = 1) {
  results$Significance <- "NS"
  results$Significance[results$FDR < fdr_threshold & results$logFC > logfc_threshold] <- "Up"
  results$Significance[results$FDR < fdr_threshold & results$logFC < -logfc_threshold] <- "Down"
  results$Significance <- factor(results$Significance, levels = c("Down", "NS", "Up"))

  p <- ggplot(results, aes(x = logCPM, y = logFC, color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "#4575B4", "NS" = "gray70", "Up" = "#D73027")) +
    geom_hline(yintercept = c(-logfc_threshold, 0, logfc_threshold),
               linetype = c("dashed", "solid", "dashed"),
               color = c("gray40", "black", "gray40")) +
    theme_bw(base_size = 12) +
    labs(title = title,
         x = "Average log2 CPM",
         y = "log2 Fold Change") +
    theme(plot.title = element_text(size = 14, face = "bold"))

  return(p)
}

cat("Creating MD plots for key contrasts...\n")

for (contrast_name in key_contrasts) {
  if (contrast_name %in% names(de_results)) {
    results <- de_results[[contrast_name]]
    p <- create_md(results, paste0("MD Plot: ", contrast_name))
    ggsave(paste0("03_results/figures/13_DE_viz/MD_", contrast_name, ".png"),
           p, width = 8, height = 6, dpi = 300)
    cat("  Saved: MD_", contrast_name, ".png\n", sep = "")
  }
}


cat("\n========================================\n")

cat("Creating DEG barplots for each line (temporal comparisons)...\n")

# Function to create DEG barplot for a line
create_deg_barplot <- function(line_name, contrast_indices, coded_matrix, filename) {
  # Extract relevant contrasts
  if (all(contrast_indices <= ncol(coded_matrix))) {
    deg_data <- t(summary(coded_matrix)[c(1,3), contrast_indices])

    # Melt for ggplot
    deg_melted <- reshape2::melt(deg_data, varnames = c("Comparison", "Direction"),
                                 value.name = "gene.no")

    # Make down values negative for visualization
    deg_melted$gene.no <- ifelse(deg_melted$Direction == "Down",
                                 -abs(deg_melted$gene.no),
                                 deg_melted$gene.no)

    # Simplify comparison names
    deg_melted$Comparison <- gsub(paste0(line_name, "_"), "", deg_melted$Comparison)

    # Create plot
    p <- ggplot(deg_melted, aes(x = Comparison, y = gene.no, fill = Direction)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
      labs(y = "Number of DE genes", x = "Comparison",
           title = paste0("DE Genes: ", line_name),
           subtitle = "FDR < 0.05, |logFC| > 1") +
      geom_text(aes(label = abs(gene.no),
                    y = ifelse(gene.no > 0, gene.no + 50, gene.no - 50)),
                color = "black", size = 3) +
      geom_hline(yintercept = 0, color = "black") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold"))

    ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
    return(TRUE)
  }
  return(FALSE)
}

# Get contrast indices for each line (temporal vs TP0)
lines <- c("PI532462A", "LD112170", "PI612713B", "PI547745")

for (line in lines) {
  # Find temporal contrasts for this line
  temporal_pattern <- paste0("^", line, "_TP[1-4]vsTP0$")
  contrast_idx <- grep(temporal_pattern, colnames(edgeR_coded))

  if (length(contrast_idx) >= 2) {
    filename <- paste0("03_results/figures/13_DE_viz/DEG_barplot_", line, ".png")
    success <- create_deg_barplot(line, contrast_idx, edgeR_coded, filename)
    if (success) {
      cat("  Saved: DEG_barplot_", line, ".png\n", sep = "")
    }
  }
}


cat("\n========================================\n")

# Get top DE genes from pooled comparison
if ("NarrowvsBroad_TP0" %in% names(de_results)) {
  pooled_results <- de_results[["NarrowvsBroad_TP0"]]
  sig_genes <- pooled_results$Gene[pooled_results$Significant == TRUE]

  if (length(sig_genes) > 0) {
    cat("Creating heatmap of top DE genes (Narrow vs Broad)...\n")

    # Select top genes (by absolute logFC)
    top_n <- min(50, length(sig_genes))
    top_genes <- pooled_results %>%
      filter(Significant == TRUE) %>%
      arrange(desc(abs(logFC))) %>%
      head(top_n) %>%
      pull(Gene)

    # Get expression data for top genes
    heat_data <- expr_mat[top_genes, ]

    # Scale rows
    heat_data_scaled <- t(scale(t(heat_data)))

    # Sample annotation
    sample_anno <- data.frame(
      Timepoint = targets$Timepoint,
      Leaf_type = targets$Leaf_type,
      Line = targets$Line,
      row.names = targets$Sample
    )

    # Get actual factor levels for colors
    tp_levels <- sort(unique(targets$Timepoint))
    leaf_levels <- sort(unique(targets$Leaf_type))
    line_levels <- sort(unique(targets$Line))

    tp_colors <- setNames(
      colorRampPalette(c("#F7FBFF", "#08306B"))(length(tp_levels)),
      tp_levels
    )
    leaf_colors <- setNames(c("#2E8B57", "#D95F02")[1:length(leaf_levels)], leaf_levels)
    line_colors <- setNames(brewer.pal(max(3, length(line_levels)), "Set1")[1:length(line_levels)], line_levels)

    anno_colors <- list(
      Timepoint = tp_colors,
      Leaf_type = leaf_colors,
      Line = line_colors
    )

    # Create heatmap
    png("03_results/figures/13_DE_viz/heatmap_top_DE_NarrowvsBroad.png",
        width = 14, height = 12, units = "in", res = 300)

    pheatmap(heat_data_scaled,
             color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
             annotation_col = sample_anno,
             annotation_colors = anno_colors,
             show_colnames = FALSE,
             fontsize_row = 6,
             main = paste0("Top ", top_n, " DE Genes: Narrow vs Broad at TP0"),
             clustering_distance_cols = "euclidean",
             clustering_method = "ward.D2")

    dev.off()
    cat("  Saved: heatmap_top_DE_NarrowvsBroad.png\n")
  }
}


cat("\n========================================\n")

# Venn diagram for Narrow vs Broad comparisons at TP0
cat("Creating Venn diagrams...\n")

# Get significant gene lists for each key contrast
sig_gene_lists <- list()
for (contrast_name in key_contrasts) {
  if (contrast_name %in% names(de_results)) {
    results <- de_results[[contrast_name]]
    sig_gene_lists[[contrast_name]] <- results$Gene[results$Significant == TRUE]
  }
}

if (length(sig_gene_lists) >= 2) {
  # Create decideTests-like matrix for Venn
  all_genes <- unique(unlist(sig_gene_lists))
  venn_matrix <- matrix(0, nrow = length(all_genes), ncol = length(sig_gene_lists))
  rownames(venn_matrix) <- all_genes
  colnames(venn_matrix) <- names(sig_gene_lists)

  for (i in 1:length(sig_gene_lists)) {
    venn_matrix[sig_gene_lists[[i]], i] <- 1
  }

  # Create Venn diagram
  png("03_results/figures/13_DE_viz/venn_NarrowvsBroad_TP0.png",
      width = 10, height = 10, units = "in", res = 300)

  vennDiagram(venn_matrix,
              include = "both",
              main = "DE Genes: Narrow vs Broad at TP0\n(All pairwise comparisons)",
              cex = c(1.2, 1, 0.8))

  dev.off()
  cat("  Saved: venn_NarrowvsBroad_TP0.png\n")

  # Also create UP-only Venn
  up_gene_lists <- list()
  for (contrast_name in key_contrasts) {
    if (contrast_name %in% names(de_results)) {
      results <- de_results[[contrast_name]]
      up_gene_lists[[contrast_name]] <- results$Gene[results$Significant == TRUE & results$Direction == "Up"]
    }
  }

  all_up_genes <- unique(unlist(up_gene_lists))
  if (length(all_up_genes) > 0) {
    venn_up_matrix <- matrix(0, nrow = length(all_up_genes), ncol = length(up_gene_lists))
    rownames(venn_up_matrix) <- all_up_genes
    colnames(venn_up_matrix) <- names(up_gene_lists)

    for (i in 1:length(up_gene_lists)) {
      venn_up_matrix[up_gene_lists[[i]], i] <- 1
    }

    png("03_results/figures/13_DE_viz/venn_NarrowvsBroad_TP0_UP.png",
        width = 10, height = 10, units = "in", res = 300)

    vennDiagram(venn_up_matrix,
                include = "both",
                main = "Genes UP in Narrow vs Broad at TP0\n(Potential JAG1 targets)",
                cex = c(1.2, 1, 0.8))

    dev.off()
    cat("  Saved: venn_NarrowvsBroad_TP0_UP.png\n")
  }
}


cat("\n========================================\n")

# Create Venn diagrams for each line (temporal comparisons)
lines <- c("PI532462A", "LD112170", "PI612713B", "PI547745")

for (line in lines) {
  # Get temporal contrasts for this line
  temporal_contrasts <- paste0(line, "_TP", 1:4, "vsTP0")
  temporal_contrasts <- temporal_contrasts[temporal_contrasts %in% names(de_results)]

  if (length(temporal_contrasts) >= 2) {
    # Get significant genes
    sig_lists <- list()
    for (tc in temporal_contrasts) {
      results <- de_results[[tc]]
      sig_lists[[tc]] <- results$Gene[results$Significant == TRUE]
    }

    # Simplify names for plot
    names(sig_lists) <- gsub(paste0(line, "_"), "", names(sig_lists))

    all_genes <- unique(unlist(sig_lists))
    if (length(all_genes) > 0) {
      venn_mat <- matrix(0, nrow = length(all_genes), ncol = length(sig_lists))
      rownames(venn_mat) <- all_genes
      colnames(venn_mat) <- names(sig_lists)

      for (i in 1:length(sig_lists)) {
        venn_mat[sig_lists[[i]], i] <- 1
      }

      png(paste0("03_results/figures/13_DE_viz/venn_temporal_", line, ".png"),
          width = 10, height = 10, units = "in", res = 300)

      vennDiagram(venn_mat,
                  include = "both",
                  main = paste0("DE Genes: ", line, " Temporal (vs TP0)"),
                  cex = c(1.2, 1, 0.8))

      dev.off()
      cat("  Saved: venn_temporal_", line, ".png\n", sep = "")
    }
  }
}


cat("\n========================================\n")

cat("Creating UpSet plots (better for >5 comparisons)...\n")

# Helper function to safely create UpSet plot
safe_upset_plot <- function(data, filename, title_main, title_sets,
                            main_color = "#4575B4", sets_color = "#D73027") {
  tryCatch({
    # Ensure data is numeric (0/1)
    data <- as.data.frame(lapply(data, as.numeric))

    # Need at least 2 sets and some genes
    if (ncol(data) < 2 || nrow(data) < 1) {
      cat("    Skipped: insufficient data for UpSet plot\n")
      return(FALSE)
    }

    # Check if there are any intersections (genes in multiple sets)
    if (max(rowSums(data)) < 1) {
      cat("    Skipped: no genes to plot\n")
      return(FALSE)
    }

    png(filename, width = 10, height = 8, units = "in", res = 300)

    print(upset(data,
                sets = colnames(data),
                order.by = "freq",
                main.bar.color = main_color,
                sets.bar.color = sets_color,
                text.scale = 1.3,
                set_size.show = TRUE))

    dev.off()
    return(TRUE)
  }, error = function(e) {
    # Close any open device
    try(dev.off(), silent = TRUE)
    cat("    Warning: UpSet plot failed -", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# UpSet plot for all temporal contrasts (all 4 lines, TP4 vs TP0)
if (exists("edgeR_coded") && ncol(edgeR_coded) > 0) {

  # Get TP4vsTP0 contrasts for all lines
  tp4_contrasts <- grep("_TP4vsTP0$", colnames(edgeR_coded), value = TRUE)

  if (length(tp4_contrasts) >= 2) {
    # Convert to binary matrix for UpSetR
    upset_data <- as.data.frame(abs(as.matrix(edgeR_coded[, tp4_contrasts])) > 0) * 1
    colnames(upset_data) <- gsub("_TP4vsTP0", "", tp4_contrasts)

    # Remove rows with no DEGs
    upset_data <- upset_data[rowSums(upset_data) > 0, , drop = FALSE]

    if (nrow(upset_data) > 0) {
      success <- safe_upset_plot(
        upset_data,
        "03_results/figures/13_DE_viz/upset_TP4vsTP0_all_lines.png",
        "Shared DE Genes", "DE Genes per Line",
        "#4575B4", "#D73027"
      )
      if (success) cat("  Saved: upset_TP4vsTP0_all_lines.png\n")
    }
  }

  # UpSet plot for Narrow vs Broad comparisons at TP0
  if (length(key_contrasts) >= 2) {
    # Get between-line contrasts that exist in edgeR_coded
    between_contrasts <- key_contrasts[key_contrasts %in% colnames(edgeR_coded)]

    if (length(between_contrasts) >= 2) {
      upset_between <- as.data.frame(abs(as.matrix(edgeR_coded[, between_contrasts])) > 0) * 1

      # Simplify names (remove special characters that might cause issues)
      simple_names <- gsub("_TP0$", "", between_contrasts)
      simple_names <- gsub("vs", "_vs_", simple_names)
      colnames(upset_between) <- simple_names

      # Remove rows with no DEGs
      upset_between <- upset_between[rowSums(upset_between) > 0, , drop = FALSE]

      if (nrow(upset_between) > 0) {
        success <- safe_upset_plot(
          upset_between,
          "03_results/figures/13_DE_viz/upset_NarrowvsBroad_TP0.png",
          "Shared DE Genes", "DE Genes per Comparison",
          "#2E8B57", "#D95F02"
        )
        if (success) cat("  Saved: upset_NarrowvsBroad_TP0.png\n")
      }
    }
  }

  # UpSet for genes UP in narrow (potential JAG1 targets)
  if (length(key_contrasts) >= 2) {
    up_contrasts <- key_contrasts[key_contrasts %in% colnames(edgeR_coded)]

    if (length(up_contrasts) >= 2) {
      # Get UP genes only (coded as 1)
      upset_up <- as.data.frame((as.matrix(edgeR_coded[, up_contrasts]) == 1) * 1)

      # Simplify names
      simple_names <- gsub("_TP0$", "", up_contrasts)
      simple_names <- gsub("vs", "_vs_", simple_names)
      colnames(upset_up) <- simple_names

      # Remove rows with no UP genes
      upset_up <- upset_up[rowSums(upset_up) > 0, , drop = FALSE]

      if (nrow(upset_up) > 0) {
        success <- safe_upset_plot(
          upset_up,
          "03_results/figures/13_DE_viz/upset_UP_in_Narrow_JAG1_targets.png",
          "Genes UP in Narrow", "UP Genes per Comparison",
          "#D73027", "#D73027"
        )
        if (success) {
          cat("  Saved: upset_UP_in_Narrow_JAG1_targets.png\n")
          cat("    (These are potential JAG1 targets - derepressed in narrow)\n")
        }
      } else {
        cat("  Note: No genes consistently UP in Narrow - skipping JAG1 targets UpSet\n")
      }
    }
  }
} else {
  cat("  Warning: edgeR_coded not found - skipping UpSet plots\n")
}


cat("\n========================================\n")
cat("SESSION INFO\n")

print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Figures generated:\n")
cat("    - Volcano plots for key contrasts\n")
cat("    - MD plots for key contrasts\n")
cat("    - DEG barplots per line (temporal)\n")
cat("    - Heatmap of top DE genes\n")
cat("    - Venn diagrams (Narrow vs Broad)\n")
cat("    - Venn diagrams (Temporal per line)\n")
cat("    - UpSet plots (multi-comparison overlaps)\n")
cat("\n")
