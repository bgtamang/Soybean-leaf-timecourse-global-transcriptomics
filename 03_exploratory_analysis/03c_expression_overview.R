
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 09: EXPRESSION OVERVIEW\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/09_expression", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "pheatmap",
  "scales",
  "readxl",
  "gplots"
)

missing <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing) > 0) {
  install.packages(missing)
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
load("03_results/checkpoints/06_validated.RData")
cat("Loaded checkpoint from script 06\n\n")

# Use primary dataset
expr_mat <- v_primary$E
targets <- targets_primary

# Ensure Sample column exists (use rownames if not)
if (!"Sample" %in% colnames(targets)) {
  targets$Sample <- rownames(targets)
}

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n\n")


# Get actual factor levels from data
tp_levels <- unique(as.character(targets$Timepoint))
leaf_levels <- unique(as.character(targets$Leaf_type))
line_levels <- unique(as.character(targets$Line))

# Timepoint colors (blue gradient)
tp_palette <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
timepoint_colors <- setNames(tp_palette[1:length(tp_levels)], sort(tp_levels))

# Leaf type colors
leaf_palette <- c("#2E8B57", "#D95F02")
leaf_colors <- setNames(leaf_palette[1:length(leaf_levels)], sort(leaf_levels))

# Line colors
line_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
line_colors <- setNames(line_palette[1:length(line_levels)], sort(line_levels))

cat("Color mappings created for:\n")
cat("  Timepoints:", paste(names(timepoint_colors), collapse = ", "), "\n")
cat("  Leaf types:", paste(names(leaf_colors), collapse = ", "), "\n")
cat("  Lines:", paste(names(line_colors), collapse = ", "), "\n\n")
# --- Plot 1: Expression density by sample ---
cat("Creating expression density plot...\n")

expr_long <- as.data.frame(expr_mat) %>%
  mutate(Gene = rownames(expr_mat)) %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(targets %>% select(Sample, Timepoint, Leaf_type, Line), by = "Sample")

p1 <- ggplot(expr_long, aes(x = Expression, color = Timepoint, group = Sample)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Expression Density by Sample",
       subtitle = "Colored by timepoint",
       x = "log2(CPM)",
       y = "Density") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  facet_wrap(~Leaf_type)

ggsave("03_results/figures/09_expression/expression_density.png", p1,
       width = 12, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/09_expression/expression_density.png\n")

# --- Plot 2: Expression boxplot by sample ---
cat("Creating expression boxplot...\n")

# Calculate median expression for ordering
sample_medians <- expr_long %>%
  group_by(Sample, Timepoint, Leaf_type) %>%
  summarize(Median = median(Expression), .groups = "drop") %>%
  arrange(Timepoint, Leaf_type)

expr_long$Sample <- factor(expr_long$Sample, levels = sample_medians$Sample)

p2 <- ggplot(expr_long, aes(x = Sample, y = Expression, fill = Timepoint)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = timepoint_colors) +
  theme_bw(base_size = 12) +
  labs(title = "Expression Distribution by Sample",
       x = "Sample",
       y = "log2(CPM)") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

ggsave("03_results/figures/09_expression/expression_boxplot.png", p2,
       width = 16, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/09_expression/expression_boxplot.png\n")

# --- Plot 3: Summary statistics ---
cat("Creating summary statistics...\n")

sample_stats <- expr_long %>%
  group_by(Sample, Timepoint, Leaf_type, Line) %>%
  summarize(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    Q25 = quantile(Expression, 0.25),
    Q75 = quantile(Expression, 0.75),
    Expressed_pct = mean(Expression > 0) * 100,
    .groups = "drop"
  )

write.csv(sample_stats, "03_results/tables/expression_summary.csv", row.names = FALSE)
cat("  Saved: 03_results/tables/expression_summary.csv\n\n")
# Find JAG genes - CORRECT GENE IDs:
# GmJAG1 = Glyma.20G116200 (controls leaf shape, D9H mutation in EAR motif in narrow genotypes)
# GmJAG2 = Glyma.10G273800 (paralog)
jag_patterns <- c("Glyma.20G116200", "Glyma.10G273800")

find_genes <- function(patterns, gene_names) {
  matches <- character()
  for (p in patterns) {
    found <- grep(p, gene_names, ignore.case = TRUE, value = TRUE)
    matches <- c(matches, found)
  }
  unique(matches)
}

jag_genes <- find_genes(jag_patterns, rownames(expr_mat))

cat("Searching for JAG genes...\n")
cat("  Found", length(jag_genes), "potential JAG-related genes\n")

if (length(jag_genes) > 0) {
  cat("  Genes:", paste(jag_genes, collapse = ", "), "\n\n")

  # Extract JAG expression
  jag_expr <- expr_mat[jag_genes, , drop = FALSE]

  # Create long format for plotting
  jag_long <- as.data.frame(t(jag_expr)) %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "Expression") %>%
    left_join(targets %>% select(Sample, Timepoint, Leaf_type, Line), by = "Sample")

  # --- Plot 4: JAG expression across timepoints ---
  cat("Creating JAG expression profile...\n")

  jag_summary <- jag_long %>%
    group_by(Gene, Timepoint, Leaf_type) %>%
    summarize(
      Mean = mean(Expression),
      SE = sd(Expression) / sqrt(n()),
      .groups = "drop"
    )

  p4 <- ggplot(jag_summary, aes(x = Timepoint, y = Mean, color = Leaf_type, group = Leaf_type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    facet_wrap(~Gene, scales = "free_y") +
    scale_color_manual(values = leaf_colors) +
    theme_bw(base_size = 14) +
    labs(title = "JAG Gene Expression Across Timepoints",
         subtitle = "Mean +/- SE by leaf type",
         x = "Timepoint",
         y = "log2(CPM)") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  ggsave("03_results/figures/09_expression/JAG_expression_profile.png", p4,
         width = 10, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/09_expression/JAG_expression_profile.png\n")

  # --- Plot 5: JAG expression by line ---
  cat("Creating JAG expression by line...\n")

  jag_by_line <- jag_long %>%
    group_by(Gene, Timepoint, Line, Leaf_type) %>%
    summarize(
      Mean = mean(Expression),
      SE = sd(Expression) / sqrt(n()),
      .groups = "drop"
    )

  p5 <- ggplot(jag_by_line, aes(x = Timepoint, y = Mean, color = Line,
                                 linetype = Leaf_type, group = interaction(Line, Leaf_type))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~Gene, scales = "free_y") +
    scale_color_manual(values = line_colors) +
    theme_bw(base_size = 14) +
    labs(title = "JAG Expression by Line and Leaf Type",
         x = "Timepoint",
         y = "log2(CPM)") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  ggsave("03_results/figures/09_expression/JAG_expression_by_line.png", p5,
         width = 12, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/09_expression/JAG_expression_by_line.png\n")

  # Save JAG expression table
  jag_table <- jag_long %>%
    group_by(Gene, Timepoint, Leaf_type) %>%
    summarize(
      Mean = round(mean(Expression), 3),
      SD = round(sd(Expression), 3),
      Min = round(min(Expression), 3),
      Max = round(max(Expression), 3),
      .groups = "drop"
    )

  write.csv(jag_table, "03_results/tables/JAG_expression.csv", row.names = FALSE)
  cat("  Saved: 03_results/tables/JAG_expression.csv\n\n")

} else {
  cat("  WARNING: No JAG genes found in dataset\n")
  cat("  Checking gene ID format...\n")
  cat("  Sample gene IDs:", head(rownames(expr_mat), 5), "\n\n")
}
# Load housekeeping genes from Machado et al. 2019
hk_file <- file.path(base_dir, "Housekeeping_genes_Machado_et_al_2019.xlsx")

if (file.exists(hk_file)) {
  cat("Loading housekeeping genes from Machado et al. 2019...\n")
  hk_ref <- read_excel(hk_file)
  cat("  Columns in file:", paste(colnames(hk_ref), collapse = ", "), "\n")

  # Try to find gene ID column (common names: Gene, GeneID, Glyma, ID, etc.)
  gene_col <- NULL
  possible_cols <- c("Gene", "GeneID", "Gene_ID", "Glyma", "ID", "gene", "gene_id")
  for (col in possible_cols) {
    if (col %in% colnames(hk_ref)) {
      gene_col <- col
      break
    }
  }

  # If no standard name found, use first column
  if (is.null(gene_col)) {
    gene_col <- colnames(hk_ref)[1]
    cat("  Using first column as gene IDs:", gene_col, "\n")
  }

  hk_gene_ids <- as.character(hk_ref[[gene_col]])
  cat("  Found", length(hk_gene_ids), "housekeeping genes in reference list\n")
  cat("  Sample IDs:", paste(head(hk_gene_ids, 3), collapse = ", "), "...\n")

  # Match with expression matrix
  hk_genes <- intersect(hk_gene_ids, rownames(expr_mat))
  cat("  Matched", length(hk_genes), "genes in expression matrix\n\n")

} else {
  cat("Housekeeping gene file not found, using pattern matching...\n")
  # Fallback to pattern matching
  hk_patterns <- c(
    "ACTIN", "ACT", "Glyma.15G080200",
    "TUBULIN", "TUB", "Glyma.03G164100",
    "ELF1", "EF1", "Glyma.05G178300",
    "GAPDH", "Glyma.06G046000",
    "UBI", "UBQ", "Glyma.20G141600"
  )
  hk_genes <- find_genes(hk_patterns, rownames(expr_mat))
  cat("Found", length(hk_genes), "potential housekeeping genes\n")
}

if (length(hk_genes) > 0) {
  # Calculate coefficient of variation for each gene
  hk_expr <- expr_mat[hk_genes, , drop = FALSE]

  hk_cv <- apply(hk_expr, 1, function(x) sd(x) / mean(x) * 100)
  hk_mean <- rowMeans(hk_expr)
  hk_sd <- apply(hk_expr, 1, sd)

  hk_stability <- data.frame(
    Gene = names(hk_cv),
    Mean_Expression = round(hk_mean, 3),
    SD = round(hk_sd, 3),
    CV_percent = round(hk_cv, 2)
  ) %>%
    arrange(CV_percent)

  cat("\nMost stable housekeeping genes (lowest CV):\n")
  print(head(hk_stability, 10))

  write.csv(hk_stability, "03_results/tables/housekeeping_stability.csv", row.names = FALSE)
  cat("\nSaved: 03_results/tables/housekeeping_stability.csv\n")

  # --- Plot 6: Housekeeping gene boxplot ---
  if (length(hk_genes) > 1) {
    cat("\nCreating housekeeping gene boxplot...\n")

    # Select top stable genes for visualization
    top_hk <- head(hk_stability$Gene, min(8, nrow(hk_stability)))

    hk_long <- as.data.frame(t(hk_expr[top_hk, , drop = FALSE])) %>%
      mutate(Sample = rownames(.)) %>%
      pivot_longer(cols = -Sample, names_to = "Gene", values_to = "Expression") %>%
      left_join(targets %>% select(Sample, Timepoint), by = "Sample")

    p6 <- ggplot(hk_long, aes(x = Gene, y = Expression, fill = Timepoint)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = timepoint_colors) +
      theme_bw(base_size = 14) +
      labs(title = "Housekeeping Gene Expression Stability",
           subtitle = "Top stable genes (lowest CV)",
           x = "Gene",
           y = "log2(CPM)") +
      theme(plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    ggsave("03_results/figures/09_expression/housekeeping_stability.png", p6,
           width = 10, height = 6, dpi = 300)
    cat("  Saved: 03_results/figures/09_expression/housekeeping_stability.png\n")
  }

  # --- Plot 7: Housekeeping gene heatmap (Publication quality - Fig S9) ---
  cat("\nCreating publication-quality housekeeping gene heatmap...\n")

  # Create supplementary figures directory
  dir.create("03_results/figures/supplementary", recursive = TRUE, showWarnings = FALSE)

  # Select top 50 most stable genes for visualization
  top50_genes <- head(hk_stability$Gene, 50)
  hk_expr_top50 <- hk_expr[top50_genes, , drop = FALSE]

  # Z-score normalize for visualization
  hk_scaled <- t(scale(t(hk_expr_top50)))
  hk_scaled[hk_scaled > 2] <- 2
  hk_scaled[hk_scaled < -2] <- -2

  # Create sample annotation (order by Leaf_type then Timepoint)
  sample_anno <- data.frame(
    Leaf_Type = targets$Leaf_type,
    Timepoint = targets$Timepoint,
    row.names = targets$Sample
  )
  sample_order <- order(sample_anno$Leaf_Type, sample_anno$Timepoint)
  hk_scaled <- hk_scaled[, sample_order]
  sample_anno <- sample_anno[sample_order, , drop = FALSE]

  # Create gene annotation (CV values)
  gene_anno <- data.frame(
    CV = hk_stability$CV_percent[match(top50_genes, hk_stability$Gene)],
    row.names = top50_genes
  )

  # Annotation colors
  anno_colors <- list(
    Leaf_Type = c("Broad" = "#2166AC", "Narrow" = "#B2182B"),
    Timepoint = c("TP0" = "#FEE5D9", "TP1" = "#FCAE91", "TP2" = "#FB6A4A",
                  "TP3" = "#DE2D26", "TP4" = "#A50F15"),
    CV = colorRampPalette(c("#F7FCF5", "#74C476", "#00441B"))(100)
  )

  # Find gap position between Broad and Narrow samples
  n_broad <- sum(sample_anno$Leaf_Type == "Broad")

  # Save as PDF (publication quality)
  pdf("03_results/figures/supplementary/FigS9_housekeeping_stability_heatmap.pdf",
      width = 14, height = 12)

  pheatmap(hk_scaled,
           annotation_col = sample_anno,
           annotation_row = gene_anno,
           annotation_colors = anno_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-2, 2, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           fontsize_row = 6,
           fontsize = 10,
           main = "Top 50 Most Stable Housekeeping Genes\n(Validation against Machado et al. 2020 reference)",
           gaps_col = n_broad,
           border_color = NA)

  dev.off()

  # Also save as PNG
  png("03_results/figures/supplementary/FigS9_housekeeping_stability_heatmap.png",
      width = 14, height = 12, units = "in", res = 300)

  pheatmap(hk_scaled,
           annotation_col = sample_anno,
           annotation_row = gene_anno,
           annotation_colors = anno_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-2, 2, length.out = 101),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           fontsize_row = 6,
           fontsize = 10,
           main = "Top 50 Most Stable Housekeeping Genes\n(Validation against Machado et al. 2020 reference)",
           gaps_col = n_broad,
           border_color = NA)

  dev.off()

  cat("  Saved: 03_results/figures/supplementary/FigS9_housekeeping_stability_heatmap.pdf\n")
  cat("  Saved: 03_results/figures/supplementary/FigS9_housekeeping_stability_heatmap.png\n")

  # Also keep original full heatmap for reference
  cat("\nCreating full housekeeping gene heatmap (all genes)...\n")

  sample_order_tp <- order(targets$Timepoint, decreasing = FALSE)
  hk_expr_ordered <- hk_expr[, sample_order_tp]
  color_scale <- gplots::colorpanel(100, "blue", "white", "red")
  tp_counts <- table(targets$Timepoint[sample_order_tp])
  col_seps <- cumsum(tp_counts)[-length(tp_counts)]

  png("03_results/figures/09_expression/housekeeping_heatmap_all.png",
      width = 10, height = 8, units = "in", res = 300)

  heatmap.2(as.matrix(hk_expr_ordered),
            col = color_scale,
            density.info = "none",
            trace = "none",
            labRow = "",
            colsep = col_seps,
            margins = c(7, 2),
            key.xlab = "log2(CPM)",
            Colv = FALSE,
            dendrogram = "row",
            main = paste("All Housekeeping Genes (n =", nrow(hk_expr), ")"),
            cexCol = 0.6)

  dev.off()
  cat("  Saved: 03_results/figures/09_expression/housekeeping_heatmap_all.png\n")

  # --- Plot 8: Histogram of housekeeping gene expression ---
  cat("Creating housekeeping expression histogram...\n")

  png("03_results/figures/09_expression/housekeeping_histogram.png",
      width = 8, height = 6, units = "in", res = 300)

  hist(as.matrix(hk_expr),
       main = "Distribution of Housekeeping Gene Expression",
       xlab = "log2(CPM)",
       col = "#6BAED6",
       breaks = 50)

  dev.off()
  cat("  Saved: 03_results/figures/09_expression/housekeeping_histogram.png\n")

} else {
  cat("  No housekeeping genes found\n")
}


cat("\n========================================\n")

# Calculate mean expression per timepoint
temporal_means <- expr_long %>%
  group_by(Gene, Timepoint) %>%
  summarize(Mean = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = Mean)

# Calculate fold changes relative to TP0
if ("TP0" %in% colnames(temporal_means)) {
  cat("Calculating temporal fold changes (relative to TP0)...\n")

  temporal_fc <- temporal_means %>%
    mutate(
      FC_TP1 = TP1 - TP0,
      FC_TP2 = TP2 - TP0,
      FC_TP3 = TP3 - TP0,
      FC_TP4 = TP4 - TP0
    )

  # Identify genes with strongest temporal changes
  temporal_fc <- temporal_fc %>%
    mutate(Max_FC = pmax(abs(FC_TP1), abs(FC_TP2), abs(FC_TP3), abs(FC_TP4), na.rm = TRUE))

  top_temporal <- temporal_fc %>%
    arrange(desc(Max_FC)) %>%
    head(100)

  cat("  Top 100 temporally dynamic genes identified\n")
  cat("  Maximum fold change range:", round(min(top_temporal$Max_FC), 2),
      "to", round(max(top_temporal$Max_FC), 2), "log2\n\n")

  # --- Plot 7: Temporal trends heatmap ---
  cat("Creating temporal trends heatmap...\n")

  top_genes <- head(top_temporal$Gene, 50)
  temporal_mat <- as.matrix(temporal_means[temporal_means$Gene %in% top_genes, -1])
  rownames(temporal_mat) <- temporal_means$Gene[temporal_means$Gene %in% top_genes]

  # Scale rows for visualization
  temporal_scaled <- t(scale(t(temporal_mat)))

  png("03_results/figures/09_expression/temporal_trends_heatmap.png",
      width = 8, height = 12, units = "in", res = 300)

  pheatmap(temporal_scaled,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           main = "Top 50 Temporally Dynamic Genes",
           fontsize = 12,
           angle_col = 0)

  dev.off()
  cat("  Saved: 03_results/figures/09_expression/temporal_trends_heatmap.png\n")

  # Save temporal dynamics table
  write.csv(top_temporal %>% select(Gene, TP0:TP4, FC_TP1:FC_TP4, Max_FC),
            "03_results/tables/temporal_dynamics.csv", row.names = FALSE)
  cat("  Saved: 03_results/tables/temporal_dynamics.csv\n")
}


cat("\n========================================\n")

# Calculate overall gene statistics
gene_stats <- data.frame(
  Gene = rownames(expr_mat),
  Mean = rowMeans(expr_mat),
  SD = apply(expr_mat, 1, sd),
  Max = apply(expr_mat, 1, max),
  Min = apply(expr_mat, 1, min),
  Median = apply(expr_mat, 1, median)
)

# Categorize expression levels
gene_stats <- gene_stats %>%
  mutate(
    Expression_Level = case_when(
      Mean < 0 ~ "Very Low",
      Mean < 3 ~ "Low",
      Mean < 6 ~ "Medium",
      Mean < 9 ~ "High",
      TRUE ~ "Very High"
    )
  )

gene_stats$Expression_Level <- factor(gene_stats$Expression_Level,
                                       levels = c("Very Low", "Low", "Medium", "High", "Very High"))

# Summary
cat("Gene expression level distribution:\n")
print(table(gene_stats$Expression_Level))

# --- Plot 8: Expression level distribution ---
cat("\nCreating expression level distribution plot...\n")

expr_level_colors <- c(
  "Very Low" = "#FEE5D9",
  "Low" = "#FCAE91",
  "Medium" = "#FB6A4A",
  "High" = "#DE2D26",
  "Very High" = "#A50F15"
)

p8 <- ggplot(gene_stats, aes(x = Mean, fill = Expression_Level)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  scale_fill_manual(values = expr_level_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Distribution of Mean Gene Expression",
       subtitle = paste(nrow(gene_stats), "genes"),
       x = "Mean log2(CPM)",
       y = "Number of Genes") +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/09_expression/expression_level_distribution.png", p8,
       width = 10, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/09_expression/expression_level_distribution.png\n")

# --- Plot 9: Mean vs CV plot ---
cat("Creating mean-variance relationship plot...\n")

gene_stats <- gene_stats %>%
  mutate(CV = SD / Mean * 100)

p9 <- ggplot(gene_stats %>% filter(Mean > 0), aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", color = "#D95F02", se = FALSE) +
  theme_bw(base_size = 14) +
  labs(title = "Mean-CV Relationship",
       subtitle = "Coefficient of variation vs mean expression",
       x = "Mean log2(CPM)",
       y = "CV (%)") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(limits = c(0, quantile(gene_stats$CV[gene_stats$Mean > 0], 0.99, na.rm = TRUE)))

ggsave("03_results/figures/09_expression/mean_cv_relationship.png", p9,
       width = 8, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/09_expression/mean_cv_relationship.png\n")

# Save gene statistics
write.csv(gene_stats, "03_results/tables/gene_expression_stats.csv", row.names = FALSE)
cat("\nSaved: 03_results/tables/gene_expression_stats.csv\n")


cat("\n========================================\n")

# Get top 50 most highly expressed genes
top_expressed <- gene_stats %>%
  arrange(desc(Mean)) %>%
  head(50)

cat("Top 20 most highly expressed genes:\n")
print(top_expressed %>% select(Gene, Mean, SD, Expression_Level) %>% head(20))

# --- Plot 10: Top expressed genes heatmap ---
cat("\nCreating top expressed genes heatmap...\n")

top_genes <- top_expressed$Gene[1:30]
top_expr <- expr_mat[top_genes, ]

# Annotation
sample_anno <- data.frame(
  Timepoint = targets$Timepoint,
  Leaf_type = targets$Leaf_type,
  row.names = targets$Sample
)

anno_colors <- list(
  Timepoint = timepoint_colors,
  Leaf_type = leaf_colors
)

png("03_results/figures/09_expression/top_expressed_heatmap.png",
    width = 14, height = 10, units = "in", res = 300)

pheatmap(top_expr,
         color = colorRampPalette(c("#FFFFCC", "#FD8D3C", "#800026"))(100),
         annotation_col = sample_anno,
         annotation_colors = anno_colors,
         show_colnames = FALSE,
         fontsize = 10,
         main = "Top 30 Most Highly Expressed Genes")

dev.off()
cat("  Saved: 03_results/figures/09_expression/top_expressed_heatmap.png\n")

# Save top expressed genes
write.csv(top_expressed, "03_results/tables/top_expressed_genes.csv", row.names = FALSE)
cat("  Saved: 03_results/tables/top_expressed_genes.csv\n")


cat("\n========================================\n")
cat("SESSION INFO\n")

print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Figures Generated:\n")
cat("    - expression_density.png\n")
cat("    - expression_boxplot.png\n")
cat("    - JAG_expression_profile.png\n")
cat("    - JAG_expression_by_line.png\n")
cat("    - housekeeping_stability.png\n")
cat("    - housekeeping_heatmap_all.png\n")
cat("    - housekeeping_histogram.png\n")
cat("    - temporal_trends_heatmap.png\n")
cat("    - expression_level_distribution.png\n")
cat("    - mean_cv_relationship.png\n")
cat("    - top_expressed_heatmap.png\n")
cat("  Supplementary Figures:\n")
cat("    - FigS9_housekeeping_stability_heatmap.pdf/png (Top 50 stable genes)\n")
cat("\n")
cat("  Tables Generated:\n")
cat("    - expression_summary.csv\n")
cat("    - JAG_expression.csv\n")
cat("    - housekeeping_stability.csv\n")
cat("    - temporal_dynamics.csv\n")
cat("    - gene_expression_stats.csv\n")
cat("    - top_expressed_genes.csv\n")
cat("\n")
