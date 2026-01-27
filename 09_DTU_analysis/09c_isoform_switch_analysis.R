
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 42: ISOFORM SWITCH ANALYSIS\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/42_isoform_switches", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "dplyr",            # Data manipulation
  "tidyr",            # Data reshaping
  "ggplot2",          # Plotting
  "patchwork",        # Combining plots
  "scales"            # Scale functions
)

# Try to load IsoformSwitchAnalyzeR if available
isar_available <- requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)

if (!isar_available) {
  cat("  Note: IsoformSwitchAnalyzeR not installed.\n")
  cat("  Will use custom isoform switching analysis.\n")
  cat("  To install: BiocManager::install('IsoformSwitchAnalyzeR')\n")
}

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
# Load DTU results
load("03_results/checkpoints/41_dtu_drimseq.RData")
cat("Loaded DTU results from script 41\n")

# Load prepared data
load("03_results/checkpoints/40_dtu_prepared.RData")
cat("Loaded prepared transcript data from script 40\n")

cat("  Transcripts:", nrow(tx_counts_filt), "\n")
cat("  Genes:", length(unique(tx2gene_filt$GENEID)), "\n")
cat("  Significant DTU genes:", sum(res_gene_final$significant_DTU, na.rm = TRUE), "\n\n")
# An isoform switch is when the dominant transcript changes between conditions
# We'll identify switches based on:
# 1. Significant DTU (from DRIMSeq)
# 2. Change in dominant isoform OR large proportion change (>= 10%)

# Calculate mean proportions by condition
cat("Calculating mean transcript proportions by condition...\n")

broad_samples <- targets$Sample[targets$Leaf_type == "Broad"]
narrow_samples <- targets$Sample[targets$Leaf_type == "Narrow"]

prop_broad <- rowMeans(tx_proportions[, broad_samples], na.rm = TRUE)
prop_narrow <- rowMeans(tx_proportions[, narrow_samples], na.rm = TRUE)

# Calculate proportion difference (delta proportion)
delta_prop <- prop_narrow - prop_broad

# Create switch analysis data frame
switch_analysis <- data.frame(
  transcript_id = names(prop_broad),
  gene_id = tx2gene_filt$GENEID[match(names(prop_broad), tx2gene_filt$TXNAME)],
  prop_broad = prop_broad,
  prop_narrow = prop_narrow,
  delta_prop = delta_prop,
  abs_delta_prop = abs(delta_prop),
  stringsAsFactors = FALSE
)

# Add DTU significance
switch_analysis$dtu_pvalue <- res_tx_final$pvalue[match(switch_analysis$transcript_id, res_tx_final$feature_id)]
switch_analysis$dtu_padj <- res_tx_final$padj_stageR[match(switch_analysis$transcript_id, res_tx_final$feature_id)]

# Identify dominant isoform in each condition
cat("Identifying dominant isoforms...\n")

dominant_broad <- switch_analysis %>%
  group_by(gene_id) %>%
  slice_max(prop_broad, n = 1, with_ties = FALSE) %>%
  select(gene_id, dominant_broad = transcript_id, prop_broad_dominant = prop_broad)

dominant_narrow <- switch_analysis %>%
  group_by(gene_id) %>%
  slice_max(prop_narrow, n = 1, with_ties = FALSE) %>%
  select(gene_id, dominant_narrow = transcript_id, prop_narrow_dominant = prop_narrow)

# Merge dominant info
switch_analysis <- switch_analysis %>%
  left_join(dominant_broad, by = "gene_id") %>%
  left_join(dominant_narrow, by = "gene_id")

# Identify genes with isoform switch (dominant isoform changes)
genes_with_switch <- switch_analysis %>%
  filter(dominant_broad != dominant_narrow) %>%
  select(gene_id) %>%
  distinct()

cat("  Genes with dominant isoform switch:", nrow(genes_with_switch), "\n")

# Mark switches
switch_analysis$is_switch_gene <- switch_analysis$gene_id %in% genes_with_switch$gene_id
# Criteria for significant isoform switch:
# 1. Gene has significant DTU (stageR OFDR < 0.05)
# 2. Either:
#    a. Dominant isoform changes between conditions, OR
#    b. Transcript proportion change >= 10% (absolute)

sig_dtu_genes <- res_gene_final$gene_id[res_gene_final$significant_DTU == TRUE]

# Define significant switches
switch_analysis$sig_dtu_gene <- switch_analysis$gene_id %in% sig_dtu_genes

significant_switches <- switch_analysis %>%
  filter(sig_dtu_gene) %>%
  filter(is_switch_gene | abs_delta_prop >= 0.10) %>%
  arrange(desc(abs_delta_prop))

cat("Significant isoform switches:\n")
cat("  Transcripts involved:", nrow(significant_switches), "\n")
cat("  Genes involved:", length(unique(significant_switches$gene_id)), "\n")
cat("  Genes with dominant switch:", sum(significant_switches$is_switch_gene), "\n\n")

# Add switch direction
significant_switches$switch_direction <- ifelse(
  significant_switches$delta_prop > 0,
  "Increased in Narrow",
  "Decreased in Narrow"
)

# Summary by direction
cat("Switch direction summary:\n")
print(table(significant_switches$switch_direction))
cat("\n")
# JAG1 gene: Glyma.20G116200
jag1_gene <- "Glyma.20G116200"
jag1_isoforms <- switch_analysis %>%
  filter(gene_id == jag1_gene) %>%
  arrange(desc(abs_delta_prop))

if (nrow(jag1_isoforms) > 0) {
  cat("JAG1 (", jag1_gene, ") isoform analysis:\n", sep = "")
  cat("  Number of isoforms:", nrow(jag1_isoforms), "\n\n")

  cat("  Isoform proportions:\n")
  for (i in 1:nrow(jag1_isoforms)) {
    row <- jag1_isoforms[i, ]
    cat("    ", row$transcript_id, ":\n", sep = "")
    cat("      Broad:  ", round(row$prop_broad * 100, 1), "%\n", sep = "")
    cat("      Narrow: ", round(row$prop_narrow * 100, 1), "%\n", sep = "")
    cat("      Change: ", ifelse(row$delta_prop > 0, "+", ""),
        round(row$delta_prop * 100, 1), " percentage points\n", sep = "")
  }

  # Check if JAG1 has dominant isoform switch
  if (jag1_isoforms$dominant_broad[1] != jag1_isoforms$dominant_narrow[1]) {
    cat("\n  ISOFORM SWITCH DETECTED!\n")
    cat("    Dominant in Broad:  ", jag1_isoforms$dominant_broad[1], "\n", sep = "")
    cat("    Dominant in Narrow: ", jag1_isoforms$dominant_narrow[1], "\n", sep = "")
  } else {
    cat("\n  No dominant isoform switch\n")
    cat("    Dominant isoform in both: ", jag1_isoforms$dominant_broad[1], "\n", sep = "")
  }

  # Check DTU significance
  jag1_dtu <- jag1_gene %in% sig_dtu_genes
  cat("  DTU significant: ", ifelse(jag1_dtu, "YES", "NO"), "\n\n", sep = "")
} else {
  cat("JAG1 not found in isoform analysis data\n\n")
}
if (!is.null(jag1_targets)) {
  # Add JAG1 target info to switch analysis
  significant_switches$is_JAG1_target <- significant_switches$gene_id %in% jag1_targets$GeneID
  significant_switches$JAG1_tier <- jag1_targets$Confidence_Tier[
    match(significant_switches$gene_id, jag1_targets$GeneID)
  ]
  significant_switches$JAG1_tier[is.na(significant_switches$JAG1_tier)] <- "Not_Target"

  # Summary by tier
  cat("Isoform switches by JAG1 target tier:\n")
  tier_switch_summary <- significant_switches %>%
    group_by(JAG1_tier) %>%
    summarise(
      n_transcripts = n(),
      n_genes = n_distinct(gene_id),
      .groups = "drop"
    )
  print(tier_switch_summary)
  cat("\n")

  # Top switches in JAG1 targets
  jag1_switches <- significant_switches %>%
    filter(is_JAG1_target) %>%
    arrange(desc(abs_delta_prop)) %>%
    select(gene_id, transcript_id, JAG1_tier, prop_broad, prop_narrow, delta_prop, is_switch_gene)

  cat("Top 15 isoform switches in JAG1 targets:\n")
  print(head(jag1_switches, 15))
  cat("\n")
}
# Analyze isoform proportions across timepoints
timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")

# Calculate mean proportions per timepoint and condition
temporal_props <- list()

for (tp in timepoints) {
  tp_samples_broad <- targets$Sample[targets$Leaf_type == "Broad" & targets$Timepoint == tp]
  tp_samples_narrow <- targets$Sample[targets$Leaf_type == "Narrow" & targets$Timepoint == tp]

  if (length(tp_samples_broad) > 0 && length(tp_samples_narrow) > 0) {
    temporal_props[[paste0(tp, "_Broad")]] <- rowMeans(tx_proportions[, tp_samples_broad, drop = FALSE], na.rm = TRUE)
    temporal_props[[paste0(tp, "_Narrow")]] <- rowMeans(tx_proportions[, tp_samples_narrow, drop = FALSE], na.rm = TRUE)
  }
}

# Convert to data frame
if (length(temporal_props) > 0) {
  temporal_df <- as.data.frame(temporal_props)
  temporal_df$transcript_id <- rownames(tx_proportions)
  temporal_df$gene_id <- tx2gene_filt$GENEID[match(temporal_df$transcript_id, tx2gene_filt$TXNAME)]

  # Reshape for plotting
  temporal_long <- temporal_df %>%
    pivot_longer(
      cols = -c(transcript_id, gene_id),
      names_to = "condition",
      values_to = "proportion"
    ) %>%
    separate(condition, into = c("timepoint", "leaf_type"), sep = "_") %>%
    mutate(timepoint = factor(timepoint, levels = timepoints))

  cat("Created temporal isoform dynamics data\n")
  cat("  Observations:", nrow(temporal_long), "\n\n")
}
# Plot 1: Distribution of proportion changes
p1 <- ggplot(switch_analysis, aes(x = delta_prop)) +
  geom_histogram(binwidth = 0.02, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "darkred", linewidth = 0.6) +
  labs(
    title = "Distribution of Transcript Proportion Changes",
    subtitle = "Narrow - Broad leaf genotypes",
    x = expression(Delta~"Proportion (Narrow - Broad)"),
    y = "Number of Transcripts"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("03_results/figures/42_isoform_switches/proportion_change_distribution.pdf", p1, width = 8, height = 6)
cat("Saved: proportion_change_distribution.pdf\n")

# Plot 2: Volcano-style plot (proportion change vs significance)
switch_analysis$neg_log10_p <- -log10(switch_analysis$dtu_pvalue + 1e-300)

p2 <- ggplot(switch_analysis, aes(x = delta_prop, y = neg_log10_p)) +
  geom_point(aes(color = abs_delta_prop >= 0.1 & sig_dtu_gene),
             alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#CCCCCC"),
                     labels = c("TRUE" = "Significant Switch", "FALSE" = "Not Significant"),
                     name = NULL) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(
    title = "Isoform Switch Volcano Plot",
    subtitle = "Transcript proportion change vs DTU significance",
    x = expression(Delta~"Proportion (Narrow - Broad)"),
    y = expression(-log[10]~"(p-value)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("03_results/figures/42_isoform_switches/switch_volcano.pdf", p2, width = 8, height = 7)
cat("Saved: switch_volcano.pdf\n")

# Plot 3: JAG1 isoform proportions (if available)
if (nrow(jag1_isoforms) > 0 && exists("temporal_long")) {
  jag1_temporal <- temporal_long %>%
    filter(gene_id == jag1_gene)

  if (nrow(jag1_temporal) > 0) {
    p3 <- ggplot(jag1_temporal, aes(x = timepoint, y = proportion,
                                     fill = transcript_id, group = transcript_id)) +
      geom_area(alpha = 0.7, position = "stack") +
      facet_wrap(~leaf_type) +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      scale_fill_brewer(palette = "Set2", name = "Transcript") +
      labs(
        title = paste("JAG1 Isoform Usage Across Development"),
        subtitle = jag1_gene,
        x = "Developmental Timepoint",
        y = "Transcript Proportion"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"
      )

    ggsave("03_results/figures/42_isoform_switches/JAG1_temporal_isoforms.pdf", p3, width = 10, height = 6)
    cat("Saved: JAG1_temporal_isoforms.pdf\n")
  }
}

# Plot 4: Top 10 genes with largest isoform switches
top_switch_genes <- significant_switches %>%
  group_by(gene_id) %>%
  summarise(max_delta = max(abs_delta_prop), .groups = "drop") %>%
  arrange(desc(max_delta)) %>%
  head(10) %>%
  pull(gene_id)

top_switches_data <- significant_switches %>%
  filter(gene_id %in% top_switch_genes) %>%
  mutate(gene_id = factor(gene_id, levels = top_switch_genes))

p4 <- ggplot(top_switches_data, aes(x = prop_broad, xend = prop_narrow,
                                      y = reorder(transcript_id, abs_delta_prop),
                                      yend = reorder(transcript_id, abs_delta_prop))) +
  geom_segment(aes(color = switch_direction), linewidth = 1.5,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_point(aes(x = prop_broad), shape = 21, fill = "#2E86AB", size = 3) +
  geom_point(aes(x = prop_narrow), shape = 21, fill = "#E94F37", size = 3) +
  scale_color_manual(values = c("Increased in Narrow" = "#E94F37",
                                "Decreased in Narrow" = "#2E86AB")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~gene_id, scales = "free_y", ncol = 2) +
  labs(
    title = "Top 10 Genes with Largest Isoform Switches",
    subtitle = "Transcript proportion change (Broad → Narrow)",
    x = "Transcript Proportion",
    y = NULL,
    color = "Direction"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

ggsave("03_results/figures/42_isoform_switches/top_isoform_switches.pdf", p4, width = 12, height = 10)
cat("Saved: top_isoform_switches.pdf\n")

cat("\n")
# Save checkpoint
save(
  switch_analysis,        # All switch analysis data
  significant_switches,   # Significant switches only
  jag1_isoforms,          # JAG1 isoform analysis
  temporal_long,          # Temporal isoform dynamics (if created)
  file = "03_results/checkpoints/42_isoform_switches.RData"
)
cat("Saved checkpoint: 42_isoform_switches.RData\n")

# Save significant switches table
write.csv(significant_switches, "03_results/tables/isoform_switches_significant.csv", row.names = FALSE)
cat("Saved: isoform_switches_significant.csv (", nrow(significant_switches), " transcripts)\n", sep = "")

# Save full switch analysis
write.csv(switch_analysis, "03_results/tables/isoform_switch_analysis_full.csv", row.names = FALSE)
cat("Saved: isoform_switch_analysis_full.csv\n")

# Save JAG1 isoform analysis
if (nrow(jag1_isoforms) > 0) {
  write.csv(jag1_isoforms, "03_results/tables/JAG1_isoform_analysis.csv", row.names = FALSE)
  cat("Saved: JAG1_isoform_analysis.csv\n")
}

cat("\n")


cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("SUMMARY:\n")
cat("  - Analyzed", nrow(switch_analysis), "transcripts for isoform switches\n")
cat("  - Found", nrow(genes_with_switch), "genes with dominant isoform switch\n")
cat("  - Identified", nrow(significant_switches), "significant switch transcripts\n")
cat("  - Genes with significant switches:", length(unique(significant_switches$gene_id)), "\n")
if (nrow(jag1_isoforms) > 0) {
  cat("  - JAG1 has", nrow(jag1_isoforms), "isoforms\n")
}
cat("\n")

