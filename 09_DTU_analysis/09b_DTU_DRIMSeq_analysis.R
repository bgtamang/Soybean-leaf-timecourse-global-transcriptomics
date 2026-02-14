# Script 41: DTU Analysis with DRIMSeq + stageR
# REVISED: TP0-focused analysis of JAG1 targets

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 41: DTU ANALYSIS (DRIMSeq + stageR)\n")
cat("  TP0-FOCUSED ANALYSIS OF JAG1 TARGETS\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/41_DTU_analysis", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "DRIMSeq",          # DTU analysis
  "stageR",           # Stage-wise testing
  "dplyr",            # Data manipulation
  "tidyr",            # Data reshaping
  "ggplot2",          # Plotting
  "BiocParallel"      # Parallel processing
)

bioc_packages <- c("DRIMSeq", "stageR", "BiocParallel")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# Set up parallel processing (Windows compatible)
if (.Platform$OS.type == "windows") {
  BPPARAM <- SnowParam(workers = 1)
} else {
  BPPARAM <- MulticoreParam(workers = 1)
}

# ===== SECTION 1: LOAD DATA =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/40_dtu_prepared.RData")
cat("Loaded checkpoint from script 40\n")

# Use UNFILTERED data to preserve TP0-specific genes
cat("\nUsing UNFILTERED transcript data to preserve TP0-specific genes\n")
cat("  Total transcripts (unfiltered):", nrow(tx_counts), "\n")
cat("  Total genes:", length(unique(tx2gene$GENEID)), "\n")
cat("  Total samples:", ncol(tx_counts), "\n\n")

# Load JAG1 targets
cat("Loading JAG1 targets...\n")
jag1_targets_file <- "03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv"
if (file.exists(jag1_targets_file)) {
  jag1_targets_df <- read.csv(jag1_targets_file, stringsAsFactors = FALSE)
  cat("  Loaded:", nrow(jag1_targets_df), "JAG1 targets\n")

  # Check column names and adapt
  if ("Confidence_Tier" %in% names(jag1_targets_df)) {
    tier_col <- "Confidence_Tier"
  } else if ("JAG1_tier" %in% names(jag1_targets_df)) {
    tier_col <- "JAG1_tier"
  } else {
    # Find tier column
    tier_col <- grep("tier|Tier", names(jag1_targets_df), value = TRUE)[1]
  }

  if (!is.null(tier_col) && !is.na(tier_col)) {
    cat("  Gold:", sum(jag1_targets_df[[tier_col]] == "Gold", na.rm = TRUE), "\n")
    cat("  Silver:", sum(jag1_targets_df[[tier_col]] == "Silver", na.rm = TRUE), "\n")
    cat("  Bronze:", sum(jag1_targets_df[[tier_col]] == "Bronze", na.rm = TRUE), "\n")
    # Standardize column name
    jag1_targets_df$Confidence_Tier <- jag1_targets_df[[tier_col]]
  }

  # Standardize GeneID column
  if (!"GeneID" %in% names(jag1_targets_df)) {
    if ("gene_id" %in% names(jag1_targets_df)) {
      jag1_targets_df$GeneID <- jag1_targets_df$gene_id
    } else {
      gene_col <- grep("gene|Gene", names(jag1_targets_df), value = TRUE)[1]
      if (!is.null(gene_col)) jag1_targets_df$GeneID <- jag1_targets_df[[gene_col]]
    }
  }
} else {
  stop("JAG1 targets file not found: ", jag1_targets_file)
}

cat("\n")

# ===== SECTION 2: SUBSET TO TP0 SAMPLES =====

cat("========================================\n")
cat("SECTION 2: SUBSET TO TP0 SAMPLES\n")
cat("========================================\n\n")

# Filter to TP0 samples only
tp0_samples <- targets$Sample[targets$Timepoint == "TP0"]
cat("TP0 samples:", length(tp0_samples), "\n")

# Check distribution
tp0_metadata <- targets[targets$Timepoint == "TP0", ]
cat("  Broad leaf:", sum(tp0_metadata$Leaf_type == "Broad"), "samples\n")
cat("  Narrow leaf:", sum(tp0_metadata$Leaf_type == "Narrow"), "samples\n")

cat("\nSample breakdown:\n")
print(table(tp0_metadata$Line, tp0_metadata$Leaf_type))

# Subset counts to TP0
tx_counts_tp0 <- tx_counts[, tp0_samples]
cat("\nTP0 count matrix:", nrow(tx_counts_tp0), "transcripts x", ncol(tx_counts_tp0), "samples\n\n")

# ===== SECTION 3: FILTER FOR TP0 EXPRESSION =====

cat("========================================\n")
cat("SECTION 3: FILTER FOR TP0 EXPRESSION\n")
cat("========================================\n\n")

# TP0-specific filtering:
# - Transcript must have >= 5 counts in >= 2 samples at TP0
# - Gene must have >= 2 expressed transcripts

min_counts_tp0 <- 5
min_samples_tp0 <- 2

cat("TP0-specific filter parameters:\n")
cat("  Minimum counts:", min_counts_tp0, "\n")
cat("  Minimum samples:", min_samples_tp0, "(out of", ncol(tx_counts_tp0), ")\n\n")

# Step 1: Filter transcripts expressed at TP0
expressed_tp0 <- rowSums(tx_counts_tp0 >= min_counts_tp0) >= min_samples_tp0
cat("Transcripts expressed at TP0:", sum(expressed_tp0), "/", nrow(tx_counts_tp0), "\n")

tx_counts_tp0_filt <- tx_counts_tp0[expressed_tp0, ]

# Update tx2gene for filtered transcripts
tx2gene_tp0 <- tx2gene[tx2gene$TXNAME %in% rownames(tx_counts_tp0_filt), ]

# Step 2: Keep genes with >= 2 transcripts
tx_per_gene <- tx2gene_tp0 %>%
  group_by(GENEID) %>%
  summarise(n_tx = n(), .groups = "drop")

genes_multi_tx <- tx_per_gene$GENEID[tx_per_gene$n_tx >= 2]
cat("Genes with >= 2 transcripts at TP0:", length(genes_multi_tx), "\n")

# Final filter
keep_tx <- tx2gene_tp0$TXNAME[tx2gene_tp0$GENEID %in% genes_multi_tx]
tx_counts_tp0_filt <- tx_counts_tp0_filt[keep_tx, ]
tx2gene_tp0_filt <- tx2gene_tp0[tx2gene_tp0$TXNAME %in% keep_tx, ]

cat("\nAfter TP0 filtering:\n")
cat("  Transcripts:", nrow(tx_counts_tp0_filt), "\n")
cat("  Genes:", length(unique(tx2gene_tp0_filt$GENEID)), "\n")

# Check JAG1 targets overlap
jag1_in_dtu <- sum(unique(tx2gene_tp0_filt$GENEID) %in% jag1_targets_df$GeneID)
cat("  JAG1 targets in dataset:", jag1_in_dtu, "\n\n")

# ===== SECTION 4: PREPARE DRIMSeq DATA =====

cat("========================================\n")
cat("SECTION 4: PREPARE DRIMSeq DATA\n")
cat("========================================\n\n")

# Prepare counts data frame for DRIMSeq
counts_df <- data.frame(
  gene_id = tx2gene_tp0_filt$GENEID[match(rownames(tx_counts_tp0_filt), tx2gene_tp0_filt$TXNAME)],
  feature_id = rownames(tx_counts_tp0_filt),
  tx_counts_tp0_filt,
  check.names = FALSE
)

cat("Counts data frame:", nrow(counts_df), "transcripts\n")

# Prepare samples data frame
samples_df <- data.frame(
  sample_id = tp0_metadata$Sample,
  condition = tp0_metadata$Leaf_type,
  line = tp0_metadata$Line,
  stringsAsFactors = FALSE
)

# Factor with Broad as reference
samples_df$condition <- factor(samples_df$condition, levels = c("Broad", "Narrow"))
samples_df$line <- factor(samples_df$line)

cat("Samples:", nrow(samples_df), "\n")
cat("  Broad:", sum(samples_df$condition == "Broad"), "\n")
cat("  Narrow:", sum(samples_df$condition == "Narrow"), "\n")
cat("\nSamples by Line and Condition:\n")
print(table(samples_df$line, samples_df$condition))
cat("\n")

# Create DRIMSeq data object
cat("Creating dmDSdata object...\n")
d <- dmDSdata(counts = counts_df, samples = samples_df)
cat("  Features:", nrow(counts(d)), "\n")
cat("  Genes:", length(unique(counts(d)$gene_id)), "\n\n")

# ===== SECTION 5: DRIMSeq FILTERING =====

cat("========================================\n")
cat("SECTION 5: DRIMSeq FILTERING\n")
cat("========================================\n\n")

# Relaxed filtering for TP0 (12 samples only)
cat("DRIMSeq filter parameters (relaxed for 12 samples):\n")
cat("  min_samps_gene_expr: 3\n")
cat("  min_samps_feature_expr: 2\n")
cat("  min_gene_expr: 10\n")
cat("  min_feature_expr: 5\n")
cat("  min_samps_feature_prop: 2\n")
cat("  min_feature_prop: 0.05\n\n")

d_filt <- dmFilter(d,
  min_samps_gene_expr = 3,
  min_samps_feature_expr = 2,
  min_gene_expr = 10,
  min_feature_expr = 5,
  min_samps_feature_prop = 2,
  min_feature_prop = 0.05
)

cat("After DRIMSeq filtering:\n")
cat("  Features:", nrow(counts(d_filt)), "\n")
cat("  Genes:", length(unique(counts(d_filt)$gene_id)), "\n")

# Check JAG1 targets after filtering
genes_after_filter <- unique(counts(d_filt)$gene_id)
jag1_after_filter <- sum(genes_after_filter %in% jag1_targets_df$GeneID)
cat("  JAG1 targets retained:", jag1_after_filter, "\n\n")

# ===== SECTION 6: FIT DIRICHLET-MULTINOMIAL MODEL =====

cat("========================================\n")
cat("SECTION 6: FIT DIRICHLET-MULTINOMIAL MODEL\n")
cat("========================================\n\n")

# Design matrix
# NOTE: Line and Condition are perfectly confounded (each line is either Broad OR Narrow)
# So we cannot include both in the model. Use ~ condition only.
# The comparison is: 2 Narrow lines (PI612713B, PI547745) vs 2 Broad lines (PI532462A, LD11-2170)
# This treats within-condition lines as biological replicates.
design_full <- model.matrix(~ condition, data = samples_df)
cat("Design matrix (~ condition):\n")
cat("  - Tests Narrow vs Broad effect\n")
cat("  - Lines within each condition treated as biological replicates\n")
cat("  - Note: Line and Condition are confounded (each line is one genotype only)\n\n")
print(design_full)
cat("\n")

# Fit model
cat("Fitting Dirichlet-multinomial model...\n")
start_time <- Sys.time()

# Estimate precision
d_filt <- dmPrecision(d_filt, design = design_full, BPPARAM = BPPARAM)
cat("  Estimated precision parameters\n")

# Fit model
d_filt <- dmFit(d_filt, design = design_full, BPPARAM = BPPARAM)
cat("  Model fitted\n")

end_time <- Sys.time()
cat("  Time elapsed:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n\n")

# ===== SECTION 7: TEST FOR DTU =====

cat("========================================\n")
cat("SECTION 7: TEST FOR DTU\n")
cat("========================================\n\n")

# Test for DTU
cat("Testing for differential transcript usage (Narrow vs Broad at TP0)...\n")
d_filt <- dmTest(d_filt, coef = "conditionNarrow", BPPARAM = BPPARAM)

# Extract results
res_gene <- DRIMSeq::results(d_filt)
res_tx <- DRIMSeq::results(d_filt, level = "feature")

cat("DTU test results:\n")
cat("  Gene-level tests:", nrow(res_gene), "\n")
cat("  Transcript-level tests:", nrow(res_tx), "\n")

# Significant genes
res_gene$padj <- p.adjust(res_gene$pvalue, method = "BH")
sig_nominal <- sum(res_gene$pvalue < 0.05, na.rm = TRUE)
sig_bh <- sum(res_gene$padj < 0.05, na.rm = TRUE)
sig_bh_10 <- sum(res_gene$padj < 0.10, na.rm = TRUE)

cat("  Significant (p < 0.05, unadjusted):", sig_nominal, "\n")
cat("  Significant (padj < 0.05, BH):", sig_bh, "\n")
cat("  Significant (padj < 0.10, BH):", sig_bh_10, "\n\n")

# ===== SECTION 8: STAGE-WISE TESTING =====

cat("========================================\n")
cat("SECTION 8: STAGE-WISE TESTING (stageR)\n")
cat("========================================\n\n")

cat("Performing stage-wise testing...\n")

# Prepare p-values
pScreen <- res_gene$pvalue
names(pScreen) <- res_gene$gene_id

# Remove NA p-values
pScreen <- pScreen[!is.na(pScreen)]

# Confirmation p-values
pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
rownames(pConfirmation) <- res_tx$feature_id

# tx2gene for stageR
tx2gene_stager <- res_tx[, c("feature_id", "gene_id")]
colnames(tx2gene_stager) <- c("transcript", "gene")

# Filter to genes with valid p-values
tx2gene_stager <- tx2gene_stager[tx2gene_stager$gene %in% names(pScreen), ]
pConfirmation <- pConfirmation[rownames(pConfirmation) %in% tx2gene_stager$transcript, , drop = FALSE]

# Create stageR object
stageRObj <- stageRTx(
  pScreen = pScreen,
  pConfirmation = pConfirmation,
  pScreenAdjusted = FALSE,
  tx2gene = tx2gene_stager
)

# Stage-wise adjustment
stageRObj <- stageWiseAdjustment(stageRObj, method = "dtu", alpha = 0.05, allowNA = TRUE)

# Get results
dtu_genes <- getSignificantGenes(stageRObj)
dtu_transcripts <- getSignificantTx(stageRObj)

cat("Stage-wise testing results (OFDR = 0.05):\n")
cat("  Genes with significant DTU:", length(dtu_genes), "\n")
cat("  Transcripts with significant DTU:", length(dtu_transcripts), "\n")

# Also check at OFDR = 0.10
stageRObj_10 <- stageWiseAdjustment(
  stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
           pScreenAdjusted = FALSE, tx2gene = tx2gene_stager),
  method = "dtu", alpha = 0.10, allowNA = TRUE
)
dtu_genes_10 <- getSignificantGenes(stageRObj_10)

cat("  Genes with significant DTU (OFDR = 0.10):", length(dtu_genes_10), "\n\n")

# Get adjusted p-values
padj_gene <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)

# ===== SECTION 9: ANNOTATE RESULTS WITH JAG1 TARGET INFO =====

cat("========================================\n")
cat("SECTION 9: ANNOTATE WITH JAG1 TARGET INFO\n")
cat("========================================\n\n")

# Add stageR p-values and JAG1 info
res_gene_final <- res_gene
res_gene_final$padj_stageR <- padj_gene[match(res_gene_final$gene_id, rownames(padj_gene)), "gene"]
res_gene_final$significant_DTU <- res_gene_final$gene_id %in% dtu_genes
res_gene_final$significant_DTU_10 <- res_gene_final$gene_id %in% dtu_genes_10

# Add JAG1 target info
res_gene_final$is_JAG1_target <- res_gene_final$gene_id %in% jag1_targets_df$GeneID
res_gene_final$JAG1_tier <- jag1_targets_df$Confidence_Tier[match(res_gene_final$gene_id, jag1_targets_df$GeneID)]
res_gene_final$JAG1_tier[is.na(res_gene_final$JAG1_tier)] <- "Not_Target"

# Summary by JAG1 tier
cat("DTU results by JAG1 target tier:\n")
cat("-" , rep("-", 50), "\n", sep = "")

for (tier in c("Gold", "Silver", "Bronze", "Not_Target")) {
  tier_genes <- res_gene_final[res_gene_final$JAG1_tier == tier, ]
  n_total <- nrow(tier_genes)
  n_sig_05 <- sum(tier_genes$significant_DTU, na.rm = TRUE)
  n_sig_10 <- sum(tier_genes$significant_DTU_10, na.rm = TRUE)
  n_nominal <- sum(tier_genes$pvalue < 0.05, na.rm = TRUE)

  cat(sprintf("  %s: %d tested, %d sig (OFDR<0.05), %d sig (OFDR<0.10), %d nominal (p<0.05)\n",
              tier, n_total, n_sig_05, n_sig_10, n_nominal))
}

cat("\n")

# Order by significance
res_gene_final <- res_gene_final[order(res_gene_final$pvalue), ]

# ===== SECTION 10: JAG1 TARGET ENRICHMENT TEST =====

cat("========================================\n")
cat("SECTION 10: JAG1 TARGET ENRICHMENT IN DTU\n")
cat("========================================\n\n")

# Test if JAG1 targets are enriched among DTU genes
# Using nominal p < 0.05 for more power (small sample size)

sig_genes_nominal <- res_gene_final$gene_id[res_gene_final$pvalue < 0.05]
all_tested <- res_gene_final$gene_id

# Enrichment test for all JAG1 targets
jag1_all <- jag1_targets_df$GeneID

a <- sum(sig_genes_nominal %in% jag1_all)  # DTU & JAG1 target
b <- length(sig_genes_nominal) - a          # DTU & Not target
c <- sum(all_tested %in% jag1_all) - a      # Not DTU & JAG1 target
d <- length(all_tested) - a - b - c         # Not DTU & Not target

cat("JAG1 Target Enrichment (using nominal p < 0.05):\n")
cat("  Contingency table:\n")
cat("                    JAG1_Target  Not_Target\n")
cat(sprintf("  DTU (p<0.05):     %d           %d\n", a, b))
cat(sprintf("  Not DTU:          %d           %d\n", c, d))

fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
cat("\n  Fisher's exact test (one-sided):\n")
cat("    Odds ratio:", round(fisher_result$estimate, 2), "\n")
cat("    p-value:", format(fisher_result$p.value, digits = 3), "\n")

# Expected vs observed
expected <- (a + c) * (a + b) / length(all_tested)
cat("    Expected:", round(expected, 1), "\n")
cat("    Observed:", a, "\n")
cat("    Fold enrichment:", round(a / expected, 2), "\n\n")

# By tier
cat("Enrichment by JAG1 tier:\n")
enrichment_results <- list()

for (tier in c("Gold", "Silver", "Bronze")) {
  tier_genes <- jag1_targets_df$GeneID[jag1_targets_df$Confidence_Tier == tier]

  a_t <- sum(sig_genes_nominal %in% tier_genes)
  b_t <- length(sig_genes_nominal) - a_t
  c_t <- sum(all_tested %in% tier_genes) - a_t
  d_t <- length(all_tested) - a_t - b_t - c_t

  test_t <- fisher.test(matrix(c(a_t, b_t, c_t, d_t), nrow = 2), alternative = "greater")
  exp_t <- (a_t + c_t) * (a_t + b_t) / length(all_tested)

  enrichment_results[[tier]] <- data.frame(
    Tier = tier,
    Observed = a_t,
    Expected = round(exp_t, 1),
    Fold_Enrichment = round(a_t / max(exp_t, 0.1), 2),
    Pvalue = test_t$p.value
  )

  cat(sprintf("  %s: %d/%d DTU (expected %.1f, fold=%.2f, p=%.3f)\n",
              tier, a_t, sum(all_tested %in% tier_genes), exp_t,
              a_t / max(exp_t, 0.1), test_t$p.value))
}

enrichment_df <- do.call(rbind, enrichment_results)

cat("\n")

# ===== SECTION 11: TOP DTU GENES =====

cat("========================================\n")
cat("SECTION 11: TOP DTU GENES\n")
cat("========================================\n\n")

# Top DTU genes overall
cat("Top 20 DTU genes (by p-value):\n")
top_dtu <- head(res_gene_final, 20)
print(top_dtu[, c("gene_id", "pvalue", "padj", "JAG1_tier")])

cat("\n")

# Top DTU genes among JAG1 targets
cat("Top DTU genes among JAG1 targets:\n")
jag1_dtu <- res_gene_final[res_gene_final$is_JAG1_target, ]
jag1_dtu <- jag1_dtu[order(jag1_dtu$pvalue), ]
print(head(jag1_dtu[, c("gene_id", "pvalue", "padj", "JAG1_tier")], 20))

cat("\n")

# ===== SECTION 12: VISUALIZATION =====

cat("========================================\n")
cat("SECTION 12: VISUALIZATION\n")
cat("========================================\n\n")

# Plot 1: P-value histogram
p1 <- ggplot(res_gene_final, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.02, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  labs(
    title = "DTU P-value Distribution at TP0",
    subtitle = paste("DRIMSeq analysis - Narrow vs Broad,", nrow(res_gene_final), "genes tested"),
    x = "P-value",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("03_results/figures/41_DTU_analysis/pvalue_histogram_TP0.pdf", p1, width = 8, height = 6)
cat("Saved: pvalue_histogram_TP0.pdf\n")

# Plot 2: DTU by JAG1 tier
tier_summary <- res_gene_final %>%
  group_by(JAG1_tier) %>%
  summarise(
    n_tested = n(),
    n_dtu_nominal = sum(pvalue < 0.05, na.rm = TRUE),
    pct_dtu = n_dtu_nominal / n_tested * 100,
    .groups = "drop"
  ) %>%
  mutate(JAG1_tier = factor(JAG1_tier, levels = c("Gold", "Silver", "Bronze", "Not_Target")))

p2 <- ggplot(tier_summary, aes(x = JAG1_tier, y = pct_dtu, fill = JAG1_tier)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0(n_dtu_nominal, "/", n_tested)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Gold" = "#D4AF37", "Silver" = "#8E8E8E",
                                "Bronze" = "#CD7F32", "Not_Target" = "#CCCCCC")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "DTU Rate by JAG1 Target Tier at TP0",
    subtitle = "Percentage of genes with nominal p < 0.05",
    x = "JAG1 Target Tier",
    y = "% Genes with DTU"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("03_results/figures/41_DTU_analysis/dtu_by_jag1_tier_TP0.pdf", p2, width = 8, height = 6)
cat("Saved: dtu_by_jag1_tier_TP0.pdf\n")

# Plot 3: Enrichment bar plot
p3 <- ggplot(enrichment_df, aes(x = Tier, y = Fold_Enrichment, fill = Tier)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = sprintf("p=%.3f", Pvalue)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Gold" = "#D4AF37", "Silver" = "#8E8E8E", "Bronze" = "#CD7F32")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "JAG1 Target Enrichment in DTU Genes",
    subtitle = "Fold enrichment at TP0 (nominal p < 0.05)",
    x = "JAG1 Target Tier",
    y = "Fold Enrichment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("03_results/figures/41_DTU_analysis/jag1_enrichment_TP0.pdf", p3, width = 8, height = 6)
cat("Saved: jag1_enrichment_TP0.pdf\n")

cat("\n")

# ===== SECTION 13: SAVE RESULTS =====

cat("========================================\n")
cat("SECTION 13: SAVE RESULTS\n")
cat("========================================\n\n")

# Prepare transcript-level results
res_tx_final <- res_tx
res_tx_final$padj <- p.adjust(res_tx_final$pvalue, method = "BH")
res_tx_final$gene_id_full <- res_tx_final$gene_id
res_tx_final$is_JAG1_target <- res_tx_final$gene_id %in% jag1_targets_df$GeneID

# Save checkpoint
save(
  d_filt,              # DRIMSeq object
  res_gene_final,      # Gene-level results with JAG1 annotation
  res_tx_final,        # Transcript-level results
  stageRObj,           # stageR object
  enrichment_df,       # Enrichment results
  tx_counts_tp0_filt,  # TP0 filtered counts
  tx2gene_tp0_filt,    # TP0 tx2gene mapping
  tp0_metadata,        # TP0 sample metadata
  file = "03_results/checkpoints/41_dtu_drimseq.RData"
)
cat("Saved checkpoint: 41_dtu_drimseq.RData\n")

# Save gene-level results
write.csv(res_gene_final, "03_results/tables/dtu_gene_level_results_TP0.csv", row.names = FALSE)
cat("Saved: dtu_gene_level_results_TP0.csv\n")

# Save transcript-level results
write.csv(res_tx_final, "03_results/tables/dtu_transcript_level_results_TP0.csv", row.names = FALSE)
cat("Saved: dtu_transcript_level_results_TP0.csv\n")

# Save significant DTU genes (nominal p < 0.05)
dtu_significant <- res_gene_final[res_gene_final$pvalue < 0.05, ]
write.csv(dtu_significant, "03_results/tables/dtu_significant_genes_TP0.csv", row.names = FALSE)
cat("Saved: dtu_significant_genes_TP0.csv (", nrow(dtu_significant), " genes)\n", sep = "")

# Save JAG1 targets with DTU
jag1_dtu_results <- res_gene_final[res_gene_final$is_JAG1_target, ]
write.csv(jag1_dtu_results, "03_results/tables/dtu_JAG1_targets_TP0.csv", row.names = FALSE)
cat("Saved: dtu_JAG1_targets_TP0.csv (", nrow(jag1_dtu_results), " genes)\n", sep = "")

# Save enrichment results
write.csv(enrichment_df, "03_results/tables/dtu_jag1_enrichment_TP0.csv", row.names = FALSE)
cat("Saved: dtu_jag1_enrichment_TP0.csv\n\n")

# ===== COMPLETION =====

cat("================================================================\n")
cat("  SCRIPT 41 COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("SUMMARY (TP0 Analysis):\n")
cat("  - Samples analyzed:", ncol(tx_counts_tp0_filt), "(TP0 only)\n")
cat("  - Design: ~ condition (Narrow vs Broad; lines treated as replicates)\n")
cat("  - Genes tested for DTU:", nrow(res_gene_final), "\n")
cat("  - Significant DTU (nominal p < 0.05):", sum(res_gene_final$pvalue < 0.05, na.rm = TRUE), "\n")
cat("  - Significant DTU (BH padj < 0.05):", sum(res_gene_final$padj < 0.05, na.rm = TRUE), "\n")
cat("  - Significant DTU (stageR OFDR < 0.05):", length(dtu_genes), "\n")
cat("  - JAG1 targets tested:", sum(res_gene_final$is_JAG1_target), "\n")
cat("  - JAG1 targets with DTU (p < 0.05):", sum(res_gene_final$is_JAG1_target & res_gene_final$pvalue < 0.05, na.rm = TRUE), "\n\n")

cat("KEY FINDINGS:\n")
cat("  - JAG1 target enrichment fold:", round(enrichment_df$Fold_Enrichment[enrichment_df$Tier == "Gold"], 2), "(Gold tier)\n")
cat("  - Enrichment p-value:", format(enrichment_df$Pvalue[enrichment_df$Tier == "Gold"], digits = 3), "(Gold tier)\n\n")

cat("NEXT: Run 42_isoform_switch_analysis.R for detailed switch analysis\n")
