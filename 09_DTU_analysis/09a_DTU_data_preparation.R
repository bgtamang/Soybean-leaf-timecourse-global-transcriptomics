# Script 40: DTU (Differential Transcript Usage) Data Preparation

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 40: DTU DATA PREPARATION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/40_DTU", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/tables", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "tximport",         # For Salmon import
  "dplyr",            # For data manipulation
  "tidyr",            # For data reshaping
  "ggplot2",          # For plotting
  "SummarizedExperiment"  # For data structures
)

bioc_packages <- c("tximport", "SummarizedExperiment")

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

# ===== SECTION 1: LOAD SALMON OUTPUT =====

cat("========================================\n")
cat("SECTION 1: LOAD SALMON OUTPUT\n")
cat("========================================\n\n")

# Load the Salmon summarized output which contains transcript-level quantification
salmon_file <- "01_data/SalmonSummarizedOutput.RData"
if (!file.exists(salmon_file)) {
  stop("Salmon output file not found: ", salmon_file)
}

load(salmon_file)
cat("Loaded Salmon summarized output\n")
cat("  Objects loaded: ", paste(ls(), collapse = ", "), "\n\n")

# Check available objects and extract transcript-level data
if (exists("tx.all")) {
  cat("Found tx.all (transcript-level summarized experiment)\n")
  cat("  Class: ", class(tx.all), "\n")

  # Extract counts and TPM
  if (inherits(tx.all, "SummarizedExperiment")) {
    tx_counts <- assay(tx.all, "counts")
    tx_abundance <- assay(tx.all, "abundance")  # TPM
    cat("  Transcripts: ", nrow(tx_counts), "\n")
    cat("  Samples: ", ncol(tx_counts), "\n")
  } else if (is.list(tx.all)) {
    tx_counts <- tx.all$counts
    tx_abundance <- tx.all$abundance
    cat("  Transcripts: ", nrow(tx_counts), "\n")
    cat("  Samples: ", ncol(tx_counts), "\n")
  }
} else {
  # Try alternative object names
  stop("tx.all object not found in Salmon output. Available objects: ", paste(ls(), collapse = ", "))
}

# Sample transcript IDs
cat("\nSample transcript IDs:\n")
head_tx <- head(rownames(tx_counts), 10)
cat("  ", paste(head_tx, collapse = "\n  "), "\n\n")

# ===== SECTION 2: LOAD SAMPLE METADATA =====

cat("========================================\n")
cat("SECTION 2: LOAD SAMPLE METADATA\n")
cat("========================================\n\n")

# Load sample metadata from checkpoint
# Try multiple possible checkpoint names
metadata_files <- c(
  "03_results/checkpoints/02_sample_metadata.RData",
  "03_results/checkpoints/02_metadata.RData",
  "03_results/checkpoints/05_normalized.RData"
)

metadata_loaded <- FALSE
for (mf in metadata_files) {
  if (file.exists(mf)) {
    load(mf)
    metadata_loaded <- TRUE
    cat("Loaded metadata from:", mf, "\n")
    break
  }
}

if (!metadata_loaded) {
  stop("Could not find metadata checkpoint file")
}

cat("  Samples in metadata:", nrow(targets), "\n")

# Check sample name formats
cat("\nSample name formats:\n")
cat("  Salmon samples: ", head(colnames(tx_counts), 3), "\n")
cat("  Metadata samples: ", head(targets$Sample, 3), "\n")

# Create mapping between Salmon names and metadata names
# Salmon format: salmo462A_T0_R1 -> PI532462A_TP0_R1
# Mapping: 462A -> PI532462A, 713B -> PI612713B, 745 -> PI547745, LD11 -> LD11-2170

salmon_samples <- colnames(tx_counts)

# Parse Salmon sample names and create metadata
salmon_metadata <- data.frame(
  Salmon_Sample = salmon_samples,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Extract Line code
    Line_code = sub("salmo([^_]+)_.*", "\\1", Salmon_Sample),
    # Extract timepoint (T0 -> TP0)
    Timepoint = paste0("TP", sub(".*_T([0-4])_.*", "\\1", Salmon_Sample)),
    # Extract replicate
    Replicate = sub(".*_R([0-9])$", "\\1", Salmon_Sample),
    # Map Line codes to full names
    Line = case_when(
      Line_code == "462A" ~ "PI532462A",
      Line_code == "713B" ~ "PI612713B",
      Line_code == "745" ~ "PI547745",
      Line_code == "LD11" ~ "LD11-2170",
      TRUE ~ Line_code
    ),
    # Assign Leaf type based on Line
    Leaf_type = case_when(
      Line %in% c("PI532462A", "LD11-2170") ~ "Broad",
      Line %in% c("PI547745", "PI612713B") ~ "Narrow",
      TRUE ~ NA_character_
    ),
    # Create matching Sample name
    Sample = Salmon_Sample  # Keep Salmon name as Sample ID
  )

# Use the Salmon-derived metadata
targets <- salmon_metadata %>%
  select(Sample, Line, Leaf_type, Timepoint, Replicate)

cat("\nCreated metadata from Salmon sample names:\n")
cat("  Samples:", nrow(targets), "\n")
cat("\nSample distribution:\n")
print(table(targets$Leaf_type, targets$Timepoint))
cat("\n")

# ===== SECTION 3: CREATE TX2GENE MAPPING =====

cat("========================================\n")
cat("SECTION 3: CREATE TX2GENE MAPPING\n")
cat("========================================\n\n")

# Extract gene IDs from transcript IDs
# Soybean transcript IDs: Glyma.20G116200.1 -> Gene: Glyma.20G116200
transcript_ids <- rownames(tx_counts)

# Create tx2gene mapping
tx2gene <- data.frame(
  TXNAME = transcript_ids,
  GENEID = sub("\\.[0-9]+$", "", transcript_ids),  # Remove last .number
  stringsAsFactors = FALSE
)

cat("Created tx2gene mapping:\n")
cat("  Transcripts:", nrow(tx2gene), "\n")
cat("  Unique genes:", length(unique(tx2gene$GENEID)), "\n")

# Check transcript counts per gene
tx_per_gene <- tx2gene %>%
  group_by(GENEID) %>%
  summarise(n_tx = n(), .groups = "drop")

cat("\nTranscripts per gene distribution:\n")
print(summary(tx_per_gene$n_tx))
cat("  Genes with 1 transcript:", sum(tx_per_gene$n_tx == 1), "\n")
cat("  Genes with 2+ transcripts:", sum(tx_per_gene$n_tx >= 2), "\n")
cat("  Genes with 5+ transcripts:", sum(tx_per_gene$n_tx >= 5), "\n\n")

# ===== SECTION 4: FILTER TRANSCRIPTS =====

cat("========================================\n")
cat("SECTION 4: FILTER TRANSCRIPTS\n")
cat("========================================\n\n")

# Filtering criteria:
# 1. Minimum expression: >= 10 counts in at least 3 samples
# 2. Focus on genes with >= 2 expressed transcripts (needed for DTU)

min_counts <- 10
min_samples <- 3

# Step 1: Filter by expression
expressed <- rowSums(tx_counts >= min_counts) >= min_samples
cat("Transcripts with >= ", min_counts, " counts in >= ", min_samples, " samples:\n", sep = "")
cat("  Passed:", sum(expressed), "/", length(expressed), "\n")

tx_counts_filt <- tx_counts[expressed, ]
tx_abundance_filt <- tx_abundance[expressed, ]

# Update tx2gene
tx2gene_filt <- tx2gene[tx2gene$TXNAME %in% rownames(tx_counts_filt), ]

# Step 2: Keep only genes with >= 2 expressed transcripts
tx_per_gene_filt <- tx2gene_filt %>%
  group_by(GENEID) %>%
  summarise(n_tx = n(), .groups = "drop")

genes_multi_tx <- tx_per_gene_filt$GENEID[tx_per_gene_filt$n_tx >= 2]
cat("\nGenes with >= 2 expressed transcripts:", length(genes_multi_tx), "\n")

# Keep only transcripts from multi-transcript genes
keep_tx <- tx2gene_filt$TXNAME[tx2gene_filt$GENEID %in% genes_multi_tx]
tx_counts_filt <- tx_counts_filt[keep_tx, ]
tx_abundance_filt <- tx_abundance_filt[keep_tx, ]
tx2gene_filt <- tx2gene_filt[tx2gene_filt$TXNAME %in% keep_tx, ]

cat("\nAfter filtering:\n")
cat("  Transcripts:", nrow(tx_counts_filt), "\n")
cat("  Genes:", length(unique(tx2gene_filt$GENEID)), "\n\n")

# ===== SECTION 5: CALCULATE TRANSCRIPT PROPORTIONS =====

cat("========================================\n")
cat("SECTION 5: CALCULATE TRANSCRIPT PROPORTIONS\n")
cat("========================================\n\n")

# Calculate transcript proportion within each gene for each sample
# Proportion = transcript_count / sum(all_transcript_counts_for_gene)

calculate_proportions <- function(counts, tx2gene) {
  # Get gene totals per sample
  gene_counts <- aggregate(counts, by = list(Gene = tx2gene$GENEID[match(rownames(counts), tx2gene$TXNAME)]), FUN = sum)
  rownames(gene_counts) <- gene_counts$Gene
  gene_counts$Gene <- NULL
  gene_counts <- as.matrix(gene_counts)

  # Calculate proportions
  proportions <- counts
  for (tx in rownames(counts)) {
    gene <- tx2gene$GENEID[tx2gene$TXNAME == tx]
    gene_total <- gene_counts[gene, ]
    proportions[tx, ] <- counts[tx, ] / gene_total
  }

  # Replace NaN with 0 (when gene total is 0)
  proportions[is.nan(proportions)] <- 0

  return(proportions)
}

tx_proportions <- calculate_proportions(tx_counts_filt, tx2gene_filt)

cat("Calculated transcript proportions\n")
cat("  Dimension:", nrow(tx_proportions), "transcripts x", ncol(tx_proportions), "samples\n")
cat("  Range:", round(min(tx_proportions, na.rm = TRUE), 3), "to",
    round(max(tx_proportions, na.rm = TRUE), 3), "\n\n")

# Summary of proportion variation
prop_sd <- apply(tx_proportions, 1, sd, na.rm = TRUE)
cat("Transcript proportion variability (SD across samples):\n")
print(summary(prop_sd))
cat("  Highly variable (SD > 0.1):", sum(prop_sd > 0.1), "transcripts\n\n")

# ===== SECTION 6: CHECK JAG1 ISOFORMS =====

cat("========================================\n")
cat("SECTION 6: CHECK JAG1 ISOFORMS\n")
cat("========================================\n\n")

# JAG1 gene: Glyma.20G116200
jag1_gene <- "Glyma.20G116200"
jag1_tx <- tx2gene_filt$TXNAME[tx2gene_filt$GENEID == jag1_gene]

cat("JAG1 gene:", jag1_gene, "\n")
cat("JAG1 transcripts in filtered data:", length(jag1_tx), "\n")

if (length(jag1_tx) > 0) {
  cat("  Transcripts:", paste(jag1_tx, collapse = ", "), "\n\n")

  # Show JAG1 expression
  cat("JAG1 transcript counts (mean across samples):\n")
  jag1_counts <- tx_counts_filt[jag1_tx, , drop = FALSE]
  jag1_mean <- rowMeans(jag1_counts)
  for (i in seq_along(jag1_tx)) {
    cat("  ", jag1_tx[i], ": ", round(jag1_mean[i], 1), "\n", sep = "")
  }

  cat("\nJAG1 transcript proportions (mean across samples):\n")
  jag1_props <- tx_proportions[jag1_tx, , drop = FALSE]
  jag1_prop_mean <- rowMeans(jag1_props)
  for (i in seq_along(jag1_tx)) {
    cat("  ", jag1_tx[i], ": ", round(jag1_prop_mean[i] * 100, 1), "%\n", sep = "")
  }
} else {
  cat("  JAG1 not found in filtered transcript data\n")
  # Check if JAG1 is in unfiltered data
  jag1_tx_all <- tx2gene$TXNAME[tx2gene$GENEID == jag1_gene]
  if (length(jag1_tx_all) > 0) {
    cat("  (JAG1 transcripts exist but were filtered out due to low expression)\n")
    cat("  Original transcripts:", paste(jag1_tx_all, collapse = ", "), "\n")
  }
}

cat("\n")

# ===== SECTION 6B: DEDICATED JAG1 ISOFORM ANALYSIS (BYPASSING FILTERS) =====

cat("========================================\n")
cat("SECTION 6B: JAG1 ISOFORM ANALYSIS (TP0 FOCUS)\n")
cat("========================================\n\n")

# JAG1 is a developmental TF expressed mainly at TP0
# Analyze isoforms using unfiltered data

jag1_tx_all <- tx2gene$TXNAME[tx2gene$GENEID == jag1_gene]

if (length(jag1_tx_all) >= 2) {
  cat("Analyzing JAG1 isoforms from unfiltered data:\n")
  cat("  Transcripts:", paste(jag1_tx_all, collapse = ", "), "\n\n")

  # Get raw counts for JAG1
  jag1_counts_raw <- tx_counts[jag1_tx_all, , drop = FALSE]
  jag1_abundance_raw <- tx_abundance[jag1_tx_all, , drop = FALSE]

  # Expression summary across all samples
  cat("JAG1 expression (counts) across ALL samples:\n")
  cat("  Mean total gene expression:", round(sum(colMeans(jag1_counts_raw)), 1), "\n")
  for (tx in jag1_tx_all) {
    counts_vec <- jag1_counts_raw[tx, ]
    cat("  ", tx, ":\n", sep = "")
    cat("    Mean:", round(mean(counts_vec), 2), "\n")
    cat("    Samples with >0 counts:", sum(counts_vec > 0), "/60\n")
    cat("    Samples with >10 counts:", sum(counts_vec >= 10), "/60\n")
  }

  # Calculate proportions for JAG1 (from unfiltered data)
  jag1_gene_total <- colSums(jag1_counts_raw)
  jag1_proportions_raw <- sweep(jag1_counts_raw, 2, jag1_gene_total, "/")
  jag1_proportions_raw[is.nan(jag1_proportions_raw)] <- 0

  # Analysis by timepoint
  cat("\n\nJAG1 isoform proportions by TIMEPOINT:\n")
  cat("=" , rep("-", 60), "=\n", sep = "")

  for (tp in c("TP0", "TP1", "TP2", "TP3", "TP4")) {
    tp_samples <- targets$Sample[targets$Timepoint == tp]
    tp_counts <- jag1_counts_raw[, tp_samples, drop = FALSE]
    tp_props <- jag1_proportions_raw[, tp_samples, drop = FALSE]

    cat("\n", tp, ":\n", sep = "")
    cat("  Total JAG1 expression (mean counts):", round(mean(colSums(tp_counts)), 1), "\n")

    for (tx in jag1_tx_all) {
      mean_prop <- mean(tp_props[tx, ], na.rm = TRUE)
      sd_prop <- sd(tp_props[tx, ], na.rm = TRUE)
      cat("  ", tx, ": ", round(mean_prop * 100, 1), "% (SD: ", round(sd_prop * 100, 1), "%)\n", sep = "")
    }
  }

  # Analysis by leaf type at TP0 (where JAG1 is most expressed)
  cat("\n\nJAG1 isoform proportions at TP0 by LEAF TYPE:\n")
  cat("=" , rep("-", 60), "=\n", sep = "")

  tp0_samples <- targets$Sample[targets$Timepoint == "TP0"]

  for (lt in c("Broad", "Narrow")) {
    lt_samples <- targets$Sample[targets$Timepoint == "TP0" & targets$Leaf_type == lt]
    lt_counts <- jag1_counts_raw[, lt_samples, drop = FALSE]
    lt_props <- jag1_proportions_raw[, lt_samples, drop = FALSE]

    cat("\n", lt, " leaf (n=", length(lt_samples), "):\n", sep = "")
    cat("  Total JAG1 expression (mean counts):", round(mean(colSums(lt_counts)), 1), "\n")

    for (tx in jag1_tx_all) {
      mean_prop <- mean(lt_props[tx, ], na.rm = TRUE)
      sd_prop <- sd(lt_props[tx, ], na.rm = TRUE)
      mean_counts <- mean(lt_counts[tx, ])
      cat("  ", tx, ": ", round(mean_prop * 100, 1), "% (SD: ", round(sd_prop * 100, 1),
          "%) [mean counts: ", round(mean_counts, 1), "]\n", sep = "")
    }
  }

  # Statistical test for isoform difference at TP0
  cat("\n\nStatistical test for JAG1 isoform usage at TP0 (Narrow vs Broad):\n")
  cat("=" , rep("-", 60), "=\n", sep = "")

  broad_tp0 <- targets$Sample[targets$Timepoint == "TP0" & targets$Leaf_type == "Broad"]
  narrow_tp0 <- targets$Sample[targets$Timepoint == "TP0" & targets$Leaf_type == "Narrow"]

  # Use the dominant isoform for testing
  for (tx in jag1_tx_all) {
    broad_props <- jag1_proportions_raw[tx, broad_tp0]
    narrow_props <- jag1_proportions_raw[tx, narrow_tp0]

    # Wilcoxon test (non-parametric)
    if (sum(!is.na(broad_props)) >= 3 && sum(!is.na(narrow_props)) >= 3) {
      test_result <- wilcox.test(broad_props, narrow_props, exact = FALSE)

      cat("\n", tx, ":\n", sep = "")
      cat("  Broad mean: ", round(mean(broad_props, na.rm = TRUE) * 100, 1), "%\n", sep = "")
      cat("  Narrow mean: ", round(mean(narrow_props, na.rm = TRUE) * 100, 1), "%\n", sep = "")
      cat("  Difference: ", round((mean(narrow_props, na.rm = TRUE) - mean(broad_props, na.rm = TRUE)) * 100, 1),
          " percentage points\n", sep = "")
      cat("  Wilcoxon p-value: ", format(test_result$p.value, digits = 3), "\n", sep = "")
    }
  }

  # Save JAG1-specific data for later use
  jag1_isoform_analysis <- list(
    transcripts = jag1_tx_all,
    counts = jag1_counts_raw,
    abundance = jag1_abundance_raw,
    proportions = jag1_proportions_raw,
    metadata = targets
  )

  cat("\n")

  # ===== SECTION 6C: SAVE JAG1 ISOFORM DETAILED RESULTS =====

  cat("========================================\n")
  cat("SECTION 6C: SAVE JAG1 ISOFORM RESULTS\n")
  cat("========================================\n\n")

  # Create comprehensive JAG1 isoform summary by timepoint and leaf type
  jag1_results_list <- list()

  for (tp in c("TP0", "TP1", "TP2", "TP3", "TP4")) {
    for (lt in c("Broad", "Narrow")) {
      samples <- targets$Sample[targets$Timepoint == tp & targets$Leaf_type == lt]

      if (length(samples) > 0) {
        for (tx in jag1_tx_all) {
          counts_vec <- jag1_counts_raw[tx, samples]
          props_vec <- jag1_proportions_raw[tx, samples]

          jag1_results_list[[length(jag1_results_list) + 1]] <- data.frame(
            Transcript = tx,
            Timepoint = tp,
            Leaf_type = lt,
            N_samples = length(samples),
            Mean_counts = round(mean(counts_vec), 2),
            SD_counts = round(sd(counts_vec), 2),
            Min_counts = round(min(counts_vec), 2),
            Max_counts = round(max(counts_vec), 2),
            Samples_gt0 = sum(counts_vec > 0),
            Samples_gt10 = sum(counts_vec >= 10),
            Mean_proportion = round(mean(props_vec, na.rm = TRUE), 4),
            SD_proportion = round(sd(props_vec, na.rm = TRUE), 4),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  jag1_detailed_results <- do.call(rbind, jag1_results_list)

  # Add total gene expression per condition
  total_expr <- jag1_detailed_results %>%
    group_by(Timepoint, Leaf_type) %>%
    summarise(Total_gene_counts = sum(Mean_counts), .groups = "drop")

  jag1_detailed_results <- jag1_detailed_results %>%
    left_join(total_expr, by = c("Timepoint", "Leaf_type"))

  # Save detailed results
  write.csv(jag1_detailed_results, "03_results/tables/JAG1_isoform_detailed_results.csv", row.names = FALSE)
  cat("Saved: JAG1_isoform_detailed_results.csv\n")

  # Create summary table: proportion by timepoint (wide format for easy reading)
  jag1_summary_wide <- jag1_detailed_results %>%
    select(Transcript, Timepoint, Leaf_type, Mean_counts, Mean_proportion, Total_gene_counts) %>%
    mutate(Condition = paste0(Timepoint, "_", Leaf_type)) %>%
    select(-Timepoint, -Leaf_type) %>%
    pivot_wider(
      names_from = Condition,
      values_from = c(Mean_counts, Mean_proportion, Total_gene_counts),
      names_sep = "_"
    )

  write.csv(jag1_summary_wide, "03_results/tables/JAG1_isoform_summary_wide.csv", row.names = FALSE)
  cat("Saved: JAG1_isoform_summary_wide.csv\n")

  # Statistical tests at each timepoint
  cat("\n=== JAG1 Isoform Statistical Tests (Narrow vs Broad) ===\n\n")

  stat_results <- list()

  for (tp in c("TP0", "TP1", "TP2", "TP3", "TP4")) {
    broad_samples <- targets$Sample[targets$Timepoint == tp & targets$Leaf_type == "Broad"]
    narrow_samples <- targets$Sample[targets$Timepoint == tp & targets$Leaf_type == "Narrow"]

    # Total expression test
    broad_total <- colSums(jag1_counts_raw[, broad_samples, drop = FALSE])
    narrow_total <- colSums(jag1_counts_raw[, narrow_samples, drop = FALSE])

    if (length(broad_total) >= 2 && length(narrow_total) >= 2) {
      expr_test <- wilcox.test(broad_total, narrow_total, exact = FALSE)

      cat(tp, ":\n", sep = "")
      cat("  Total JAG1 expression:\n")
      cat("    Broad mean: ", round(mean(broad_total), 1), " counts\n", sep = "")
      cat("    Narrow mean: ", round(mean(narrow_total), 1), " counts\n", sep = "")
      cat("    Wilcoxon p-value: ", format(expr_test$p.value, digits = 3), "\n", sep = "")

      # Per-isoform proportion tests
      for (tx in jag1_tx_all) {
        broad_props <- jag1_proportions_raw[tx, broad_samples]
        narrow_props <- jag1_proportions_raw[tx, narrow_samples]

        # Remove NaN values
        broad_props <- broad_props[!is.nan(broad_props)]
        narrow_props <- narrow_props[!is.nan(narrow_props)]

        if (length(broad_props) >= 2 && length(narrow_props) >= 2) {
          prop_test <- wilcox.test(broad_props, narrow_props, exact = FALSE)

          stat_results[[length(stat_results) + 1]] <- data.frame(
            Timepoint = tp,
            Transcript = tx,
            Broad_mean_prop = round(mean(broad_props), 4),
            Narrow_mean_prop = round(mean(narrow_props), 4),
            Delta_prop = round(mean(narrow_props) - mean(broad_props), 4),
            Wilcox_pvalue = prop_test$p.value,
            stringsAsFactors = FALSE
          )

          cat("  ", tx, " proportion:\n", sep = "")
          cat("    Broad: ", round(mean(broad_props) * 100, 1), "%\n", sep = "")
          cat("    Narrow: ", round(mean(narrow_props) * 100, 1), "%\n", sep = "")
          cat("    Delta: ", round((mean(narrow_props) - mean(broad_props)) * 100, 1), " pp\n", sep = "")
          cat("    p-value: ", format(prop_test$p.value, digits = 3), "\n", sep = "")
        }
      }
      cat("\n")
    }
  }

  # Save statistical test results
  if (length(stat_results) > 0) {
    jag1_stat_tests <- do.call(rbind, stat_results)
    jag1_stat_tests$Significant <- jag1_stat_tests$Wilcox_pvalue < 0.05
    write.csv(jag1_stat_tests, "03_results/tables/JAG1_isoform_statistical_tests.csv", row.names = FALSE)
    cat("Saved: JAG1_isoform_statistical_tests.csv\n\n")
  }

} else {
  cat("JAG1 has fewer than 2 transcripts, skipping isoform analysis\n\n")
  jag1_isoform_analysis <- NULL
}

# ===== SECTION 7: LOAD JAG1 TARGETS =====

cat("========================================\n")
cat("SECTION 7: LOAD JAG1 TARGETS\n")
cat("========================================\n\n")

# Load JAG1 targets for enrichment analysis
jag1_targets_file <- "03_results/tables/JAG1_unified_targets_comprehensive.csv"
if (file.exists(jag1_targets_file)) {
  jag1_targets <- read.csv(jag1_targets_file, stringsAsFactors = FALSE)
  cat("Loaded JAG1 targets:", nrow(jag1_targets), "genes\n")

  # Check overlap with DTU genes
  dtu_genes <- unique(tx2gene_filt$GENEID)
  jag1_in_dtu <- jag1_targets$GeneID[jag1_targets$GeneID %in% dtu_genes]
  cat("JAG1 targets in DTU dataset:", length(jag1_in_dtu), "\n")

  # By tier
  cat("\nJAG1 targets by tier in DTU dataset:\n")
  for (tier in c("Gold", "Silver", "Bronze")) {
    tier_genes <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == tier]
    tier_in_dtu <- sum(tier_genes %in% dtu_genes)
    cat("  ", tier, ": ", tier_in_dtu, "/", length(tier_genes), "\n", sep = "")
  }
} else {
  cat("JAG1 targets file not found, skipping...\n")
  jag1_targets <- NULL
}

cat("\n")

# ===== SECTION 8: VISUALIZATION =====

cat("========================================\n")
cat("SECTION 8: VISUALIZATION\n")
cat("========================================\n\n")

# Plot 1: Transcripts per gene distribution
tx_per_gene_final <- tx2gene_filt %>%
  group_by(GENEID) %>%
  summarise(n_tx = n(), .groups = "drop")

p1 <- ggplot(tx_per_gene_final, aes(x = n_tx)) +
  geom_histogram(binwidth = 1, fill = "#2166AC", color = "white", alpha = 0.8) +
  scale_x_continuous(breaks = seq(2, max(tx_per_gene_final$n_tx), by = 2)) +
  labs(
    title = "Distribution of Expressed Transcripts per Gene",
    subtitle = paste("n =", nrow(tx_per_gene_final), "genes with 2+ transcripts"),
    x = "Number of Transcripts per Gene",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("03_results/figures/40_DTU/transcripts_per_gene.pdf", p1, width = 8, height = 6)
cat("Saved: transcripts_per_gene.pdf\n")

# Plot 2: Proportion variability distribution
prop_var_df <- data.frame(
  SD = prop_sd,
  Transcript = names(prop_sd)
)

p2 <- ggplot(prop_var_df, aes(x = SD)) +
  geom_histogram(binwidth = 0.02, fill = "#B2182B", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  annotate("text", x = 0.15, y = Inf, vjust = 2, hjust = 0,
           label = paste("Highly variable:", sum(prop_sd > 0.1)), color = "darkred") +
  labs(
    title = "Distribution of Transcript Proportion Variability",
    subtitle = "SD of transcript proportions across samples",
    x = "Standard Deviation of Proportion",
    y = "Number of Transcripts"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("03_results/figures/40_DTU/proportion_variability.pdf", p2, width = 8, height = 6)
cat("Saved: proportion_variability.pdf\n")

# Plot 3: JAG1 isoform proportions by condition (using unfiltered data)
if (!is.null(jag1_isoform_analysis) && length(jag1_isoform_analysis$transcripts) >= 2) {
  jag1_tx_plot <- jag1_isoform_analysis$transcripts
  jag1_props_plot <- jag1_isoform_analysis$proportions

  jag1_prop_df <- data.frame(
    Transcript = rep(jag1_tx_plot, ncol(jag1_props_plot)),
    Proportion = as.vector(t(jag1_props_plot)),
    Sample = rep(colnames(jag1_props_plot), each = length(jag1_tx_plot))
  )

  # Add metadata
  jag1_prop_df <- jag1_prop_df %>%
    left_join(targets[, c("Sample", "Leaf_type", "Timepoint", "Line")], by = "Sample")

  # Plot 3a: Across all timepoints
  p3a <- ggplot(jag1_prop_df, aes(x = Timepoint, y = Proportion, fill = Transcript)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21) +
    facet_wrap(~Leaf_type) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "JAG1 Isoform Proportions Across Development",
      subtitle = "Glyma.20G116200 transcript usage (from unfiltered data)",
      x = "Timepoint",
      y = "Transcript Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave("03_results/figures/40_DTU/JAG1_isoform_proportions_all.pdf", p3a, width = 10, height = 6)
  cat("Saved: JAG1_isoform_proportions_all.pdf\n")

  # Plot 3b: Focus on TP0 (where JAG1 is most expressed)
  jag1_tp0_df <- jag1_prop_df %>%
    filter(Timepoint == "TP0")

  p3b <- ggplot(jag1_tp0_df, aes(x = Leaf_type, y = Proportion, fill = Transcript)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21, width = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5, size = 2) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = "JAG1 Isoform Usage at TP0",
      subtitle = "Comparing Broad vs Narrow leaf genotypes",
      x = "Leaf Type",
      y = "Transcript Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave("03_results/figures/40_DTU/JAG1_isoform_TP0_comparison.pdf", p3b, width = 8, height = 6)
  cat("Saved: JAG1_isoform_TP0_comparison.pdf\n")

  # Plot 3c: JAG1 expression trajectory
  jag1_counts_plot <- jag1_isoform_analysis$counts
  jag1_expr_df <- data.frame(
    Sample = colnames(jag1_counts_plot),
    Total_Counts = colSums(jag1_counts_plot)
  ) %>%
    left_join(targets, by = "Sample")

  p3c <- ggplot(jag1_expr_df, aes(x = Timepoint, y = Total_Counts, color = Leaf_type, group = Leaf_type)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 4) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8) +
    scale_color_manual(values = c("Broad" = "#2E86AB", "Narrow" = "#E94F37")) +
    labs(
      title = "JAG1 Expression Across Development",
      subtitle = "Total counts (both isoforms combined)",
      x = "Timepoint",
      y = "Total Counts",
      color = "Leaf Type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave("03_results/figures/40_DTU/JAG1_expression_trajectory.pdf", p3c, width = 8, height = 6)
  cat("Saved: JAG1_expression_trajectory.pdf\n")
}

cat("\n")

# ===== SECTION 9: SAVE DATA =====

cat("========================================\n")
cat("SECTION 9: SAVE DATA\n")
cat("========================================\n\n")

# Save checkpoint for next script
save(
  tx_counts_filt,        # Filtered transcript counts
  tx_abundance_filt,     # Filtered transcript TPM
  tx_proportions,        # Transcript proportions
  tx2gene_filt,          # tx2gene mapping (filtered)
  targets,               # Sample metadata
  jag1_targets,          # JAG1 target genes (if loaded)
  jag1_isoform_analysis, # JAG1 isoform-specific analysis (bypasses filters)
  tx_counts,             # Unfiltered counts (for JAG1 analysis)
  tx_abundance,          # Unfiltered abundance
  tx2gene,               # Full tx2gene mapping
  file = "03_results/checkpoints/40_dtu_prepared.RData"
)
cat("Saved checkpoint: 40_dtu_prepared.RData\n")

# Save summary table
dtu_summary <- data.frame(
  Metric = c(
    "Total transcripts (raw)",
    "Total transcripts (filtered)",
    "Total genes (filtered)",
    "Min transcripts per gene",
    "Max transcripts per gene",
    "Median transcripts per gene",
    "Highly variable transcripts (SD > 0.1)",
    "JAG1 transcripts"
  ),
  Value = c(
    nrow(tx_counts),
    nrow(tx_counts_filt),
    length(unique(tx2gene_filt$GENEID)),
    min(tx_per_gene_final$n_tx),
    max(tx_per_gene_final$n_tx),
    median(tx_per_gene_final$n_tx),
    sum(prop_sd > 0.1),
    length(jag1_tx)
  )
)

write.csv(dtu_summary, "03_results/tables/dtu_data_summary.csv", row.names = FALSE)
cat("Saved: dtu_data_summary.csv\n\n")

# ===== COMPLETION =====

cat("================================================================\n")
cat("  SCRIPT 40 COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("SUMMARY:\n")
cat("  - Loaded transcript-level quantification from Salmon\n")
cat("  - Created tx2gene mapping\n")
cat("  - Filtered to", nrow(tx_counts_filt), "transcripts from",
    length(unique(tx2gene_filt$GENEID)), "genes\n")
cat("  - Calculated transcript proportions\n")
cat("  - Identified", sum(prop_sd > 0.1), "highly variable transcripts\n")
cat("  - JAG1 has", length(jag1_tx), "transcripts in dataset\n\n")

cat("NEXT: Run 41_DTU_DRIMSeq_analysis.R for statistical DTU testing\n")
