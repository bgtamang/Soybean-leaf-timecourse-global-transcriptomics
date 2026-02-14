# Script 31: Phenotype Integration

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 31: PHENOTYPE INTEGRATION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/phenotype", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/31_phenotype", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "minpack.lm",  # For nonlinear least squares (growth curves)
  "RColorBrewer"
)

# Check and load packages
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    cat("  Loaded:", pkg, "\n")
  } else {
    cat("  Missing:", pkg, "- install with install.packages('", pkg, "')\n")
  }
}

cat("\n")

# ===== CONFIGURATION =====

cat("========================================\n")
cat("SECTION 1: CONFIGURATION\n")
cat("========================================\n\n")

#--- PHENOTYPE DATA FILE ---#
# Your actual phenotype data file
phenotype_file <- "01_data/soybean_leaf_data_long_format.csv"

# Expected format (long format):
# Leaf,Rep,Day,Leaf_Type,Measurement,value
# Broad,1,1,VC1,Length,2.4
# Broad,1,1,VC1,Width,2
# Broad,1,1,VC1,Ratio,1.2
# ...

#--- GENOTYPE INFORMATION ---#
# RNA-seq lines (for mapping)
broad_lines <- c("PI532462A", "LD112170")
narrow_lines <- c("PI612713B", "PI547745")
all_lines <- c(broad_lines, narrow_lines)

#--- LEAF TYPE INFORMATION ---#
# VC = cotyledon leaves (no phenotypic difference expected)
# V1, V2 = true leaves where Broad vs Narrow differences occur
cotyledon_leaves <- c("VC1", "VC2")
true_leaves <- c("V1", "V2")

cat("Configuration set:\n")
cat("  Phenotype file:", phenotype_file, "\n")
cat("  Cotyledon leaves (no difference):", paste(cotyledon_leaves, collapse = ", "), "\n")
cat("  True leaves (phenotype differs):", paste(true_leaves, collapse = ", "), "\n\n")

# ===== LOAD RNA-SEQ DATA =====

cat("========================================\n")
cat("SECTION 2: LOAD RNA-SEQ DATA\n")
cat("========================================\n\n")

# Load expression data
if (file.exists("03_results/checkpoints/06_validated.RData")) {
  load("03_results/checkpoints/06_validated.RData")
  cat("Loaded: validated expression data\n")
  cat("  Samples:", ncol(v_primary$E), "\n")
  cat("  Genes:", nrow(v_primary$E), "\n")
} else {
  stop("Required checkpoint 06_validated.RData not found!")
}

# Load JAG1 targets
if (file.exists("03_results/checkpoints/14_JAG1_targets.RData")) {
  load("03_results/checkpoints/14_JAG1_targets.RData")
  jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
  cat("Loaded: JAG1 targets (", nrow(jag1_targets), " genes)\n")
} else {
  cat("Warning: JAG1 targets not found, continuing without\n")
  jag1_targets <- NULL
}

# Load WGCNA if available
wgcna_available <- FALSE
if (file.exists("03_results/checkpoints/21_WGCNA_JAG1.RData")) {
  load("03_results/checkpoints/21_WGCNA_JAG1.RData")
  wgcna_available <- TRUE
  cat("Loaded: WGCNA results\n")
}

cat("\n")

# ===== LOAD PHENOTYPE DATA =====

cat("========================================\n")
cat("SECTION 3: LOAD PHENOTYPE DATA\n")
cat("========================================\n\n")

if (file.exists(phenotype_file)) {
  # Load actual phenotype data (long format)
  pheno_long <- read.csv(phenotype_file, stringsAsFactors = FALSE)
  cat("Loaded phenotype data:", nrow(pheno_long), "observations\n")
  cat("  Columns:", paste(colnames(pheno_long), collapse = ", "), "\n")
  cat("  Leaf types:", paste(unique(pheno_long$Leaf), collapse = ", "), "\n")
  cat("  Leaf positions:", paste(unique(pheno_long$Leaf_Type), collapse = ", "), "\n")
  cat("  Day range:", min(pheno_long$Day), "-", max(pheno_long$Day), "\n")
  cat("  Replicates:", length(unique(pheno_long$Rep)), "\n")
  cat("  Measurements:", paste(unique(pheno_long$Measurement), collapse = ", "), "\n\n")

  # Convert from long to wide format
  pheno_wide <- reshape(pheno_long,
                        idvar = c("Leaf", "Rep", "Day", "Leaf_Type"),
                        timevar = "Measurement",
                        direction = "wide")
  colnames(pheno_wide) <- gsub("value.", "", colnames(pheno_wide))

  cat("Converted to wide format:", nrow(pheno_wide), "observations\n")

  # Preview data
  cat("\nPhenotype data preview (wide format):\n")
  print(head(pheno_wide, 10))

  phenotype_available <- TRUE

} else {
  cat("ERROR: Phenotype file not found:", phenotype_file, "\n")
  cat("Please ensure the file exists at the specified path.\n\n")
  phenotype_available <- FALSE
  stop("Phenotype file required. Stopping.")
}

# ===== CALCULATE DERIVED TRAITS =====

cat("\n========================================\n")
cat("SECTION 4: ANALYZE PHENOTYPE DATA\n")
cat("========================================\n\n")

# Work with the wide format data
pheno_data <- pheno_wide

# Ensure Ratio is calculated if not present
if (!"Ratio" %in% colnames(pheno_data) && "Length" %in% colnames(pheno_data) && "Width" %in% colnames(pheno_data)) {
  pheno_data$Ratio <- pheno_data$Length / pheno_data$Width
}

# Focus on V1 and V2 leaves where phenotypic differences occur
pheno_v1v2 <- pheno_data[pheno_data$Leaf_Type %in% true_leaves, ]

cat("V1/V2 leaf measurements:", nrow(pheno_v1v2), "observations\n")
cat("  V1 observations:", sum(pheno_v1v2$Leaf_Type == "V1"), "\n")
cat("  V2 observations:", sum(pheno_v1v2$Leaf_Type == "V2"), "\n\n")

# Calculate summary statistics for each Leaf type (Broad vs Narrow) and leaf position (V1, V2)
pheno_summary <- pheno_v1v2 %>%
  filter(!is.na(Ratio)) %>%
  group_by(Leaf, Leaf_Type) %>%
  summarize(
    N_observations = n(),
    N_reps = length(unique(Rep)),
    Days_measured = length(unique(Day)),
    Mean_Length = mean(Length, na.rm = TRUE),
    SD_Length = sd(Length, na.rm = TRUE),
    Mean_Width = mean(Width, na.rm = TRUE),
    SD_Width = sd(Width, na.rm = TRUE),
    Mean_Ratio = mean(Ratio, na.rm = TRUE),
    SD_Ratio = sd(Ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("Phenotype summary (V1/V2 leaves):\n")
print(as.data.frame(pheno_summary))
cat("\n")

# Calculate mature leaf measurements (using later days when leaves are fully developed)
# Get the last measurement for each rep (mature leaf)
mature_leaves <- pheno_v1v2 %>%
  filter(!is.na(Ratio)) %>%
  group_by(Leaf, Rep, Leaf_Type) %>%
  filter(Day == max(Day)) %>%
  ungroup()

mature_summary <- mature_leaves %>%
  group_by(Leaf, Leaf_Type) %>%
  summarize(
    N_plants = n(),
    Final_Length = mean(Length, na.rm = TRUE),
    Final_Width = mean(Width, na.rm = TRUE),
    Final_Ratio = mean(Ratio, na.rm = TRUE),
    Ratio_SD = sd(Ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("Mature leaf summary:\n")
print(as.data.frame(mature_summary))
cat("\n")

# KEY COMPARISON: Broad vs Narrow at maturity
cat("========================================\n")
cat("KEY PHENOTYPE COMPARISON\n")
cat("========================================\n\n")

for (lt in c("V1", "V2")) {
  broad_ratio <- mature_summary$Final_Ratio[mature_summary$Leaf == "Broad" & mature_summary$Leaf_Type == lt]
  narrow_ratio <- mature_summary$Final_Ratio[mature_summary$Leaf == "Narrow" & mature_summary$Leaf_Type == lt]

  if (length(broad_ratio) > 0 && length(narrow_ratio) > 0) {
    cat(lt, "Leaf L:W Ratio:\n")
    cat("  Broad:", round(broad_ratio, 3), "\n")
    cat("  Narrow:", round(narrow_ratio, 3), "\n")
    cat("  Difference:", round(narrow_ratio - broad_ratio, 3), "\n")
    cat("  Fold change:", round(narrow_ratio / broad_ratio, 3), "\n\n")
  }
}

# Save summaries
write.csv(pheno_summary, "03_results/tables/phenotype/phenotype_summary.csv", row.names = FALSE)
write.csv(mature_summary, "03_results/tables/phenotype/mature_leaf_summary.csv", row.names = FALSE)

# ===== GROWTH TRAJECTORY ANALYSIS =====

cat("\n========================================\n")
cat("SECTION 5: GROWTH TRAJECTORY ANALYSIS\n")
cat("========================================\n\n")

# Analyze how L:W ratio changes over time for V1 leaves
v1_trajectory <- pheno_v1v2 %>%
  filter(Leaf_Type == "V1", !is.na(Ratio)) %>%
  group_by(Leaf, Day) %>%
  summarize(
    Mean_Length = mean(Length, na.rm = TRUE),
    Mean_Width = mean(Width, na.rm = TRUE),
    Mean_Ratio = mean(Ratio, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

cat("V1 Leaf growth trajectory (by day):\n")
cat("\nBroad leaf V1:\n")
broad_v1 <- v1_trajectory[v1_trajectory$Leaf == "Broad", ]
print(head(broad_v1[order(broad_v1$Day), ], 10))

cat("\nNarrow leaf V1:\n")
narrow_v1 <- v1_trajectory[v1_trajectory$Leaf == "Narrow", ]
print(head(narrow_v1[order(narrow_v1$Day), ], 10))

# Calculate when ratio divergence becomes significant
cat("\n========================================\n")
cat("L:W RATIO DIVERGENCE ANALYSIS\n")
cat("========================================\n\n")

# For each day, test if Broad vs Narrow ratios are different
ratio_comparison <- pheno_v1v2 %>%
  filter(Leaf_Type == "V1", !is.na(Ratio)) %>%
  group_by(Day) %>%
  summarize(
    Broad_Mean = mean(Ratio[Leaf == "Broad"], na.rm = TRUE),
    Narrow_Mean = mean(Ratio[Leaf == "Narrow"], na.rm = TRUE),
    Difference = Narrow_Mean - Broad_Mean,
    N_Broad = sum(Leaf == "Broad"),
    N_Narrow = sum(Leaf == "Narrow"),
    .groups = "drop"
  ) %>%
  filter(N_Broad > 0, N_Narrow > 0)

# Run t-test for each day with sufficient samples
ratio_comparison$P_value <- NA
for (i in 1:nrow(ratio_comparison)) {
  day_data <- pheno_v1v2[pheno_v1v2$Leaf_Type == "V1" &
                          pheno_v1v2$Day == ratio_comparison$Day[i] &
                          !is.na(pheno_v1v2$Ratio), ]
  if (sum(day_data$Leaf == "Broad") >= 3 && sum(day_data$Leaf == "Narrow") >= 3) {
    tt <- t.test(Ratio ~ Leaf, data = day_data)
    ratio_comparison$P_value[i] <- tt$p.value
  }
}

ratio_comparison$Significant <- ratio_comparison$P_value < 0.05

cat("Ratio comparison by day (V1 leaves):\n")
print(ratio_comparison)

# Find day when difference becomes significant
first_sig_day <- min(ratio_comparison$Day[ratio_comparison$Significant], na.rm = TRUE)
cat("\nFirst day with significant B vs N difference:", first_sig_day, "\n")

write.csv(ratio_comparison, "03_results/tables/phenotype/ratio_divergence_by_day.csv", row.names = FALSE)

# ===== CREATE PHENOTYPE TRAITS FOR RNA-SEQ SAMPLES =====

cat("\n========================================\n")
cat("SECTION 6: CREATE PHENOTYPE TRAITS\n")
cat("========================================\n\n")

# Since phenotype data is collected by Broad/Narrow (not by specific genotype),
# we assign phenotype values to RNA-seq samples based on their leaf type

# Get mature V1 ratio values for Broad and Narrow
v1_mature <- mature_summary[mature_summary$Leaf_Type == "V1", ]
broad_v1_ratio <- v1_mature$Final_Ratio[v1_mature$Leaf == "Broad"]
narrow_v1_ratio <- v1_mature$Final_Ratio[v1_mature$Leaf == "Narrow"]

cat("Mature V1 L:W Ratios:\n")
cat("  Broad:", round(broad_v1_ratio, 3), "\n")
cat("  Narrow:", round(narrow_v1_ratio, 3), "\n\n")

# Get V2 values if available
v2_mature <- mature_summary[mature_summary$Leaf_Type == "V2", ]
if (nrow(v2_mature) > 0) {
  broad_v2_ratio <- v2_mature$Final_Ratio[v2_mature$Leaf == "Broad"]
  narrow_v2_ratio <- v2_mature$Final_Ratio[v2_mature$Leaf == "Narrow"]
  cat("Mature V2 L:W Ratios:\n")
  cat("  Broad:", round(broad_v2_ratio, 3), "\n")
  cat("  Narrow:", round(narrow_v2_ratio, 3), "\n\n")
} else {
  broad_v2_ratio <- NA
  narrow_v2_ratio <- NA
}

# Store phenotype characteristics
pheno_traits <- data.frame(
  Leaf_Type = c("Broad", "Narrow"),
  V1_Ratio = c(broad_v1_ratio, narrow_v1_ratio),
  V2_Ratio = c(broad_v2_ratio, narrow_v2_ratio),
  V1_Length = c(v1_mature$Final_Length[v1_mature$Leaf == "Broad"],
                v1_mature$Final_Length[v1_mature$Leaf == "Narrow"]),
  V1_Width = c(v1_mature$Final_Width[v1_mature$Leaf == "Broad"],
               v1_mature$Final_Width[v1_mature$Leaf == "Narrow"])
)

cat("Phenotype trait values to assign to RNA-seq samples:\n")
print(pheno_traits)

write.csv(pheno_traits, "03_results/tables/phenotype/phenotype_traits.csv", row.names = FALSE)

# ===== CREATE TRAIT MATRIX FOR RNA-SEQ SAMPLES =====

cat("\n========================================\n")
cat("SECTION 7: CREATE TRAIT MATRIX\n")
cat("========================================\n\n")

# Create trait matrix aligned with RNA-seq samples
# Each row = sample, columns = phenotypic traits

if (exists("targets_primary")) {

  # Get sample names that match expression matrix column names
  # Expression matrix uses rownames of targets_primary (or Sample_name column)
  expr_samples <- colnames(v_primary$E)

  # Determine which identifier matches expression matrix
  if (all(expr_samples %in% rownames(targets_primary))) {
    sample_ids <- rownames(targets_primary)
    cat("Using rownames of targets_primary for sample matching\n")
  } else if ("Sample_name" %in% colnames(targets_primary) &&
             all(expr_samples %in% targets_primary$Sample_name)) {
    sample_ids <- targets_primary$Sample_name
    cat("Using Sample_name column for sample matching\n")
  } else if (all(expr_samples %in% targets_primary$Sample)) {
    sample_ids <- targets_primary$Sample
    cat("Using Sample column for sample matching\n")
  } else {
    # Find the column that best matches
    cat("WARNING: Sample name mismatch detected!\n")
    cat("Expression matrix samples:", head(expr_samples, 3), "...\n")
    cat("targets_primary rownames:", head(rownames(targets_primary), 3), "...\n")
    if ("Sample_name" %in% colnames(targets_primary)) {
      cat("targets_primary$Sample_name:", head(targets_primary$Sample_name, 3), "...\n")
    }
    cat("targets_primary$Sample:", head(targets_primary$Sample, 3), "...\n")
    sample_ids <- rownames(targets_primary)  # Default fallback
  }

  # Initialize trait matrix with proper sample identifiers
  trait_matrix <- data.frame(
    Sample = sample_ids,
    Line = targets_primary$Line,
    Leaf_type = as.character(targets_primary$Leaf_type),  # Ensure character, not factor
    Timepoint = as.character(targets_primary$Timepoint),
    stringsAsFactors = FALSE
  )

  cat("Leaf_type values in trait_matrix:", paste(unique(trait_matrix$Leaf_type), collapse = ", "), "\n")
  cat("Broad samples:", sum(trait_matrix$Leaf_type == "Broad"), "\n")
  cat("Narrow samples:", sum(trait_matrix$Leaf_type == "Narrow"), "\n\n")

  # ===== LINE-SPECIFIC L:W RATIOS (ACTUAL MEASURED PHENOTYPES) =====
  # These are the actual measured leaf L:W ratios for each soybean line
  # Using line-specific values gives 4 distinct phenotype values instead of 2
  line_lw_ratios <- c(
    "PI612713B" = 2.36,  # Narrow line 1
    "PI547745" = 2.16,   # Narrow line 2
    "LD112170" = 1.59,   # Broad line 1
    "PI532462A" = 1.55   # Broad line 2
  )

  cat("Using LINE-SPECIFIC L:W ratios (4 values, not 2):\n")
  cat("  PI612713B (Narrow): 2.36\n")
  cat("  PI547745 (Narrow): 2.16\n")
  cat("  LD112170 (Broad): 1.59\n")
  cat("  PI532462A (Broad): 1.55\n\n")

  # Assign L:W ratio by LINE (not by Leaf_type)
  trait_matrix$Line_LW_Ratio <- line_lw_ratios[as.character(trait_matrix$Line)]

  # Verify assignment
  cat("Line-specific L:W ratio assignment:\n")
  for (line in names(line_lw_ratios)) {
    n_samples <- sum(trait_matrix$Line == line)
    assigned_val <- unique(trait_matrix$Line_LW_Ratio[trait_matrix$Line == line])
    cat("  ", line, ":", n_samples, "samples -> L:W =", assigned_val, "\n")
  }
  cat("\n")

  # Also keep the old 2-value approach for comparison (renamed)
  trait_matrix$V1_Ratio_Binary <- ifelse(trait_matrix$Leaf_type == "Broad",
                                          broad_v1_ratio,
                                          narrow_v1_ratio)

  # Set the primary phenotype column to use line-specific values
  trait_matrix$V1_Ratio <- trait_matrix$Line_LW_Ratio

  trait_matrix$V1_Length <- ifelse(trait_matrix$Leaf_type == "Broad",
                                    pheno_traits$V1_Length[pheno_traits$Leaf_Type == "Broad"],
                                    pheno_traits$V1_Length[pheno_traits$Leaf_Type == "Narrow"])

  trait_matrix$V1_Width <- ifelse(trait_matrix$Leaf_type == "Broad",
                                   pheno_traits$V1_Width[pheno_traits$Leaf_Type == "Broad"],
                                   pheno_traits$V1_Width[pheno_traits$Leaf_Type == "Narrow"])

  # Assign V2 ratio if available
  if (!is.na(broad_v2_ratio) && !is.na(narrow_v2_ratio)) {
    trait_matrix$V2_Ratio <- ifelse(trait_matrix$Leaf_type == "Broad",
                                     broad_v2_ratio,
                                     narrow_v2_ratio)
  }

  # Create timepoint-specific phenotype traits
  # V1 phenotype is relevant at TP1+ (when V1 starts developing)
  trait_matrix$V1_Ratio_TP1plus <- trait_matrix$V1_Ratio *
    as.numeric(trait_matrix$Timepoint %in% c("TP1", "TP2", "TP3", "TP4"))

  # V1 phenotype at mature stages (TP2+)
  trait_matrix$V1_Ratio_TP2plus <- trait_matrix$V1_Ratio *
    as.numeric(trait_matrix$Timepoint %in% c("TP2", "TP3", "TP4"))

  # Binary leaf type for correlation
  trait_matrix$Is_Narrow <- as.numeric(trait_matrix$Leaf_type == "Narrow")

  # Reorder to match expression data (use the identified sample column)
  match_idx <- match(colnames(v_primary$E), trait_matrix$Sample)
  if (any(is.na(match_idx))) {
    cat("WARNING: Some samples did not match! Non-matching:", sum(is.na(match_idx)), "\n")
  }
  trait_matrix <- trait_matrix[match_idx, ]

  cat("Trait matrix created:\n")
  cat("  Samples:", nrow(trait_matrix), "\n")
  cat("  Trait columns:\n")
  trait_cols <- c("V1_Ratio", "V1_Length", "V1_Width", "V2_Ratio",
                  "V1_Ratio_TP1plus", "V1_Ratio_TP2plus", "Is_Narrow")
  for (col in trait_cols) {
    if (col %in% colnames(trait_matrix)) {
      vals <- trait_matrix[[col]]
      cat("    ", col, ": range =", round(min(vals, na.rm = TRUE), 3), "-",
          round(max(vals, na.rm = TRUE), 3), "\n")
    }
  }

  write.csv(trait_matrix, "03_results/tables/phenotype/trait_matrix.csv", row.names = FALSE)

} else {
  cat("targets_primary not found. Skipping trait matrix creation.\n")
  cat("Load validated expression data (checkpoint 06) first.\n")
}

# ===== GENE-PHENOTYPE CORRELATIONS =====

cat("\n========================================\n")
cat("SECTION 8: GENE-PHENOTYPE CORRELATIONS\n")
cat("========================================\n\n")

# Calculate correlation between each gene and phenotypic traits
# Focus on V1 L:W ratio as the key trait

if (exists("trait_matrix") && exists("v_primary")) {

  # Use V1_Ratio as the primary phenotype trait
  valid_samples <- !is.na(trait_matrix$V1_Ratio)
  cat("Samples with V1 phenotype data:", sum(valid_samples), "of", length(valid_samples), "\n")

  if (sum(valid_samples) >= 10) {

    expr_data <- v_primary$E[, valid_samples]
    v1_ratio <- trait_matrix$V1_Ratio[valid_samples]
    is_narrow <- trait_matrix$Is_Narrow[valid_samples]

    # Calculate correlations for all genes
    cat("Calculating gene-phenotype correlations...\n")

    gene_pheno_cor <- data.frame(
      GeneID = rownames(expr_data),
      stringsAsFactors = FALSE
    )

    # Correlation with V1 L:W ratio
    gene_pheno_cor$Cor_LW_Ratio <- apply(expr_data, 1, function(x) {
      cor(x, v1_ratio, method = "spearman", use = "complete.obs")
    })

    # P-value for correlation
    gene_pheno_cor$Pval_LW_Ratio <- apply(expr_data, 1, function(x) {
      tryCatch({
        cor.test(x, v1_ratio, method = "spearman")$p.value
      }, error = function(e) NA)
    })

    # FDR correction
    gene_pheno_cor$FDR_LW_Ratio <- p.adjust(gene_pheno_cor$Pval_LW_Ratio, method = "BH")

    # Correlation with V1 width (negative = associated with narrow)
    gene_pheno_cor$Cor_Width <- apply(expr_data, 1, function(x) {
      cor(x, trait_matrix$V1_Width[valid_samples], method = "spearman", use = "complete.obs")
    })

    # Classify genes based on phenotype correlation
    gene_pheno_cor$Pheno_Association <- "None"
    gene_pheno_cor$Pheno_Association[gene_pheno_cor$Cor_LW_Ratio > 0.3 & gene_pheno_cor$FDR_LW_Ratio < 0.05] <- "High_LW_Ratio"
    gene_pheno_cor$Pheno_Association[gene_pheno_cor$Cor_LW_Ratio < -0.3 & gene_pheno_cor$FDR_LW_Ratio < 0.05] <- "Low_LW_Ratio"

    # Summary
    cat("\nGene-phenotype correlation summary:\n")
    cat("  Positive correlation with V1 L:W ratio (r > 0.3, FDR < 0.05):",
        sum(gene_pheno_cor$Pheno_Association == "High_LW_Ratio"), "\n")
    cat("  Negative correlation with V1 L:W ratio (r < -0.3, FDR < 0.05):",
        sum(gene_pheno_cor$Pheno_Association == "Low_LW_Ratio"), "\n")

    # Add JAG1 target status
    if (!is.null(jag1_targets)) {
      gene_pheno_cor$Is_JAG1_Target <- gene_pheno_cor$GeneID %in% jag1_targets$GeneID
      gene_pheno_cor$Target_Tier <- jag1_targets$Confidence_Tier[
        match(gene_pheno_cor$GeneID, jag1_targets$GeneID)
      ]

      # How many JAG1 targets are correlated with phenotype?
      target_pheno <- gene_pheno_cor[gene_pheno_cor$Is_JAG1_Target, ]
      cat("\nJAG1 targets with phenotype correlation:\n")
      cat("  High L:W ratio (narrow phenotype):", sum(target_pheno$Pheno_Association == "High_LW_Ratio", na.rm = TRUE), "\n")
      cat("  Low L:W ratio (broad phenotype):", sum(target_pheno$Pheno_Association == "Low_LW_Ratio", na.rm = TRUE), "\n")

      # Summary by tier
      cat("\nJAG1 target phenotype association by tier:\n")
      tier_pheno <- table(target_pheno$Target_Tier, target_pheno$Pheno_Association)
      print(tier_pheno)
    }

    # Save results
    write.csv(gene_pheno_cor, "03_results/tables/phenotype/gene_phenotype_correlations.csv",
              row.names = FALSE)

    # Top genes correlated with L:W ratio
    top_lw_genes <- gene_pheno_cor %>%
      filter(FDR_LW_Ratio < 0.05) %>%
      arrange(desc(abs(Cor_LW_Ratio))) %>%
      head(100)

    write.csv(top_lw_genes, "03_results/tables/phenotype/top_LW_ratio_genes.csv",
              row.names = FALSE)
    cat("\nSaved: gene_phenotype_correlations.csv\n")
    cat("Saved: top_LW_ratio_genes.csv (top 100)\n")
  }
} else {
  cat("Expression data or trait matrix not available.\n")
  cat("Skipping gene-phenotype correlation analysis.\n")
  gene_pheno_cor <- NULL
}

# ===== WGCNA MODULE-PHENOTYPE CORRELATION =====

cat("\n========================================\n")
cat("SECTION 9: MODULE-PHENOTYPE CORRELATION\n")
cat("========================================\n\n")

if (wgcna_available && exists("MEs") && exists("trait_matrix")) {

  # Prepare numeric trait matrix for WGCNA correlation
  valid_samples <- !is.na(trait_matrix$V1_Ratio)

  # Select numeric phenotype columns
  pheno_cols <- c("V1_Ratio", "V1_Length", "V1_Width", "Is_Narrow")
  if ("V2_Ratio" %in% colnames(trait_matrix)) {
    pheno_cols <- c(pheno_cols, "V2_Ratio")
  }

  numeric_traits <- trait_matrix[valid_samples, pheno_cols, drop = FALSE]

  # Subset MEs to match samples
  ME_samples <- rownames(MEs)
  trait_samples <- trait_matrix$Sample[valid_samples]
  common_samples <- intersect(ME_samples, trait_samples)

  if (length(common_samples) >= 10) {
    ME_subset <- MEs[common_samples, ]
    numeric_traits_subset <- numeric_traits[match(common_samples, trait_samples), ]

    # Calculate module-trait correlations
    moduleTraitCor_pheno <- cor(ME_subset, numeric_traits_subset, use = "pairwise.complete.obs")
    moduleTraitPval_pheno <- matrix(NA, nrow = ncol(ME_subset), ncol = ncol(numeric_traits_subset))

    for (i in 1:ncol(ME_subset)) {
      for (j in 1:ncol(numeric_traits_subset)) {
        test <- cor.test(ME_subset[, i], numeric_traits_subset[, j])
        moduleTraitPval_pheno[i, j] <- test$p.value
      }
    }

    rownames(moduleTraitCor_pheno) <- colnames(ME_subset)
    rownames(moduleTraitPval_pheno) <- colnames(ME_subset)
    colnames(moduleTraitPval_pheno) <- colnames(numeric_traits_subset)

    cat("Module-phenotype correlations (V1_Ratio):\n")
    v1_cor <- moduleTraitCor_pheno[, "V1_Ratio"]
    v1_cor_sorted <- sort(abs(v1_cor), decreasing = TRUE)
    print(head(round(v1_cor[names(v1_cor_sorted)], 3), 10))

    # Find modules most correlated with V1 L:W ratio
    top_modules <- names(v1_cor_sorted)[1:5]

    cat("\nTop 5 modules correlated with V1 L:W ratio:\n")
    for (mod in top_modules) {
      pval <- moduleTraitPval_pheno[mod, "V1_Ratio"]
      cat("  ", mod, ": r =", round(v1_cor[mod], 3), ", p =", format(pval, digits = 3), "\n")
    }

    # Save
    write.csv(moduleTraitCor_pheno, "03_results/tables/phenotype/module_phenotype_correlations.csv")
    cat("\nSaved: module_phenotype_correlations.csv\n")

  } else {
    cat("Not enough common samples between WGCNA and trait matrix.\n")
  }

} else {
  cat("WGCNA results not available. Skipping module-phenotype analysis.\n")
  cat("Run WGCNA scripts (18-21) first, then re-run this script.\n")
}

# ===== VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 10: VISUALIZATIONS\n")
cat("========================================\n\n")

# Plot 1: V1 Leaf growth trajectories (Broad vs Narrow)
png("03_results/figures/31_phenotype/v1_growth_trajectories.png",
    width = 12, height = 8, units = "in", res = 150)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# V1 Length over time
plot(NULL, xlim = c(min(v1_trajectory$Day), max(v1_trajectory$Day)),
     ylim = c(0, max(v1_trajectory$Mean_Length, na.rm = TRUE) * 1.1),
     xlab = "Day", ylab = "V1 Leaf Length (cm)", main = "A. V1 Leaf Length Growth")

lines(broad_v1$Day, broad_v1$Mean_Length, col = "#E69F00", lwd = 2)
points(broad_v1$Day, broad_v1$Mean_Length, col = "#E69F00", pch = 16)
lines(narrow_v1$Day, narrow_v1$Mean_Length, col = "#56B4E9", lwd = 2)
points(narrow_v1$Day, narrow_v1$Mean_Length, col = "#56B4E9", pch = 17)
legend("bottomright", legend = c("Broad", "Narrow"),
       col = c("#E69F00", "#56B4E9"), pch = c(16, 17), lwd = 2)

# V1 Width over time
plot(NULL, xlim = c(min(v1_trajectory$Day), max(v1_trajectory$Day)),
     ylim = c(0, max(v1_trajectory$Mean_Width, na.rm = TRUE) * 1.1),
     xlab = "Day", ylab = "V1 Leaf Width (cm)", main = "B. V1 Leaf Width Growth")

lines(broad_v1$Day, broad_v1$Mean_Width, col = "#E69F00", lwd = 2)
points(broad_v1$Day, broad_v1$Mean_Width, col = "#E69F00", pch = 16)
lines(narrow_v1$Day, narrow_v1$Mean_Width, col = "#56B4E9", lwd = 2)
points(narrow_v1$Day, narrow_v1$Mean_Width, col = "#56B4E9", pch = 17)

# V1 L:W Ratio over time
plot(NULL, xlim = c(min(v1_trajectory$Day), max(v1_trajectory$Day)),
     ylim = c(min(v1_trajectory$Mean_Ratio, na.rm = TRUE) * 0.9,
              max(v1_trajectory$Mean_Ratio, na.rm = TRUE) * 1.1),
     xlab = "Day", ylab = "V1 L:W Ratio", main = "C. V1 L:W Ratio Trajectory")

lines(broad_v1$Day, broad_v1$Mean_Ratio, col = "#E69F00", lwd = 2)
points(broad_v1$Day, broad_v1$Mean_Ratio, col = "#E69F00", pch = 16)
lines(narrow_v1$Day, narrow_v1$Mean_Ratio, col = "#56B4E9", lwd = 2)
points(narrow_v1$Day, narrow_v1$Mean_Ratio, col = "#56B4E9", pch = 17)

# Mark significant difference days
if (exists("ratio_comparison") && any(ratio_comparison$Significant, na.rm = TRUE)) {
  sig_days <- ratio_comparison$Day[ratio_comparison$Significant]
  abline(v = sig_days, col = "red", lty = 3, lwd = 0.5)
}

# Final L:W ratio comparison (mature leaves)
final_ratios <- mature_summary[mature_summary$Leaf_Type == "V1", ]
barplot(final_ratios$Final_Ratio,
        names.arg = final_ratios$Leaf,
        col = c("#E69F00", "#56B4E9"),
        main = "D. Mature V1 L:W Ratio",
        ylab = "Length:Width Ratio")

# Add error bars
if ("Ratio_SD" %in% colnames(final_ratios)) {
  x_pos <- barplot(final_ratios$Final_Ratio, plot = FALSE)
  arrows(x_pos, final_ratios$Final_Ratio - final_ratios$Ratio_SD,
         x_pos, final_ratios$Final_Ratio + final_ratios$Ratio_SD,
         angle = 90, code = 3, length = 0.1)
}

dev.off()
cat("Saved: v1_growth_trajectories.png\n")

# Plot 2: Gene-phenotype correlation distribution
if (exists("gene_pheno_cor") && !is.null(gene_pheno_cor)) {
  png("03_results/figures/31_phenotype/gene_pheno_correlation.png",
      width = 10, height = 6, units = "in", res = 150)

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Distribution of correlations
  hist(gene_pheno_cor$Cor_LW_Ratio, breaks = 50, col = "steelblue",
       main = "Gene-V1 L:W Ratio Correlations",
       xlab = "Spearman Correlation", ylab = "Frequency")
  abline(v = c(-0.3, 0.3), lty = 2, col = "red")

  # Volcano-style plot
  plot(gene_pheno_cor$Cor_LW_Ratio, -log10(gene_pheno_cor$Pval_LW_Ratio + 1e-300),
       pch = 16, cex = 0.5, col = "gray",
       xlab = "Correlation with V1 L:W Ratio",
       ylab = "-log10(p-value)",
       main = "Gene-Phenotype Association")

  # Highlight significant genes
  sig_genes <- gene_pheno_cor$FDR_LW_Ratio < 0.05 & abs(gene_pheno_cor$Cor_LW_Ratio) > 0.3
  if (sum(sig_genes, na.rm = TRUE) > 0) {
    points(gene_pheno_cor$Cor_LW_Ratio[sig_genes],
           -log10(gene_pheno_cor$Pval_LW_Ratio[sig_genes] + 1e-300),
           pch = 16, cex = 0.7, col = "red")
  }

  # Highlight JAG1 targets
  if ("Is_JAG1_Target" %in% colnames(gene_pheno_cor)) {
    jag1_idx <- gene_pheno_cor$Is_JAG1_Target
    if (sum(jag1_idx, na.rm = TRUE) > 0) {
      points(gene_pheno_cor$Cor_LW_Ratio[jag1_idx],
             -log10(gene_pheno_cor$Pval_LW_Ratio[jag1_idx] + 1e-300),
             pch = 1, cex = 1, col = "darkgreen", lwd = 2)
      legend("topright", legend = c("Significant", "JAG1 Target"),
             col = c("red", "darkgreen"), pch = c(16, 1), cex = 0.8)
    }
  }

  dev.off()
  cat("Saved: gene_pheno_correlation.png\n")
}

# Plot 3: Ratio divergence over time
if (exists("ratio_comparison")) {
  png("03_results/figures/31_phenotype/ratio_divergence.png",
      width = 10, height = 6, units = "in", res = 150)

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Plot ratio difference over time
  plot(ratio_comparison$Day, ratio_comparison$Difference,
       type = "b", pch = 16, col = "darkblue", lwd = 2,
       xlab = "Day", ylab = "Narrow - Broad Ratio Difference",
       main = "V1 L:W Ratio Divergence Over Time")
  abline(h = 0, lty = 2, col = "gray")

  # Mark significant days
  if (any(ratio_comparison$Significant, na.rm = TRUE)) {
    sig_idx <- ratio_comparison$Significant
    points(ratio_comparison$Day[sig_idx], ratio_comparison$Difference[sig_idx],
           pch = 16, cex = 1.5, col = "red")
    legend("topleft", legend = c("Non-significant", "Significant (p<0.05)"),
           col = c("darkblue", "red"), pch = 16)
  }

  # Plot p-values over time
  plot(ratio_comparison$Day, -log10(ratio_comparison$P_value + 1e-10),
       type = "b", pch = 16, col = "darkblue", lwd = 2,
       xlab = "Day", ylab = "-log10(p-value)",
       main = "Statistical Significance Over Time")
  abline(h = -log10(0.05), lty = 2, col = "red")
  text(max(ratio_comparison$Day), -log10(0.05), "p=0.05", pos = 3, col = "red")

  dev.off()
  cat("Saved: ratio_divergence.png\n")
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 11: SAVE CHECKPOINT\n")
cat("========================================\n\n")

phenotype_results <- list(
  pheno_data = pheno_data,
  pheno_long = pheno_long,
  pheno_summary = pheno_summary,
  mature_summary = mature_summary,
  v1_trajectory = v1_trajectory,
  ratio_comparison = ratio_comparison,
  pheno_traits = pheno_traits,
  trait_matrix = if (exists("trait_matrix")) trait_matrix else NULL,
  gene_pheno_cor = if (exists("gene_pheno_cor")) gene_pheno_cor else NULL
)

save(
  phenotype_results,
  pheno_data,
  pheno_long,
  pheno_summary,
  mature_summary,
  pheno_traits,
  file = "03_results/checkpoints/31_phenotype_integration.RData"
)

cat("Checkpoint saved: 31_phenotype_integration.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 31 COMPLETE: PHENOTYPE INTEGRATION\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("PHENOTYPE ANALYSIS SUMMARY:\n")
cat("  Leaf types analyzed:", paste(unique(pheno_data$Leaf), collapse = ", "), "\n")
cat("  Leaf positions (true leaves):", paste(true_leaves, collapse = ", "), "\n")
cat("  Days of measurement:", min(pheno_data$Day), "-", max(pheno_data$Day), "\n")
cat("  Replicates per leaf type:", length(unique(pheno_long$Rep)), "\n\n")

cat("MATURE LEAF L:W RATIOS:\n")
v1_mat <- mature_summary[mature_summary$Leaf_Type == "V1", ]
if (nrow(v1_mat) > 0) {
  cat("  V1 Broad:", round(v1_mat$Final_Ratio[v1_mat$Leaf == "Broad"], 3), "\n")
  cat("  V1 Narrow:", round(v1_mat$Final_Ratio[v1_mat$Leaf == "Narrow"], 3), "\n")
}

if (exists("gene_pheno_cor") && !is.null(gene_pheno_cor)) {
  cat("\nGENE-PHENOTYPE CORRELATIONS:\n")
  cat("  Genes positively correlated with L:W ratio:",
      sum(gene_pheno_cor$Pheno_Association == "High_LW_Ratio", na.rm = TRUE), "\n")
  cat("  Genes negatively correlated with L:W ratio:",
      sum(gene_pheno_cor$Pheno_Association == "Low_LW_Ratio", na.rm = TRUE), "\n")

  if ("Is_JAG1_Target" %in% colnames(gene_pheno_cor)) {
    cat("\nJAG1 TARGET PHENOTYPE OVERLAP:\n")
    target_cor <- gene_pheno_cor[gene_pheno_cor$Is_JAG1_Target, ]
    cat("  JAG1 targets with high L:W correlation:",
        sum(target_cor$Pheno_Association == "High_LW_Ratio", na.rm = TRUE), "\n")
    cat("  JAG1 targets with low L:W correlation:",
        sum(target_cor$Pheno_Association == "Low_LW_Ratio", na.rm = TRUE), "\n")
  }
}

cat("\nOUTPUT FILES:\n")
cat("  - 03_results/tables/phenotype/phenotype_summary.csv\n")
cat("  - 03_results/tables/phenotype/mature_leaf_summary.csv\n")
cat("  - 03_results/tables/phenotype/ratio_divergence_by_day.csv\n")
cat("  - 03_results/tables/phenotype/phenotype_traits.csv\n")
cat("  - 03_results/tables/phenotype/trait_matrix.csv\n")
cat("  - 03_results/tables/phenotype/gene_phenotype_correlations.csv\n")
cat("  - 03_results/tables/phenotype/top_LW_ratio_genes.csv\n")
cat("  - 03_results/figures/31_phenotype/*.png\n\n")

cat("NEXT: Script 32 - Phenotype Validation\n\n")
