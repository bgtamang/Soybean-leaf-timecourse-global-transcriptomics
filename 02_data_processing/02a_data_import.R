
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 01: DATA IMPORT\n")


# Define base directory
# MODIFY THIS PATH if running on a different machine
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"

# Set working directory to Phase 2
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/checkpoints", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/01_import", recursive = TRUE, showWarnings = FALSE)
# List of required packages
required_packages <- c(
  "tximport",       # For Salmon import (if needed)
  "edgeR",          # For DGEList
  "ggplot2",        # For plotting
  "dplyr",          # For data manipulation
  "tidyr"           # For data reshaping
)

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_packages)
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))
cat("  All packages loaded successfully\n\n")


# Define key parameters for the analysis
PARAMS <- list(
  # Gene of interest
  JAG1 = "Glyma.20G116200",
  JAG2 = "Glyma.10G273800",

  # Minimum count threshold (genes must have at least this many counts in some samples)
  min_count = 10,

  # Minimum samples (genes must pass min_count in at least this many samples)
  min_samples = 3
)

cat("Analysis Parameters:\n")
cat("  JAG1 gene ID:", PARAMS$JAG1, "\n")
cat("  JAG2 gene ID:", PARAMS$JAG2, "\n")
cat("  Minimum count threshold:", PARAMS$min_count, "\n")
cat("  Minimum samples:", PARAMS$min_samples, "\n\n")
# All input data is stored in 01_data/ for self-contained analysis
data_dir <- "01_data"

# Primary data file (Salmon summarized output)
salmon_file <- file.path(data_dir, "SalmonSummarizedOutput.RData")

# Backup: original large file (if salmon summary not available)
backup_data_file <- file.path(base_dir, "Soybean_RNASeq.RData")

# Try to load the Salmon summarized output first (smaller, faster)
if (file.exists(salmon_file)) {
  cat("Loading SalmonSummarizedOutput.RData from 01_data/...\n")
  load(salmon_file)
  cat("  File loaded successfully\n")

  # List objects loaded
  cat("  Objects in file:\n")
  # The loaded objects should now be in the environment

} else if (file.exists(backup_data_file)) {
  cat("Loading Soybean_RNASeq.RData (this may take a minute)...\n")
  load(backup_data_file)
  cat("  File loaded successfully\n")

} else {
  stop("ERROR: No data files found! Expected:\n",
       "  ", salmon_file, "\n",
       "  OR ", backup_data_file)
}


cat("\nExamining loaded objects...\n")

# List all objects in environment (excluding functions and parameters)
loaded_objects <- ls()
loaded_objects <- loaded_objects[!loaded_objects %in% c("base_dir", "salmon_file",
                                                         "main_data_file", "required_packages",
                                                         "missing_packages", "PARAMS")]

cat("  Loaded objects:", paste(loaded_objects, collapse = ", "), "\n\n")
# The Salmon output should contain counts - find them
# Common object names from tximport: txi, tx.all, salmon_data, counts, etc.

# Check for common tximport output structures
if (exists("tx.all")) {
  cat("Found tximport object 'tx.all'\n")

  # tximport returns a list with counts, abundance, length
  if (is.list(tx.all) && "counts" %in% names(tx.all)) {
    raw_counts <- tx.all$counts
    cat("  Extracted counts matrix from tx.all$counts\n")
  } else if (is.list(tx.all) && "abundance" %in% names(tx.all)) {
    # If counts not available, use abundance (TPM)
    raw_counts <- tx.all$abundance
    cat("  Extracted abundance matrix from tx.all$abundance\n")
    cat("  WARNING: Using TPM values, not raw counts\n")
  } else {
    stop("tx.all exists but doesn't have expected structure")
  }

} else if (exists("txi")) {
  cat("Found tximport object 'txi'\n")
  raw_counts <- txi$counts
  cat("  Extracted counts matrix\n")

} else if (exists("salmon_gene_counts")) {
  cat("Found 'salmon_gene_counts'\n")
  raw_counts <- salmon_gene_counts

} else if (exists("counts")) {
  cat("Found 'counts' object\n")
  raw_counts <- counts

} else {
  # List all objects and their types
  cat("Looking for count matrix in loaded objects...\n")
  for (obj_name in loaded_objects) {
    obj <- get(obj_name)
    if (is.matrix(obj) || is.data.frame(obj)) {
      dims <- dim(obj)
      cat("  ", obj_name, ": ", dims[1], " rows x ", dims[2], " cols (matrix/df)\n", sep = "")
    } else if (is.list(obj)) {
      cat("  ", obj_name, ": list with ", length(obj), " elements\n", sep = "")
      if (length(obj) > 0) {
        cat("    Names:", paste(head(names(obj), 5), collapse = ", "), "\n")
      }
    }
  }

  # Try to identify the count matrix (should have ~55000 genes x 60 samples)
  for (obj_name in loaded_objects) {
    obj <- get(obj_name)
    if ((is.matrix(obj) || is.data.frame(obj)) && nrow(obj) > 50000 && ncol(obj) == 60) {
      cat("\n  Identified likely count matrix:", obj_name, "\n")
      raw_counts <- as.matrix(obj)
      break
    }
  }
}

# Verify we have the count matrix
if (!exists("raw_counts")) {
  cat("\nERROR: Could not automatically identify count matrix.\n")
  cat("Please inspect the loaded objects and modify this script.\n")
  cat("Available objects:\n")
  print(ls())
  stop("Count matrix not found")
}

# Ensure it's a matrix
raw_counts <- as.matrix(raw_counts)

# Sample names have "salmo" prefix that needs to be removed to match metadata

cat("\nFixing sample names...\n")
original_names <- colnames(raw_counts)
# Remove "salmo" prefix
fixed_names <- gsub("^salmo", "", original_names)
colnames(raw_counts) <- fixed_names
cat("  Original:", head(original_names, 3), "...\n")
cat("  Fixed:", head(fixed_names, 3), "...\n")

# Data is at transcript level (Glyma.20G116200.1, .3, etc.)
# Need to sum counts across transcripts for each gene

cat("\nSummarizing transcripts to gene level...\n")
cat("  Transcripts before:", nrow(raw_counts), "\n")

# Extract gene ID by removing transcript suffix (last .N)
transcript_ids <- rownames(raw_counts)
gene_ids <- sub("\\.[0-9]+$", "", transcript_ids)

cat("  Example: ", transcript_ids[1], " -> ", gene_ids[1], "\n", sep = "")

# Sum counts by gene
unique_genes <- unique(gene_ids)
cat("  Unique genes:", length(unique_genes), "\n")

# Create gene-level count matrix by summing transcripts
gene_counts <- matrix(0, nrow = length(unique_genes), ncol = ncol(raw_counts))
rownames(gene_counts) <- unique_genes
colnames(gene_counts) <- colnames(raw_counts)

# Use rowsum for efficient aggregation
gene_counts <- rowsum(raw_counts, gene_ids)

cat("  Genes after summarization:", nrow(gene_counts), "\n")

# Replace transcript-level with gene-level
raw_counts <- gene_counts


cat("\nCount matrix dimensions:\n")
cat("  Genes:", nrow(raw_counts), "\n")
cat("  Samples:", ncol(raw_counts), "\n")

# Check sample names
cat("\nSample names (first 10):\n")
print(head(colnames(raw_counts), 10))

# Check gene names
cat("\nGene names (first 10):\n")
print(head(rownames(raw_counts), 10))

# Check if JAG1 is present
cat("\nChecking for JAG1 (", PARAMS$JAG1, "):\n", sep = "")
if (PARAMS$JAG1 %in% rownames(raw_counts)) {
  cat("  JAG1 FOUND in count matrix\n")
  cat("  JAG1 counts (first 5 samples):", head(raw_counts[PARAMS$JAG1, ], 5), "\n")
} else {
  cat("  WARNING: JAG1 not found! Check gene ID format.\n")
  # Try to find similar IDs
  similar <- grep("20G116", rownames(raw_counts), value = TRUE)
  if (length(similar) > 0) {
    cat("  Similar gene IDs found:", paste(similar, collapse = ", "), "\n")
  }
}


cat("\n========================================\n")

# Total counts per sample
lib_sizes <- colSums(raw_counts)

cat("Library sizes (total counts per sample):\n")
cat("  Mean:", format(mean(lib_sizes), big.mark = ","), "\n")
cat("  Median:", format(median(lib_sizes), big.mark = ","), "\n")
cat("  Min:", format(min(lib_sizes), big.mark = ","), "(", names(which.min(lib_sizes)), ")\n")
cat("  Max:", format(max(lib_sizes), big.mark = ","), "(", names(which.max(lib_sizes)), ")\n")

# Genes detected
genes_detected <- rowSums(raw_counts > 0)
cat("\nGenes detected:\n")
cat("  Total genes in matrix:", nrow(raw_counts), "\n")
cat("  Genes detected in at least 1 sample:", sum(genes_detected > 0), "\n")
cat("  Genes detected in all 60 samples:", sum(genes_detected == 60), "\n")
cat("  Genes with zero counts in all samples:", sum(genes_detected == 0), "\n")

# Count distribution
cat("\nCount distribution (non-zero values):\n")
nonzero_counts <- raw_counts[raw_counts > 0]
cat("  Median:", median(nonzero_counts), "\n")
cat("  Mean:", round(mean(nonzero_counts), 2), "\n")
cat("  Max:", max(nonzero_counts), "\n")


# Summary by sample
sample_summary <- data.frame(
  Sample = colnames(raw_counts),
  Total_Counts = lib_sizes,
  Genes_Detected = colSums(raw_counts > 0),
  Genes_GT10 = colSums(raw_counts > 10),
  Median_Count = apply(raw_counts, 2, function(x) median(x[x > 0])),
  stringsAsFactors = FALSE
)

# Save summary
write.csv(sample_summary,
          file = "03_results/tables/raw_counts_summary.csv",
          row.names = FALSE)
cat("\nSaved: 03_results/tables/raw_counts_summary.csv\n")


cat("\n========================================\n")

# Plot 1: Library sizes
cat("Creating library size plot...\n")

png("03_results/figures/01_import/library_sizes.png",
    width = 12, height = 6, units = "in", res = 300)

# Bar plot of library sizes
par(mar = c(10, 5, 4, 2))
barplot(lib_sizes / 1e6,
        names.arg = names(lib_sizes),
        las = 2,
        col = "steelblue",
        main = "Library Sizes (Raw Counts)",
        ylab = "Total Counts (millions)",
        cex.names = 0.7,
        cex.main = 1.2,
        cex.lab = 1.1)
abline(h = mean(lib_sizes) / 1e6, col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Mean:", round(mean(lib_sizes)/1e6, 1), "M"),
       col = "red", lty = 2, bty = "n", cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/01_import/library_sizes.png\n")

# Plot 2: Genes detected per sample
png("03_results/figures/01_import/genes_detected.png",
    width = 12, height = 6, units = "in", res = 300)

par(mar = c(10, 5, 4, 2))
barplot(colSums(raw_counts > 0) / 1000,
        names.arg = colnames(raw_counts),
        las = 2,
        col = "darkgreen",
        main = "Genes Detected per Sample",
        ylab = "Genes Detected (thousands)",
        cex.names = 0.7,
        cex.main = 1.2,
        cex.lab = 1.1)
abline(h = mean(colSums(raw_counts > 0)) / 1000, col = "red", lty = 2, lwd = 2)

dev.off()
cat("  Saved: 03_results/figures/01_import/genes_detected.png\n")

# Plot 3: Count distribution (log scale)
png("03_results/figures/01_import/count_distribution.png",
    width = 10, height = 6, units = "in", res = 300)

# Log-transform for visualization (add 1 to avoid log(0))
log_counts <- log2(raw_counts + 1)

# Boxplot of log counts
par(mar = c(10, 5, 4, 2))
boxplot(log_counts,
        las = 2,
        col = "lightblue",
        main = "Count Distribution (log2 scale)",
        ylab = "log2(count + 1)",
        cex.axis = 0.7,
        cex.main = 1.2,
        cex.lab = 1.1,
        outline = FALSE)  # Don't show outliers for clarity

dev.off()
cat("  Saved: 03_results/figures/01_import/count_distribution.png\n")


cat("\n========================================\n")

# Save checkpoint with all important objects
save(
  raw_counts,           # The main count matrix
  sample_summary,       # Summary statistics per sample
  PARAMS,               # Analysis parameters
  file = "03_results/checkpoints/01_data_imported.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/01_data_imported.RData\n")
cat("  Contains: raw_counts, sample_summary, PARAMS\n")


cat("\n========================================\n")
cat("SESSION INFO\n")

# Print session info for reproducibility
print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Genes:", nrow(raw_counts), "\n")
cat("    - Samples:", ncol(raw_counts), "\n")
cat("    - Mean library size:", format(round(mean(lib_sizes)), big.mark = ","), "\n")
cat("    - JAG1 present:", PARAMS$JAG1 %in% rownames(raw_counts), "\n")
cat("\n")
