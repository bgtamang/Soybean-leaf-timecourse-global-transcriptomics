
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 02: SAMPLE METADATA\n")


# Define base directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/02_metadata", recursive = TRUE, showWarnings = FALSE)
required_packages <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer")

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")


cat("Loading previous checkpoint...\n")
load("03_results/checkpoints/01_data_imported.RData")
cat("  Loaded: raw_counts (", nrow(raw_counts), " genes x ", ncol(raw_counts), " samples)\n\n", sep = "")
# All input data is stored in 01_data/ for self-contained analysis
data_dir <- "01_data"

# Primary metadata file location
metadata_file <- file.path(data_dir, "Targets_Final.csv")

# Backup locations (if not in 01_data)
backup_paths <- c(
  file.path(base_dir, "Phase1-Exploratory/Targets_Final.csv"),
  file.path(base_dir, "Targets_Final.csv"),
  file.path(base_dir, "Targets_Final.txt")
)

# Try primary location first
targets <- NULL
if (file.exists(metadata_file)) {
  cat("Loading metadata from 01_data/Targets_Final.csv...\n")
  targets <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
} else {
  # Try backup locations
  for (path in backup_paths) {
    if (file.exists(path)) {
      cat("Found metadata file:", path, "\n")
      if (grepl("\\.csv$", path)) {
        targets <- read.csv(path, row.names = 1, stringsAsFactors = FALSE)
      } else {
        targets <- read.delim(path, row.names = 1, stringsAsFactors = FALSE)
      }
      break
    }
  }
}

if (is.null(targets)) {
  stop("ERROR: Could not find Targets_Final.csv\n",
       "  Expected: ", metadata_file)
}

# Keep only the 60 samples (in case there are extra rows)
targets <- targets[1:60, ]

cat("\nMetadata loaded:\n")
cat("  Samples:", nrow(targets), "\n")
cat("  Variables:", ncol(targets), "\n")
cat("  Column names:", paste(colnames(targets), collapse = ", "), "\n\n")
# Get sample names from count matrix
count_samples <- colnames(raw_counts)
meta_samples <- rownames(targets)

cat("Sample name comparison:\n")
cat("  Samples in count matrix:", length(count_samples), "\n")
cat("  Samples in metadata:", length(meta_samples), "\n")

# Check if they match
if (all(count_samples %in% meta_samples) && all(meta_samples %in% count_samples)) {
  cat("  STATUS: All samples match!\n\n")
} else {
  cat("  WARNING: Sample mismatch detected!\n")
  cat("  In counts but not metadata:", paste(setdiff(count_samples, meta_samples), collapse = ", "), "\n")
  cat("  In metadata but not counts:", paste(setdiff(meta_samples, count_samples), collapse = ", "), "\n\n")
}

# Reorder targets to match count matrix
targets <- targets[count_samples, ]
cat("Reordered metadata to match count matrix column order\n\n")
# Ensure key columns exist and are properly formatted
# Line
if ("Line" %in% colnames(targets)) {
  targets$Line <- factor(targets$Line,
                         levels = c("PI532462A", "LD112170", "PI612713B", "PI547745"))
  cat("Line factor levels:", paste(levels(targets$Line), collapse = ", "), "\n")
}

# Leaf type
if ("Leaf_type" %in% colnames(targets)) {
  targets$Leaf_type <- factor(targets$Leaf_type, levels = c("Broad", "Narrow"))
  cat("Leaf_type factor levels:", paste(levels(targets$Leaf_type), collapse = ", "), "\n")
}

# Timepoint
if ("Timepoint" %in% colnames(targets)) {
  targets$Timepoint <- factor(targets$Timepoint,
                              levels = c("TP0", "TP1", "TP2", "TP3", "TP4"))
  cat("Timepoint factor levels:", paste(levels(targets$Timepoint), collapse = ", "), "\n")
}

# Batch
if ("Batch" %in% colnames(targets)) {
  targets$Batch <- factor(targets$Batch)
  cat("Batch factor levels:", paste(levels(targets$Batch), collapse = ", "), "\n")
}

# Replicate
if ("Rep" %in% colnames(targets)) {
  targets$Rep <- factor(targets$Rep)
  cat("Rep factor levels:", paste(levels(targets$Rep), collapse = ", "), "\n")
}


cat("\nCreating derived variables...\n")

# Group: Line_Timepoint combination
targets$Group <- paste(targets$Line, targets$Timepoint, sep = "_")
targets$Group <- factor(targets$Group)
cat("  Created Group variable (", length(levels(targets$Group)), " levels)\n", sep = "")

# LeafLine: Leaf_type_Line combination
targets$LeafLine <- paste(targets$Leaf_type, targets$Line, sep = "_")
targets$LeafLine <- factor(targets$LeafLine)

# Full group: Line_Timepoint_Rep
targets$FullGroup <- paste(targets$Line, targets$Timepoint, targets$Rep, sep = "_")

# Numeric timepoint (for modeling)
targets$Timepoint_numeric <- as.numeric(gsub("TP", "", targets$Timepoint))

# Genotype (Broad vs Narrow - alternative name for Leaf_type)
targets$Genotype <- targets$Leaf_type

cat("  Created LeafLine, FullGroup, Timepoint_numeric, Genotype\n\n")
cat("Experimental Design:\n")
cat("  Total samples: 60\n")
cat("  Lines: 4 (2 Broad, 2 Narrow)\n")
cat("  Timepoints: 5 (TP0-TP4)\n")
cat("  Replicates: 3 per condition\n")
cat("  Batches: 2 (2021, 2022)\n\n")

# Design table
design_table <- table(targets$Line, targets$Timepoint)
cat("Samples per Line x Timepoint:\n")
print(design_table)

# Batch distribution
cat("\nBatch distribution:\n")
batch_table <- table(targets$Batch, targets$Timepoint)
print(batch_table)

cat("\nBatch distribution by Line:\n")
batch_line_table <- table(targets$Batch, targets$Line)
print(batch_line_table)

# Note about batch confounding
cat("\nIMPORTANT - Batch Structure:\n")
cat("  Batch 2022: TP0 samples (meristem tissue)\n")
cat("  Batch 2021: TP1-TP4 samples (leaf development)\n")
cat("  NOTE: Batch is partially confounded with timepoint!\n")
cat("        This will be addressed in batch correction step.\n\n")
# Check if QC columns exist in targets
qc_cols <- c("NumReads", "Assigned_to_gene", "Multimapped", "Unmapped")
available_qc <- qc_cols[qc_cols %in% colnames(targets)]

if (length(available_qc) > 0) {
  cat("Available QC metrics:", paste(available_qc, collapse = ", "), "\n\n")

  # Calculate mapping rate if possible
  if ("NumReads" %in% colnames(targets) && "Assigned_to_gene" %in% colnames(targets)) {
    targets$Mapping_rate <- targets$Assigned_to_gene / targets$NumReads * 100

    cat("Mapping rates:\n")
    cat("  Mean:", round(mean(targets$Mapping_rate), 2), "%\n")
    cat("  Min:", round(min(targets$Mapping_rate), 2), "% (",
        rownames(targets)[which.min(targets$Mapping_rate)], ")\n")
    cat("  Max:", round(max(targets$Mapping_rate), 2), "% (",
        rownames(targets)[which.max(targets$Mapping_rate)], ")\n\n")
  }
} else {
  cat("No QC metrics found in metadata\n\n")
}
# Initialize flag column
targets$QC_flag <- "OK"

# Flag sample 745_T2_R2 (identified in Phase 1 as potentially mislabeled)
problematic_sample <- "745_T2_R2"
if (problematic_sample %in% rownames(targets)) {
  targets[problematic_sample, "QC_flag"] <- "CHECK_TIMEPOINT"
  cat("Flagged sample:", problematic_sample, "\n")
  cat("  Reason: Clusters with TP3 samples, not TP2 (identified in Phase 1)\n")
  cat("  Action: Will be evaluated in QC script\n\n")
}

# Summary of flags
flag_summary <- table(targets$QC_flag)
cat("Sample QC flags:\n")
print(flag_summary)


cat("\n========================================\n")

# Plot 1: Experimental design layout
cat("Creating experimental design visualization...\n")

# Design by Line and Timepoint
design_df <- as.data.frame(table(targets$Line, targets$Timepoint))
colnames(design_df) <- c("Line", "Timepoint", "Count")

plot1 <- ggplot(design_df, aes(x = Timepoint, y = Line, fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(title = "Samples per Line x Timepoint",
       x = "Timepoint", y = "Line") +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/02_metadata/experimental_design.png", plot1,
       width = 10, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/02_metadata/experimental_design.png\n")

# Batch distribution
batch_df <- as.data.frame(table(targets$Batch, targets$Timepoint))
colnames(batch_df) <- c("Batch", "Timepoint", "Count")

plot2 <- ggplot(batch_df, aes(x = Timepoint, y = Batch, fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), size = 6) +
  scale_fill_gradient(low = "white", high = "coral") +
  theme_minimal(base_size = 14) +
  labs(title = "Batch Distribution",
       subtitle = "Note: TP0 in 2022, TP1-4 in 2021",
       x = "Timepoint", y = "Batch") +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/02_metadata/batch_distribution.png", plot2,
       width = 10, height = 5, dpi = 300)
cat("  Saved: 03_results/figures/02_metadata/batch_distribution.png\n")

# Plot 2: Sample information overview
# Create a visual sample map
sample_df <- targets %>%
  select(Line, Leaf_type, Timepoint, Batch, Rep) %>%
  mutate(Sample = rownames(targets))

# Color by leaf type
sample_df$Color <- ifelse(sample_df$Leaf_type == "Broad", "Broad (Green)", "Narrow (Purple)")

plot3 <- ggplot(sample_df, aes(x = Timepoint, y = interaction(Line, Rep),
                                fill = Leaf_type, alpha = Batch)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("Broad" = "darkgreen", "Narrow" = "purple")) +
  scale_alpha_manual(values = c("2021" = 0.6, "2022" = 1.0)) +
  theme_minimal(base_size = 14) +
  labs(title = "Complete Sample Layout",
       subtitle = "Alpha indicates batch (lighter = 2021, darker = 2022)",
       x = "Timepoint",
       y = "Line.Replicate") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/02_metadata/sample_overview.png", plot3,
       width = 12, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/02_metadata/sample_overview.png\n")

# Plot 3: Library sizes by group (if available from previous script)
if ("NumReads" %in% colnames(targets)) {
  plot4 <- ggplot(targets, aes(x = Timepoint, y = NumReads/1e6, fill = Leaf_type)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Broad" = "darkgreen", "Narrow" = "purple")) +
    facet_wrap(~Line, nrow = 1) +
    theme_bw(base_size = 14) +
    labs(title = "Sequencing Depth by Group",
         x = "Timepoint",
         y = "Total Reads (millions)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 16, face = "bold"))

  ggsave("03_results/figures/02_metadata/reads_by_group.png", plot4,
         width = 12, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/02_metadata/reads_by_group.png\n")
}


cat("\nSaving experimental design table...\n")

# Select key columns for export
export_cols <- c("Sample_name", "Line", "Leaf_type", "Timepoint", "Batch", "Rep",
                 "Group", "Genotype", "Timepoint_numeric", "QC_flag")
export_cols <- export_cols[export_cols %in% colnames(targets)]

design_export <- targets[, export_cols]
design_export$Sample <- rownames(targets)
design_export <- design_export[, c("Sample", export_cols)]

write.csv(design_export,
          file = "03_results/tables/experimental_design.csv",
          row.names = FALSE)
cat("  Saved: 03_results/tables/experimental_design.csv\n")


cat("\n========================================\n")

# Save checkpoint with all important objects
save(
  raw_counts,           # Count matrix (from script 01)
  targets,              # Complete metadata
  PARAMS,               # Analysis parameters
  file = "03_results/checkpoints/02_sample_metadata.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/02_sample_metadata.RData\n")
cat("  Contains: raw_counts, targets, PARAMS\n")


cat("\n========================================\n")
cat("SESSION INFO\n")

print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Experimental Design Summary:\n")
cat("    - 4 Lines: PI532462A, LD112170 (Broad); PI612713B, PI547745 (Narrow)\n")
cat("    - 5 Timepoints: TP0 (meristem) to TP4 (mature leaf)\n")
cat("    - 3 Replicates per condition\n")
cat("    - 2 Batches: 2022 (TP0), 2021 (TP1-TP4)\n")
cat("    - Flagged samples:", sum(targets$QC_flag != "OK"), "\n")
cat("\n")
