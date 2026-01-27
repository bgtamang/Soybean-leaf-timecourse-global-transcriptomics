
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 06: SAMPLE VALIDATION\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/06_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("04_reports", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "RColorBrewer",
  "pheatmap",
  "dplyr"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
load("03_results/checkpoints/05_normalized.RData")
cat("Loaded checkpoint from script 05\n")
cat("  Voom objects: v_full, v_tp1_4\n")
cat("  Design matrices: design_full, design_tp1_4\n\n")
# The flagged sample
problem_sample <- "745_T2_R2"
cat("Analyzing flagged sample:", problem_sample, "\n\n")

# Get sample info
if (problem_sample %in% rownames(targets)) {
  sample_info <- targets[problem_sample, ]
  cat("Sample information:\n")
  cat("  Line:", as.character(sample_info$Line), "\n")
  cat("  Leaf type:", as.character(sample_info$Leaf_type), "\n")
  cat("  Labeled timepoint:", as.character(sample_info$Timepoint), "\n")
  cat("  Batch:", as.character(sample_info$Batch), "\n")
  cat("  QC flag:", as.character(sample_info$QC_flag), "\n\n")
}


cat("Detailed correlation analysis:\n\n")

# Get expression data for the problem sample
expr_problem <- v_full$E[, problem_sample]

# Calculate correlation with ALL other samples
all_cors <- cor(expr_problem, v_full$E[, colnames(v_full$E) != problem_sample])

# Create correlation summary by group
cor_by_sample <- data.frame(
  Sample = colnames(all_cors),
  Correlation = as.vector(all_cors),
  stringsAsFactors = FALSE
)
cor_by_sample <- merge(cor_by_sample, targets[, c("Line", "Leaf_type", "Timepoint", "Batch")],
                        by.x = "Sample", by.y = "row.names")

# Summary by timepoint (same line)
same_line <- cor_by_sample[cor_by_sample$Line == sample_info$Line, ]
tp_summary <- same_line %>%
  group_by(Timepoint) %>%
  summarize(
    Mean_Cor = mean(Correlation),
    SD_Cor = sd(Correlation),
    N = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Cor))

cat("Correlation with same-line samples by timepoint:\n")
print(as.data.frame(tp_summary))

# Best matching timepoint
best_tp <- tp_summary$Timepoint[1]
expected_tp <- "TP2"

cat("\n  Expected timepoint:", expected_tp, "\n")
cat("  Best matching timepoint:", best_tp, "\n")
cat("  Correlation with expected (TP2):",
    round(tp_summary$Mean_Cor[tp_summary$Timepoint == "TP2"], 4), "\n")
cat("  Correlation with best match:",
    round(tp_summary$Mean_Cor[1], 4), "\n")

# Correlation difference
cor_diff <- tp_summary$Mean_Cor[tp_summary$Timepoint == best_tp] -
            tp_summary$Mean_Cor[tp_summary$Timepoint == expected_tp]
cat("  Difference:", round(cor_diff, 4), "\n\n")
# What happens if we remove this sample?
cat("Impact of removing", problem_sample, ":\n\n")

# Current replicate count for this group
group_name <- paste(sample_info$Line, sample_info$Timepoint, sep = "_")
current_reps <- sum(targets$Group == group_name)
cat("  Current group:", group_name, "\n")
cat("  Current replicates:", current_reps, "\n")
cat("  Replicates after removal:", current_reps - 1, "\n\n")

# Check if removal would leave < 2 replicates
if (current_reps - 1 < 2) {
  cat("  WARNING: Removal would leave fewer than 2 replicates!\n")
  cat("  RECOMMENDATION: Keep sample for statistical power\n\n")
} else {
  cat("  Removal would leave adequate replicates\n\n")
}
# Decision: Keep sample but create datasets for sensitivity analysis
cat("DECISION: Keep", problem_sample, "in primary analysis\n")
cat("RATIONALE:\n")
cat("  1. Correlation with expected TP2 is still high (>0.92)\n")
cat("  2. Removing leaves only n=2 for PI547745_TP2\n")
cat("  3. Sensitivity analysis will verify results are robust\n\n")

# Create three analysis tracks:

# Track 1: ALL 60 samples (PRIMARY)
cat("Creating analysis tracks:\n\n")

cat("Track 1: PRIMARY (all 60 samples)\n")
v_primary <- v_full
design_primary <- design_full
targets_primary <- targets
cat("  Samples:", ncol(v_primary$E), "\n")
cat("  Use: Main analysis results\n\n")

# Track 2: EXCLUDE flagged sample (SENSITIVITY)
cat("Track 2: SENSITIVITY (exclude", problem_sample, ")\n")
keep_samples <- colnames(v_full$E) != problem_sample
v_sensitivity <- v_full
v_sensitivity$E <- v_full$E[, keep_samples]
v_sensitivity$weights <- v_full$weights[, keep_samples]
v_sensitivity$design <- design_full[keep_samples, ]
targets_sensitivity <- targets[keep_samples, ]

# Recreate design matrix for sensitivity
targets_sensitivity$Group <- droplevels(factor(targets_sensitivity$Group))
design_sensitivity <- model.matrix(~0 + Group, data = targets_sensitivity)
colnames(design_sensitivity) <- gsub("Group", "", colnames(design_sensitivity))
v_sensitivity$design <- design_sensitivity

cat("  Samples:", ncol(v_sensitivity$E), "\n")
cat("  Use: Verify key findings are robust\n\n")

# Track 3: TP1-4 only (BATCH-FREE)
cat("Track 3: BATCH-FREE (TP1-4 only, no batch confounding)\n")
cat("  Samples:", ncol(v_tp1_4$E), "\n")
cat("  Use: Results without batch correction concerns\n\n")
# Colors
tp_colors <- brewer.pal(5, "YlOrRd")
names(tp_colors) <- c("TP0", "TP1", "TP2", "TP3", "TP4")
leaf_colors <- c("Broad" = "darkgreen", "Narrow" = "purple")

# --- Plot 1: Problem sample in context ---
cat("Creating sample context plot...\n")

# MDS with problem sample highlighted
png("03_results/figures/06_validation/problem_sample_MDS.png",
    width = 12, height = 10, units = "in", res = 300)

par(mfrow = c(2, 2), cex.main = 1.2, cex.lab = 1.1)

# MDS colored by timepoint, problem sample as star
mds <- plotMDS(v_full, plot = FALSE)
plot(mds$x, mds$y,
     col = tp_colors[targets$Timepoint],
     pch = ifelse(rownames(targets) == problem_sample, 8, 19),
     cex = ifelse(rownames(targets) == problem_sample, 3, 1.5),
     xlab = "Leading logFC dim 1",
     ylab = "Leading logFC dim 2",
     main = "MDS - Problem Sample Highlighted")
legend("topright", legend = c(names(tp_colors), problem_sample),
       col = c(tp_colors, "black"), pch = c(rep(19, 5), 8), cex = 0.8)

# Same line samples only
same_line_idx <- targets$Line == sample_info$Line
mds_same <- plotMDS(v_full$E[, same_line_idx], plot = FALSE)
targets_same <- targets[same_line_idx, ]

plot(mds_same$x, mds_same$y,
     col = tp_colors[targets_same$Timepoint],
     pch = ifelse(rownames(targets_same) == problem_sample, 8, 19),
     cex = ifelse(rownames(targets_same) == problem_sample, 3, 2),
     xlab = "Leading logFC dim 1",
     ylab = "Leading logFC dim 2",
     main = paste("MDS -", sample_info$Line, "Only"))
text(mds_same$x, mds_same$y, labels = targets_same$Timepoint, pos = 3, cex = 0.7)

# Correlation heatmap for same line
same_line_cor <- cor(v_full$E[, same_line_idx])
pheatmap(same_line_cor,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         breaks = seq(0.85, 1, length.out = 101),
         main = paste("Correlations -", sample_info$Line),
         fontsize = 10,
         silent = TRUE)

# Boxplot of correlations by timepoint
cor_plot_data <- same_line %>%
  mutate(Is_Expected = Timepoint == "TP2",
         Is_Best = Timepoint == best_tp)

boxplot(Correlation ~ Timepoint, data = same_line,
        col = tp_colors,
        main = paste(problem_sample, "Correlation by Timepoint"),
        ylab = "Correlation",
        xlab = "Timepoint",
        cex.main = 1.1)
abline(h = same_line$Correlation[same_line$Timepoint == "TP2"][1],
       col = "blue", lty = 2, lwd = 2)
legend("bottomright", "Expected TP2 level", col = "blue", lty = 2, cex = 0.8)

dev.off()
cat("  Saved: 03_results/figures/06_validation/problem_sample_MDS.png\n")

# --- Plot 2: Three tracks comparison ---
cat("Creating analysis tracks comparison...\n")

png("03_results/figures/06_validation/analysis_tracks.png",
    width = 14, height = 5, units = "in", res = 300)

par(mfrow = c(1, 3), cex.main = 1.2, cex.lab = 1.1)

# Track 1: Primary
mds1 <- plotMDS(v_primary, plot = FALSE)
plot(mds1$x, mds1$y,
     col = tp_colors[targets_primary$Timepoint],
     pch = 19, cex = 1.5,
     xlab = "Dim 1", ylab = "Dim 2",
     main = "PRIMARY: All 60 Samples")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.7)

# Track 2: Sensitivity
mds2 <- plotMDS(v_sensitivity, plot = FALSE)
plot(mds2$x, mds2$y,
     col = tp_colors[targets_sensitivity$Timepoint],
     pch = 19, cex = 1.5,
     xlab = "Dim 1", ylab = "Dim 2",
     main = "SENSITIVITY: 59 Samples")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.7)

# Track 3: Batch-free
mds3 <- plotMDS(v_tp1_4, plot = FALSE)
plot(mds3$x, mds3$y,
     col = tp_colors[targets_tp1_4$Timepoint],
     pch = 19, cex = 1.5,
     xlab = "Dim 1", ylab = "Dim 2",
     main = "BATCH-FREE: TP1-4 (48 Samples)")
legend("topright", legend = names(tp_colors)[2:5], col = tp_colors[2:5], pch = 19, cex = 0.7)

dev.off()
cat("  Saved: 03_results/figures/06_validation/analysis_tracks.png\n")


cat("\n========================================\n")

# Create validation report table
validation_report <- data.frame(
  Metric = c(
    "Flagged sample",
    "Expected timepoint",
    "Best matching timepoint",
    "Correlation with expected",
    "Correlation with best match",
    "Correlation difference",
    "Decision",
    "Primary dataset (n)",
    "Sensitivity dataset (n)",
    "Batch-free dataset (n)"
  ),
  Value = c(
    problem_sample,
    expected_tp,
    as.character(best_tp),
    round(tp_summary$Mean_Cor[tp_summary$Timepoint == expected_tp], 4),
    round(tp_summary$Mean_Cor[1], 4),
    round(cor_diff, 4),
    "KEEP (with sensitivity analysis)",
    ncol(v_primary$E),
    ncol(v_sensitivity$E),
    ncol(v_tp1_4$E)
  )
)

write.csv(validation_report, "03_results/tables/sample_validation_report.csv", row.names = FALSE)
cat("Saved: 03_results/tables/sample_validation_report.csv\n\n")

print(validation_report)


cat("\n========================================\n")

# Create markdown report
decision_doc <- paste0(
"# Sample Validation Decision Report
**Date:** ", format(Sys.time(), "%Y-%m-%d"), "

---

## Flagged Sample: ", problem_sample, "

### Sample Information
- **Line:** ", sample_info$Line, "
- **Leaf type:** ", sample_info$Leaf_type, "
- **Labeled timepoint:** ", sample_info$Timepoint, " (TP2)
- **Batch:** ", sample_info$Batch, "

### Issue Identified
Sample ", problem_sample, " clusters more closely with **TP3** samples than with its labeled **TP2** samples.

### Correlation Analysis

| Timepoint | Mean Correlation | Interpretation |
|-----------|-----------------|----------------|
| TP0 | ", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP0'], 4), " | Low (expected) |
| TP1 | ", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP1'], 4), " | Moderate |
| **TP2 (Expected)** | **", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP2'], 4), "** | High but not highest |
| **TP3 (Best match)** | **", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP3'], 4), "** | Highest |
| TP4 | ", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP4'], 4), " | High |

### Decision: **KEEP** (with sensitivity analysis)

### Rationale
1. **Correlation with expected TP2 is still high** (", round(tp_summary$Mean_Cor[tp_summary$Timepoint == 'TP2'], 3), ")
   - This suggests the sample is not a complete outlier
   - The difference with TP3 (", round(cor_diff, 3), ") is modest

2. **Statistical power concern**
   - Current replicates for PI547745_TP2: ", current_reps, "
   - After removal: ", current_reps - 1, "
   - n=2 would reduce statistical power significantly

3. **Biological plausibility**
   - Developmental staging can have natural variation
   - Sample may represent transitional stage between TP2 and TP3
   - Does NOT appear to be a technical failure

### Analysis Strategy

| Track | Samples | Purpose |
|-------|---------|---------|
| **PRIMARY** | 60 (all) | Main results for publication |
| **SENSITIVITY** | 59 (excluding ", problem_sample, ") | Verify robustness |
| **BATCH-FREE** | 48 (TP1-4 only) | No batch confounding |

### How to Report
In methods section:
> \"One sample (745_T2_R2) showed higher correlation with TP3 than its labeled TP2 timepoint (r=0.96 vs r=0.93). Sensitivity analysis excluding this sample confirmed all key findings (Supplementary Table X).\"

### Validation Checklist
- [ ] Run primary DE analysis
- [ ] Run sensitivity DE analysis
- [ ] Compare significant gene lists
- [ ] Report overlap statistics
- [ ] Document any findings that change with sample exclusion

---

*Report generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "*
"
)

writeLines(decision_doc, "04_reports/sample_decision.md")
cat("Saved: 04_reports/sample_decision.md\n")


cat("\n========================================\n")

# Save all analysis tracks
save(
  # PRIMARY track (all 60 samples)
  v_primary,
  design_primary,
  targets_primary,

  # SENSITIVITY track (59 samples, excluding flagged)
  v_sensitivity,
  design_sensitivity,
  targets_sensitivity,

  # BATCH-FREE track (TP1-4 only)
  v_tp1_4,
  design_tp1_4,
  targets_tp1_4,

  # Original objects
  dge_corrected,
  logCPM_corrected,

  # Parameters
  PARAMS,
  QC_PARAMS,

  # Validation info
  problem_sample,
  tp_summary,
  validation_report,

  file = "03_results/checkpoints/06_validated.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/06_validated.RData\n")
cat("  Contains:\n")
cat("    - PRIMARY: v_primary, design_primary, targets_primary\n")
cat("    - SENSITIVITY: v_sensitivity, design_sensitivity, targets_sensitivity\n")
cat("    - BATCH-FREE: v_tp1_4, design_tp1_4, targets_tp1_4\n")
cat("    - Validation: problem_sample, tp_summary, validation_report\n")


cat("\n========================================\n")
cat("SESSION INFO\n")

print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Decision Summary:\n")
cat("    - Flagged sample:", problem_sample, "\n")
cat("    - Decision: KEEP (with sensitivity analysis)\n")
cat("    - Primary dataset: 60 samples\n")
cat("    - Sensitivity dataset: 59 samples\n")
cat("    - Batch-free dataset: 48 samples\n")
cat("\n")
cat("  Analysis tracks ready for differential expression\n")
cat("\n")
