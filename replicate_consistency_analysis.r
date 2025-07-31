#-------------------------------------------------------#
#     Replicate Consistency and Sample Identity Analysis#
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# This script analyzes replicate consistency and identifies samples
# that may be more similar to different timepoints than their labeled timepoint

# Load required libraries
library(ggplot2)
library(pheatmap)
library(dendextend)
library(corrplot)
library(reshape2)
library(gridExtra)

# Load data if not already in environment
if (!exists("logCPM") || !exists("targets")) {
  cat("Loading required data files...\n")
  
  # Load targets file
  if (file.exists("Targets_Final.csv")) {
    targets <- read.csv("Targets_Final.csv", row.names = 1)
    targets <- targets[1:60,]  # First 60 samples
  } else {
    stop("Targets file not found. Please ensure Targets_Final.csv is in working directory.")
  }
  
  # Add grouping variables
  targets$group <- paste(targets$Line, targets$Timepoint, sep = "_")
  targets$group_rep <- paste(targets$Line, targets$Timepoint, targets$Rep, sep = "_")
  
  # Load expression data
  if (file.exists("results/tables/logCPM_all_samples.csv")) {
    logCPM <- read.csv("results/tables/logCPM_all_samples.csv", row.names = 1)
  } else if (file.exists("results/checkpoints/04_preprocessed_data.RData")) {
    load("results/checkpoints/04_preprocessed_data.RData")
  } else {
    stop("Expression data not found. Please run preprocessing scripts first.")
  }
  
  cat("Data loaded successfully.\n")
  cat("Samples:", ncol(logCPM), "\n")
  cat("Genes:", nrow(logCPM), "\n")
}

# Create output directory
dir.create("results/replicate_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/replicate_analysis/figures", showWarnings = FALSE)

#=============================================================
# SECTION 1: WITHIN-GROUP REPLICATE CORRELATION
#=============================================================

cat("\n=== WITHIN-GROUP REPLICATE CORRELATION ANALYSIS ===\n")

# Calculate correlations between all samples
sample_cors <- cor(logCPM, method = "pearson")

# Extract within-group replicate correlations
replicate_cors <- data.frame(
  Group = character(),
  Rep1_vs_Rep2 = numeric(),
  Rep1_vs_Rep3 = numeric(),
  Rep2_vs_Rep3 = numeric(),
  Mean_Correlation = numeric(),
  Min_Correlation = numeric(),
  stringsAsFactors = FALSE
)

for (grp in unique(targets$group)) {
  grp_samples <- which(targets$group == grp)
  
  if (length(grp_samples) == 3) {
    # Get pairwise correlations
    cors <- c(
      sample_cors[grp_samples[1], grp_samples[2]],
      sample_cors[grp_samples[1], grp_samples[3]],
      sample_cors[grp_samples[2], grp_samples[3]]
    )
    
    replicate_cors <- rbind(replicate_cors, data.frame(
      Group = grp,
      Rep1_vs_Rep2 = cors[1],
      Rep1_vs_Rep3 = cors[2],
      Rep2_vs_Rep3 = cors[3],
      Mean_Correlation = mean(cors),
      Min_Correlation = min(cors)
    ))
  }
}

# Sort by minimum correlation to identify problematic groups
replicate_cors <- replicate_cors[order(replicate_cors$Min_Correlation), ]

cat("\nGroups with lowest replicate correlations:\n")
print(head(replicate_cors, 10))

# Save replicate correlation summary
write.csv(replicate_cors, 
          file = "results/replicate_analysis/replicate_correlations.csv",
          row.names = FALSE)

# Plot distribution of replicate correlations
pdf("results/replicate_analysis/figures/01_replicate_correlation_distribution.pdf", 
    width = 10, height = 6)

par(mfrow = c(1, 2))

# Histogram
hist(replicate_cors$Mean_Correlation, 
     breaks = 20,
     col = "lightblue",
     main = "Distribution of Mean Replicate Correlations",
     xlab = "Mean Correlation",
     xlim = c(0.8, 1))
abline(v = mean(replicate_cors$Mean_Correlation), col = "red", lwd = 2)
abline(v = 0.95, col = "green", lty = 2)

# Boxplot by line
replicate_cors$Line <- sub("_TP[0-9]", "", replicate_cors$Group)
boxplot(Mean_Correlation ~ Line, data = replicate_cors,
        col = rainbow(4),
        main = "Replicate Correlations by Line",
        ylab = "Mean Correlation",
        las = 2)
abline(h = 0.95, col = "green", lty = 2)

dev.off()

#=============================================================
# SECTION 2: SAMPLE-TO-TIMEPOINT SIMILARITY ANALYSIS
#=============================================================

cat("\n=== SAMPLE-TO-TIMEPOINT SIMILARITY ANALYSIS ===\n")

# For each sample, calculate its correlation to all timepoint groups
sample_timepoint_similarity <- data.frame()

for (i in 1:nrow(targets)) {
  sample_name <- rownames(targets)[i]
  sample_line <- targets$Line[i]
  sample_tp <- targets$Timepoint[i]
  sample_rep <- targets$Rep[i]
  
  # Calculate correlation to each timepoint group in the same line
  tp_cors <- data.frame(
    Sample = sample_name,
    Line = sample_line,
    Labeled_Timepoint = sample_tp,
    Rep = sample_rep
  )
  
  for (tp in unique(targets$Timepoint)) {
    # Get all samples from this timepoint in the same line
    tp_samples <- which(targets$Line == sample_line & targets$Timepoint == tp)
    
    if (length(tp_samples) > 0) {
      # Calculate mean correlation to this timepoint group
      # Exclude self-correlation
      cors_to_tp <- sample_cors[i, tp_samples]
      cors_to_tp <- cors_to_tp[tp_samples != i]
      
      if (length(cors_to_tp) > 0) {
        tp_cors[[paste0("Cor_to_", tp)]] <- mean(cors_to_tp)
      } else {
        tp_cors[[paste0("Cor_to_", tp)]] <- NA
      }
    }
  }
  
  sample_timepoint_similarity <- rbind(sample_timepoint_similarity, tp_cors)
}

# Identify best matching timepoint for each sample
cor_cols <- grep("^Cor_to_", colnames(sample_timepoint_similarity))
sample_timepoint_similarity$Best_Match_TP <- apply(
  sample_timepoint_similarity[, cor_cols], 1, 
  function(x) {
    if (all(is.na(x))) return(NA)
    sub("Cor_to_", "", names(which.max(x)))
  }
)

# Flag mismatches
sample_timepoint_similarity$Mismatch <- 
  sample_timepoint_similarity$Labeled_Timepoint != sample_timepoint_similarity$Best_Match_TP

# Calculate mismatch score (difference between labeled TP correlation and best TP correlation)
sample_timepoint_similarity$Mismatch_Score <- apply(
  sample_timepoint_similarity[, cor_cols], 1,
  function(x) {
    if (all(is.na(x))) return(NA)
    max(x, na.rm = TRUE) - x[paste0("Cor_to_", sample_timepoint_similarity$Labeled_Timepoint[1])]
  }
)

# Identify problematic samples
problematic_samples <- sample_timepoint_similarity[
  sample_timepoint_similarity$Mismatch & 
  !is.na(sample_timepoint_similarity$Mismatch), 
]

cat("\nSamples with potential timepoint mismatches:\n")
if (nrow(problematic_samples) > 0) {
  print(problematic_samples[, c("Sample", "Line", "Labeled_Timepoint", 
                                "Best_Match_TP", "Mismatch_Score")])
} else {
  cat("No clear timepoint mismatches detected.\n")
}

# Save full analysis
write.csv(sample_timepoint_similarity, 
          file = "results/replicate_analysis/sample_timepoint_similarity.csv",
          row.names = FALSE)

#=============================================================
# SECTION 3: VISUALIZATION OF SAMPLE RELATIONSHIPS
#=============================================================

cat("\n=== CREATING SAMPLE RELATIONSHIP VISUALIZATIONS ===\n")

# Create heatmaps for each line showing sample-to-sample correlations
pdf("results/replicate_analysis/figures/02_sample_correlation_heatmaps.pdf", 
    width = 14, height = 12)

par(mfrow = c(2, 2))

for (line in unique(targets$Line)) {
  line_samples <- which(targets$Line == line)
  line_cors <- sample_cors[line_samples, line_samples]
  
  # Order by timepoint and rep
  sample_order <- order(targets$Timepoint[line_samples], targets$Rep[line_samples])
  line_cors <- line_cors[sample_order, sample_order]
  
  # Create labels
  labels <- paste(targets$Timepoint[line_samples], 
                  "R", targets$Rep[line_samples], sep = "")[sample_order]
  
  # Create heatmap
  heatmap(line_cors,
          Rowv = NA, Colv = NA,
          scale = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = paste("Sample Correlations -", line),
          labRow = labels,
          labCol = labels,
          margins = c(6, 6))
}

dev.off()

# Create detailed heatmap for samples with mismatches
if (nrow(problematic_samples) > 0) {
  pdf("results/replicate_analysis/figures/03_problematic_samples_detail.pdf", 
      width = 12, height = 10)
  
  for (i in 1:min(nrow(problematic_samples), 4)) {
    sample_idx <- which(rownames(targets) == problematic_samples$Sample[i])
    sample_line <- problematic_samples$Line[i]
    
    # Get all samples from this line
    line_samples <- which(targets$Line == sample_line)
    
    # Create correlation profile
    cors_to_line <- sample_cors[sample_idx, line_samples]
    
    # Organize by timepoint
    tp_means <- sapply(unique(targets$Timepoint), function(tp) {
      tp_idx <- line_samples[targets$Timepoint[line_samples] == tp]
      mean(cors_to_line[tp_idx])
    })
    
    # Plot
    par(mfrow = c(2, 1))
    
    # Bar plot of correlations to each timepoint
    barplot(tp_means,
            names.arg = unique(targets$Timepoint),
            col = ifelse(names(tp_means) == problematic_samples$Labeled_Timepoint[i], 
                        "red", "gray"),
            main = paste("Sample:", problematic_samples$Sample[i],
                        "\nLabeled as", problematic_samples$Labeled_Timepoint[i],
                        "but most similar to", problematic_samples$Best_Match_TP[i]),
            ylab = "Mean Correlation to Timepoint",
            ylim = c(0.8, 1))
    
    # Individual sample correlations
    plot(cors_to_line,
         pch = 19,
         col = as.numeric(factor(targets$Timepoint[line_samples])),
         xlab = "Sample Index",
         ylab = "Correlation",
         main = "Correlation to Individual Samples")
    legend("bottomleft", 
           legend = unique(targets$Timepoint),
           col = 1:5,
           pch = 19)
  }
  
  dev.off()
}

#=============================================================
# SECTION 4: REPLICATE OUTLIER ANALYSIS
#=============================================================

cat("\n=== REPLICATE OUTLIER ANALYSIS ===\n")

# For each group, identify outlier replicates
outlier_analysis <- data.frame()

for (grp in unique(targets$group)) {
  grp_samples <- which(targets$group == grp)
  
  if (length(grp_samples) == 3) {
    # Calculate mean correlation for each replicate to the other two
    rep_cors <- numeric(3)
    
    for (i in 1:3) {
      other_reps <- grp_samples[-i]
      rep_cors[i] <- mean(sample_cors[grp_samples[i], other_reps])
    }
    
    # Identify potential outlier (lowest correlation)
    outlier_idx <- which.min(rep_cors)
    outlier_score <- mean(rep_cors[-outlier_idx]) - rep_cors[outlier_idx]
    
    outlier_analysis <- rbind(outlier_analysis, data.frame(
      Group = grp,
      Outlier_Rep = targets$Rep[grp_samples[outlier_idx]],
      Outlier_Sample = rownames(targets)[grp_samples[outlier_idx]],
      Outlier_Correlation = rep_cors[outlier_idx],
      Other_Reps_Mean_Cor = mean(rep_cors[-outlier_idx]),
      Outlier_Score = outlier_score
    ))
  }
}

# Sort by outlier score
outlier_analysis <- outlier_analysis[order(outlier_analysis$Outlier_Score, 
                                           decreasing = TRUE), ]

cat("\nTop replicate outliers:\n")
print(head(outlier_analysis[outlier_analysis$Outlier_Score > 0.02, ], 10))

write.csv(outlier_analysis, 
          file = "results/replicate_analysis/replicate_outliers.csv",
          row.names = FALSE)

#=============================================================
# SECTION 5: DEVELOPMENTAL ASYNCHRONY ANALYSIS
#=============================================================

cat("\n=== DEVELOPMENTAL ASYNCHRONY ANALYSIS ===\n")

# For each replicate, calculate its "developmental position" 
# relative to the timepoint progression
developmental_position <- data.frame()

for (line in unique(targets$Line)) {
  line_samples <- which(targets$Line == line)
  
  # Calculate mean expression profile for each timepoint
  tp_profiles <- sapply(paste0("TP", 0:4), function(tp) {
    tp_samples <- line_samples[targets$Timepoint[line_samples] == tp]
    if (length(tp_samples) > 0) {
      rowMeans(logCPM[, tp_samples, drop = FALSE])
    } else {
      rep(NA, nrow(logCPM))
    }
  })
  
  # For each sample, calculate correlation to each timepoint profile
  for (i in line_samples) {
    sample_profile <- logCPM[, i]
    
    # Correlate with each timepoint profile
    tp_cors <- cor(sample_profile, tp_profiles, use = "complete.obs")
    
    # Convert to vector if it's a matrix
    if (is.matrix(tp_cors)) {
      tp_cors <- as.vector(tp_cors)
    }
    
    # Fit a curve to find "developmental position"
    tp_numeric <- 0:4
    
    # Skip if not enough valid correlations
    if (sum(!is.na(tp_cors)) >= 3) {
      # Create data frame for fitting
      fit_data <- data.frame(
        tp_numeric = tp_numeric[!is.na(tp_cors)],
        tp_cors = tp_cors[!is.na(tp_cors)]
      )
      
      # Fit polynomial to find peak
      fit <- lm(tp_cors ~ poly(tp_numeric, 2), data = fit_data)
      
      # Predict on fine grid to find maximum
      pred_data <- data.frame(tp_numeric = seq(0, 4, 0.1))
      pred_cors <- predict(fit, newdata = pred_data)
      dev_position <- seq(0, 4, 0.1)[which.max(pred_cors)]
      
      developmental_position <- rbind(developmental_position, data.frame(
        Sample = rownames(targets)[i],
        Line = line,
        Timepoint = targets$Timepoint[i],
        Rep = targets$Rep[i],
        Timepoint_Numeric = as.numeric(sub("TP", "", targets$Timepoint[i])),
        Developmental_Position = dev_position,
        Position_Difference = dev_position - as.numeric(sub("TP", "", targets$Timepoint[i])),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Identify samples with significant developmental asynchrony
if (nrow(developmental_position) > 0) {
  async_samples <- developmental_position[
    abs(developmental_position$Position_Difference) > 0.5, 
  ]
} else {
  async_samples <- data.frame()
}

cat("\nSamples with potential developmental asynchrony:\n")
if (nrow(async_samples) > 0) {
  print(async_samples)
} else {
  cat("No significant developmental asynchrony detected.\n")
}

# Plot developmental positions
if (nrow(developmental_position) > 0) {
  pdf("results/replicate_analysis/figures/04_developmental_positions.pdf", 
      width = 12, height = 10)
  
  p <- ggplot(developmental_position, aes(x = Timepoint_Numeric, y = Developmental_Position)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(aes(color = factor(Rep)), size = 3, alpha = 0.7) +
    facet_wrap(~ Line, scales = "free") +
    scale_color_manual(values = c("1" = "red", "2" = "green", "3" = "blue")) +
    labs(title = "Developmental Position Analysis",
         x = "Labeled Timepoint",
         y = "Estimated Developmental Position",
         color = "Replicate") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  dev.off()
} else {
  cat("No developmental position data to plot.\n")
}

#=============================================================
# SECTION 6: RECOMMENDATIONS AND SUMMARY
#=============================================================

cat("\n=== GENERATING RECOMMENDATIONS ===\n")

# Compile problematic samples
problematic_summary <- data.frame(
  Sample = character(),
  Issue = character(),
  Severity = character(),
  Recommendation = character(),
  stringsAsFactors = FALSE
)

# Add timepoint mismatches
if (nrow(problematic_samples) > 0) {
  for (i in 1:nrow(problematic_samples)) {
    if (problematic_samples$Mismatch_Score[i] > 0.05) {
      severity <- "High"
      rec <- "Consider removing or re-labeling"
    } else if (problematic_samples$Mismatch_Score[i] > 0.02) {
      severity <- "Medium"
      rec <- "Flag for careful interpretation"
    } else {
      severity <- "Low"
      rec <- "Monitor in downstream analyses"
    }
    
    problematic_summary <- rbind(problematic_summary, data.frame(
      Sample = problematic_samples$Sample[i],
      Issue = paste("Best matches", problematic_samples$Best_Match_TP[i], 
                    "not", problematic_samples$Labeled_Timepoint[i]),
      Severity = severity,
      Recommendation = rec
    ))
  }
}

# Add replicate outliers
high_outliers <- outlier_analysis[outlier_analysis$Outlier_Score > 0.05, ]
if (nrow(high_outliers) > 0) {
  for (i in 1:nrow(high_outliers)) {
    problematic_summary <- rbind(problematic_summary, data.frame(
      Sample = high_outliers$Outlier_Sample[i],
      Issue = paste("Low replicate correlation (", 
                    round(high_outliers$Outlier_Correlation[i], 3), ")"),
      Severity = "High",
      Recommendation = "Consider removing from analysis"
    ))
  }
}

# Add developmental asynchrony
if (nrow(async_samples) > 0) {
  for (i in 1:nrow(async_samples)) {
    problematic_summary <- rbind(problematic_summary, data.frame(
      Sample = async_samples$Sample[i],
      Issue = paste("Developmental asynchrony (", 
                    round(async_samples$Position_Difference[i], 2), 
                    " timepoints off)"),
      Severity = ifelse(abs(async_samples$Position_Difference[i]) > 1, 
                       "High", "Medium"),
      Recommendation = "Flag for biological interpretation"
    ))
  }
}

# Remove duplicates and save
problematic_summary <- unique(problematic_summary)
write.csv(problematic_summary, 
          file = "results/replicate_analysis/problematic_samples_summary.csv",
          row.names = FALSE)

cat("\nProblematic samples summary:\n")
print(table(problematic_summary$Severity))

# Create final summary report
summary_text <- paste0(
  "REPLICATE CONSISTENCY ANALYSIS SUMMARY\n",
  "======================================\n\n",
  
  "1. OVERALL REPLICATE QUALITY:\n",
  "   - Mean replicate correlation: ", 
  round(mean(replicate_cors$Mean_Correlation), 3), "\n",
  "   - Groups with correlation < 0.95: ", 
  sum(replicate_cors$Mean_Correlation < 0.95), " out of ", 
  nrow(replicate_cors), "\n\n",
  
  "2. TIMEPOINT MISMATCH ANALYSIS:\n",
  "   - Samples with wrong best match: ", 
  sum(sample_timepoint_similarity$Mismatch, na.rm = TRUE), " out of ", 
  nrow(sample_timepoint_similarity), "\n",
  "   - Samples with high mismatch score (>0.05): ",
  sum(sample_timepoint_similarity$Mismatch_Score > 0.05, na.rm = TRUE), "\n\n",
  
  "3. REPLICATE OUTLIERS:\n",
  "   - High outlier score (>0.05): ", 
  sum(outlier_analysis$Outlier_Score > 0.05), "\n",
  "   - Medium outlier score (0.02-0.05): ",
  sum(outlier_analysis$Outlier_Score > 0.02 & 
      outlier_analysis$Outlier_Score <= 0.05), "\n\n",
  
  "4. DEVELOPMENTAL ASYNCHRONY:\n",
  "   - Samples >0.5 timepoints off: ",
  sum(abs(developmental_position$Position_Difference) > 0.5), "\n",
  "   - Samples >1 timepoint off: ",
  sum(abs(developmental_position$Position_Difference) > 1), "\n\n",
  
  "5. RECOMMENDATIONS:\n"
)

if (nrow(problematic_summary) > 0) {
  summary_text <- paste0(summary_text,
    "   - Total problematic samples: ", 
    length(unique(problematic_summary$Sample)), "\n",
    "   - High severity issues: ",
    sum(problematic_summary$Severity == "High"), "\n",
    "   - Recommended for removal: ",
    sum(grepl("removing", problematic_summary$Recommendation)), "\n\n",
    
    "   SPECIFIC ACTIONS:\n"
  )
  
  # Add specific sample recommendations
  high_severity <- problematic_summary[problematic_summary$Severity == "High", ]
  if (nrow(high_severity) > 0) {
    summary_text <- paste0(summary_text,
      "   High priority samples to review:\n"
    )
    for (i in 1:min(nrow(high_severity), 5)) {
      summary_text <- paste0(summary_text,
        "   - ", high_severity$Sample[i], ": ",
        high_severity$Issue[i], "\n"
      )
    }
  }
} else {
  summary_text <- paste0(summary_text,
    "   - No significant issues detected\n",
    "   - All replicates show good consistency\n"
  )
}

writeLines(summary_text, "results/replicate_analysis/replicate_analysis_summary.txt")
cat("\n", summary_text)

# Create a CSV of samples recommended for removal
remove_samples <- problematic_summary[
  grepl("removing", problematic_summary$Recommendation), "Sample"
]

if (length(remove_samples) > 0) {
  writeLines(remove_samples, 
             "results/replicate_analysis/samples_to_remove.txt")
  cat("\nSamples recommended for removal saved to: samples_to_remove.txt\n")
}

cat("\n=== REPLICATE CONSISTENCY ANALYSIS COMPLETE ===\n")
cat("Results saved to: results/replicate_analysis/\n")

# Session info
sink("results/replicate_analysis/session_info.txt")
sessionInfo()
sink()
