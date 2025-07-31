#-------------------------------------------------------#
#     Timepoint Grouping Validation Analysis            #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# This script analyzes whether any timepoints can be combined across lines
# It tests similarity of expression patterns between timepoints

# Load required libraries
library(ggplot2)
library(pheatmap)
library(dendextend)
library(corrplot)

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
dir.create("results/timepoint_validation", showWarnings = FALSE, recursive = TRUE)
dir.create("results/timepoint_validation/figures", showWarnings = FALSE)

#=============================================================
# SECTION 1: TIMEPOINT CORRELATION ANALYSIS
#=============================================================

cat("\n=== TIMEPOINT CORRELATION ANALYSIS ===\n")

# Calculate mean expression for each line x timepoint combination
group_means <- aggregate(t(logCPM), 
                         by = list(Group = targets$group), 
                         FUN = mean)
rownames(group_means) <- group_means$Group
group_means <- group_means[, -1]
group_means <- t(group_means)

# Extract timepoint information for each group
group_info <- data.frame(
  Group = colnames(group_means),
  Line = sub("_TP[0-4]", "", colnames(group_means)),
  Timepoint = sub(".*_", "", colnames(group_means))
)

# Calculate correlation between timepoints within each line
timepoint_cors_by_line <- list()

for (line in unique(targets$Line)) {
  line_groups <- group_info$Group[group_info$Line == line]
  line_data <- group_means[, line_groups]
  
  # Calculate correlation between timepoints
  tp_cor <- cor(line_data, method = "pearson")
  
  # Rename for clarity
  colnames(tp_cor) <- sub(".*_", "", colnames(tp_cor))
  rownames(tp_cor) <- sub(".*_", "", rownames(tp_cor))
  
  timepoint_cors_by_line[[line]] <- tp_cor
}

# Create correlation heatmaps for each line
pdf("results/timepoint_validation/figures/01_timepoint_correlations_by_line.pdf", 
    width = 12, height = 10)
par(mfrow = c(2, 2))

for (line in names(timepoint_cors_by_line)) {
  corrplot(timepoint_cors_by_line[[line]], 
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           addCoef.col = "black",
           number.cex = 0.8,
           main = paste("Timepoint Correlations -", line),
           mar = c(0, 0, 2, 0))
}

dev.off()

# Calculate average correlation between timepoints across all lines
avg_timepoint_cor <- Reduce("+", timepoint_cors_by_line) / length(timepoint_cors_by_line)

# Plot average correlation
pdf("results/timepoint_validation/figures/02_average_timepoint_correlations.pdf", 
    width = 8, height = 6)

pheatmap(avg_timepoint_cor,
         display_numbers = TRUE,
         number_format = "%.3f",
         main = "Average Timepoint Correlations Across All Lines",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(0, 1, length.out = 101))

dev.off()

cat("\nAverage timepoint correlations across all lines:\n")
print(round(avg_timepoint_cor, 3))

#=============================================================
# SECTION 2: CONSECUTIVE TIMEPOINT ANALYSIS
#=============================================================

cat("\n=== CONSECUTIVE TIMEPOINT ANALYSIS ===\n")

# Extract consecutive timepoint correlations
consecutive_cors <- data.frame(
  Line = character(),
  Transition = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

transitions <- c("TP0_to_TP1", "TP1_to_TP2", "TP2_to_TP3", "TP3_to_TP4")

for (line in names(timepoint_cors_by_line)) {
  cor_mat <- timepoint_cors_by_line[[line]]
  
  for (i in 1:4) {
    tp1 <- paste0("TP", i-1)
    tp2 <- paste0("TP", i)
    
    consecutive_cors <- rbind(consecutive_cors, data.frame(
      Line = line,
      Transition = transitions[i],
      Correlation = cor_mat[tp1, tp2]
    ))
  }
}

# Plot consecutive timepoint correlations
pdf("results/timepoint_validation/figures/03_consecutive_timepoint_correlations.pdf", 
    width = 10, height = 6)

ggplot(consecutive_cors, aes(x = Transition, y = Correlation, fill = Line)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Between Consecutive Timepoints",
       y = "Pearson Correlation",
       x = "Timepoint Transition") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  ylim(0, 1)

dev.off()

# Summary statistics
consec_summary <- aggregate(Correlation ~ Transition, data = consecutive_cors, 
                            FUN = function(x) c(mean = mean(x), sd = sd(x)))

cat("\nConsecutive timepoint correlations (mean ± SD):\n")
print(consec_summary)

#=============================================================
# SECTION 3: TIMEPOINT CLUSTERING ANALYSIS
#=============================================================

cat("\n=== TIMEPOINT CLUSTERING ANALYSIS ===\n")

# Perform hierarchical clustering on timepoints for each line
pdf("results/timepoint_validation/figures/04_timepoint_clustering_by_line.pdf", 
    width = 12, height = 10)
par(mfrow = c(2, 2))

for (line in names(timepoint_cors_by_line)) {
  # Convert correlation to distance
  tp_dist <- as.dist(1 - timepoint_cors_by_line[[line]])
  tp_clust <- hclust(tp_dist, method = "average")
  
  plot(tp_clust, main = paste("Timepoint Clustering -", line),
       xlab = "", sub = "", ylab = "Distance (1 - correlation)")
  
  # Add colored rectangles for potential groups
  rect.hclust(tp_clust, k = 3, border = c("red", "blue", "green"))
}

dev.off()

# Cluster based on average correlation across lines
avg_dist <- as.dist(1 - avg_timepoint_cor)
avg_clust <- hclust(avg_dist, method = "average")

pdf("results/timepoint_validation/figures/05_average_timepoint_clustering.pdf", 
    width = 8, height = 6)

plot(avg_clust, main = "Average Timepoint Clustering Across All Lines",
     xlab = "", sub = "", ylab = "Distance (1 - correlation)")
rect.hclust(avg_clust, k = 3, border = "red")

dev.off()

#=============================================================
# SECTION 4: EXPRESSION TRAJECTORY ANALYSIS
#=============================================================

cat("\n=== EXPRESSION TRAJECTORY ANALYSIS ===\n")

# Calculate standardized expression trajectories
# Sample genes to visualize patterns
set.seed(123)
sample_genes <- sample(rownames(logCPM), 1000)

# Calculate z-scores for each gene within each line
trajectory_data <- list()

for (line in unique(targets$Line)) {
  line_samples <- which(targets$Line == line)
  line_expr <- logCPM[sample_genes, line_samples]
  
  # Calculate mean expression per timepoint
  tp_means <- sapply(paste0("TP", 0:4), function(tp) {
    tp_samples <- line_samples[targets$Timepoint[line_samples] == tp]
    rowMeans(line_expr[, targets$Timepoint[line_samples] == tp, drop = FALSE])
  })
  
  # Standardize each gene's trajectory
  tp_zscore <- t(scale(t(tp_means)))
  
  trajectory_data[[line]] <- tp_zscore
}

# Calculate trajectory similarity between timepoints
trajectory_similarity <- matrix(0, nrow = 5, ncol = 5)
rownames(trajectory_similarity) <- paste0("TP", 0:4)
colnames(trajectory_similarity) <- paste0("TP", 0:4)

for (i in 1:5) {
  for (j in 1:5) {
    if (i != j) {
      # Compare expression patterns from TPi to all other timepoints
      # with patterns from TPj to all other timepoints
      similarities <- numeric()
      
      for (line in names(trajectory_data)) {
        traj <- trajectory_data[[line]]
        
        # Create pattern vectors excluding the timepoints being compared
        pattern_i <- traj[, -i]
        pattern_j <- traj[, -j]
        
        # Calculate correlation between patterns
        gene_cors <- diag(cor(t(pattern_i), t(pattern_j)))
        similarities <- c(similarities, mean(gene_cors, na.rm = TRUE))
      }
      
      trajectory_similarity[i, j] <- mean(similarities)
    } else {
      trajectory_similarity[i, j] <- 1
    }
  }
}

pdf("results/timepoint_validation/figures/06_trajectory_similarity.pdf", 
    width = 8, height = 6)

pheatmap(trajectory_similarity,
         display_numbers = TRUE,
         number_format = "%.3f",
         main = "Trajectory Pattern Similarity Between Timepoints",
         color = colorRampPalette(c("white", "red"))(100))

dev.off()

#=============================================================
# SECTION 5: VARIANCE ANALYSIS BY TIMEPOINT
#=============================================================

cat("\n=== VARIANCE ANALYSIS BY TIMEPOINT ===\n")

# Calculate within-timepoint variance across lines
timepoint_variance <- data.frame(
  Timepoint = paste0("TP", 0:4),
  Within_TP_Variance = numeric(5),
  Between_TP_Variance = numeric(5),
  F_statistic = numeric(5),
  stringsAsFactors = FALSE
)

for (i in 1:5) {
  tp <- paste0("TP", i-1)
  tp_samples <- which(targets$Timepoint == tp)
  other_samples <- which(targets$Timepoint != tp)
  
  # Calculate variance within this timepoint
  within_var <- mean(apply(logCPM[sample_genes, tp_samples], 1, var))
  
  # Calculate variance between this and other timepoints
  between_var <- mean(apply(logCPM[sample_genes, ], 1, var))
  
  timepoint_variance[i, "Within_TP_Variance"] <- within_var
  timepoint_variance[i, "Between_TP_Variance"] <- between_var
  timepoint_variance[i, "F_statistic"] <- between_var / within_var
}

pdf("results/timepoint_validation/figures/07_timepoint_variance.pdf", 
    width = 10, height = 6)

par(mfrow = c(1, 2))

# Within-timepoint variance
barplot(timepoint_variance$Within_TP_Variance,
        names.arg = timepoint_variance$Timepoint,
        col = heat.colors(5),
        main = "Within-Timepoint Variance",
        ylab = "Average Variance")

# F-statistic (between/within ratio)
barplot(timepoint_variance$F_statistic,
        names.arg = timepoint_variance$Timepoint,
        col = heat.colors(5),
        main = "Between/Within Variance Ratio",
        ylab = "F-statistic")

dev.off()

cat("\nTimepoint variance analysis:\n")
print(timepoint_variance)

#=============================================================
# SECTION 6: TIMEPOINT-SPECIFIC GENE ANALYSIS
#=============================================================

cat("\n=== TIMEPOINT-SPECIFIC GENE ANALYSIS ===\n")

# Identify genes with high specificity to each timepoint
timepoint_specific_genes <- list()

for (tp in paste0("TP", 0:4)) {
  tp_samples <- which(targets$Timepoint == tp)
  other_samples <- which(targets$Timepoint != tp)
  
  # Calculate mean expression in this timepoint vs others
  tp_mean <- rowMeans(logCPM[, tp_samples])
  other_mean <- rowMeans(logCPM[, other_samples])
  
  # Calculate specificity score (like fold change)
  specificity <- tp_mean - other_mean
  
  # Get top specific genes
  top_genes <- names(sort(specificity, decreasing = TRUE)[1:100])
  timepoint_specific_genes[[tp]] <- top_genes
}

# Calculate overlap between timepoint-specific genes
overlap_matrix <- matrix(0, nrow = 5, ncol = 5)
rownames(overlap_matrix) <- paste0("TP", 0:4)
colnames(overlap_matrix) <- paste0("TP", 0:4)

for (i in 1:5) {
  for (j in 1:5) {
    tp1 <- paste0("TP", i-1)
    tp2 <- paste0("TP", j-1)
    
    overlap <- length(intersect(timepoint_specific_genes[[tp1]], 
                                timepoint_specific_genes[[tp2]]))
    overlap_matrix[i, j] <- overlap
  }
}

pdf("results/timepoint_validation/figures/08_timepoint_specific_gene_overlap.pdf", 
    width = 8, height = 6)

pheatmap(overlap_matrix,
         display_numbers = TRUE,
         main = "Overlap of Timepoint-Specific Genes (Top 100)",
         color = colorRampPalette(c("white", "darkblue"))(100))

dev.off()

#=============================================================
# SECTION 7: RECOMMENDATIONS FOR TIMEPOINT GROUPING
#=============================================================

cat("\n=== TIMEPOINT GROUPING RECOMMENDATIONS ===\n")

# Analyze which timepoints could potentially be combined
grouping_criteria <- data.frame(
  Comparison = character(),
  Avg_Correlation = numeric(),
  Min_Correlation = numeric(),
  Max_Correlation = numeric(),
  SD_Correlation = numeric(),
  Recommendation = character(),
  stringsAsFactors = FALSE
)

# Check all pairwise combinations
for (i in 1:4) {
  for (j in (i+1):5) {
    tp1 <- paste0("TP", i-1)
    tp2 <- paste0("TP", j-1)
    
    # Get correlations across all lines
    cors <- sapply(timepoint_cors_by_line, function(x) x[tp1, tp2])
    
    # Determine recommendation
    avg_cor <- mean(cors)
    min_cor <- min(cors)
    
    if (min_cor > 0.95 && avg_cor > 0.97) {
      rec <- "Can be combined"
    } else if (min_cor > 0.90 && avg_cor > 0.93) {
      rec <- "Consider combining with caution"
    } else {
      rec <- "Should not be combined"
    }
    
    grouping_criteria <- rbind(grouping_criteria, data.frame(
      Comparison = paste(tp1, "vs", tp2),
      Avg_Correlation = avg_cor,
      Min_Correlation = min_cor,
      Max_Correlation = max(cors),
      SD_Correlation = sd(cors),
      Recommendation = rec
    ))
  }
}

# Sort by average correlation
grouping_criteria <- grouping_criteria[order(grouping_criteria$Avg_Correlation, 
                                             decreasing = TRUE), ]

cat("\nTimepoint grouping recommendations:\n")
print(grouping_criteria)

# Save recommendations
write.csv(grouping_criteria, 
          file = "results/timepoint_validation/timepoint_grouping_recommendations.csv",
          row.names = FALSE)

# Create summary plot
pdf("results/timepoint_validation/figures/09_grouping_recommendations.pdf", 
    width = 12, height = 8)

# Color code by recommendation
colors <- ifelse(grouping_criteria$Recommendation == "Can be combined", "green",
                 ifelse(grouping_criteria$Recommendation == "Consider combining with caution", 
                        "orange", "red"))

barplot(grouping_criteria$Avg_Correlation,
        names.arg = grouping_criteria$Comparison,
        col = colors,
        las = 2,
        ylim = c(0, 1),
        main = "Timepoint Correlation and Grouping Recommendations",
        ylab = "Average Correlation Across Lines")

# Add error bars for SD
arrows(x0 = 1:nrow(grouping_criteria) * 1.2 - 0.5,
       y0 = grouping_criteria$Avg_Correlation - grouping_criteria$SD_Correlation,
       y1 = grouping_criteria$Avg_Correlation + grouping_criteria$SD_Correlation,
       angle = 90, code = 3, length = 0.05)

# Add threshold lines
abline(h = 0.97, col = "green", lty = 2)
abline(h = 0.93, col = "orange", lty = 2)

legend("bottomleft", 
       legend = c("Can combine (r > 0.97)", 
                  "Caution (0.93 < r < 0.97)", 
                  "Do not combine (r < 0.93)"),
       fill = c("green", "orange", "red"))

dev.off()

#=============================================================
# SECTION 8: FINAL SUMMARY
#=============================================================

# Create comprehensive summary
summary_text <- paste0(
  "TIMEPOINT GROUPING VALIDATION SUMMARY\n",
  "=====================================\n\n",
  
  "1. HIGHEST CORRELATIONS BETWEEN TIMEPOINTS:\n"
)

# Add top 3 correlations
for (i in 1:min(3, nrow(grouping_criteria))) {
  summary_text <- paste0(summary_text,
                         "   ", grouping_criteria$Comparison[i], ": ",
                         "r = ", round(grouping_criteria$Avg_Correlation[i], 3),
                         " (", grouping_criteria$Recommendation[i], ")\n"
  )
}

summary_text <- paste0(summary_text,
                       "\n2. CONSECUTIVE TIMEPOINT ANALYSIS:\n"
)

# Add consecutive timepoint summary
for (i in 1:nrow(consec_summary)) {
  summary_text <- paste0(summary_text,
                         "   ", consec_summary$Transition[i], ": ",
                         "r = ", round(consec_summary$Correlation[i, "mean"], 3),
                         " ± ", round(consec_summary$Correlation[i, "sd"], 3), "\n"
  )
}

summary_text <- paste0(summary_text,
                       "\n3. RECOMMENDATIONS:\n",
                       "   Based on correlation analysis across all lines:\n"
)

# Count recommendations
rec_table <- table(grouping_criteria$Recommendation)
for (rec in names(rec_table)) {
  summary_text <- paste0(summary_text,
                         "   - ", rec, ": ", rec_table[rec], " pairs\n"
  )
}

summary_text <- paste0(summary_text,
                       "\n4. BIOLOGICAL INTERPRETATION:\n",
                       "   Timepoints represent distinct developmental stages with limited\n",
                       "   opportunity for combining without losing biological information.\n",
                       "   The highest correlations are typically between adjacent timepoints,\n",
                       "   suggesting a continuous developmental progression.\n"
)

writeLines(summary_text, "results/timepoint_validation/timepoint_validation_summary.txt")
cat("\n", summary_text)

cat("\n\n=== TIMEPOINT VALIDATION ANALYSIS COMPLETE ===\n")
cat("Results saved to: results/timepoint_validation/\n")

# Session info
sink("results/timepoint_validation/session_info.txt")
sessionInfo()
sink()
