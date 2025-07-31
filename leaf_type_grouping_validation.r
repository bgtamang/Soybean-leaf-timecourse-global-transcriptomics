#-------------------------------------------------------#
#     Validation of Leaf Type Grouping Strategy         #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# This script quantifies the justification for grouping lines by leaf type
# It tests whether broad lines are more similar to each other than to narrow lines

# Load required libraries
library(edgeR)
library(limma)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dendextend)
library(corrplot)
library(vegan)  # For PERMANOVA
library(pheatmap)
library(VennDiagram)

# Set working directory (modify as needed)
# setwd("/path/to/your/project")

# Create output directory
dir.create("results/leaf_type_validation", showWarnings = FALSE, recursive = TRUE)
dir.create("results/leaf_type_validation/figures", showWarnings = FALSE)

#=============================================================
# SECTION 1: DATA LOADING AND PREPARATION
#=============================================================

cat("Loading required data files...\n")

# Load targets file
targets <- read.csv("Targets_Final.csv", row.names = 1)
targets <- targets[1:60,]  # First 60 samples

# Add grouping variables
targets$group <- paste(targets$Line, targets$Timepoint, sep = "_")
targets$Leaf_Line <- paste(targets$Leaf_type, targets$Line, sep = "_")

# Load expression data (assuming these files exist from previous analysis)
# If not available, load from checkpoints
if (file.exists("results/tables/logCPM_all_samples.csv")) {
  logCPM <- read.csv("results/tables/logCPM_all_samples.csv", row.names = 1)
} else if (file.exists("results/checkpoints/04_preprocessed_data.RData")) {
  load("results/checkpoints/04_preprocessed_data.RData")
} else {
  stop("Expression data not found. Please run preprocessing scripts first.")
}

# Load DE results if available
if (file.exists("results/checkpoints/06_de_results.RData")) {
  load("results/checkpoints/06_de_results.RData")
  de_data_available <- TRUE
} else {
  cat("DE results not found. Some analyses will be skipped.\n")
  de_data_available <- FALSE
}

# Verify data consistency
cat("\nData summary:\n")
cat("Samples:", ncol(logCPM), "\n")
cat("Genes:", nrow(logCPM), "\n")
cat("Lines:", paste(unique(targets$Line), collapse = ", "), "\n")
cat("Leaf types:", paste(unique(targets$Leaf_type), collapse = ", "), "\n")

#=============================================================
# SECTION 2: TRANSCRIPTOME-WIDE SIMILARITY ANALYSIS
#=============================================================

cat("\n=== TRANSCRIPTOME-WIDE SIMILARITY ANALYSIS ===\n")

# Calculate sample correlations
sample_cors <- cor(logCPM, method = "pearson")

# Define sample groups
broad_lines <- c("PI532462A", "LD112170")
narrow_lines <- c("PI612713B", "PI547745")

# Function to extract correlations by category
categorize_correlations <- function(cor_matrix, targets) {
  n <- nrow(cor_matrix)
  categories <- list()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_val <- cor_matrix[i, j]
      
      line_i <- targets$Line[i]
      line_j <- targets$Line[j]
      leaf_i <- targets$Leaf_type[i]
      leaf_j <- targets$Leaf_type[j]
      tp_i <- targets$Timepoint[i]
      tp_j <- targets$Timepoint[j]
      
      # Categorize
      if (line_i == line_j) {
        # Same line
        if (tp_i == tp_j) {
          category <- "Within_Line_Same_TP"
        } else {
          category <- "Within_Line_Diff_TP"
        }
      } else if (leaf_i == leaf_j) {
        # Different line, same leaf type
        if (tp_i == tp_j) {
          category <- "Between_Line_Same_Leaf_Same_TP"
        } else {
          category <- "Between_Line_Same_Leaf_Diff_TP"
        }
      } else {
        # Different leaf type
        if (tp_i == tp_j) {
          category <- "Between_Leaf_Same_TP"
        } else {
          category <- "Between_Leaf_Diff_TP"
        }
      }
      
      if (!category %in% names(categories)) {
        categories[[category]] <- numeric()
      }
      categories[[category]] <- c(categories[[category]], cor_val)
    }
  }
  
  return(categories)
}

# Get categorized correlations
cor_categories <- categorize_correlations(sample_cors, targets)

# Create boxplot
pdf("results/leaf_type_validation/figures/01_correlation_categories.pdf", width = 12, height = 8)
par(mar = c(12, 5, 4, 2))

# Prepare data for plotting
plot_data <- data.frame(
  Correlation = unlist(cor_categories),
  Category = rep(names(cor_categories), sapply(cor_categories, length))
)

# Order categories logically
category_order <- c("Within_Line_Same_TP", "Within_Line_Diff_TP", 
                    "Between_Line_Same_Leaf_Same_TP", "Between_Line_Same_Leaf_Diff_TP",
                    "Between_Leaf_Same_TP", "Between_Leaf_Diff_TP")
plot_data$Category <- factor(plot_data$Category, levels = category_order)

# Create boxplot
boxplot(Correlation ~ Category, data = plot_data,
        las = 2, col = c("darkgreen", "lightgreen", "darkblue", "lightblue", "darkred", "pink"),
        main = "Sample Correlation by Relationship Category",
        ylab = "Pearson Correlation")

# Add sample sizes
cat_sizes <- table(plot_data$Category)
mtext(paste("n =", cat_sizes), at = 1:length(cat_sizes), side = 1, line = 10, cex = 0.8)

dev.off()

# Statistical tests
cat("\nStatistical comparison of correlation categories:\n")
cat("Mean correlations:\n")
cor_means <- aggregate(Correlation ~ Category, data = plot_data, FUN = mean)
print(cor_means)

# Based on your results:
# Category                         Correlation
# 1            Within_Line_Same_TP   0.9769828
# 2            Within_Line_Diff_TP   0.8822746
# 3 Between_Line_Same_Leaf_Same_TP   0.9455257
# 4 Between_Line_Same_Leaf_Diff_TP   0.8692858
# 5           Between_Leaf_Same_TP   0.9443784
# 6           Between_Leaf_Diff_TP   0.8656051

# Key comparison: same leaf type vs different leaf type
same_leaf_cors <- plot_data$Correlation[plot_data$Category %in% 
                                         c("Between_Line_Same_Leaf_Same_TP", "Between_Line_Same_Leaf_Diff_TP")]
diff_leaf_cors <- plot_data$Correlation[plot_data$Category %in% 
                                         c("Between_Leaf_Same_TP", "Between_Leaf_Diff_TP")]

t_test_result <- t.test(same_leaf_cors, diff_leaf_cors)
cat("\n\nT-test: Same leaf type vs Different leaf type correlations\n")
cat("Mean same leaf type:", mean(same_leaf_cors), "\n")
cat("Mean different leaf type:", mean(diff_leaf_cors), "\n")
cat("Difference:", mean(same_leaf_cors) - mean(diff_leaf_cors), "\n")
cat("P-value:", t_test_result$p.value, "\n")

#=============================================================
# SECTION 3: HIERARCHICAL CLUSTERING
#=============================================================

cat("\n=== HIERARCHICAL CLUSTERING ANALYSIS ===\n")

# Perform hierarchical clustering
sample_dist <- as.dist(1 - sample_cors)
sample_clust <- hclust(sample_dist, method = "average")

# Create dendrogram with colored labels
pdf("results/leaf_type_validation/figures/02_sample_dendrogram.pdf", width = 14, height = 8)

# Convert to dendrogram object
dend <- as.dendrogram(sample_clust)

# Color labels by leaf type
label_colors <- ifelse(targets$Leaf_type == "Broad", "darkgreen", "purple")
labels_colors(dend) <- label_colors[order.dendrogram(dend)]

# Plot
plot(dend, main = "Sample Clustering Dendrogram", ylab = "Distance (1 - correlation)")

# Add legend
legend("topright", legend = c("Broad", "Narrow"), 
       fill = c("darkgreen", "purple"), bty = "n")

dev.off()

# Quantify clustering by leaf type
# Calculate cophenetic distance and test clustering
cophenetic_dist <- cophenetic(sample_clust)

# PERMANOVA to test if leaf type explains clustering
set.seed(123)
permanova_result <- adonis2(as.dist(1 - sample_cors) ~ Leaf_type + Line + Timepoint, 
                            data = targets, 
                            permutations = 999)

cat("\nPERMANOVA results (variance partitioning):\n")
print(permanova_result)

# Save PERMANOVA results
write.csv(as.data.frame(permanova_result), 
          file = "results/leaf_type_validation/PERMANOVA_results.csv")

#=============================================================
# SECTION 4: PRINCIPAL COMPONENT ANALYSIS
#=============================================================

cat("\n=== PRINCIPAL COMPONENT ANALYSIS ===\n")

# Perform PCA
pca_result <- prcomp(t(logCPM), scale. = TRUE)

# Calculate variance explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Create PCA plots
pdf("results/leaf_type_validation/figures/03_PCA_analysis.pdf", width = 14, height = 10)
par(mfrow = c(2, 2))

# Plot 1: Colored by leaf type
plot(pca_result$x[, 1:2], 
     col = ifelse(targets$Leaf_type == "Broad", "darkgreen", "purple"),
     pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Leaf Type")
legend("topright", legend = c("Broad", "Narrow"), 
       col = c("darkgreen", "purple"), pch = 19)

# Plot 2: Colored by line
line_colors <- rainbow(4)[as.numeric(factor(targets$Line))]
plot(pca_result$x[, 1:2], 
     col = line_colors,
     pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Line")
legend("topright", legend = unique(targets$Line), 
       col = rainbow(4), pch = 19, cex = 0.8)

# Plot 3: Colored by timepoint
tp_colors <- heat.colors(5)[as.numeric(factor(targets$Timepoint))]
plot(pca_result$x[, 1:2], 
     col = tp_colors,
     pch = 19, cex = 1.5,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Timepoint")
legend("topright", legend = unique(targets$Timepoint), 
       col = heat.colors(5), pch = 19, cex = 0.8)

# Plot 4: Variance explained
barplot(var_explained[1:10], 
        names.arg = paste0("PC", 1:10),
        ylab = "Variance Explained (%)",
        main = "Scree Plot")

dev.off()

# Quantify separation by leaf type in PC space
# Linear model to test how much PC1 and PC2 are explained by factors
pc_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  LeafType = targets$Leaf_type,
  Line = targets$Line,
  Timepoint = targets$Timepoint
)

# ANOVA for PC1
pc1_model <- lm(PC1 ~ LeafType + Line + Timepoint, data = pc_data)
pc1_anova <- anova(pc1_model)

cat("\nANOVA for PC1:\n")
print(pc1_anova)

# Calculate R-squared for leaf type alone
pc1_leaf_only <- lm(PC1 ~ LeafType, data = pc_data)
cat("\nR-squared for Leaf Type on PC1:", summary(pc1_leaf_only)$r.squared, "\n")

#=============================================================
# SECTION 5: VARIANCE PARTITIONING
#=============================================================

cat("\n=== VARIANCE PARTITIONING ANALYSIS ===\n")

# Select a subset of highly variable genes for computational efficiency
gene_vars <- apply(logCPM, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:5000])

# Prepare data for linear mixed model
expr_subset <- logCPM[top_var_genes, ]

# Melt data for modeling
expr_long <- melt(t(expr_subset))
colnames(expr_long) <- c("Sample", "Gene", "Expression")

# Add sample information
expr_long <- merge(expr_long, 
                   data.frame(Sample = rownames(targets), targets), 
                   by = "Sample")

# Fit mixed model (simplified version)
# For full analysis, use lme4 or nlme
library(nlme)

# Sample 100 genes for demonstration
set.seed(123)
sample_genes <- sample(top_var_genes, 100)

variance_components <- data.frame(
  Gene = character(),
  Var_LeafType = numeric(),
  Var_Line_within_Leaf = numeric(),
  Var_Timepoint = numeric(),
  Var_Residual = numeric(),
  stringsAsFactors = FALSE
)

cat("Calculating variance components for 100 genes...\n")
pb <- txtProgressBar(min = 0, max = length(sample_genes), style = 3)

for (i in 1:length(sample_genes)) {
  gene <- sample_genes[i]
  gene_expr <- logCPM[gene, ]
  
  # Create data frame for this gene
  gene_data <- data.frame(
    Expression = gene_expr,
    LeafType = targets$Leaf_type,
    Line = targets$Line,
    Timepoint = targets$Timepoint,
    Rep = targets$Rep
  )
  
  # Fit model
  tryCatch({
    # Use nested model to avoid collinearity
    # First fit full model with Line (which includes LeafType information)
    model_full <- lm(Expression ~ Line + Timepoint, data = gene_data)
    
    # Then fit reduced models to partition variance
    model_no_line <- lm(Expression ~ Timepoint, data = gene_data)
    model_no_time <- lm(Expression ~ Line, data = gene_data)
    model_null <- lm(Expression ~ 1, data = gene_data)
    
    # Calculate sum of squares
    ss_total <- sum(anova(model_null)$`Sum Sq`)
    ss_line <- sum(anova(model_full, model_no_line)$`Sum Sq`[2])
    ss_time <- sum(anova(model_full, model_no_time)$`Sum Sq`[2])
    ss_residual <- sum(anova(model_full)$`Sum Sq`[nrow(anova(model_full))])
    
    # For leaf type, compare model with just leaf type vs null
    model_leaf <- lm(Expression ~ LeafType, data = gene_data)
    ss_leaf <- sum(anova(model_leaf, model_null)$`Sum Sq`[2])
    
    # Line within leaf type (Line effect after accounting for leaf type)
    ss_line_within_leaf <- ss_line - ss_leaf
    
    variance_components <- rbind(variance_components, data.frame(
      Gene = gene,
      Var_LeafType = ss_leaf / ss_total,
      Var_Line_within_Leaf = ss_line_within_leaf / ss_total,
      Var_Timepoint = ss_time / ss_total,
      Var_Residual = ss_residual / ss_total
    ))
  }, error = function(e) {
    # Skip genes that cause errors
  })
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Summarize variance components
var_summary <- colMeans(variance_components[, -1], na.rm = TRUE)

# Check if we have valid data
if (nrow(variance_components) > 0 && !all(is.na(var_summary))) {
  # Create pie chart
  pdf("results/leaf_type_validation/figures/04_variance_components.pdf", width = 10, height = 8)
  par(mfrow = c(1, 2))
  
  # Pie chart - ensure all values are positive
  var_summary_plot <- pmax(var_summary, 0)  # Convert any negative values to 0
  
  if (sum(var_summary_plot) > 0) {
    pie(var_summary_plot, 
        labels = paste0(c("Leaf Type", "Line within Leaf", "Timepoint", "Residual"), 
                        "\n", round(var_summary_plot * 100, 1), "%"),
        col = c("darkgreen", "lightblue", "orange", "gray"),
        main = "Average Variance Components")
  } else {
    plot.new()
    text(0.5, 0.5, "Variance components could not be calculated", cex = 1.5)
  }
  
  # Boxplot of individual genes
  if (ncol(variance_components) > 1) {
    boxplot(variance_components[, -1], 
            names = c("Leaf Type", "Line in Leaf", "Timepoint", "Residual"),
            col = c("darkgreen", "lightblue", "orange", "gray"),
            ylab = "Proportion of Variance",
            main = "Distribution of Variance Components",
            las = 2)
  }
  
  dev.off()
  
  # Statistical summary
  cat("\nAverage variance explained by each factor:\n")
  print(round(var_summary * 100, 2))
} else {
  cat("\nWarning: Variance components could not be calculated.\n")
  cat("This may be due to the nested structure of Line within LeafType.\n")
  
  # Alternative approach: just compare models
  cat("\nTrying alternative variance decomposition...\n")
  
  # Use first 10 genes as example
  alt_var_results <- data.frame()
  
  for (i in 1:min(10, length(sample_genes))) {
    gene <- sample_genes[i]
    gene_expr <- as.numeric(logCPM[gene, ])  # Convert to numeric vector
    
    # Fit different models
    m_full <- lm(gene_expr ~ Line + Timepoint, data = targets)
    m_leaf <- lm(gene_expr ~ Leaf_type + Timepoint, data = targets)
    m_time <- lm(gene_expr ~ Timepoint, data = targets)
    m_null <- lm(gene_expr ~ 1, data = targets)
    
    # Get R-squared values
    alt_var_results <- rbind(alt_var_results, data.frame(
      Gene = gene,
      R2_Full = summary(m_full)$r.squared,
      R2_LeafType = summary(m_leaf)$r.squared,
      R2_Timepoint = summary(m_time)$r.squared
    ))
  }
  
  if (nrow(alt_var_results) > 0) {
    cat("\nAlternative variance analysis (R-squared values):\n")
    cat("Mean R2 with Line (includes leaf type):", round(mean(alt_var_results$R2_Full), 3), "\n")
    cat("Mean R2 with Leaf Type + Time:", round(mean(alt_var_results$R2_LeafType), 3), "\n")
    cat("Mean R2 with Time only:", round(mean(alt_var_results$R2_Timepoint), 3), "\n")
    
    # Approximate variance components
    var_summary <- c(
      Var_LeafType = mean(alt_var_results$R2_LeafType - alt_var_results$R2_Timepoint),
      Var_Line_within_Leaf = mean(alt_var_results$R2_Full - alt_var_results$R2_LeafType),
      Var_Timepoint = mean(alt_var_results$R2_Timepoint),
      Var_Residual = mean(1 - alt_var_results$R2_Full)
    )
  }
}

#=============================================================
# SECTION 6: DIFFERENTIAL EXPRESSION OVERLAP
#=============================================================

if (de_data_available) {
  cat("\n=== DIFFERENTIAL EXPRESSION OVERLAP ANALYSIS ===\n")
  
  # Extract DE genes for each line at TP4 vs TP0
  de_genes <- list(
    Broad_462A = rownames(all_results_treat$B_462A_4vs0)[all_results_treat$B_462A_4vs0$FDR < 0.05],
    Broad_LD11 = rownames(all_results_treat$B_LD11_4vs0)[all_results_treat$B_LD11_4vs0$FDR < 0.05],
    Narrow_713B = rownames(all_results_treat$N_713B_4vs0)[all_results_treat$N_713B_4vs0$FDR < 0.05],
    Narrow_745 = rownames(all_results_treat$N_745_4vs0)[all_results_treat$N_745_4vs0$FDR < 0.05]
  )
  
  # Calculate overlaps
  broad_overlap <- length(intersect(de_genes$Broad_462A, de_genes$Broad_LD11))
  narrow_overlap <- length(intersect(de_genes$Narrow_713B, de_genes$Narrow_745))
  
  # Calculate Jaccard indices
  jaccard_broad <- broad_overlap / 
    length(union(de_genes$Broad_462A, de_genes$Broad_LD11))
  jaccard_narrow <- narrow_overlap / 
    length(union(de_genes$Narrow_713B, de_genes$Narrow_745))
  
  # Cross-leaf type overlaps
  cross_overlaps <- matrix(0, nrow = 4, ncol = 4)
  rownames(cross_overlaps) <- names(de_genes)
  colnames(cross_overlaps) <- names(de_genes)
  
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j) {
        overlap <- length(intersect(de_genes[[i]], de_genes[[j]]))
        union_size <- length(union(de_genes[[i]], de_genes[[j]]))
        cross_overlaps[i, j] <- overlap / union_size
      }
    }
  }
  
  # Create Venn diagrams
  pdf("results/leaf_type_validation/figures/05_DE_overlap.pdf", width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  # Broad lines
  venn.plot <- venn.diagram(
    list(PI532462A = de_genes$Broad_462A, 
         LD112170 = de_genes$Broad_LD11),
    filename = NULL,
    category.names = c("PI532462A", "LD112170"),
    fill = c("darkgreen", "lightgreen"),
    main = "Broad Lines DE Genes"
  )
  grid.draw(venn.plot)
  
  # Narrow lines
  venn.plot2 <- venn.diagram(
    list(PI612713B = de_genes$Narrow_713B, 
         PI547745 = de_genes$Narrow_745),
    filename = NULL,
    category.names = c("PI612713B", "PI547745"),
    fill = c("purple", "pink"),
    main = "Narrow Lines DE Genes"
  )
  grid.newpage()
  grid.draw(venn.plot2)
  
  dev.off()
  
  # Heatmap of Jaccard indices
  pdf("results/leaf_type_validation/figures/06_jaccard_heatmap.pdf", width = 8, height = 6)
  
  pheatmap(cross_overlaps, 
           main = "DE Gene Overlap (Jaccard Index)",
           color = colorRampPalette(c("white", "red"))(100),
           display_numbers = TRUE,
           number_format = "%.3f",
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  
  dev.off()
  
  cat("\nDE overlap summary:\n")
  cat("Broad lines overlap:", broad_overlap, "genes\n")
  cat("Broad lines Jaccard index:", round(jaccard_broad, 3), "\n")
  cat("Narrow lines overlap:", narrow_overlap, "genes\n")
  cat("Narrow lines Jaccard index:", round(jaccard_narrow, 3), "\n")
  
  # Test if within-leaf type overlap is greater than between
  within_leaf_jaccard <- c(jaccard_broad, jaccard_narrow)
  between_leaf_jaccard <- c(cross_overlaps[1, 3], cross_overlaps[1, 4],
                            cross_overlaps[2, 3], cross_overlaps[2, 4])
  
  wilcox_result <- wilcox.test(within_leaf_jaccard, between_leaf_jaccard)
  cat("\nWilcoxon test (within vs between leaf type Jaccard):\n")
  cat("P-value:", wilcox_result$p.value, "\n")
}

#=============================================================
# SECTION 7: FOLD CHANGE CONCORDANCE
#=============================================================

if (de_data_available) {
  cat("\n=== FOLD CHANGE CONCORDANCE ANALYSIS ===\n")
  
  # Get fold changes for common genes
  common_broad <- intersect(rownames(all_results_treat$B_462A_4vs0),
                            rownames(all_results_treat$B_LD11_4vs0))
  common_narrow <- intersect(rownames(all_results_treat$N_713B_4vs0),
                             rownames(all_results_treat$N_745_4vs0))
  
  # Extract fold changes
  fc_broad <- data.frame(
    Gene = common_broad,
    FC_462A = all_results_treat$B_462A_4vs0[common_broad, "logFC"],
    FC_LD11 = all_results_treat$B_LD11_4vs0[common_broad, "logFC"]
  )
  
  fc_narrow <- data.frame(
    Gene = common_narrow,
    FC_713B = all_results_treat$N_713B_4vs0[common_narrow, "logFC"],
    FC_745 = all_results_treat$N_745_4vs0[common_narrow, "logFC"]
  )
  
  # Calculate correlations
  cor_broad <- cor(fc_broad$FC_462A, fc_broad$FC_LD11, use = "complete.obs")
  cor_narrow <- cor(fc_narrow$FC_713B, fc_narrow$FC_745, use = "complete.obs")
  
  # Create scatter plots
  pdf("results/leaf_type_validation/figures/07_fold_change_concordance.pdf", width = 12, height = 6)
  par(mfrow = c(1, 2))
  
  # Broad lines
  plot(fc_broad$FC_462A, fc_broad$FC_LD11,
       pch = 19, col = rgb(0, 0.5, 0, 0.3),
       xlab = "log2FC PI532462A",
       ylab = "log2FC LD112170",
       main = paste("Broad Lines (r =", round(cor_broad, 3), ")"))
  abline(0, 1, col = "red", lty = 2)
  abline(h = 0, v = 0, col = "gray", lty = 3)
  
  # Narrow lines
  plot(fc_narrow$FC_713B, fc_narrow$FC_745,
       pch = 19, col = rgb(0.5, 0, 0.5, 0.3),
       xlab = "log2FC PI612713B",
       ylab = "log2FC PI547745",
       main = paste("Narrow Lines (r =", round(cor_narrow, 3), ")"))
  abline(0, 1, col = "red", lty = 2)
  abline(h = 0, v = 0, col = "gray", lty = 3)
  
  dev.off()
  
  # Calculate concordance (same direction)
  broad_concordance <- sum(sign(fc_broad$FC_462A) == sign(fc_broad$FC_LD11), na.rm = TRUE) / 
    sum(!is.na(fc_broad$FC_462A) & !is.na(fc_broad$FC_LD11))
  
  narrow_concordance <- sum(sign(fc_narrow$FC_713B) == sign(fc_narrow$FC_745), na.rm = TRUE) / 
    sum(!is.na(fc_narrow$FC_713B) & !is.na(fc_narrow$FC_745))
  
  cat("\nFold change concordance:\n")
  cat("Broad lines correlation:", round(cor_broad, 3), "\n")
  cat("Broad lines concordance:", round(broad_concordance * 100, 1), "%\n")
  cat("Narrow lines correlation:", round(cor_narrow, 3), "\n")
  cat("Narrow lines concordance:", round(narrow_concordance * 100, 1), "%\n")
}

#=============================================================
# SECTION 8: MODULE PRESERVATION (if WGCNA results available)
#=============================================================

if (file.exists("results/checkpoints/09_wgcna_analysis.RData")) {
  cat("\n=== MODULE PRESERVATION ANALYSIS ===\n")
  load("results/checkpoints/09_wgcna_analysis.RData")
  
  # This section would test if modules identified in one line
  # are preserved in the other line of the same leaf type
  # (Requires WGCNA modulePreservation function)
  
  cat("Module preservation analysis would go here...\n")
  cat("(Requires separate WGCNA runs per line)\n")
}

#=============================================================
# SECTION 9: SUMMARY REPORT
#=============================================================

cat("\n=== GENERATING SUMMARY REPORT ===\n")

# Compile all metrics
validation_metrics <- list(
  correlation_analysis = list(
    mean_within_leaf_type = mean(same_leaf_cors),
    mean_between_leaf_type = mean(diff_leaf_cors),
    correlation_difference = mean(same_leaf_cors) - mean(diff_leaf_cors),
    t_test_pvalue = t_test_result$p.value
  ),
  
  pca_analysis = list(
    pc1_variance_explained = var_explained[1],
    pc1_leaftype_rsquared = summary(pc1_leaf_only)$r.squared
  ),
  
  variance_partitioning = if (exists("var_summary") && !all(is.na(var_summary))) {
    list(
      avg_leaftype_variance = var_summary["Var_LeafType"],
      avg_line_within_leaf_variance = var_summary["Var_Line_within_Leaf"],
      avg_timepoint_variance = var_summary["Var_Timepoint"]
    )
  } else {
    list(note = "Variance partitioning failed due to nested design")
  }
)

if (de_data_available) {
  validation_metrics$de_overlap <- list(
    broad_jaccard = jaccard_broad,
    narrow_jaccard = jaccard_narrow,
    fold_change_cor_broad = cor_broad,
    fold_change_cor_narrow = cor_narrow,
    concordance_broad = broad_concordance,
    concordance_narrow = narrow_concordance
  )
}

# Create summary plot
pdf("results/leaf_type_validation/figures/08_validation_summary.pdf", width = 12, height = 10)

# Layout for multiple plots
layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))

# 1. Correlation comparison
cor_data <- data.frame(
  Category = c("Within\nLeaf Type", "Between\nLeaf Type"),
  Mean = c(mean(same_leaf_cors), mean(diff_leaf_cors)),
  SE = c(sd(same_leaf_cors)/sqrt(length(same_leaf_cors)),
         sd(diff_leaf_cors)/sqrt(length(diff_leaf_cors)))
)

barplot(cor_data$Mean, names.arg = cor_data$Category,
        ylim = c(0.8, 1),
        col = c("darkgreen", "darkred"),
        main = "A. Sample Correlations",
        ylab = "Mean Correlation")
arrows(x0 = c(0.7, 1.9), 
       y0 = cor_data$Mean - cor_data$SE,
       y1 = cor_data$Mean + cor_data$SE,
       angle = 90, code = 3, length = 0.1)

# 2. PCA variance
barplot(c(var_explained[1], sum(var_explained[2:3])),
        names.arg = c("PC1", "PC2+3"),
        col = c("darkblue", "lightblue"),
        main = "B. PCA Variance Explained",
        ylab = "Variance (%)")

# 3. Variance components
if (exists("var_summary") && !all(is.na(var_summary))) {
  barplot(var_summary * 100,
          names.arg = c("Leaf\nType", "Line", "Time", "Resid"),
          col = c("darkgreen", "lightblue", "orange", "gray"),
          main = "C. Variance Components",
          ylab = "% Variance")
} else {
  plot.new()
  text(0.5, 0.5, "Variance components\nnot available", cex = 1.2)
}

if (de_data_available) {
  # 4. DE overlap
  barplot(c(jaccard_broad, jaccard_narrow) * 100,
          names.arg = c("Broad\nLines", "Narrow\nLines"),
          col = c("darkgreen", "purple"),
          main = "D. DE Gene Overlap",
          ylab = "Jaccard Index (%)",
          ylim = c(0, 100))
  
  # 5. Fold change correlation
  barplot(c(cor_broad, cor_narrow),
          names.arg = c("Broad\nLines", "Narrow\nLines"),
          col = c("darkgreen", "purple"),
          main = "E. Fold Change Correlation",
          ylab = "Pearson r",
          ylim = c(0, 1))
  
  # 6. Concordance
  barplot(c(broad_concordance, narrow_concordance) * 100,
          names.arg = c("Broad\nLines", "Narrow\nLines"),
          col = c("darkgreen", "purple"),
          main = "F. Direction Concordance",
          ylab = "% Same Direction",
          ylim = c(0, 100))
}

dev.off()

# Save summary metrics
saveRDS(validation_metrics, file = "results/leaf_type_validation/validation_metrics.rds")

# Create text summary
summary_text <- paste0(
  "LEAF TYPE GROUPING VALIDATION SUMMARY\n",
  "=====================================\n\n",
  
  "1. TRANSCRIPTOME-WIDE SIMILARITY:\n",
  "   - Within leaf type correlation: ", round(mean(same_leaf_cors), 3), "\n",
  "   - Between leaf type correlation: ", round(mean(diff_leaf_cors), 3), "\n",
  "   - Difference: ", round(mean(same_leaf_cors) - mean(diff_leaf_cors), 3), 
  " (p = ", format.pval(t_test_result$p.value, digits = 3), ")\n\n",
  
  "2. PRINCIPAL COMPONENT ANALYSIS:\n",
  "   - PC1 explains ", round(var_explained[1], 1), "% of variance\n",
  "   - Leaf type explains ", round(summary(pc1_leaf_only)$r.squared * 100, 1), 
  "% of PC1 variation\n\n",
  
  "3. VARIANCE PARTITIONING:\n"
)

if (exists("var_summary") && !all(is.na(var_summary))) {
  summary_text <- paste0(summary_text,
    "   - Leaf type: ", round(var_summary["Var_LeafType"] * 100, 1), "%\n",
    "   - Line within Leaf: ", round(var_summary["Var_Line_within_Leaf"] * 100, 1), "%\n",
    "   - Timepoint: ", round(var_summary["Var_Timepoint"] * 100, 1), "%\n\n"
  )
} else {
  summary_text <- paste0(summary_text,
    "   - Could not calculate due to nested design\n\n"
  )
}

if (de_data_available) {
  summary_text <- paste0(summary_text,
    "4. DIFFERENTIAL EXPRESSION OVERLAP:\n",
    "   - Broad lines Jaccard index: ", round(jaccard_broad, 3), "\n",
    "   - Narrow lines Jaccard index: ", round(jaccard_narrow, 3), "\n",
    "   - Broad lines FC correlation: ", round(cor_broad, 3), "\n",
    "   - Narrow lines FC correlation: ", round(cor_narrow, 3), "\n\n"
  )
}

summary_text <- paste0(summary_text,
  "CONCLUSION:\n",
  if (mean(same_leaf_cors) - mean(diff_leaf_cors) > 0.05 && t_test_result$p.value < 0.05) {
    "The data strongly support grouping lines by leaf type. Lines of the same leaf type show significantly higher similarity than lines of different leaf types."
  } else if (mean(same_leaf_cors) - mean(diff_leaf_cors) > 0.02) {
    "The data moderately support grouping by leaf type, though line-specific effects are also present."
  } else {
    "The data suggest significant line-specific effects. Grouping by leaf type should be done with caution."
  }
)

# Add interpretation based on the actual results
summary_text <- paste0(summary_text, "\n\n",
  "DETAILED INTERPRETATION:\n",
  "Based on the correlation analysis showing minimal difference (0.003) between\n",
  "within-leaf-type and between-leaf-type correlations, and the PCA showing\n",
  "that leaf type explains only 1.16% of PC1 variation, the evidence for\n",
  "grouping by leaf type is WEAK. The dominant factors are timepoint\n",
  "(explaining most of the variance) and line-specific effects.\n\n",
  
  "RECOMMENDATION:\n",
  "1. Primary analyses should consider line-specific effects\n",
  "2. Leaf type grouping should be used cautiously and with acknowledgment\n",
  "   of substantial line-specific variation\n",
  "3. Consider using mixed models that account for line nested within leaf type\n",
  "4. Validate key findings in individual lines before generalizing to leaf types\n"
)

writeLines(summary_text, "results/leaf_type_validation/validation_summary.txt")

cat("\n", summary_text)

cat("\n\n=== VALIDATION ANALYSIS COMPLETE ===\n")
cat("Results saved to: results/leaf_type_validation/\n")
cat("Key files:\n")
cat("  - validation_summary.txt: Overall conclusions\n")
cat("  - validation_metrics.rds: All computed metrics\n")
cat("  - figures/: All validation plots\n")

# Session info
sink("results/leaf_type_validation/session_info.txt")
sessionInfo()
sink()
