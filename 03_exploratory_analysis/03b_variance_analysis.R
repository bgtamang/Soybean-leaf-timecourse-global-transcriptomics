# Script 08: Variance Analysis

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 08: VARIANCE ANALYSIS\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/08_variance", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "RColorBrewer",
  "vegan",       # For PERMANOVA
  "dplyr",
  "tidyr"
)

missing <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing) > 0) {
  install.packages(missing)
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded checkpoint from script 06\n\n")

# Use primary dataset
expr_mat <- v_primary$E
targets <- targets_primary

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n\n")

# ===== DEFINE COLORS =====

factor_colors <- c(
  "Leaf_type" = "#2E8B57",
  "Line" = "#D95F02",
  "Timepoint" = "#E7298A",
  "Batch" = "#377EB8",
  "Residual" = "#999999"
)

# ===== SECTION 2: VARIANCE PARTITIONING (Per-Gene) =====

cat("========================================\n")
cat("SECTION 2: PER-GENE VARIANCE PARTITIONING\n")
cat("========================================\n\n")

cat("Calculating variance explained by each factor for all genes...\n")
cat("  This may take a minute...\n\n")

# Function to calculate R-squared for each factor
calc_variance_per_gene <- function(expr_vec, metadata) {
  # Total variance
  total_var <- var(expr_vec)

  # Leaf type only
  fit_leaf <- lm(expr_vec ~ Leaf_type, data = metadata)
  rsq_leaf <- summary(fit_leaf)$r.squared

  # Line only
  fit_line <- lm(expr_vec ~ Line, data = metadata)
  rsq_line <- summary(fit_line)$r.squared

  # Timepoint only
  fit_tp <- lm(expr_vec ~ Timepoint, data = metadata)
  rsq_tp <- summary(fit_tp)$r.squared

  # Batch only
  fit_batch <- lm(expr_vec ~ Batch, data = metadata)
  rsq_batch <- summary(fit_batch)$r.squared

  # Full model (all factors)
  fit_full <- lm(expr_vec ~ Leaf_type + Line + Timepoint + Batch, data = metadata)
  rsq_full <- summary(fit_full)$r.squared

  return(c(
    Leaf_type = rsq_leaf,
    Line = rsq_line,
    Timepoint = rsq_tp,
    Batch = rsq_batch,
    Full_model = rsq_full,
    Residual = 1 - rsq_full
  ))
}

# Apply to all genes (sample for speed if needed)
n_genes <- nrow(expr_mat)
cat("Analyzing", n_genes, "genes...\n")

variance_results <- t(apply(expr_mat, 1, calc_variance_per_gene, metadata = targets))
variance_df <- as.data.frame(variance_results)
variance_df$Gene <- rownames(variance_df)

cat("  Complete!\n\n")

# Summary statistics
cat("Variance Explained Summary (mean across all genes):\n")
cat("  Leaf_type:", round(mean(variance_df$Leaf_type) * 100, 2), "%\n")
cat("  Line:", round(mean(variance_df$Line) * 100, 2), "%\n")
cat("  Timepoint:", round(mean(variance_df$Timepoint) * 100, 2), "%\n")
cat("  Batch:", round(mean(variance_df$Batch) * 100, 2), "%\n")
cat("  Full model:", round(mean(variance_df$Full_model) * 100, 2), "%\n")
cat("  Residual:", round(mean(variance_df$Residual) * 100, 2), "%\n\n")

# ===== SECTION 3: VISUALIZE VARIANCE COMPONENTS =====

cat("========================================\n")
cat("SECTION 3: VARIANCE VISUALIZATIONS\n")
cat("========================================\n\n")

# --- Plot 1: Bar plot of mean variance explained ---
cat("Creating variance bar plot...\n")

var_summary <- data.frame(
  Factor = c("Leaf_type", "Line", "Timepoint", "Batch", "Residual"),
  Mean_Variance = c(
    mean(variance_df$Leaf_type),
    mean(variance_df$Line),
    mean(variance_df$Timepoint),
    mean(variance_df$Batch),
    mean(variance_df$Residual)
  ) * 100
)
var_summary$Factor <- factor(var_summary$Factor,
                              levels = c("Timepoint", "Line", "Batch", "Leaf_type", "Residual"))

p1 <- ggplot(var_summary, aes(x = Factor, y = Mean_Variance, fill = Factor)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = factor_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Mean Variance Explained by Experimental Factors",
       subtitle = paste("Across", n_genes, "genes"),
       x = "Factor",
       y = "Mean % Variance Explained") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = paste0(round(Mean_Variance, 1), "%")),
            vjust = -0.5, size = 4)

ggsave("03_results/figures/08_variance/variance_barplot.png", p1,
       width = 8, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/08_variance/variance_barplot.png\n")

# --- Plot 2: Distribution of variance explained ---
cat("Creating variance distribution plot...\n")

var_long <- variance_df %>%
  select(Gene, Leaf_type, Line, Timepoint, Batch) %>%
  pivot_longer(cols = -Gene, names_to = "Factor", values_to = "Variance") %>%
  mutate(Variance_Pct = Variance * 100)

var_long$Factor <- factor(var_long$Factor,
                          levels = c("Timepoint", "Line", "Batch", "Leaf_type"))

p2 <- ggplot(var_long, aes(x = Factor, y = Variance_Pct, fill = Factor)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = factor_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Distribution of Variance Explained",
       subtitle = "Per-gene variance partitioning",
       x = "Factor",
       y = "% Variance Explained") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none")

ggsave("03_results/figures/08_variance/variance_distribution.png", p2,
       width = 8, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/08_variance/variance_distribution.png\n")

# --- Plot 3: Pie chart of average variance ---
cat("Creating variance pie chart...\n")

png("03_results/figures/08_variance/variance_pie.png",
    width = 8, height = 8, units = "in", res = 300)

pie_data <- var_summary$Mean_Variance
names(pie_data) <- paste0(var_summary$Factor, "\n(", round(pie_data, 1), "%)")

pie(pie_data,
    col = factor_colors[as.character(var_summary$Factor)],
    main = "Average Variance Partition",
    cex.main = 1.3)

dev.off()
cat("  Saved: 03_results/figures/08_variance/variance_pie.png\n")

# ===== SECTION 4: PERMANOVA =====

cat("\n========================================\n")
cat("SECTION 4: PERMANOVA ANALYSIS\n")
cat("========================================\n\n")

cat("Running PERMANOVA to test significance of factors...\n")

# Calculate distance matrix
sample_dist <- vegdist(t(expr_mat), method = "euclidean")

# PERMANOVA with all factors
set.seed(42)
permanova_full <- adonis2(sample_dist ~ Leaf_type + Line + Timepoint + Batch,
                          data = targets,
                          permutations = 999)

cat("\nPERMANOVA Results (Full Model):\n")
print(permanova_full)

# Individual factor tests
cat("\nIndividual Factor Tests:\n")

permanova_leaf <- adonis2(sample_dist ~ Leaf_type, data = targets, permutations = 999)
cat("\nLeaf_type alone:\n")
print(permanova_leaf)

permanova_line <- adonis2(sample_dist ~ Line, data = targets, permutations = 999)
cat("\nLine alone:\n")
print(permanova_line)

permanova_tp <- adonis2(sample_dist ~ Timepoint, data = targets, permutations = 999)
cat("\nTimepoint alone:\n")
print(permanova_tp)

permanova_batch <- adonis2(sample_dist ~ Batch, data = targets, permutations = 999)
cat("\nBatch alone:\n")
print(permanova_batch)

# Save PERMANOVA results
permanova_summary <- data.frame(
  Factor = c("Leaf_type", "Line", "Timepoint", "Batch"),
  R2_Individual = c(
    permanova_leaf$R2[1],
    permanova_line$R2[1],
    permanova_tp$R2[1],
    permanova_batch$R2[1]
  ),
  P_value = c(
    permanova_leaf$`Pr(>F)`[1],
    permanova_line$`Pr(>F)`[1],
    permanova_tp$`Pr(>F)`[1],
    permanova_batch$`Pr(>F)`[1]
  ),
  R2_FullModel = permanova_full$R2[1:4],
  P_value_Full = permanova_full$`Pr(>F)`[1:4]
)

write.csv(permanova_summary, "03_results/tables/PERMANOVA_results.csv", row.names = FALSE)
cat("\nSaved: 03_results/tables/PERMANOVA_results.csv\n")

# ===== SECTION 5: LEAF-TYPE VALIDATION =====

cat("\n========================================\n")
cat("SECTION 5: LEAF-TYPE GROUPING VALIDATION\n")
cat("========================================\n\n")

cat("Evaluating whether grouping by Leaf_type is appropriate:\n\n")

leaf_var <- mean(variance_df$Leaf_type) * 100
line_var <- mean(variance_df$Line) * 100

cat("Mean variance explained:\n")
cat("  Leaf_type:", round(leaf_var, 2), "%\n")
cat("  Line:", round(line_var, 2), "%\n")
cat("  Ratio (Line/Leaf_type):", round(line_var/leaf_var, 2), "\n\n")

# PERMANOVA significance
cat("PERMANOVA p-values:\n")
cat("  Leaf_type:", permanova_leaf$`Pr(>F)`[1], "\n")
cat("  Line:", permanova_line$`Pr(>F)`[1], "\n\n")

# Assessment
cat("ASSESSMENT:\n")
if (leaf_var < 5) {
  cat("  WARNING: Leaf_type explains very little variance (<5%)\n")
  cat("  Line-specific effects may dominate\n")
  cat("  Consider line-specific analysis in addition to grouped analysis\n")
} else if (line_var > leaf_var * 2) {
  cat("  CAUTION: Line explains >2x more variance than Leaf_type\n")
  cat("  Consider both grouped and line-specific analyses\n")
} else {
  cat("  Leaf_type grouping appears reasonable\n")
}

# --- Plot 4: Leaf vs Line variance comparison ---
cat("\nCreating Leaf vs Line comparison plot...\n")

var_comparison <- variance_df %>%
  mutate(Leaf_pct = Leaf_type * 100,
         Line_pct = Line * 100,
         Higher = ifelse(Leaf_pct > Line_pct, "Leaf_type", "Line"))

p4 <- ggplot(var_comparison, aes(x = Line_pct, y = Leaf_pct, color = Higher)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Leaf_type" = "#2E8B57", "Line" = "#D95F02")) +
  theme_bw(base_size = 14) +
  labs(title = "Leaf_type vs Line Variance per Gene",
       subtitle = paste("Points above line: Leaf_type > Line (",
                        sum(var_comparison$Higher == "Leaf_type"), "genes)"),
       x = "% Variance Explained by Line",
       y = "% Variance Explained by Leaf_type") +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/08_variance/leaf_vs_line_variance.png", p4,
       width = 8, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/08_variance/leaf_vs_line_variance.png\n")

# ===== SECTION 6: SAVE RESULTS =====

cat("\n========================================\n")
cat("SECTION 6: SAVE RESULTS\n")
cat("========================================\n\n")

# Save variance components
write.csv(variance_df, "03_results/tables/variance_components.csv", row.names = FALSE)
cat("Saved: 03_results/tables/variance_components.csv\n")

# Save summary
variance_summary <- data.frame(
  Factor = c("Leaf_type", "Line", "Timepoint", "Batch", "Full_model", "Residual"),
  Mean_Variance_Pct = c(
    mean(variance_df$Leaf_type),
    mean(variance_df$Line),
    mean(variance_df$Timepoint),
    mean(variance_df$Batch),
    mean(variance_df$Full_model),
    mean(variance_df$Residual)
  ) * 100,
  Median_Variance_Pct = c(
    median(variance_df$Leaf_type),
    median(variance_df$Line),
    median(variance_df$Timepoint),
    median(variance_df$Batch),
    median(variance_df$Full_model),
    median(variance_df$Residual)
  ) * 100
)

write.csv(variance_summary, "03_results/tables/variance_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/variance_summary.csv\n\n")

print(variance_summary)

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 08: VARIANCE ANALYSIS - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Key Findings:\n")
cat("    - Timepoint variance:", round(mean(variance_df$Timepoint) * 100, 1), "%\n")
cat("    - Line variance:", round(mean(variance_df$Line) * 100, 1), "%\n")
cat("    - Leaf_type variance:", round(mean(variance_df$Leaf_type) * 100, 1), "%\n")
cat("    - Batch variance:", round(mean(variance_df$Batch) * 100, 1), "%\n")
cat("\n")
cat("  Next: Run 09_expression_overview.R\n")
cat("================================================================\n")
