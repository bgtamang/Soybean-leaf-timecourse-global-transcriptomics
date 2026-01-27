
# Set base directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase2-Refined-Analysis"
setwd(base_dir)

library(tidyverse)
library(edgeR)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)

# 1. Load Required Data

# Load DE results by timepoint
load(file.path(base_dir, "03_results/checkpoints/12_DE_organized.RData"))

# Load JAG1 targets
jag1_targets <- read_csv(file.path(base_dir, "03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv"),
                         show_col_types = FALSE)

# Load normalized expression
load(file.path(base_dir, "03_results/checkpoints/05_normalized.RData"))

# Load binding data if available (updated to use new binding integration file)
binding_file <- file.path(base_dir, "03_results/tables/binding_integration/integrated_binding_analysis.csv")
if(file.exists(binding_file)) {
  binding_data <- read_csv(binding_file, show_col_types = FALSE)
  # Rename columns to match expected names
  if("Has_DAPseq" %in% names(binding_data)) {
    binding_data <- binding_data %>%
      rename(has_wang = Has_DAPseq, has_huang = Has_ChIPseq)
  }
  has_binding <- TRUE
} else {
  has_binding <- FALSE
}

cat("Loaded", nrow(jag1_targets), "JAG1 targets for temporal classification\n\n")

# 2. Extract DE Status by Timepoint


# Get DE status for each timepoint comparison (Narrow vs Broad)
# Need to check what contrasts are available in DE results

# List DE result files (individual comparison files, not comprehensive)
de_files <- list.files(file.path(base_dir, "03_results/tables/DE/"),
                       pattern = "^[A-Z].*\\.csv$", full.names = TRUE)
# Exclude comprehensive files
de_files <- de_files[!grepl("comprehensive|top100", basename(de_files))]
cat("Found", length(de_files), "DE result files\n")

# Read and combine DE results
de_combined <- data.frame()

for(f in de_files) {
  # Extract comparison info from filename
  fname <- basename(f)

  # Read file
  de_tmp <- read_csv(f, show_col_types = FALSE)

  # Check for Gene or GeneID column and standardize
  if(nrow(de_tmp) > 0) {
    # Rename Gene to gene_id if needed
    if("Gene" %in% names(de_tmp) && !"gene_id" %in% names(de_tmp)) {
      de_tmp <- de_tmp %>% rename(gene_id = Gene)
    }
    if("GeneID" %in% names(de_tmp) && !"gene_id" %in% names(de_tmp)) {
      de_tmp <- de_tmp %>% rename(gene_id = GeneID)
    }
    if("gene_id" %in% names(de_tmp) && "logFC" %in% names(de_tmp) && "FDR" %in% names(de_tmp)) {
      de_tmp$source_file <- fname
      de_combined <- bind_rows(de_combined, de_tmp %>% select(gene_id, logFC, FDR, source_file))
    }
  }
}

cat("Combined DE results:", nrow(de_combined), "rows\n")
cat("Files used:\n")
print(unique(de_combined$source_file))

# 3. Define Temporal Categories

# Define timepoint groups
early_timepoints <- c("TP0", "TP1")  # SAM and primordia - JAG1 most active
mid_timepoints <- c("TP2")           # Young leaf - transitional
late_timepoints <- c("TP3", "TP4")   # Expanding/mature leaf - consequences

# Function to classify temporal pattern
classify_temporal <- function(de_data, gene_list) {

  # Initialize classification dataframe
  classification <- data.frame(
    gene_id = gene_list,
    DE_TP0 = FALSE,
    DE_TP1 = FALSE,
    DE_TP2 = FALSE,
    DE_TP3 = FALSE,
    DE_TP4 = FALSE,
    stringsAsFactors = FALSE
  )

  # Check each timepoint
  for(tp in c("TP0", "TP1", "TP2", "TP3", "TP4")) {
    # Find DE genes at this timepoint
    # Looking for Narrow vs Broad comparisons at each timepoint
    tp_de <- de_data %>%
      filter(grepl(tp, source_file, ignore.case = TRUE)) %>%
      filter(abs(logFC) > 1 & FDR < 0.05) %>%
      pull(gene_id) %>%
      unique()

    col_name <- paste0("DE_", tp)
    classification[[col_name]] <- classification$gene_id %in% tp_de
  }

  # Calculate derived categories
  classification <- classification %>%
    mutate(
      # Count timepoints with DE
      n_timepoints_DE = DE_TP0 + DE_TP1 + DE_TP2 + DE_TP3 + DE_TP4,

      # Early DE (TP0 or TP1)
      early_DE = DE_TP0 | DE_TP1,

      # Late only (TP3 or TP4 but not early)
      late_only = (DE_TP3 | DE_TP4) & !early_DE,

      # Persistent (early AND late)
      persistent = early_DE & (DE_TP3 | DE_TP4),

      # Transient (only middle timepoints)
      transient = DE_TP2 & !DE_TP0 & !DE_TP1 & !DE_TP3 & !DE_TP4,

      # Classify
      temporal_category = case_when(
        early_DE & persistent ~ "A_Persistent_Direct",    # Direct + maintained
        early_DE & !persistent ~ "B_Early_Transient",     # Direct, resolves
        late_only ~ "C_Late_Indirect",                    # Indirect effects
        transient ~ "D_Mid_Transient",                    # Transitional
        n_timepoints_DE == 0 ~ "E_Not_DE_Individual",     # Only DE in pooled
        TRUE ~ "F_Other"
      )
    )

  return(classification)
}

# Apply classification to JAG1 targets
temporal_class <- classify_temporal(de_combined, jag1_targets$GeneID)

# Merge with JAG1 target info
temporal_class <- temporal_class %>%
  left_join(jag1_targets %>% select(GeneID, Confidence_Tier, Mean_logFC_Pairwise, Pattern, N_Pairwise_UP),
            by = c("gene_id" = "GeneID"))

cat("\n=== Temporal Classification Results ===\n")
print(table(temporal_class$temporal_category))

# 4. Direct Target Evidence Scoring

cat("\n=== Direct Target Evidence Scoring ===\n")

# Add binding information if available
if(has_binding) {
  # Check column names in binding data
  if("GeneID" %in% names(binding_data) && !"gene_id" %in% names(binding_data)) {
    binding_data <- binding_data %>% rename(gene_id = GeneID)
  }

  temporal_class <- temporal_class %>%
    left_join(binding_data %>% select(any_of(c("gene_id", "has_wang", "has_huang", "binding_source"))),
              by = "gene_id")

  temporal_class <- temporal_class %>%
    mutate(
      has_wang = ifelse(is.na(has_wang), FALSE, has_wang),
      has_huang = ifelse(is.na(has_huang), FALSE, has_huang)
    )
} else {
  temporal_class$has_wang <- FALSE
  temporal_class$has_huang <- FALSE
}

# Calculate direct target score
temporal_class <- temporal_class %>%
  mutate(
    # Direct target evidence
    direct_score = 0,

    # Add points for early DE
    direct_score = direct_score + ifelse(DE_TP0, 3, 0),
    direct_score = direct_score + ifelse(DE_TP1, 2, 0),

    # Add points for tier
    direct_score = direct_score + case_when(
      Confidence_Tier == "Gold" ~ 3,
      Confidence_Tier == "Silver" ~ 2,
      Confidence_Tier == "Bronze" ~ 1,
      TRUE ~ 0
    ),

    # Add points for binding (if available)
    direct_score = direct_score + ifelse(has_wang, 2, 0),
    direct_score = direct_score + ifelse(has_huang, 2, 0),

    # Classify by direct evidence
    direct_confidence = case_when(
      direct_score >= 8 ~ "Very_High",
      direct_score >= 5 ~ "High",
      direct_score >= 3 ~ "Moderate",
      TRUE ~ "Low"
    )
  )

cat("\nDirect Target Confidence Distribution:\n")
print(table(temporal_class$direct_confidence))

# 5. Top Direct Target Candidates

cat("\n=== Top Direct Target Candidates ===\n")

top_direct <- temporal_class %>%
  filter(early_DE == TRUE) %>%
  arrange(desc(direct_score), desc(abs(Mean_logFC_Pairwise))) %>%
  head(50)

cat("\nTop 20 Direct Target Candidates:\n")
print(top_direct %>%
        select(gene_id, Confidence_Tier, temporal_category, direct_score,
               direct_confidence, Mean_logFC_Pairwise, n_timepoints_DE) %>%
        head(20))

# 6. Create Visualizations

cat("\n=== Creating Visualizations ===\n")

# Output directory
dir.create(file.path(base_dir, "03_results/figures/35_temporal"), showWarnings = FALSE, recursive = TRUE)

# 6a. UpSet plot for temporal patterns
upset_data <- temporal_class %>%
  select(gene_id, DE_TP0, DE_TP1, DE_TP2, DE_TP3, DE_TP4) %>%
  mutate(across(starts_with("DE_"), as.integer))

pdf(file.path(base_dir, "03_results/figures/35_temporal/temporal_upset_plot.pdf"), width = 10, height = 6)
upset(as.data.frame(upset_data),
      sets = c("DE_TP0", "DE_TP1", "DE_TP2", "DE_TP3", "DE_TP4"),
      order.by = "freq",
      mainbar.y.label = "Number of Genes",
      sets.x.label = "Genes DE at Timepoint",
      main.bar.color = "#4292C6",
      sets.bar.color = "#2171B5")
dev.off()
cat("Saved upset plot\n")

# 6b. Bar chart of temporal categories
category_counts <- temporal_class %>%
  count(temporal_category) %>%
  mutate(temporal_category = factor(temporal_category,
                                     levels = c("A_Persistent_Direct", "B_Early_Transient",
                                               "C_Late_Indirect", "D_Mid_Transient",
                                               "E_Not_DE_Individual", "F_Other")))

p_categories <- ggplot(category_counts, aes(x = temporal_category, y = n, fill = temporal_category)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Temporal Classification of JAG1 Targets",
       subtitle = "Categories A-B are likely direct targets; C-D are indirect effects",
       x = "Temporal Category",
       y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_text(aes(label = n), vjust = -0.5)

ggsave(file.path(base_dir, "03_results/figures/35_temporal/temporal_category_barplot.pdf"),
       p_categories, width = 10, height = 6)
cat("Saved category bar plot\n")

# 6c. Direct confidence by tier
tier_confidence <- temporal_class %>%
  filter(!is.na(Confidence_Tier)) %>%
  count(Confidence_Tier, direct_confidence) %>%
  mutate(Confidence_Tier = factor(Confidence_Tier, levels = c("Gold", "Silver", "Bronze")),
         direct_confidence = factor(direct_confidence,
                                    levels = c("Very_High", "High", "Moderate", "Low")))

p_tier_conf <- ggplot(tier_confidence, aes(x = Confidence_Tier, y = n, fill = direct_confidence)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Very_High" = "#238B45", "High" = "#74C476",
                               "Moderate" = "#FDD49E", "Low" = "#EF6548")) +
  theme_bw() +
  labs(title = "Direct Target Evidence by DE Tier",
       x = "DE Tier",
       y = "Proportion",
       fill = "Direct\nConfidence") +
  scale_y_continuous(labels = scales::percent)

ggsave(file.path(base_dir, "03_results/figures/35_temporal/tier_direct_confidence.pdf"),
       p_tier_conf, width = 8, height = 5)
cat("Saved tier confidence plot\n")

# 6d. Heatmap of early vs late expression
# Get expression for a subset of interesting genes
early_genes <- temporal_class %>%
  filter(temporal_category == "A_Persistent_Direct") %>%
  arrange(desc(direct_score)) %>%
  head(30) %>%
  pull(gene_id)

late_genes <- temporal_class %>%
  filter(temporal_category == "C_Late_Indirect") %>%
  arrange(desc(abs(Mean_logFC_Pairwise))) %>%
  head(30) %>%
  pull(gene_id)

genes_for_heatmap <- c(early_genes, late_genes)
genes_for_heatmap <- genes_for_heatmap[genes_for_heatmap %in% rownames(v_full$E)]

if(length(genes_for_heatmap) > 5) {
  expr_heatmap <- v_full$E[genes_for_heatmap, ]
  expr_scaled <- t(scale(t(expr_heatmap)))

  # Row annotation
  row_anno <- temporal_class %>%
    filter(gene_id %in% genes_for_heatmap) %>%
    select(gene_id, temporal_category, Confidence_Tier) %>%
    column_to_rownames("gene_id")

  # Column annotation
  col_anno <- targets %>%
    select(Leaf_type, Timepoint) %>%
    as.data.frame()
  rownames(col_anno) <- colnames(expr_scaled)

  anno_colors <- list(
    Leaf_type = c("Broad" = "#2166AC", "Narrow" = "#B2182B"),
    Timepoint = c("TP0" = "#FEE5D9", "TP1" = "#FCAE91", "TP2" = "#FB6A4A",
                  "TP3" = "#DE2D26", "TP4" = "#A50F15"),
    temporal_category = c("A_Persistent_Direct" = "#1B9E77", "C_Late_Indirect" = "#D95F02"),
    Confidence_Tier = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")
  )

  pdf(file.path(base_dir, "03_results/figures/35_temporal/early_vs_late_heatmap.pdf"), width = 14, height = 10)
  pheatmap(expr_scaled,
           annotation_row = row_anno,
           annotation_col = col_anno,
           annotation_colors = anno_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           main = "Early (Direct) vs Late (Indirect) JAG1 Targets")
  dev.off()
  cat("Saved heatmap\n")
}

# 7. Summary Statistics

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n=== TEMPORAL CLASSIFICATION SUMMARY ===\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n")

# Count by category
cat("\nTemporal Categories:\n")
cat("- A. Persistent Direct (early DE + maintained):", sum(temporal_class$temporal_category == "A_Persistent_Direct"), "\n")
cat("- B. Early Transient (early DE, resolves):", sum(temporal_class$temporal_category == "B_Early_Transient"), "\n")
cat("- C. Late Indirect (late DE only):", sum(temporal_class$temporal_category == "C_Late_Indirect"), "\n")
cat("- D. Mid Transient (TP2 only):", sum(temporal_class$temporal_category == "D_Mid_Transient"), "\n")
cat("- E. Not DE Individually:", sum(temporal_class$temporal_category == "E_Not_DE_Individual"), "\n")

# Direct target counts
cat("\nDirect Target Evidence:\n")
cat("- Very High confidence:", sum(temporal_class$direct_confidence == "Very_High"), "\n")
cat("- High confidence:", sum(temporal_class$direct_confidence == "High"), "\n")
cat("- Moderate confidence:", sum(temporal_class$direct_confidence == "Moderate"), "\n")
cat("- Low confidence:", sum(temporal_class$direct_confidence == "Low"), "\n")

# Overlap with Gold tier
gold_early <- temporal_class %>%
  filter(Confidence_Tier == "Gold" & early_DE == TRUE)

cat("\nGold Tier + Early DE (highest confidence direct targets):", nrow(gold_early), "\n")

# 8. Save Results

dir.create(file.path(base_dir, "03_results/tables/temporal"), showWarnings = FALSE, recursive = TRUE)

# Save full classification
write_csv(temporal_class, file.path(base_dir, "03_results/tables/temporal/temporal_classification_all.csv"))

# Save top direct targets
write_csv(top_direct, file.path(base_dir, "03_results/tables/temporal/top_direct_targets.csv"))

# Save category summaries
category_summary <- temporal_class %>%
  group_by(temporal_category, Confidence_Tier) %>%
  summarize(
    n_genes = n(),
    mean_logFC = mean(abs(Mean_logFC_Pairwise), na.rm = TRUE),
    mean_direct_score = mean(direct_score, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(category_summary, file.path(base_dir, "03_results/tables/temporal/temporal_category_summary.csv"))

# Save direct target list for publication
direct_targets_pub <- temporal_class %>%
  filter(direct_confidence %in% c("Very_High", "High")) %>%
  select(gene_id, Confidence_Tier, temporal_category, direct_score, direct_confidence,
         Mean_logFC_Pairwise, Pattern, n_timepoints_DE, DE_TP0, DE_TP1) %>%
  arrange(desc(direct_score))

write_csv(direct_targets_pub, file.path(base_dir, "03_results/tables/temporal/high_confidence_direct_targets.csv"))

cat("\nResults saved to 03_results/tables/temporal/\n")

# 9. Save Checkpoint

save(temporal_class, top_direct, category_summary,
     file = file.path(base_dir, "03_results/checkpoints/35_temporal_classification.RData"))

cat("\nCheckpoint saved: 35_temporal_classification.RData\n")
cat("\n=== Script 35 Complete ===\n")
