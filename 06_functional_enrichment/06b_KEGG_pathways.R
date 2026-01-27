# Classifies JAG1 targets based on Arabidopsis ortholog annotation keywords

rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 23: FUNCTIONAL CATEGORY ANALYSIS\n")
cat("  Keyword-Based Classification\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/23_functional_categories", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/tables/functional", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "RColorBrewer"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
# Keywords to search in Arabidopsis ortholog deflines
functional_categories <- list(

  "Transcription Factor" = c("transcription factor", "DNA binding", "DNA-binding",
                              "zinc finger", "bHLH", "MYB", "WRKY", "NAC", "AP2",
                              "ERF", "bZIP", "homeodomain", "MADS", "TCP"),

  "Protein Kinase" = c("kinase", "phosphatase", "receptor kinase", "MAP kinase",
                       "MAPK", "phosphorylation"),

  "Hormone Pathway" = c("auxin", "cytokinin", "gibberellin", "ethylene",
                        "abscisic", "brassinosteroid", "jasmonic", "salicylic",
                        "ABA", "GA", "IAA", "hormone"),

  "Cell Wall" = c("cell wall", "cellulose", "pectin", "xylan", "lignin",
                  "expansin", "xyloglucan", "extensin", "pectinesterase",
                  "glycosyl hydrolase"),

  "Transporter" = c("transporter", "channel", "permease", "pump", "ABC transporter",
                    "carrier", "antiporter", "symporter"),

  "Signaling" = c("signal", "receptor", "G protein", "calmodulin", "calcium",
                  "phospholipase", "GTPase"),

  "Cell Cycle" = c("cyclin", "cell cycle", "cell division", "CDK", "mitotic",
                   "proliferation"),

  "Metabolism" = c("synthase", "reductase", "oxidase", "transferase",
                   "dehydrogenase", "carboxylase", "hydrolase", "isomerase"),

  "Stress/Defense" = c("stress", "defense", "resistance", "pathogen", "disease",
                       "heat shock", "drought", "cold", "salt", "NBS-LRR", "PR protein"),

  "Development" = c("development", "morphogenesis", "differentiation", "meristem",
                    "flower", "leaf", "root", "embryo", "seed")
)

cat("Defined", length(functional_categories), "functional categories\n")
cat("Categories:", paste(names(functional_categories), collapse = ", "), "\n\n")
# Load JAG1 target checkpoint
load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n")

# Get targets
jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("JAG1 targets:", nrow(jag1_targets), "\n")

# Load annotation file for Arabidopsis deflines
annotation_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (file.exists(annotation_file)) {
  annot <- read.delim(annotation_file, header = TRUE, stringsAsFactors = FALSE,
                      comment.char = "", quote = "")
  cat("Loaded annotation file:", nrow(annot), "entries\n")

  # Check for Arabidopsis columns
  if ("Best.hit.arabi.defline" %in% colnames(annot)) {
    cat("Arabidopsis defline column found\n\n")
  } else {
    stop("Arabidopsis defline column not found!")
  }
} else {
  stop("Annotation file not found!")
}
# Get Arabidopsis deflines for JAG1 targets
annot_genes <- annot %>%
  select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline) %>%
  filter(locusName %in% jag1_targets$GeneID) %>%
  distinct() %>%
  group_by(locusName) %>%
  slice(1) %>%
  ungroup()

cat("Extracted annotations for", nrow(annot_genes), "JAG1 targets\n")

# Merge with JAG1 targets
jag1_with_annot <- jag1_targets %>%
  left_join(annot_genes, by = c("GeneID" = "locusName"))

n_with_annot <- sum(!is.na(jag1_with_annot$Best.hit.arabi.defline) &
                     jag1_with_annot$Best.hit.arabi.defline != "")
cat("Targets with Arabidopsis annotation:", n_with_annot,
    "(", round(n_with_annot/nrow(jag1_with_annot)*100, 1), "%)\n\n")
# Function to assign category based on defline keywords
assign_category <- function(defline) {
  if (is.na(defline) || defline == "") {
    return("Unknown/No Annotation")
  }

  defline_lower <- tolower(defline)

  # Check each category
  matched_categories <- c()

  for (cat_name in names(functional_categories)) {
    keywords <- functional_categories[[cat_name]]
    # Check if any keyword matches
    if (any(sapply(keywords, function(kw) grepl(tolower(kw), defline_lower, fixed = FALSE)))) {
      matched_categories <- c(matched_categories, cat_name)
    }
  }

  if (length(matched_categories) == 0) {
    return("Other Annotated")
  } else if (length(matched_categories) == 1) {
    return(matched_categories[1])
  } else {
    # If multiple categories, prioritize
    priority_order <- c("Transcription Factor", "Protein Kinase", "Hormone Pathway",
                        "Cell Cycle", "Signaling", "Cell Wall", "Transporter",
                        "Development", "Stress/Defense", "Metabolism")
    for (cat in priority_order) {
      if (cat %in% matched_categories) {
        return(cat)
      }
    }
    return(matched_categories[1])
  }
}

# Apply categorization
cat("Assigning functional categories...\n")
jag1_with_annot$Functional_Category <- sapply(jag1_with_annot$Best.hit.arabi.defline, assign_category)


cat("\n========================================\n")

# Overall summary
category_summary <- jag1_with_annot %>%
  group_by(Functional_Category) %>%
  summarise(
    Count = n(),
    Pct = round(n() / nrow(jag1_with_annot) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

cat("Functional category distribution:\n")
print(as.data.frame(category_summary))

# Summary by tier
category_by_tier <- jag1_with_annot %>%
  group_by(Confidence_Tier, Functional_Category) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  )


cat("\n========================================\n")

# For each category, list example genes
detailed_results <- list()

for (cat_name in unique(category_summary$Functional_Category)) {
  cat_genes <- jag1_with_annot %>%
    filter(Functional_Category == cat_name) %>%
    select(GeneID, Confidence_Tier, Best.hit.arabi.name, Best.hit.arabi.defline)

  if (nrow(cat_genes) > 0) {
    detailed_results[[cat_name]] <- cat_genes

    # Print top genes in category
    if (cat_name != "Unknown/No Annotation" && cat_name != "Other Annotated") {
      cat(cat_name, "- Top 3 genes:\n")
      print(head(cat_genes[, c("GeneID", "Best.hit.arabi.name")], 3))
      cat("\n")
    }
  }
}
# Search for leaf development keywords in defline
leaf_keywords <- c(
  "leaf", "blade", "lamina", "margin",
  "cell cycle", "cyclin", "cell division",
  "polarity", "adaxial", "abaxial",
  "boundary", "meristem", "primordium",
  "auxin", "cytokinin", "gibberellin",
  "TCP", "CUC", "NAC", "KNOX", "WOX", "YABBY"
)

leaf_results <- data.frame()

for (kw in leaf_keywords) {
  matches <- grepl(kw, jag1_with_annot$Best.hit.arabi.defline, ignore.case = TRUE)
  n_matches <- sum(matches, na.rm = TRUE)

  if (n_matches > 0) {
    leaf_results <- rbind(leaf_results, data.frame(
      Keyword = kw,
      N_Genes = n_matches,
      Example_Genes = paste(head(jag1_with_annot$GeneID[matches], 3), collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(leaf_results) > 0) {
  leaf_results <- leaf_results[order(-leaf_results$N_Genes), ]
  cat("Leaf development-related genes in JAG1 targets:\n")
  print(leaf_results)
}


cat("\n========================================\n")

# Publication theme
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Color palette
category_colors <- c(
  "Transcription Factor" = "#E41A1C",
  "Protein Kinase" = "#377EB8",
  "Hormone Pathway" = "#FF7F00",
  "Cell Cycle" = "#984EA3",
  "Signaling" = "#4DAF4A",
  "Cell Wall" = "#A65628",
  "Transporter" = "#F781BF",
  "Development" = "#999999",
  "Stress/Defense" = "#66C2A5",
  "Metabolism" = "#FFFF33",
  "Other Annotated" = "#B3B3B3",
  "Unknown/No Annotation" = "#D9D9D9"
)

# Plot 1: Bar chart
cat("Creating bar chart...\n")

plot_data <- category_summary %>%
  filter(Count >= 5) %>%
  mutate(Functional_Category = factor(Functional_Category,
                                       levels = rev(Functional_Category)))

p_bar <- ggplot(plot_data, aes(x = Functional_Category, y = Count, fill = Functional_Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "",
    y = "Number of JAG1 Target Genes",
    title = "Functional Classification of JAG1 Targets",
    subtitle = paste0("Based on Arabidopsis ortholog annotations (n = ", nrow(jag1_with_annot), " targets)")
  ) +
  theme_publication(base_size = 12) +
  theme(axis.text.y = element_text(size = 11))

ggsave("03_results/figures/23_functional_categories/functional_categories_bar.png",
       p_bar, width = 10, height = 7, dpi = 300)
ggsave("03_results/figures/23_functional_categories/functional_categories_bar.pdf",
       p_bar, width = 10, height = 7, dpi = 300)
cat("  Saved bar chart\n")

# Plot 2: Stacked bar by tier (normalized)
cat("Creating tier comparison...\n")

tier_plot_data <- category_by_tier %>%
  filter(Functional_Category %in% plot_data$Functional_Category) %>%
  group_by(Functional_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() %>%
  mutate(
    Confidence_Tier = factor(Confidence_Tier, levels = c("Gold", "Silver", "Bronze")),
    Functional_Category = factor(Functional_Category, levels = levels(plot_data$Functional_Category))
  )

p_tier <- ggplot(tier_plot_data, aes(x = Functional_Category, y = Percentage, fill = Confidence_Tier)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32"),
                    name = "Confidence Tier") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "",
    y = "Proportion of Genes",
    title = "Tier Composition by Functional Category"
  ) +
  theme_publication(base_size = 11) +
  theme(legend.position = "top")

ggsave("03_results/figures/23_functional_categories/category_by_tier.png",
       p_tier, width = 10, height = 7, dpi = 300)
ggsave("03_results/figures/23_functional_categories/category_by_tier.pdf",
       p_tier, width = 10, height = 7, dpi = 300)
cat("  Saved tier comparison plot\n")


cat("\n========================================\n")

# Overall summary
write.csv(category_summary, "03_results/tables/functional/functional_categories.csv",
          row.names = FALSE)
cat("Saved: functional_categories.csv\n")

# By tier
write.csv(category_by_tier, "03_results/tables/functional/category_by_tier.csv",
          row.names = FALSE)
cat("Saved: category_by_tier.csv\n")

# Full target list with categories
jag1_export <- jag1_with_annot %>%
  select(GeneID, Confidence_Tier, Functional_Category, Best.hit.arabi.name, Best.hit.arabi.defline)
write.csv(jag1_export, "03_results/tables/functional/genes_by_category.csv",
          row.names = FALSE)
cat("Saved: genes_by_category.csv\n")

# Leaf development genes
if (nrow(leaf_results) > 0) {
  write.csv(leaf_results, "03_results/tables/functional/leaf_development_genes.csv",
            row.names = FALSE)
  cat("Saved: leaf_development_genes.csv\n")
}


cat("\n========================================\n")

functional_analysis <- list(
  category_summary = category_summary,
  category_by_tier = category_by_tier,
  detailed_results = detailed_results,
  leaf_results = if(exists("leaf_results") && nrow(leaf_results) > 0) leaf_results else NULL,
  functional_categories = functional_categories,
  jag1_with_categories = jag1_with_annot
)

save(
  functional_analysis,
  jag1_targets,
  file = "03_results/checkpoints/23_functional_categories.RData"
)

cat("Checkpoint saved: 23_functional_categories.RData\n")


cat("\n================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("KEY RESULTS:\n")
cat("  Total JAG1 targets:", nrow(jag1_with_annot), "\n")
cat("  With Arabidopsis annotation:", n_with_annot, "(", round(n_with_annot/nrow(jag1_with_annot)*100, 1), "%)\n\n")

cat("  Top functional categories:\n")
for (i in 1:min(5, nrow(category_summary))) {
  cat("    ", category_summary$Functional_Category[i], ": ",
      category_summary$Count[i], " genes (", category_summary$Pct[i], "%)\n", sep = "")
}

cat("\nMETHOD: Keyword-based classification\n")
cat("  - ", length(functional_categories), " functional categories defined\n", sep = "")
cat("  - Searches Arabidopsis ortholog deflines for keywords\n")
cat("  - Priority ordering for multi-category genes\n\n")

cat("OUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/23_functional_categories.RData\n")
cat("  - Tables: 03_results/tables/functional/\n")
cat("    - functional_categories.csv\n")
cat("    - category_by_tier.csv\n")
cat("    - genes_by_category.csv\n")
cat("    - leaf_development_genes.csv\n")
cat("  - Figures: 03_results/figures/23_functional_categories/\n")
cat("    - functional_categories_bar.png/pdf\n")
cat("    - category_by_tier.png/pdf\n\n")

