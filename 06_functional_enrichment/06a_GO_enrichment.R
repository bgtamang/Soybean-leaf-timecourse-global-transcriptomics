# Author: Bishal Tamang
# Date: 2026-01-19
#          with proper TERM2GENE mapping from Phytozome annotations

rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 22: GO ENRICHMENT ANALYSIS (clusterProfiler)\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/functional", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/22_GO_enrichment", recursive = TRUE, showWarnings = FALSE)
# Install if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("clusterProfiler", "enrichplot", "GO.db", "AnnotationDbi")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

cat("  All packages loaded\n\n")
# Load JAG1 targets
load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n")

# Load expression data for background
load("03_results/checkpoints/06_validated.RData")
cat("Loaded expression data\n")

# Get all expressed genes as background
background_genes <- rownames(v_primary$E)
cat("Background genes:", length(background_genes), "\n\n")

# Filter to actual targets
jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("JAG1 targets:", nrow(jag1_targets), "\n")


cat("\n========================================\n")

# Load Phytozome annotation file
annotation_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

cat("Loading annotation file...\n")
annot <- read.delim(annotation_file, header = TRUE, stringsAsFactors = FALSE,
                    comment.char = "", quote = "")
cat("  Loaded", nrow(annot), "transcript annotations\n")

# Extract GO column - column 10 is "GO"
# GO terms are space-separated
go_col <- "GO"
gene_col <- "locusName"

# Create TERM2GENE mapping (GO term -> Gene)
cat("Creating TERM2GENE mapping...\n")

term2gene_list <- list()

for (i in 1:nrow(annot)) {
  gene <- annot[[gene_col]][i]
  go_terms <- annot[[go_col]][i]

  if (!is.na(go_terms) && go_terms != "" && go_terms != ".") {
    # Parse space-separated GO terms
    terms <- unlist(strsplit(as.character(go_terms), " "))
    terms <- trimws(terms)
    terms <- terms[grepl("^GO:", terms)]

    if (length(terms) > 0) {
      for (term in terms) {
        if (is.null(term2gene_list[[term]])) {
          term2gene_list[[term]] <- character()
        }
        term2gene_list[[term]] <- c(term2gene_list[[term]], gene)
      }
    }
  }
}

# Convert to data frame (TERM2GENE format for clusterProfiler)
# Format: term | gene
term2gene <- do.call(rbind, lapply(names(term2gene_list), function(term) {
  data.frame(
    GO_ID = term,
    GeneID = unique(term2gene_list[[term]]),
    stringsAsFactors = FALSE
  )
}))

cat("  Total GO terms:", length(unique(term2gene$GO_ID)), "\n")
cat("  Genes with GO annotations:", length(unique(term2gene$GeneID)), "\n")

# Get unique genes per locus (remove transcript variants)
term2gene$GeneID <- gsub("\\.\\d+$", "", term2gene$GeneID)
term2gene <- unique(term2gene)

cat("  After deduplication:", nrow(term2gene), "gene-GO associations\n")


cat("\nRetrieving GO term names from GO.db...\n")

# Get all GO terms from our mapping
all_go_terms <- unique(term2gene$GO_ID)

# Query GO.db for term names and ontology
go_info <- AnnotationDbi::select(GO.db,
                                  keys = all_go_terms,
                                  columns = c("GOID", "TERM", "ONTOLOGY"),
                                  keytype = "GOID")

# Handle missing terms
go_info <- go_info[!is.na(go_info$TERM), ]

cat("  Retrieved info for", nrow(go_info), "GO terms\n")

# Create TERM2NAME mapping
term2name <- go_info %>%
  select(GO_ID = GOID, Name = TERM) %>%
  distinct()

# Separate by ontology
bp_terms <- go_info$GOID[go_info$ONTOLOGY == "BP"]
mf_terms <- go_info$GOID[go_info$ONTOLOGY == "MF"]
cc_terms <- go_info$GOID[go_info$ONTOLOGY == "CC"]

cat("  BP terms:", length(bp_terms), "\n")
cat("  MF terms:", length(mf_terms), "\n")
cat("  CC terms:", length(cc_terms), "\n")

# Create ontology-specific TERM2GENE mappings
term2gene_bp <- term2gene[term2gene$GO_ID %in% bp_terms, ]
term2gene_mf <- term2gene[term2gene$GO_ID %in% mf_terms, ]
term2gene_cc <- term2gene[term2gene$GO_ID %in% cc_terms, ]


cat("\n========================================\n")

run_go_enrichment <- function(gene_list, universe, term2gene, term2name,
                               ontology_name = "GO", pval_cutoff = 0.05,
                               qval_cutoff = 0.1, min_gs = 5, max_gs = 500) {

  cat("  Running enrichment for", ontology_name, "...\n")
  cat("    Input genes:", length(gene_list), "\n")

  # Filter to genes in universe
  gene_list <- intersect(gene_list, universe)
  cat("    Genes in universe:", length(gene_list), "\n")

  # Run enricher
  result <- tryCatch({
    enricher(
      gene = gene_list,
      universe = universe,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = min_gs,
      maxGSSize = max_gs
    )
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(result) && nrow(result@result) > 0) {
    cat("    Significant terms (q <", qval_cutoff, "):", sum(result@result$qvalue < qval_cutoff), "\n")
  } else {
    cat("    No significant terms found\n")
  }

  return(result)
}


# Define gene lists by tier
gold_genes <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Gold"]
silver_genes <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Silver"]
bronze_genes <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Bronze"]
all_target_genes <- jag1_targets$GeneID

# Use genes with GO annotations as universe
universe_genes <- unique(term2gene$GeneID)
universe_genes <- intersect(universe_genes, background_genes)

cat("Universe (expressed genes with GO):", length(universe_genes), "\n\n")

# Store all results
go_results <- list()

cat("Total genes:", length(all_target_genes), "\n\n")

go_results$all_BP <- run_go_enrichment(all_target_genes, universe_genes,
                                        term2gene_bp, term2name, "BP")
go_results$all_MF <- run_go_enrichment(all_target_genes, universe_genes,
                                        term2gene_mf, term2name, "MF")
go_results$all_CC <- run_go_enrichment(all_target_genes, universe_genes,
                                        term2gene_cc, term2name, "CC")

cat("\n=== GOLD TIER ===\n")
cat("Total genes:", length(gold_genes), "\n\n")

go_results$gold_BP <- run_go_enrichment(gold_genes, universe_genes,
                                         term2gene_bp, term2name, "BP")
go_results$gold_MF <- run_go_enrichment(gold_genes, universe_genes,
                                         term2gene_mf, term2name, "MF")
go_results$gold_CC <- run_go_enrichment(gold_genes, universe_genes,
                                         term2gene_cc, term2name, "CC")

cat("\n=== SILVER TIER ===\n")
cat("Total genes:", length(silver_genes), "\n\n")

go_results$silver_BP <- run_go_enrichment(silver_genes, universe_genes,
                                           term2gene_bp, term2name, "BP")
go_results$silver_MF <- run_go_enrichment(silver_genes, universe_genes,
                                           term2gene_mf, term2name, "MF")
go_results$silver_CC <- run_go_enrichment(silver_genes, universe_genes,
                                           term2gene_cc, term2name, "CC")

cat("\n=== BRONZE TIER ===\n")
cat("Total genes:", length(bronze_genes), "\n\n")

go_results$bronze_BP <- run_go_enrichment(bronze_genes, universe_genes,
                                           term2gene_bp, term2name, "BP")
go_results$bronze_MF <- run_go_enrichment(bronze_genes, universe_genes,
                                           term2gene_mf, term2name, "MF")
go_results$bronze_CC <- run_go_enrichment(bronze_genes, universe_genes,
                                           term2gene_cc, term2name, "CC")


cat("\n========================================\n")

# Function to save enrichment result
save_enrichment <- function(result, filename) {
  if (!is.null(result) && nrow(result@result) > 0) {
    df <- as.data.frame(result)
    write.csv(df, filename, row.names = FALSE)
    cat("  Saved:", basename(filename), "(", nrow(df), "terms )\n")
    return(df)
  } else {
    cat("  Skipped:", basename(filename), "(no results)\n")
    return(NULL)
  }
}

# Save all results
save_enrichment(go_results$all_BP, "03_results/tables/functional/GO_all_targets_BP.csv")
save_enrichment(go_results$all_MF, "03_results/tables/functional/GO_all_targets_MF.csv")
save_enrichment(go_results$all_CC, "03_results/tables/functional/GO_all_targets_CC.csv")

save_enrichment(go_results$gold_BP, "03_results/tables/functional/GO_gold_tier_BP.csv")
save_enrichment(go_results$gold_MF, "03_results/tables/functional/GO_gold_tier_MF.csv")
save_enrichment(go_results$gold_CC, "03_results/tables/functional/GO_gold_tier_CC.csv")

save_enrichment(go_results$silver_BP, "03_results/tables/functional/GO_silver_tier_BP.csv")
save_enrichment(go_results$silver_MF, "03_results/tables/functional/GO_silver_tier_MF.csv")
save_enrichment(go_results$silver_CC, "03_results/tables/functional/GO_silver_tier_CC.csv")

save_enrichment(go_results$bronze_BP, "03_results/tables/functional/GO_bronze_tier_BP.csv")
save_enrichment(go_results$bronze_MF, "03_results/tables/functional/GO_bronze_tier_MF.csv")
save_enrichment(go_results$bronze_CC, "03_results/tables/functional/GO_bronze_tier_CC.csv")

# Create combined summary for all targets (backward compatibility)
all_combined <- bind_rows(
  if (!is.null(go_results$all_BP)) as.data.frame(go_results$all_BP) %>% mutate(Ontology = "BP") else NULL,
  if (!is.null(go_results$all_MF)) as.data.frame(go_results$all_MF) %>% mutate(Ontology = "MF") else NULL,
  if (!is.null(go_results$all_CC)) as.data.frame(go_results$all_CC) %>% mutate(Ontology = "CC") else NULL
)

if (nrow(all_combined) > 0) {
  write.csv(all_combined, "03_results/tables/functional/GO_all_targets_combined.csv", row.names = FALSE)
  cat("  Saved: GO_all_targets_combined.csv (", nrow(all_combined), "terms )\n")
}


cat("\n========================================\n")

# Publication theme
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.title = element_text(size = base_size + 1, face = "bold"),
      axis.text = element_text(size = base_size, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}


if (!is.null(go_results$all_BP) && nrow(go_results$all_BP@result) > 0) {
  cat("Creating BP dotplot...\n")

  p_bp <- dotplot(go_results$all_BP, showCategory = 20,
                  title = "GO Biological Process - All JAG1 Targets") +
    theme_publication(base_size = 10)

  ggsave("03_results/figures/22_GO_enrichment/GO_BP_all_targets_dotplot.png",
         p_bp, width = 10, height = 8, dpi = 300)
  ggsave("03_results/figures/22_GO_enrichment/GO_BP_all_targets_dotplot.pdf",
         p_bp, width = 10, height = 8, device = cairo_pdf)

  cat("  Saved BP dotplot\n")
}

if (!is.null(go_results$all_MF) && nrow(go_results$all_MF@result) > 0) {
  cat("Creating MF dotplot...\n")

  p_mf <- dotplot(go_results$all_MF, showCategory = 20,
                  title = "GO Molecular Function - All JAG1 Targets") +
    theme_publication(base_size = 10)

  ggsave("03_results/figures/22_GO_enrichment/GO_MF_all_targets_dotplot.png",
         p_mf, width = 10, height = 8, dpi = 300)
  ggsave("03_results/figures/22_GO_enrichment/GO_MF_all_targets_dotplot.pdf",
         p_mf, width = 10, height = 8, device = cairo_pdf)

  cat("  Saved MF dotplot\n")
}

if (!is.null(go_results$all_CC) && nrow(go_results$all_CC@result) > 0) {
  cat("Creating CC dotplot...\n")

  p_cc <- dotplot(go_results$all_CC, showCategory = 20,
                  title = "GO Cellular Component - All JAG1 Targets") +
    theme_publication(base_size = 10)

  ggsave("03_results/figures/22_GO_enrichment/GO_CC_all_targets_dotplot.png",
         p_cc, width = 10, height = 8, dpi = 300)
  ggsave("03_results/figures/22_GO_enrichment/GO_CC_all_targets_dotplot.pdf",
         p_cc, width = 10, height = 8, device = cairo_pdf)

  cat("  Saved CC dotplot\n")
}


if (!is.null(go_results$all_BP) && nrow(go_results$all_BP@result) > 0) {
  cat("Creating BP barplot...\n")

  p_bar <- barplot(go_results$all_BP, showCategory = 15,
                   title = "GO Biological Process Enrichment") +
    theme_publication(base_size = 10)

  ggsave("03_results/figures/22_GO_enrichment/GO_BP_all_targets_barplot.png",
         p_bar, width = 10, height = 7, dpi = 300)

  cat("  Saved BP barplot\n")
}


if (!is.null(go_results$all_BP) && nrow(go_results$all_BP@result) >= 5) {
  cat("Creating concept network plot...\n")

  tryCatch({
    p_cnet <- cnetplot(go_results$all_BP, showCategory = 5,
                       categorySize = "pvalue",
                       foldChange = NULL) +
      ggtitle("Gene-Concept Network: Top BP Terms")

    ggsave("03_results/figures/22_GO_enrichment/GO_BP_cnetplot.png",
           p_cnet, width = 12, height = 10, dpi = 300)

    cat("  Saved concept network plot\n")
  }, error = function(e) {
    cat("  Could not create cnetplot:", e$message, "\n")
  })
}


if (!is.null(go_results$all_BP) && nrow(go_results$all_BP@result) >= 3) {
  cat("Creating upset plot...\n")

  tryCatch({
    p_upset <- upsetplot(go_results$all_BP, n = 10)

    ggsave("03_results/figures/22_GO_enrichment/GO_BP_upset.png",
           p_upset, width = 10, height = 6, dpi = 300)

    cat("  Saved upset plot\n")
  }, error = function(e) {
    cat("  Could not create upset plot:", e$message, "\n")
  })
}


cat("\nCreating tier comparison...\n")

# Combine BP results from all tiers
tier_comparison <- bind_rows(
  if (!is.null(go_results$gold_BP)) as.data.frame(go_results$gold_BP) %>% mutate(Tier = "Gold") else NULL,
  if (!is.null(go_results$silver_BP)) as.data.frame(go_results$silver_BP) %>% mutate(Tier = "Silver") else NULL,
  if (!is.null(go_results$bronze_BP)) as.data.frame(go_results$bronze_BP) %>% mutate(Tier = "Bronze") else NULL
)

if (nrow(tier_comparison) > 0) {
  # Get top terms from each tier
  top_terms <- tier_comparison %>%
    group_by(Tier) %>%
    slice_min(qvalue, n = 5) %>%
    ungroup()

  p_tier <- ggplot(top_terms, aes(x = reorder(Description, -qvalue), y = -log10(qvalue),
                                   fill = Tier)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")) +
    labs(x = "", y = "-log10(q-value)",
         title = "Top GO BP Terms by Confidence Tier") +
    theme_publication(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave("03_results/figures/22_GO_enrichment/GO_BP_tier_comparison.png",
         p_tier, width = 12, height = 8, dpi = 300)
  ggsave("03_results/figures/22_GO_enrichment/GO_BP_tier_comparison.pdf",
         p_tier, width = 12, height = 8, device = cairo_pdf)

  cat("  Saved tier comparison plot\n")
}


cat("\n========================================\n")

# Count significant terms per analysis
summary_stats <- data.frame(
  Analysis = character(),
  Ontology = character(),
  Input_Genes = integer(),
  Genes_With_GO = integer(),
  Sig_Terms_q01 = integer(),
  Sig_Terms_q05 = integer(),
  stringsAsFactors = FALSE
)

count_sig <- function(result, analysis, ontology, n_input) {
  # Handle NULL or empty results
  if (is.null(result)) {
    return(data.frame(
      Analysis = analysis,
      Ontology = ontology,
      Input_Genes = n_input,
      Genes_With_GO = 0,
      Sig_Terms_q01 = 0,
      Sig_Terms_q05 = 0,
      stringsAsFactors = FALSE
    ))
  }

  # Try to convert to data frame
  df <- tryCatch({
    as.data.frame(result)
  }, error = function(e) {
    return(NULL)
  })

  # Handle empty data frame
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      Analysis = analysis,
      Ontology = ontology,
      Input_Genes = n_input,
      Genes_With_GO = 0,
      Sig_Terms_q01 = 0,
      Sig_Terms_q05 = 0,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Analysis = analysis,
    Ontology = ontology,
    Input_Genes = n_input,
    Genes_With_GO = length(unique(unlist(strsplit(paste(df$geneID, collapse = "/"), "/")))),
    Sig_Terms_q01 = sum(df$qvalue < 0.01, na.rm = TRUE),
    Sig_Terms_q05 = sum(df$qvalue < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

summary_stats <- bind_rows(
  count_sig(go_results$all_BP, "All Targets", "BP", length(all_target_genes)),
  count_sig(go_results$all_MF, "All Targets", "MF", length(all_target_genes)),
  count_sig(go_results$all_CC, "All Targets", "CC", length(all_target_genes)),
  count_sig(go_results$gold_BP, "Gold", "BP", length(gold_genes)),
  count_sig(go_results$gold_MF, "Gold", "MF", length(gold_genes)),
  count_sig(go_results$gold_CC, "Gold", "CC", length(gold_genes)),
  count_sig(go_results$silver_BP, "Silver", "BP", length(silver_genes)),
  count_sig(go_results$silver_MF, "Silver", "MF", length(silver_genes)),
  count_sig(go_results$silver_CC, "Silver", "CC", length(silver_genes)),
  count_sig(go_results$bronze_BP, "Bronze", "BP", length(bronze_genes)),
  count_sig(go_results$bronze_MF, "Bronze", "MF", length(bronze_genes)),
  count_sig(go_results$bronze_CC, "Bronze", "CC", length(bronze_genes))
)

print(summary_stats)

write.csv(summary_stats, "03_results/tables/functional/GO_enrichment_summary.csv", row.names = FALSE)
cat("\nSaved: GO_enrichment_summary.csv\n")


cat("\n========================================\n")

go_analysis <- list(
  go_results = go_results,
  term2gene = term2gene,
  term2name = term2name,
  universe_genes = universe_genes,
  summary_stats = summary_stats
)

save(
  go_analysis,
  go_results,
  term2gene,
  term2name,
  jag1_targets,
  file = "03_results/checkpoints/22_GO_enrichment.RData"
)

cat("Checkpoint saved: 22_GO_enrichment.RData\n")


cat("\n================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("KEY RESULTS:\n")
cat("  - Total GO terms in database:", length(unique(term2gene$GO_ID)), "\n")
cat("  - Genes with GO annotations:", length(unique(term2gene$GeneID)), "\n")
cat("  - Universe size:", length(universe_genes), "\n\n")

cat("ENRICHMENT RESULTS (q < 0.05):\n")
for (i in 1:nrow(summary_stats)) {
  if (summary_stats$Sig_Terms_q05[i] > 0) {
    cat("  -", summary_stats$Analysis[i], summary_stats$Ontology[i], ":",
        summary_stats$Sig_Terms_q05[i], "terms\n")
  }
}

cat("\nOUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/22_GO_enrichment.RData\n")
cat("  - Tables: 03_results/tables/functional/GO_*.csv\n")
cat("  - Figures: 03_results/figures/22_GO_enrichment/\n\n")

