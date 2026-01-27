
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 43: DTU FUNCTIONAL ANNOTATION\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/43_DTU_annotation", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "clusterProfiler",  # GO enrichment
  "org.At.tair.db",   # Arabidopsis annotation
  "dplyr",            # Data manipulation
  "tidyr",            # Data reshaping
  "ggplot2",          # Plotting
  "VennDiagram",      # Venn diagrams
  "RColorBrewer"      # Color palettes
)

bioc_packages <- c("clusterProfiler", "org.At.tair.db")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load packages with error handling
for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
  }, error = function(e) {
    cat("  Warning: Could not load", pkg, "\n")
  })
}

cat("  Packages loaded\n\n")
# Load DTU results
load("03_results/checkpoints/41_dtu_drimseq.RData")
cat("Loaded DTU results\n")

# Load isoform switch results
load("03_results/checkpoints/42_isoform_switches.RData")
cat("Loaded isoform switch results\n")

# Load prepared data (for background)
load("03_results/checkpoints/40_dtu_prepared.RData")
cat("Loaded prepared transcript data\n")

# Load JAG1 targets
if (!is.null(jag1_targets)) {
  cat("JAG1 targets available:", nrow(jag1_targets), "genes\n")
}

# Load gene annotation (Arabidopsis homologs)
anno_file <- "03_results/tables/JAG1_unified_targets_comprehensive.csv"
if (file.exists(anno_file)) {
  gene_anno <- read.csv(anno_file, stringsAsFactors = FALSE)
  cat("Loaded gene annotations:", nrow(gene_anno), "genes\n")
} else {
  gene_anno <- NULL
  cat("Gene annotation file not found\n")
}

cat("\n")
# DTU genes
dtu_genes <- res_gene_final$gene_id[res_gene_final$significant_DTU == TRUE]
cat("DTU genes:", length(dtu_genes), "\n")

# Genes with significant isoform switches
switch_genes <- unique(significant_switches$gene_id)
cat("Isoform switch genes:", length(switch_genes), "\n")

# Background genes (all genes tested for DTU)
background_genes <- unique(tx2gene_filt$GENEID)
cat("Background genes:", length(background_genes), "\n")

# JAG1 targets (Gold tier)
jag1_gold <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Gold"]
jag1_silver <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Silver"]
jag1_bronze <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Bronze"]
cat("JAG1 Gold targets:", length(jag1_gold), "\n")
cat("JAG1 Silver targets:", length(jag1_silver), "\n")
cat("JAG1 Bronze targets:", length(jag1_bronze), "\n\n")
# Fisher's exact test for JAG1 target enrichment in DTU genes

run_enrichment_test <- function(target_set, dtu_set, background_set, label) {
  # Contingency table
  a <- sum(dtu_set %in% target_set)                    # DTU & Target
  b <- sum(dtu_set %in% background_set) - a            # DTU & Not Target
  c <- sum(target_set %in% background_set) - a         # Not DTU & Target
  d <- length(background_set) - a - b - c              # Not DTU & Not Target

  # Fisher's exact test
  mat <- matrix(c(a, b, c, d), nrow = 2)
  test <- fisher.test(mat, alternative = "greater")

  cat(label, ":\n")
  cat("  DTU genes in target set:", a, "/", length(dtu_set), "\n")
  cat("  Expected:", round(length(target_set) * length(dtu_set) / length(background_set), 1), "\n")
  cat("  Fold enrichment:", round(a / (length(target_set) * length(dtu_set) / length(background_set)), 2), "\n")
  cat("  Fisher p-value:", format(test$p.value, digits = 3), "\n\n")

  return(list(
    observed = a,
    expected = length(target_set) * length(dtu_set) / length(background_set),
    fold_enrichment = a / (length(target_set) * length(dtu_set) / length(background_set)),
    pvalue = test$p.value
  ))
}

# Test enrichment for each tier
enrich_gold <- run_enrichment_test(jag1_gold, dtu_genes, background_genes, "Gold tier")
enrich_silver <- run_enrichment_test(jag1_silver, dtu_genes, background_genes, "Silver tier")
enrich_bronze <- run_enrichment_test(jag1_bronze, dtu_genes, background_genes, "Bronze tier")
enrich_all <- run_enrichment_test(jag1_targets$GeneID, dtu_genes, background_genes, "All JAG1 targets")

# Save enrichment results
enrichment_results <- data.frame(
  Tier = c("Gold", "Silver", "Bronze", "All"),
  DTU_in_Target = c(enrich_gold$observed, enrich_silver$observed, enrich_bronze$observed, enrich_all$observed),
  Expected = c(enrich_gold$expected, enrich_silver$expected, enrich_bronze$expected, enrich_all$expected),
  Fold_Enrichment = c(enrich_gold$fold_enrichment, enrich_silver$fold_enrichment, enrich_bronze$fold_enrichment, enrich_all$fold_enrichment),
  Pvalue = c(enrich_gold$pvalue, enrich_silver$pvalue, enrich_bronze$pvalue, enrich_all$pvalue)
)
# Map soybean genes to Arabidopsis for GO annotation
if (!is.null(gene_anno) && "Best_hit_arabi_name" %in% names(gene_anno)) {
  # Create mapping from soybean to Arabidopsis
  gene_to_arabi <- gene_anno %>%
    select(GeneID, Best_hit_arabi_name) %>%
    filter(!is.na(Best_hit_arabi_name) & Best_hit_arabi_name != "") %>%
    mutate(TAIR_ID = gsub("\\..*", "", Best_hit_arabi_name))  # Remove version

  cat("Soybean-Arabidopsis mapping available for", nrow(gene_to_arabi), "genes\n\n")

  # Map DTU genes to Arabidopsis
  dtu_arabi <- gene_to_arabi$TAIR_ID[gene_to_arabi$GeneID %in% dtu_genes]
  dtu_arabi <- unique(dtu_arabi[dtu_arabi != ""])
  cat("DTU genes with Arabidopsis mapping:", length(dtu_arabi), "\n")

  # Map background to Arabidopsis
  bg_arabi <- gene_to_arabi$TAIR_ID[gene_to_arabi$GeneID %in% background_genes]
  bg_arabi <- unique(bg_arabi[bg_arabi != ""])
  cat("Background genes with Arabidopsis mapping:", length(bg_arabi), "\n\n")

  # Run GO enrichment
  if (length(dtu_arabi) >= 10 && requireNamespace("clusterProfiler", quietly = TRUE)) {
    cat("Running GO enrichment analysis...\n")

    # Biological Process
    go_bp <- tryCatch({
      enrichGO(
        gene = dtu_arabi,
        universe = bg_arabi,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      cat("  Error in BP enrichment:", e$message, "\n")
      NULL
    })

    if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
      cat("  Biological Process: ", nrow(as.data.frame(go_bp)), " enriched terms\n", sep = "")
    } else {
      cat("  Biological Process: No significant terms\n")
    }

    # Molecular Function
    go_mf <- tryCatch({
      enrichGO(
        gene = dtu_arabi,
        universe = bg_arabi,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      cat("  Error in MF enrichment:", e$message, "\n")
      NULL
    })

    if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
      cat("  Molecular Function: ", nrow(as.data.frame(go_mf)), " enriched terms\n", sep = "")
    } else {
      cat("  Molecular Function: No significant terms\n")
    }

    # Cellular Component
    go_cc <- tryCatch({
      enrichGO(
        gene = dtu_arabi,
        universe = bg_arabi,
        OrgDb = org.At.tair.db,
        keyType = "TAIR",
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      cat("  Error in CC enrichment:", e$message, "\n")
      NULL
    })

    if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
      cat("  Cellular Component: ", nrow(as.data.frame(go_cc)), " enriched terms\n", sep = "")
    } else {
      cat("  Cellular Component: No significant terms\n")
    }
  } else {
    cat("Not enough mapped genes for GO enrichment\n")
    go_bp <- NULL
    go_mf <- NULL
    go_cc <- NULL
  }
} else {
  cat("Gene annotation not available for GO enrichment\n")
  go_bp <- NULL
  go_mf <- NULL
  go_cc <- NULL
}

cat("\n")
# Search for splicing-related terms in Arabidopsis annotations
splicing_keywords <- c("splic", "spliceosome", "RNA binding", "snRNP",
                        "serine/arginine", "SR protein", "hnRNP",
                        "U1", "U2", "U4", "U5", "U6")

if (!is.null(gene_anno)) {
  # Search in Arabidopsis gene names and descriptions
  splicing_genes <- gene_anno %>%
    filter(
      grepl(paste(splicing_keywords, collapse = "|"),
            Best_hit_arabi_name, ignore.case = TRUE)
    )

  cat("Genes with splicing-related annotations:", nrow(splicing_genes), "\n")

  if (nrow(splicing_genes) > 0) {
    # Check which are DTU genes
    splicing_dtu <- splicing_genes %>%
      filter(GeneID %in% dtu_genes)

    cat("  Among DTU genes:", nrow(splicing_dtu), "\n")

    if (nrow(splicing_dtu) > 0) {
      cat("\nSplicing factors with DTU:\n")
      print(splicing_dtu[, c("GeneID", "Best_hit_arabi_name", "Confidence_Tier")])
    }
  }
} else {
  splicing_genes <- NULL
  splicing_dtu <- NULL
}

cat("\n")
# Calculate overlaps between DTU, JAG1 targets, and DE genes

# Load DE results if available
de_file <- "03_results/tables/JAG1_unified_targets_comprehensive.csv"
if (file.exists(de_file)) {
  de_genes <- read.csv(de_file, stringsAsFactors = FALSE)$GeneID
  cat("DE genes (JAG1 targets):", length(de_genes), "\n")
} else {
  de_genes <- jag1_targets$GeneID
  cat("Using JAG1 targets as DE genes\n")
}

# Calculate overlaps
overlap_dtu_jag1 <- length(intersect(dtu_genes, jag1_targets$GeneID))
overlap_dtu_de <- length(intersect(dtu_genes, de_genes))
overlap_switch_jag1 <- length(intersect(switch_genes, jag1_targets$GeneID))

cat("\nOverlap summary:\n")
cat("  DTU + JAG1 targets:", overlap_dtu_jag1, "\n")
cat("  DTU + DE genes:", overlap_dtu_de, "\n")
cat("  Isoform switches + JAG1 targets:", overlap_switch_jag1, "\n\n")

# Unique DTU genes (not in JAG1 targets)
unique_dtu <- setdiff(dtu_genes, jag1_targets$GeneID)
cat("Unique DTU genes (not JAG1 targets):", length(unique_dtu), "\n")
# Plot 1: JAG1 target enrichment in DTU
enrichment_results$Tier <- factor(enrichment_results$Tier, levels = c("Gold", "Silver", "Bronze", "All"))

p1 <- ggplot(enrichment_results, aes(x = Tier, y = Fold_Enrichment, fill = Tier)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = paste0("p=", format(Pvalue, digits = 2))), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Gold" = "#D4AF37", "Silver" = "#8E8E8E",
                                "Bronze" = "#CD7F32", "All" = "#5DADE2")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "JAG1 Target Enrichment in DTU Genes",
    subtitle = "Fold enrichment vs random expectation",
    x = "JAG1 Target Tier",
    y = "Fold Enrichment"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("03_results/figures/43_DTU_annotation/jag1_enrichment_in_dtu.pdf", p1, width = 8, height = 6)
cat("Saved: jag1_enrichment_in_dtu.pdf\n")

# Plot 2: GO enrichment dot plot (if available)
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p2 <- dotplot(go_bp, showCategory = 15) +
    labs(title = "GO Biological Process Enrichment",
         subtitle = "DTU genes vs background") +
    theme(plot.title = element_text(face = "bold"))

  ggsave("03_results/figures/43_DTU_annotation/go_bp_dotplot.pdf", p2, width = 10, height = 8)
  cat("Saved: go_bp_dotplot.pdf\n")
}

# Plot 3: Venn diagram of overlaps
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  # Prepare sets for Venn
  venn_list <- list(
    "DTU genes" = dtu_genes,
    "JAG1 targets" = jag1_targets$GeneID,
    "Isoform switches" = switch_genes
  )

  # Generate Venn diagram
  pdf("03_results/figures/43_DTU_annotation/dtu_jag1_venn.pdf", width = 8, height = 8)
  venn.plot <- venn.diagram(
    x = venn_list,
    filename = NULL,
    fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.2,
    cat.fontface = "bold",
    main = "Overlap: DTU, JAG1 Targets, and Isoform Switches",
    main.cex = 1.5
  )
  grid.draw(venn.plot)
  dev.off()
  cat("Saved: dtu_jag1_venn.pdf\n")
}

cat("\n")
# Create comprehensive DTU gene annotation table
dtu_annotation <- res_gene_final %>%
  filter(significant_DTU == TRUE) %>%
  select(gene_id, pvalue, padj_stageR, is_JAG1_target, JAG1_tier)

# Add annotation from gene_anno
if (!is.null(gene_anno)) {
  dtu_annotation <- dtu_annotation %>%
    left_join(
      gene_anno %>% select(GeneID, Best_hit_arabi_name) %>% distinct(),
      by = c("gene_id" = "GeneID")
    )
}

# Add isoform switch info
switch_summary <- significant_switches %>%
  group_by(gene_id) %>%
  summarise(
    n_switching_transcripts = n(),
    max_delta_proportion = max(abs_delta_prop),
    has_dominant_switch = any(is_switch_gene),
    .groups = "drop"
  )

dtu_annotation <- dtu_annotation %>%
  left_join(switch_summary, by = "gene_id")

# Order by significance
dtu_annotation <- dtu_annotation %>%
  arrange(pvalue)

cat("Created comprehensive DTU annotation for", nrow(dtu_annotation), "genes\n\n")

# Display top DTU genes
cat("Top 20 DTU genes with annotation:\n")
print(head(dtu_annotation, 20))
cat("\n")
# Save checkpoint
save(
  enrichment_results,    # JAG1 enrichment test results
  go_bp, go_mf, go_cc,   # GO enrichment results
  splicing_genes,        # Splicing-related genes
  splicing_dtu,          # Splicing factors with DTU
  dtu_annotation,        # Comprehensive DTU annotation
  file = "03_results/checkpoints/43_dtu_annotation.RData"
)
cat("Saved checkpoint: 43_dtu_annotation.RData\n")

# Save tables
write.csv(enrichment_results, "03_results/tables/dtu_jag1_enrichment.csv", row.names = FALSE)
cat("Saved: dtu_jag1_enrichment.csv\n")

write.csv(dtu_annotation, "03_results/tables/dtu_genes_annotated.csv", row.names = FALSE)
cat("Saved: dtu_genes_annotated.csv\n")

# Save GO results
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  write.csv(as.data.frame(go_bp), "03_results/tables/dtu_go_bp.csv", row.names = FALSE)
  cat("Saved: dtu_go_bp.csv\n")
}

if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) > 0) {
  write.csv(as.data.frame(go_mf), "03_results/tables/dtu_go_mf.csv", row.names = FALSE)
  cat("Saved: dtu_go_mf.csv\n")
}

if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) > 0) {
  write.csv(as.data.frame(go_cc), "03_results/tables/dtu_go_cc.csv", row.names = FALSE)
  cat("Saved: dtu_go_cc.csv\n")
}

# Save overlap summary
overlap_summary <- data.frame(
  Comparison = c("DTU + JAG1 targets", "DTU + DE genes", "Isoform switches + JAG1 targets",
                 "Unique DTU (not JAG1)"),
  N_genes = c(overlap_dtu_jag1, overlap_dtu_de, overlap_switch_jag1, length(unique_dtu))
)
write.csv(overlap_summary, "03_results/tables/dtu_overlap_summary.csv", row.names = FALSE)
cat("Saved: dtu_overlap_summary.csv\n\n")


cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("SUMMARY:\n")
cat("  - JAG1 target enrichment tested (Gold/Silver/Bronze)\n")
if (!is.null(go_bp)) {
  cat("  - GO BP enrichment:", nrow(as.data.frame(go_bp)), "terms\n")
}
if (!is.null(splicing_dtu) && nrow(splicing_dtu) > 0) {
  cat("  - Splicing factors with DTU:", nrow(splicing_dtu), "\n")
}
cat("  - DTU genes annotated:", nrow(dtu_annotation), "\n")
cat("  - DTU genes overlapping JAG1 targets:", overlap_dtu_jag1, "\n\n")

