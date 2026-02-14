# Script 22.1: WGCNA Module GO Enrichment

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 22.1: WGCNA MODULE GO ENRICHMENT\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/functional/modules", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/22.1_module_GO", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "WGCNA",
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "pheatmap"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINTS =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

# Load WGCNA results
load("03_results/checkpoints/21_WGCNA_JAG1.RData")
cat("Loaded WGCNA JAG1 analysis data\n")

# Load expression data for background
load("03_results/checkpoints/06_validated.RData")
cat("Loaded expression data\n")

# Get all expressed genes as background
background_genes <- rownames(v_primary$E)
cat("Background genes:", length(background_genes), "\n\n")

# Get module assignments
genes_in_wgcna <- colnames(datExpr)
module_assignments <- net$colors
names(module_assignments) <- genes_in_wgcna
module_colors_all <- labels2colors(module_assignments)
names(module_colors_all) <- genes_in_wgcna

cat("Genes in WGCNA:", length(genes_in_wgcna), "\n")
cat("Modules:", length(unique(module_assignments)), "\n\n")

# ===== LOAD GO ANNOTATIONS =====

cat("========================================\n")
cat("SECTION 2: LOAD GO ANNOTATIONS\n")
cat("========================================\n\n")

# Load annotation file
annotation_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (file.exists(annotation_file)) {
  cat("Loading annotation file...\n")
  annot <- read.delim(annotation_file, header = TRUE, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(annot), "gene annotations\n")

  # Identify gene and GO columns
  gene_col <- NULL
  for (col in c("locusName", "transcript", "gene_id", "Gene", "GeneID")) {
    if (col %in% colnames(annot)) {
      gene_col <- col
      break
    }
  }

  go_col <- NULL
  for (col in c("GO", "GO.terms", "Gene.Ontology", "GO_terms", "go_id")) {
    if (col %in% colnames(annot)) {
      go_col <- col
      break
    }
  }

  # Create GO mapping
  if (!is.null(gene_col) && !is.null(go_col)) {
    cat("Using gene column:", gene_col, "\n")
    cat("Using GO column:", go_col, "\n")

    go_mapping <- list()

    for (i in 1:nrow(annot)) {
      gene <- annot[[gene_col]][i]
      go_terms <- annot[[go_col]][i]

      if (!is.na(go_terms) && go_terms != "" && go_terms != ".") {
        terms <- unlist(strsplit(as.character(go_terms), "[,;|]"))
        terms <- trimws(terms)
        terms <- terms[grepl("^GO:", terms)]

        if (length(terms) > 0) {
          go_mapping[[gene]] <- terms
        }
      }
    }

    cat("Genes with GO annotations:", length(go_mapping), "\n")
    has_go_mapping <- TRUE

  } else {
    cat("Could not identify GO columns in annotation file\n")
    has_go_mapping <- FALSE
  }
} else {
  cat("Annotation file not found\n")
  has_go_mapping <- FALSE
}

# ===== DEFINE GO ENRICHMENT FUNCTION =====

cat("\n========================================\n")
cat("SECTION 3: GO ENRICHMENT FUNCTION\n")
cat("========================================\n\n")

perform_go_enrichment <- function(gene_list, background, go_mapping,
                                  pval_cutoff = 0.05, min_genes = 5) {

  # Filter to genes with GO annotations
  genes_with_go <- intersect(gene_list, names(go_mapping))
  bg_with_go <- intersect(background, names(go_mapping))

  if (length(genes_with_go) < min_genes) {
    return(NULL)
  }

  # Get GO terms for target genes
  go_counts <- table(unlist(go_mapping[genes_with_go]))
  bg_go_counts <- table(unlist(go_mapping[bg_with_go]))

  # Perform hypergeometric test for each term
  results <- data.frame()

  for (term in names(go_counts)) {
    k <- go_counts[term]  # Module genes with this term
    n <- length(genes_with_go)  # Total module genes with GO
    K <- bg_go_counts[term]  # Background with this term
    N <- length(bg_with_go)  # Total background

    if (is.na(K)) K <- 0
    if (K < 3) next  # Skip rare terms

    # Hypergeometric test
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

    # Calculate fold enrichment
    expected <- n * K / N
    fold_enrichment <- k / expected

    results <- rbind(results, data.frame(
      GO_ID = term,
      Count = as.numeric(k),
      Background_Count = as.numeric(K),
      List_Size = n,
      Background_Size = N,
      Expected = round(expected, 2),
      Fold_Enrichment = round(fold_enrichment, 2),
      P_value = pval,
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(results) == 0) return(NULL)

  # Adjust p-values
  results$FDR <- p.adjust(results$P_value, method = "BH")

  # Filter significant
  results <- results[results$FDR < pval_cutoff, ]
  results <- results[order(results$FDR), ]

  return(results)
}

cat("GO enrichment function defined\n")

# ===== IDENTIFY KEY MODULES =====

cat("\n========================================\n")
cat("SECTION 4: IDENTIFY KEY MODULES\n")
cat("========================================\n\n")

# Key modules based on our previous analysis:
# 1. Modules enriched for JAG1 targets
# 2. Modules correlated with Narrow phenotype
# 3. JAG1's own module

# Get module info from enrichment results
key_modules <- data.frame(
  Module_Num = c(18, 16, 15, 6, 1, 4),
  Module_Color = c("lightgreen", "lightcyan", "midnightblue", "red", "turquoise", "yellow"),
  Reason = c(
    "r=0.955 with Narrow, 10.34x target enrichment",
    "15.99x JAG1 target enrichment (highest)",
    "9.2x JAG1 target enrichment",
    "Contains JAG1 (transcription factor)",
    "Largest module, 717 JAG1 targets",
    "1.24x enriched, correlates with phenotype"
  ),
  stringsAsFactors = FALSE
)

# Add module sizes
for (i in 1:nrow(key_modules)) {
  mod_num <- key_modules$Module_Num[i]
  key_modules$Size[i] <- sum(module_assignments == mod_num)
  key_modules$JAG1_Targets[i] <- sum(module_assignments == mod_num &
                                       names(module_assignments) %in% jag1_target_genes)
}

cat("Key modules to analyze:\n")
print(key_modules)
cat("\n")

# ===== RUN GO ENRICHMENT ON KEY MODULES =====

cat("========================================\n")
cat("SECTION 5: MODULE GO ENRICHMENT\n")
cat("========================================\n\n")

module_go_results <- list()

if (has_go_mapping) {

  for (i in 1:nrow(key_modules)) {
    mod_num <- key_modules$Module_Num[i]
    mod_color <- key_modules$Module_Color[i]

    cat("Processing module", mod_num, "(", mod_color, ")...\n")

    # Get genes in this module
    module_genes <- names(module_assignments)[module_assignments == mod_num]

    cat("  Module genes:", length(module_genes), "\n")

    # Run GO enrichment
    go_result <- perform_go_enrichment(module_genes, background_genes, go_mapping)

    if (!is.null(go_result) && nrow(go_result) > 0) {
      go_result$Module <- mod_color
      go_result$Module_Num <- mod_num
      module_go_results[[mod_color]] <- go_result

      cat("  Significant GO terms:", nrow(go_result), "\n")

      # Save individual module results
      write.csv(go_result,
                paste0("03_results/tables/functional/modules/GO_module_", mod_color, ".csv"),
                row.names = FALSE)

      # Print top 5
      cat("  Top 5 terms:\n")
      print(head(go_result[, c("GO_ID", "Count", "Fold_Enrichment", "FDR")], 5))

    } else {
      cat("  No significant GO terms found\n")
    }

    cat("\n")
  }

  # Also analyze all other modules with decent size
  cat("Analyzing additional modules...\n\n")

  all_modules <- unique(module_assignments)
  all_modules <- all_modules[all_modules != 0]  # Exclude grey
  additional_modules <- setdiff(all_modules, key_modules$Module_Num)

  for (mod_num in additional_modules) {
    mod_color <- labels2colors(mod_num)
    module_genes <- names(module_assignments)[module_assignments == mod_num]

    if (length(module_genes) >= 50) {  # Only analyze larger modules
      go_result <- perform_go_enrichment(module_genes, background_genes, go_mapping)

      if (!is.null(go_result) && nrow(go_result) > 0) {
        go_result$Module <- mod_color
        go_result$Module_Num <- mod_num
        module_go_results[[mod_color]] <- go_result

        cat("  Module", mod_num, "(", mod_color, "):", nrow(go_result), "significant terms\n")

        write.csv(go_result,
                  paste0("03_results/tables/functional/modules/GO_module_", mod_color, ".csv"),
                  row.names = FALSE)
      }
    }
  }

} else {
  cat("No GO mapping available. Cannot perform module enrichment.\n")
}

# ===== CREATE COMBINED RESULTS TABLE =====

cat("\n========================================\n")
cat("SECTION 6: COMBINED RESULTS\n")
cat("========================================\n\n")

if (length(module_go_results) > 0) {

  # Combine all module GO results
  combined_go <- do.call(rbind, module_go_results)

  cat("Combined GO results:\n")
  cat("  Total modules with GO enrichment:", length(module_go_results), "\n")
  cat("  Total significant terms:", nrow(combined_go), "\n")

  write.csv(combined_go,
            "03_results/tables/functional/modules/GO_all_modules_combined.csv",
            row.names = FALSE)
  cat("  Saved: GO_all_modules_combined.csv\n\n")

  # Summarize top terms per module
  module_summary <- data.frame()

  for (mod in names(module_go_results)) {
    go_data <- module_go_results[[mod]]

    if (nrow(go_data) > 0) {
      top_terms <- head(go_data$GO_ID, 5)

      module_summary <- rbind(module_summary, data.frame(
        Module = mod,
        N_Significant_Terms = nrow(go_data),
        Top_Terms = paste(top_terms, collapse = "; "),
        Min_FDR = min(go_data$FDR),
        stringsAsFactors = FALSE
      ))
    }
  }

  module_summary <- module_summary[order(module_summary$Min_FDR), ]

  cat("Module GO summary:\n")
  print(module_summary[, c("Module", "N_Significant_Terms", "Min_FDR")])

  write.csv(module_summary,
            "03_results/tables/functional/modules/GO_module_summary.csv",
            row.names = FALSE)
}

# ===== VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 7: VISUALIZATIONS\n")
cat("========================================\n\n")

if (length(module_go_results) > 0) {

  # Plot 1: Top GO terms for key modules
  for (mod_color in c("lightgreen", "lightcyan", "midnightblue")) {

    if (mod_color %in% names(module_go_results)) {
      go_data <- module_go_results[[mod_color]]
      top_go <- head(go_data, 15)

      if (nrow(top_go) > 0) {
        png(paste0("03_results/figures/22.1_module_GO/GO_", mod_color, "_top15.png"),
            width = 1000, height = 700, res = 120)

        par(mar = c(5, 20, 4, 2))
        barplot(rev(-log10(top_go$FDR)),
                names.arg = rev(substr(top_go$GO_ID, 1, 50)),
                horiz = TRUE, las = 1,
                main = paste0("GO Enrichment - ", mod_color, " module"),
                xlab = "-log10(FDR)",
                col = mod_color,
                border = "white")
        abline(v = -log10(0.05), lty = 2, col = "red")

        dev.off()
        cat("Saved: GO_", mod_color, "_top15.png\n")
      }
    }
  }

  # Plot 2: Comparison across key modules
  # Get top 5 terms from each key module
  top_terms_list <- list()

  for (mod in names(module_go_results)[names(module_go_results) %in% key_modules$Module_Color]) {
    top_terms_list[[mod]] <- head(module_go_results[[mod]]$GO_ID, 5)
  }

  all_top_terms <- unique(unlist(top_terms_list))

  if (length(all_top_terms) > 0 && length(all_top_terms) <= 50) {

    # Create matrix of enrichment values
    enrich_matrix <- matrix(0,
                            nrow = length(all_top_terms),
                            ncol = length(top_terms_list),
                            dimnames = list(all_top_terms, names(top_terms_list)))

    for (mod in names(top_terms_list)) {
      go_data <- module_go_results[[mod]]
      for (term in all_top_terms) {
        if (term %in% go_data$GO_ID) {
          enrich_matrix[term, mod] <- go_data$Fold_Enrichment[go_data$GO_ID == term]
        }
      }
    }

    # Heatmap
    png("03_results/figures/22.1_module_GO/GO_module_comparison_heatmap.png",
        width = 900, height = 800, res = 120)

    pheatmap(enrich_matrix,
             main = "GO Term Enrichment Across Key Modules",
             color = colorRampPalette(c("white", "orange", "red"))(50),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             fontsize_row = 7,
             fontsize_col = 10,
             display_numbers = TRUE,
             number_format = "%.1f",
             number_color = "black")

    dev.off()
    cat("Saved: GO_module_comparison_heatmap.png\n")
  }

  # Plot 3: Summary bar plot - number of significant terms per module
  if (nrow(module_summary) > 0) {
    png("03_results/figures/22.1_module_GO/GO_terms_per_module.png",
        width = 900, height = 600, res = 120)

    par(mar = c(7, 5, 4, 2))

    # Sort by number of terms
    plot_data <- module_summary[order(-module_summary$N_Significant_Terms), ]
    plot_data <- head(plot_data, 15)

    # Color bars by module color
    bar_colors <- plot_data$Module

    barplot(plot_data$N_Significant_Terms,
            names.arg = plot_data$Module,
            col = bar_colors,
            las = 2,
            main = "Number of Significant GO Terms per Module",
            ylab = "N Significant GO Terms (FDR < 0.05)",
            border = "white")

    dev.off()
    cat("Saved: GO_terms_per_module.png\n")
  }

  # Plot 4: Dotplot for lightgreen module (phenotype-correlated)
  if ("lightgreen" %in% names(module_go_results)) {
    go_data <- module_go_results[["lightgreen"]]
    top_go <- head(go_data, 20)
    top_go$GeneRatio <- top_go$Count / top_go$List_Size

    p <- ggplot(top_go, aes(x = GeneRatio, y = reorder(GO_ID, GeneRatio))) +
      geom_point(aes(size = Count, color = -log10(FDR))) +
      scale_color_gradient(low = "lightgreen", high = "darkgreen") +
      theme_bw() +
      labs(title = "GO Enrichment - Lightgreen Module\n(r=0.955 correlation with Narrow phenotype)",
           x = "Gene Ratio",
           y = "GO Term",
           size = "Gene Count",
           color = "-log10(FDR)") +
      theme(axis.text.y = element_text(size = 7))

    ggsave("03_results/figures/22.1_module_GO/GO_lightgreen_dotplot.png",
           p, width = 10, height = 8, dpi = 120)
    cat("Saved: GO_lightgreen_dotplot.png\n")
  }
}

# ===== BIOLOGICAL INTERPRETATION =====

cat("\n========================================\n")
cat("SECTION 8: BIOLOGICAL INTERPRETATION\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n\n")

if (length(module_go_results) > 0) {

  # Lightgreen module (phenotype-correlated)
  if ("lightgreen" %in% names(module_go_results)) {
    cat("LIGHTGREEN MODULE (r=0.955 with Narrow phenotype):\n")
    cat("  - This module almost perfectly tracks leaf shape\n")
    cat("  - Top enriched functions:\n")
    top_lg <- head(module_go_results[["lightgreen"]], 5)
    for (i in 1:min(5, nrow(top_lg))) {
      cat("    ", i, ". ", top_lg$GO_ID[i], " (FDR=",
          formatC(top_lg$FDR[i], format = "e", digits = 2), ")\n", sep = "")
    }
    cat("\n")
  }

  # Lightcyan module (highest target enrichment)
  if ("lightcyan" %in% names(module_go_results)) {
    cat("LIGHTCYAN MODULE (15.99x JAG1 target enrichment):\n")
    cat("  - Contains highest concentration of JAG1 targets\n")
    cat("  - Top enriched functions:\n")
    top_lc <- head(module_go_results[["lightcyan"]], 5)
    for (i in 1:min(5, nrow(top_lc))) {
      cat("    ", i, ". ", top_lc$GO_ID[i], " (FDR=",
          formatC(top_lc$FDR[i], format = "e", digits = 2), ")\n", sep = "")
    }
    cat("\n")
  }

  # Midnightblue module
  if ("midnightblue" %in% names(module_go_results)) {
    cat("MIDNIGHTBLUE MODULE (9.2x JAG1 target enrichment):\n")
    cat("  - Second-highest concentration of JAG1 targets\n")
    cat("  - Top enriched functions:\n")
    top_mb <- head(module_go_results[["midnightblue"]], 5)
    for (i in 1:min(5, nrow(top_mb))) {
      cat("    ", i, ". ", top_mb$GO_ID[i], " (FDR=",
          formatC(top_mb$FDR[i], format = "e", digits = 2), ")\n", sep = "")
    }
    cat("\n")
  }
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 9: SAVE CHECKPOINT\n")
cat("========================================\n\n")

module_go_analysis <- list(
  module_go_results = module_go_results,
  key_modules = key_modules,
  has_go_mapping = has_go_mapping,
  n_modules_analyzed = length(module_go_results)
)

if (exists("module_summary")) {
  module_go_analysis$module_summary <- module_summary
}

save(
  module_go_analysis,
  module_go_results,
  key_modules,
  file = "03_results/checkpoints/22.1_module_GO.RData"
)

cat("Checkpoint saved: 22.1_module_GO.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 22.1 COMPLETE: WGCNA MODULE GO ENRICHMENT\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat("  - Modules analyzed:", length(module_go_results), "\n")
if (length(module_go_results) > 0) {
  cat("  - Total significant GO terms:", sum(sapply(module_go_results, nrow)), "\n")
  cat("\n  Modules with enrichment:\n")
  for (mod in names(module_go_results)) {
    cat("    - ", mod, ": ", nrow(module_go_results[[mod]]), " terms\n", sep = "")
  }
}

cat("\nOUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/22.1_module_GO.RData\n")
cat("  - Tables: 03_results/tables/functional/modules/\n")
cat("    - GO_module_*.csv (individual modules)\n")
cat("    - GO_all_modules_combined.csv\n")
cat("    - GO_module_summary.csv\n")
cat("  - Figures: 03_results/figures/22.1_module_GO/\n\n")

cat("INTERPRETATION:\n")
cat("  - Lightgreen module (r=0.955): Functions directly linked to leaf shape\n")
cat("  - Lightcyan module (15.99x): Core JAG1-regulated functions\n")
cat("  - Different modules = different regulatory pathways\n\n")

cat("NEXT STEP: Continue to Script 23 (KEGG pathways) or integrate results\n\n")
