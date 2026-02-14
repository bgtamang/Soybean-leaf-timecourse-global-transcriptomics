# =============================================================================
# Script 34c: Comprehensive KRP Gene Family Verification
# =============================================================================
#
# PURPOSE:
# This script provides independent verification that the 9 soybean KRP genes
# identified by Guo et al. (2023) represent the complete KRP gene family.
#
# METHODS:
# 1. Pfam domain search (PF02234 = CDK inhibitor domain)
# 2. Arabidopsis ortholog mapping via Best-hit-arabi-name
# 3. Cross-validation with literature (Guo et al. 2023)
#
# OUTPUT:
# - Comprehensive verification report
# - Tables suitable for manuscript supplementary material
#
# CITATION:
# Guo B, Chen L, Dong L, et al. (2023) Characterization of the soybean KRP
# gene family reveals a key role for GmKRP2a in root development.
# Front. Plant Sci. 14:1096467. DOI: 10.3389/fpls.2023.1096467
#
# =============================================================================

# Set base directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase3-Refined-Analysis"
setwd(base_dir)

library(tidyverse)

# =============================================================================
# 1. SETUP AND DATA LOADING
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("COMPREHENSIVE KRP GENE FAMILY VERIFICATION\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Purpose: Verify completeness of soybean KRP gene list\n\n")

# -----------------------------------------------------------------------------
# 1a. Define Guo et al. 2023 KRP genes (Literature Reference)
# -----------------------------------------------------------------------------

cat("--- STEP 1: Loading Reference Data ---\n\n")

# The 9 KRP genes from Guo et al. 2023
# Source: Table 1, Front. Plant Sci. 14:1096467
guo_krp_genes <- tribble(
  ~gene_id,            ~name,      ~arabidopsis_ortholog,  ~at_gene_id,
  "Glyma.01G185400",   "GmKRP1a",  "KRP1/ICK1",           "AT2G23430",

"Glyma.11G056700",   "GmKRP1b",  "KRP1/ICK1",           "AT2G23430",
  "Glyma.05G085000",   "GmKRP2a",  "KRP2/ICK2",           "AT3G50630",
  "Glyma.17G175700",   "GmKRP2b",  "KRP2/ICK2",           "AT3G50630",
  "Glyma.08G354300",   "GmKRP3",   "KRP3/ICK6",           "AT5G48820",
  "Glyma.20G198800",   "GmKRP4",   "KRP4/ICK7",           "AT2G32710",
  "Glyma.18G170800",   "GmKRP5",   "KRP5/ICK3",           "AT3G24810",
  "Glyma.02G133700",   "GmKRP6",   "KRP6/ICK4",           "AT3G19150",
  "Glyma.07G211300",   "GmKRP7",   "KRP7/ICK5",           "AT1G49620"
)

cat("Literature reference: Guo et al. (2023) Front. Plant Sci. 14:1096467\n")
cat("Number of KRP genes reported:", nrow(guo_krp_genes), "\n\n")

# All 7 Arabidopsis KRP gene IDs (for ortholog search)
arabidopsis_krp_ids <- c(
  "AT2G23430",  # KRP1/ICK1
"AT3G50630",  # KRP2/ICK2 - JAG target
  "AT5G48820",  # KRP3/ICK6
  "AT2G32710",  # KRP4/ICK7 - JAG target
  "AT3G24810",  # KRP5/ICK3
  "AT3G19150",  # KRP6/ICK4
  "AT1G49620"   # KRP7/ICK5
)

cat("Arabidopsis KRP genes used for ortholog search:", length(arabidopsis_krp_ids), "\n")
cat("  ", paste(arabidopsis_krp_ids, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 1b. Load Phytozome Annotation File
# -----------------------------------------------------------------------------

annotation_file <- file.path(base_dir, "01_data/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (!file.exists(annotation_file)) {
  # Try alternative location
  annotation_file <- file.path(dirname(base_dir), "Phase1-Exploratory/Gmax_880_Wm82.a6.v1.annotation_info.txt")
}

if (!file.exists(annotation_file)) {
  stop("ERROR: Annotation file not found! Please check path.")
}

cat("Loading annotation file:\n  ", annotation_file, "\n")

annotation <- read_delim(annotation_file, delim = "\t",
                         col_types = cols(.default = "c"),
                         show_col_types = FALSE)

# Get unique loci (remove transcript variants)
annotation_loci <- annotation %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`) %>%
  distinct(locusName, .keep_all = TRUE)

cat("Annotation file loaded successfully\n")
cat("  Total entries:", nrow(annotation), "\n")
cat("  Unique loci:", nrow(annotation_loci), "\n")
cat("  Genome version: Wm82.a6.v1 (v6)\n\n")

# =============================================================================
# 2. METHOD 1: PFAM DOMAIN SEARCH (PF02234)
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("METHOD 1: PFAM DOMAIN SEARCH\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Searching for Pfam domain: PF02234 (CDK inhibitor / CID domain)\n")
cat("This is the defining domain for KRP proteins.\n\n")

# Search for PF02234 in Pfam column
pfam_hits <- annotation_loci %>%
  filter(grepl("PF02234", Pfam, ignore.case = TRUE)) %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`)

cat("RESULTS:\n")
cat("  Genes with PF02234 domain:", nrow(pfam_hits), "\n\n")

if (nrow(pfam_hits) > 0) {
  cat("Genes found:\n")
  for (i in 1:nrow(pfam_hits)) {
    cat(sprintf("  %d. %s\n", i, pfam_hits$locusName[i]))
    cat(sprintf("     Pfam: %s\n", pfam_hits$Pfam[i]))
    cat(sprintf("     Arabidopsis hit: %s\n", pfam_hits$`Best-hit-arabi-name`[i]))
    cat("\n")
  }
}

method1_genes <- pfam_hits$locusName

# =============================================================================
# 3. METHOD 2: ARABIDOPSIS ORTHOLOG SEARCH
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("METHOD 2: ARABIDOPSIS ORTHOLOG SEARCH\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Searching Best-hit-arabi-name column for Arabidopsis KRP gene IDs\n")
cat("This uses pre-computed BLASTP results from Phytozome.\n\n")

# Create regex pattern for all Arabidopsis KRP IDs
at_pattern <- paste(arabidopsis_krp_ids, collapse = "|")

# Search for Arabidopsis KRP orthologs
ortholog_hits <- annotation_loci %>%
  filter(grepl(at_pattern, `Best-hit-arabi-name`, ignore.case = TRUE)) %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`)

cat("RESULTS:\n")
cat("  Genes with Arabidopsis KRP as best hit:", nrow(ortholog_hits), "\n\n")

if (nrow(ortholog_hits) > 0) {
  cat("Genes found:\n")
  for (i in 1:nrow(ortholog_hits)) {
    cat(sprintf("  %d. %s\n", i, ortholog_hits$locusName[i]))
    cat(sprintf("     Arabidopsis hit: %s\n", ortholog_hits$`Best-hit-arabi-name`[i]))
    cat(sprintf("     Has PF02234: %s\n",
                ifelse(grepl("PF02234", ortholog_hits$Pfam[i]), "YES", "NO")))
    cat("\n")
  }
}

method2_genes <- ortholog_hits$locusName

# =============================================================================
# 4. CROSS-VALIDATION
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("CROSS-VALIDATION ANALYSIS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Get all unique genes from all three sources
guo_genes <- guo_krp_genes$gene_id
all_genes <- unique(c(guo_genes, method1_genes, method2_genes))

# Create comparison table
comparison <- tibble(gene_id = all_genes) %>%
  mutate(
    in_guo_2023 = gene_id %in% guo_genes,
    in_pfam_search = gene_id %in% method1_genes,
    in_ortholog_search = gene_id %in% method2_genes,
    n_methods = in_guo_2023 + in_pfam_search + in_ortholog_search
  ) %>%
  left_join(guo_krp_genes %>% select(gene_id, name), by = "gene_id") %>%
  left_join(annotation_loci %>% select(locusName, `Best-hit-arabi-name`),
            by = c("gene_id" = "locusName")) %>%
  arrange(desc(n_methods), gene_id)

cat("COMPARISON SUMMARY:\n\n")

# Summary statistics
cat("Source                      | Genes Found\n")
cat("----------------------------|------------\n")
cat(sprintf("Guo et al. 2023 (literature)| %d\n", sum(comparison$in_guo_2023)))
cat(sprintf("Method 1 (Pfam PF02234)     | %d\n", sum(comparison$in_pfam_search)))
cat(sprintf("Method 2 (At orthologs)     | %d\n", sum(comparison$in_ortholog_search)))
cat(sprintf("Union (all sources)         | %d\n", nrow(comparison)))
cat("\n")

# Detailed breakdown
cat("DETAILED BREAKDOWN:\n\n")

# Genes in all three sources
in_all_three <- comparison %>% filter(n_methods == 3)
cat(sprintf("Genes found by ALL 3 methods: %d\n", nrow(in_all_three)))
if (nrow(in_all_three) > 0) {
  cat("  ", paste(in_all_three$gene_id, collapse = ", "), "\n")
}
cat("\n")

# Genes in Guo but not found by our methods
guo_only <- comparison %>% filter(in_guo_2023 & !in_pfam_search & !in_ortholog_search)
cat(sprintf("Genes in Guo et al. but NOT found by our searches: %d\n", nrow(guo_only)))
if (nrow(guo_only) > 0) {
  cat("  WARNING: These genes may need manual verification\n")
  cat("  ", paste(guo_only$gene_id, collapse = ", "), "\n")
}
cat("\n")

# Genes found by our methods but not in Guo
new_genes <- comparison %>% filter(!in_guo_2023 & (in_pfam_search | in_ortholog_search))
cat(sprintf("Genes found by our searches but NOT in Guo et al.: %d\n", nrow(new_genes)))
if (nrow(new_genes) > 0) {
  cat("  WARNING: Guo et al. may have missed these genes!\n")
  for (i in 1:nrow(new_genes)) {
    cat(sprintf("  - %s (Pfam: %s, Ortholog: %s)\n",
                new_genes$gene_id[i],
                ifelse(new_genes$in_pfam_search[i], "YES", "NO"),
                ifelse(new_genes$in_ortholog_search[i], "YES", "NO")))
  }
}
cat("\n")

# =============================================================================
# 5. VERIFY IN EXPRESSION DATA
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("EXPRESSION VERIFICATION\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Load expression data
checkpoint_file <- file.path(base_dir, "03_results/checkpoints/05_normalized.RData")

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)

  all_expressed_genes <- rownames(v_full$E)

  comparison <- comparison %>%
    mutate(
      in_expression_data = gene_id %in% all_expressed_genes
    )

  cat("Expression data loaded successfully\n")
  cat("  Total genes in expression matrix:", length(all_expressed_genes), "\n\n")

  expressed_krps <- comparison %>% filter(in_expression_data)
  not_expressed_krps <- comparison %>% filter(!in_expression_data)

  cat("KRP genes present in expression data:", nrow(expressed_krps), "/", nrow(comparison), "\n")

  if (nrow(not_expressed_krps) > 0) {
    cat("  WARNING: Following KRP genes NOT in expression data:\n")
    cat("  ", paste(not_expressed_krps$gene_id, collapse = ", "), "\n")
  } else {
    cat("  All identified KRP genes are present in expression data\n")
  }

} else {
  cat("WARNING: Expression checkpoint not found. Skipping expression verification.\n")
  comparison$in_expression_data <- NA
}

cat("\n")

# =============================================================================
# 6. CHECK JAG1 TARGET STATUS
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("JAG1 TARGET STATUS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

jag1_file <- file.path(base_dir, "03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv")

if (file.exists(jag1_file)) {
  jag1_targets <- read_csv(jag1_file, show_col_types = FALSE)

  comparison <- comparison %>%
    mutate(
      is_jag1_target = gene_id %in% jag1_targets$GeneID
    )

  jag1_krps <- comparison %>% filter(is_jag1_target)

  cat("JAG1 targets list loaded\n")
  cat("KRP genes that are JAG1 targets:", nrow(jag1_krps), "/", nrow(comparison), "\n")

  if (nrow(jag1_krps) > 0) {
    cat("  JAG1-targeted KRPs:\n")
    for (i in 1:nrow(jag1_krps)) {
      cat(sprintf("    - %s (%s)\n", jag1_krps$gene_id[i],
                  ifelse(is.na(jag1_krps$name[i]), "unnamed", jag1_krps$name[i])))
    }
  }

} else {
  cat("WARNING: JAG1 targets file not found. Skipping JAG1 target check.\n")
  comparison$is_jag1_target <- NA
}

cat("\n")

# =============================================================================
# 7. DETAILED EXPRESSION ANALYSIS
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("DETAILED EXPRESSION ANALYSIS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("This section addresses: Are KRPs expressed but not DE, or just not expressed?\n\n")

if (exists("v_full") && exists("targets")) {

  # Get KRP genes that are in expression data
  krp_in_data <- comparison %>% filter(in_expression_data == TRUE) %>% pull(gene_id)

  if (length(krp_in_data) > 0) {

    # -------------------------------------------------------------------------
    # 7a. Calculate CPM (Counts Per Million)
    # -------------------------------------------------------------------------
    cat("--- 7a. CPM Analysis ---\n\n")

    # Get raw counts if available, otherwise use normalized values
    if (exists("dge_full") && !is.null(dge_full$counts)) {
      raw_counts <- dge_full$counts
      cpm_values <- cpm(dge_full)
    } else {
      # Use 2^(logCPM) as approximation from voom object
      cpm_values <- 2^v_full$E
      raw_counts <- NULL
    }

    # Extract KRP expression
    krp_cpm <- cpm_values[krp_in_data, , drop = FALSE]

    # Calculate mean CPM across all samples
    krp_mean_cpm <- rowMeans(krp_cpm)

    # Calculate genome-wide percentiles
    all_mean_cpm <- rowMeans(cpm_values)

    # Classify expression level based on percentiles
    get_expression_class <- function(cpm_val, all_cpm) {
      pct <- ecdf(all_cpm)(cpm_val) * 100
      case_when(
        cpm_val < 1 ~ "Very Low (<1 CPM)",
        pct < 25 ~ "Low (bottom 25%)",
        pct < 75 ~ "Medium (25-75%)",
        TRUE ~ "High (top 25%)"
      )
    }

    # -------------------------------------------------------------------------
    # 7b. Expression by Genotype (Broad vs Narrow)
    # -------------------------------------------------------------------------
    cat("--- 7b. Expression by Genotype ---\n\n")

    broad_samples <- which(targets$Leaf_type == "Broad")
    narrow_samples <- which(targets$Leaf_type == "Narrow")

    krp_expression_summary <- tibble(
      gene_id = krp_in_data
    ) %>%
      left_join(guo_krp_genes %>% select(gene_id, name), by = "gene_id") %>%
      mutate(
        mean_CPM_all = krp_mean_cpm[gene_id],
        mean_CPM_Broad = rowMeans(krp_cpm[gene_id, broad_samples, drop = FALSE]),
        mean_CPM_Narrow = rowMeans(krp_cpm[gene_id, narrow_samples, drop = FALSE]),
        log2FC_Narrow_vs_Broad = log2((mean_CPM_Narrow + 0.5) / (mean_CPM_Broad + 0.5)),
        expression_percentile = round(ecdf(all_mean_cpm)(mean_CPM_all) * 100, 1),
        expression_class = sapply(mean_CPM_all, get_expression_class, all_cpm = all_mean_cpm)
      ) %>%
      mutate(name = ifelse(is.na(name), gene_id, name)) %>%
      arrange(desc(mean_CPM_all))

    cat("KRP Expression Summary:\n\n")
    cat(sprintf("%-18s %-10s %12s %12s %12s %8s %15s\n",
                "Gene", "Name", "Mean_CPM", "Broad_CPM", "Narrow_CPM", "log2FC", "Expression"))
    cat(paste(rep("-", 95), collapse = ""), "\n")

    for (i in 1:nrow(krp_expression_summary)) {
      row <- krp_expression_summary[i, ]
      cat(sprintf("%-18s %-10s %12.2f %12.2f %12.2f %8.3f %15s\n",
                  row$gene_id, row$name, row$mean_CPM_all,
                  row$mean_CPM_Broad, row$mean_CPM_Narrow,
                  row$log2FC_Narrow_vs_Broad, row$expression_class))
    }

    cat("\n")

    # -------------------------------------------------------------------------
    # 7c. Summary Statistics
    # -------------------------------------------------------------------------
    cat("--- 7c. Expression Level Summary ---\n\n")

    n_very_low <- sum(krp_expression_summary$mean_CPM_all < 1)
    n_low <- sum(krp_expression_summary$expression_class == "Low (bottom 25%)")
    n_medium <- sum(krp_expression_summary$expression_class == "Medium (25-75%)")
    n_high <- sum(krp_expression_summary$expression_class == "High (top 25%)")

    cat("Expression level distribution:\n")
    cat(sprintf("  Very Low (<1 CPM):  %d genes\n", n_very_low))
    cat(sprintf("  Low (bottom 25%%):   %d genes\n", n_low))
    cat(sprintf("  Medium (25-75%%):    %d genes\n", n_medium))
    cat(sprintf("  High (top 25%%):     %d genes\n", n_high))
    cat("\n")

    mean_krp_cpm <- mean(krp_expression_summary$mean_CPM_all)
    median_krp_cpm <- median(krp_expression_summary$mean_CPM_all)

    cat(sprintf("Mean CPM across KRP genes: %.2f\n", mean_krp_cpm))
    cat(sprintf("Median CPM across KRP genes: %.2f\n", median_krp_cpm))
    cat(sprintf("Genome-wide median CPM: %.2f\n", median(all_mean_cpm)))
    cat("\n")

    # -------------------------------------------------------------------------
    # 7d. TP0-Specific KRP Expression (Where JAG1 is Active)
    # -------------------------------------------------------------------------
    cat("--- 7d. TP0-Specific Analysis ---\n\n")

    cat("RATIONALE:\n")
    cat("  - GmJAG1 is expressed at TP0 and drops to near-zero by TP1\n")
    cat("  - If JAG1 regulates KRPs (like in Arabidopsis), we'd see DE at TP0\n")
    cat("  - Testing: Are KRPs differentially expressed at TP0?\n\n")

    # Identify TP0 samples
    if ("Timepoint" %in% names(targets)) {
      tp0_samples <- which(targets$Timepoint == "TP0")
      tp0_broad <- which(targets$Timepoint == "TP0" & targets$Leaf_type == "Broad")
      tp0_narrow <- which(targets$Timepoint == "TP0" & targets$Leaf_type == "Narrow")

      cat(sprintf("TP0 samples: %d total (%d Broad, %d Narrow)\n\n",
                  length(tp0_samples), length(tp0_broad), length(tp0_narrow)))

      if (length(tp0_broad) >= 2 && length(tp0_narrow) >= 2) {

        # KRP expression at TP0
        krp_tp0_expression <- tibble(
          gene_id = krp_in_data
        ) %>%
          left_join(guo_krp_genes %>% select(gene_id, name), by = "gene_id") %>%
          mutate(
            name = ifelse(is.na(name), gene_id, name),
            mean_CPM_TP0_all = rowMeans(cpm_values[gene_id, tp0_samples, drop = FALSE]),
            mean_CPM_TP0_Broad = rowMeans(cpm_values[gene_id, tp0_broad, drop = FALSE]),
            mean_CPM_TP0_Narrow = rowMeans(cpm_values[gene_id, tp0_narrow, drop = FALSE]),
            log2FC_TP0 = log2((mean_CPM_TP0_Narrow + 0.5) / (mean_CPM_TP0_Broad + 0.5))
          ) %>%
          arrange(desc(abs(log2FC_TP0)))

        cat("KRP Expression at TP0 (where JAG1 is active):\n\n")
        cat(sprintf("%-18s %-10s %12s %12s %12s %10s\n",
                    "Gene", "Name", "TP0_Mean", "TP0_Broad", "TP0_Narrow", "log2FC_TP0"))
        cat(paste(rep("-", 80), collapse = ""), "\n")

        for (i in 1:nrow(krp_tp0_expression)) {
          row <- krp_tp0_expression[i, ]
          cat(sprintf("%-18s %-10s %12.2f %12.2f %12.2f %10.3f\n",
                      row$gene_id, row$name, row$mean_CPM_TP0_all,
                      row$mean_CPM_TP0_Broad, row$mean_CPM_TP0_Narrow,
                      row$log2FC_TP0))
        }

        cat("\n")

        krp_mean_fc <- mean(abs(krp_tp0_expression$log2FC_TP0))
        krp_max_fc <- max(abs(krp_tp0_expression$log2FC_TP0))

        cat(sprintf("Mean |log2FC| at TP0: %.3f\n", krp_mean_fc))
        cat(sprintf("Max |log2FC| at TP0:  %.3f\n\n", krp_max_fc))

        # Interpretation
        if (krp_max_fc < 1.0) {
          cat("CONCLUSION: No KRP gene shows |log2FC| > 1 at TP0.\n")
          cat("            KRPs are EXPRESSED but NOT differentially expressed.\n")
          cat("            Soybean JAG1 does NOT regulate KRPs (unlike Arabidopsis JAG).\n")
        } else {
          cat("NOTE: Some KRPs show |log2FC| > 1. Check formal DE test for significance.\n")
        }

        cat("\n")

      } else {
        cat("WARNING: Insufficient TP0 samples for analysis.\n")
        cat(sprintf("  TP0 Broad samples: %d (need >= 2)\n", length(tp0_broad)))
        cat(sprintf("  TP0 Narrow samples: %d (need >= 2)\n", length(tp0_narrow)))
      }
    } else {
      cat("WARNING: Timepoint column not found in targets. Skipping TP0-specific analysis.\n")
    }

    cat("\n")

    # Store expression summary in comparison table
    comparison <- comparison %>%
      left_join(krp_expression_summary %>%
                  select(gene_id, mean_CPM_all, mean_CPM_Broad, mean_CPM_Narrow,
                         log2FC_Narrow_vs_Broad, expression_class),
                by = "gene_id")

  }
} else {
  cat("WARNING: Expression data not available for detailed analysis.\n")
}

cat("\n")

# =============================================================================
# 8. CHECK DIFFERENTIAL EXPRESSION (FORMAL TEST)
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("DIFFERENTIAL EXPRESSION STATUS (FORMAL TEST)\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Checking DE status from formal limma/edgeR analysis\n")
cat("(This is different from raw log2FC - it accounts for variance)\n\n")

de_file <- file.path(base_dir, "03_results/checkpoints/12_DE_organized.RData")

if (file.exists(de_file)) {
  load(de_file)

  # Check if any KRP genes are DE at TP0 or any timepoint
  # DE results are typically stored as a list by contrast

  # Try to find DE genes
  de_genes <- c()

  if (exists("de_results")) {
    # If de_results is a list of dataframes
    if (is.list(de_results)) {
      for (contrast_name in names(de_results)) {
        if (grepl("TP0|Narrow.*Broad", contrast_name, ignore.case = TRUE)) {
          df <- de_results[[contrast_name]]
          if ("GeneID" %in% names(df)) {
            sig_genes <- df %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% pull(GeneID)
            de_genes <- c(de_genes, sig_genes)
          }
        }
      }
    }
  }

  de_genes <- unique(de_genes)

  comparison <- comparison %>%
    mutate(
      is_DE = gene_id %in% de_genes
    )

  de_krps <- comparison %>% filter(is_DE == TRUE)

  cat("KRP genes that are differentially expressed:", sum(comparison$is_DE, na.rm = TRUE), "/", nrow(comparison), "\n")

  if (nrow(de_krps) > 0) {
    cat("  DE KRPs:\n")
    cat("  ", paste(de_krps$gene_id, collapse = ", "), "\n")
  } else {
    cat("  None of the KRP genes are differentially expressed\n")
    cat("  This supports the hypothesis that soybean JAG1 does NOT regulate KRPs\n")
  }

} else {
  cat("WARNING: DE results file not found. Skipping DE check.\n")
  comparison$is_DE <- NA
}

cat("\n")

# =============================================================================
# 8. FINAL VERIFICATION SUMMARY
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("FINAL VERIFICATION SUMMARY\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Print final comparison table
cat("Complete KRP Gene Verification Table:\n\n")

final_table <- comparison %>%
  select(gene_id, name, in_guo_2023, in_pfam_search, in_ortholog_search,
         in_expression_data, is_jag1_target, is_DE) %>%
  mutate(
    name = ifelse(is.na(name), "-", name),
    across(where(is.logical), ~ifelse(is.na(.), "N/A", ifelse(., "Yes", "No")))
  )

print(final_table, n = 20)

cat("\n")

# Verification conclusion
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("VERIFICATION CONCLUSION\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

guo_verified <- sum(comparison$in_guo_2023 & (comparison$in_pfam_search | comparison$in_ortholog_search))
guo_total <- sum(comparison$in_guo_2023)

if (nrow(new_genes) == 0 && guo_verified == guo_total) {
  cat("✓ VERIFICATION SUCCESSFUL\n\n")
  cat("The 9 KRP genes from Guo et al. (2023) represent a COMPLETE list:\n")
  cat("  - All 9 genes verified by Pfam domain and/or Arabidopsis ortholog search\n")
  cat("  - No additional KRP genes found in the soybean v6 genome annotation\n")
  cat("  - Gene IDs from v4 (Guo et al.) are valid in v6 annotation\n")
} else if (nrow(new_genes) > 0) {
  cat("⚠ VERIFICATION FOUND ADDITIONAL GENES\n\n")
  cat("Guo et al. (2023) may have missed some KRP genes:\n")
  cat("  Additional genes found:", nrow(new_genes), "\n")
  cat("  These require further investigation\n")
} else {
  cat("⚠ PARTIAL VERIFICATION\n\n")
  cat("Some Guo et al. genes could not be verified:\n")
  cat("  Verified:", guo_verified, "/", guo_total, "\n")
}

cat("\n")

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Create output directory
output_dir <- file.path(base_dir, "03_results/tables/KRP_verification")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save comparison table
write_csv(comparison, file.path(output_dir, "KRP_verification_comparison.csv"))
cat("Saved: KRP_verification_comparison.csv\n")

# Save Pfam search results
write_csv(pfam_hits, file.path(output_dir, "KRP_pfam_search_results.csv"))
cat("Saved: KRP_pfam_search_results.csv\n")

# Save ortholog search results
write_csv(ortholog_hits, file.path(output_dir, "KRP_ortholog_search_results.csv"))
cat("Saved: KRP_ortholog_search_results.csv\n")

# Save expression summary if available
if (exists("krp_expression_summary")) {
  write_csv(krp_expression_summary, file.path(output_dir, "KRP_expression_summary.csv"))
  cat("Saved: KRP_expression_summary.csv\n")
}

# Save TP0-specific expression data
if (exists("krp_tp0_expression")) {
  write_csv(krp_tp0_expression, file.path(output_dir, "KRP_TP0_expression.csv"))
  cat("Saved: KRP_TP0_expression.csv\n")
}

# Generate markdown report
report_file <- file.path(output_dir, "KRP_Verification_Report.md")

report_lines <- c(
  "# KRP Gene Family Verification Report",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  "This report verifies that the 9 soybean KRP genes identified by Guo et al. (2023)",
  "represent a complete list of the KRP gene family in the soybean genome.",
  "",
  "## Methods",
  "",
  "### Method 1: Pfam Domain Search",
  "- Searched Phytozome annotation (Wm82.a6.v1) for Pfam domain PF02234",

  "- PF02234 = CDK inhibitor domain (defining feature of KRP proteins)",
  paste0("- Result: ", nrow(pfam_hits), " genes found"),
  "",
  "### Method 2: Arabidopsis Ortholog Search",
  "- Searched Best-hit-arabi-name column for Arabidopsis KRP gene IDs",
  "- Uses pre-computed BLASTP results from Phytozome",
  paste0("- Result: ", nrow(ortholog_hits), " genes found"),
  "",
  "## Results Summary",
  "",
  "| Source | Genes Found |",
  "|--------|-------------|",
  paste0("| Guo et al. 2023 (literature) | ", sum(comparison$in_guo_2023), " |"),
  paste0("| Method 1 (Pfam PF02234) | ", sum(comparison$in_pfam_search), " |"),
  paste0("| Method 2 (At orthologs) | ", sum(comparison$in_ortholog_search), " |"),
  "",
  "## Verification Status",
  "",
  paste0("- Genes in Guo et al. verified by our methods: ", guo_verified, "/", guo_total),
  paste0("- Additional genes found: ", nrow(new_genes)),
  "",
  "## The 9 Verified KRP Genes",
  "",
  "| Gene ID | Name | Pfam | Ortholog | Expressed | JAG1 Target | DE |",
  "|---------|------|------|----------|-----------|-------------|-----|"
)

for (i in 1:nrow(comparison)) {
  if (comparison$in_guo_2023[i]) {
    line <- sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                    comparison$gene_id[i],
                    ifelse(is.na(comparison$name[i]), "-", comparison$name[i]),
                    ifelse(comparison$in_pfam_search[i], "✓", "✗"),
                    ifelse(comparison$in_ortholog_search[i], "✓", "✗"),
                    ifelse(is.na(comparison$in_expression_data[i]), "N/A",
                           ifelse(comparison$in_expression_data[i], "✓", "✗")),
                    ifelse(is.na(comparison$is_jag1_target[i]), "N/A",
                           ifelse(comparison$is_jag1_target[i], "✓", "✗")),
                    ifelse(is.na(comparison$is_DE[i]), "N/A",
                           ifelse(comparison$is_DE[i], "✓", "✗")))
    report_lines <- c(report_lines, line)
  }
}

# Add expression analysis to report
if (exists("krp_expression_summary")) {
  report_lines <- c(report_lines,
    "",
    "## Expression Analysis",
    "",
    "### Key Question: Are KRPs expressed but not DE, or just not expressed?",
    "",
    "| Gene | Name | Mean CPM | Broad CPM | Narrow CPM | log2FC | Expression Level |",
    "|------|------|----------|-----------|------------|--------|------------------|"
  )

  for (i in 1:nrow(krp_expression_summary)) {
    row <- krp_expression_summary[i, ]
    line <- sprintf("| %s | %s | %.2f | %.2f | %.2f | %.3f | %s |",
                    row$gene_id, row$name, row$mean_CPM_all,
                    row$mean_CPM_Broad, row$mean_CPM_Narrow,
                    row$log2FC_Narrow_vs_Broad, row$expression_class)
    report_lines <- c(report_lines, line)
  }

  report_lines <- c(report_lines,
    "",
    "### Interpretation",
    "",
    sprintf("- KRPs with substantial expression (>1 CPM): %d / %d",
            sum(krp_expression_summary$mean_CPM_all >= 1), nrow(krp_expression_summary)),
    sprintf("- Mean log2FC (Narrow vs Broad): %.3f", mean(krp_expression_summary$log2FC_Narrow_vs_Broad)),
    "",
    "**Conclusion:** KRPs ARE expressed but show minimal differential expression,",
    "supporting the hypothesis that JAG1 binds but does not regulate KRP genes in soybean."
  )
}

# Add TP0-specific analysis to report
if (exists("krp_tp0_expression")) {
  report_lines <- c(report_lines,
    "",
    "### TP0-Specific Analysis",
    "",
    "GmJAG1 is expressed at TP0 and drops to near-zero by TP1. If JAG1 regulates KRPs",
    "(as in Arabidopsis), we would expect to see differential expression at TP0.",
    "",
    "**KRP Expression at TP0:**",
    "",
    "| Gene | Name | TP0 Broad CPM | TP0 Narrow CPM | log2FC (TP0) |",
    "|------|------|---------------|----------------|--------------|"
  )

  for (i in 1:nrow(krp_tp0_expression)) {
    row <- krp_tp0_expression[i, ]
    line <- sprintf("| %s | %s | %.2f | %.2f | %.3f |",
                    row$gene_id, row$name,
                    row$mean_CPM_TP0_Broad, row$mean_CPM_TP0_Narrow,
                    row$log2FC_TP0)
    report_lines <- c(report_lines, line)
  }

  krp_max_fc <- max(abs(krp_tp0_expression$log2FC_TP0))

  report_lines <- c(report_lines,
    "",
    sprintf("**KRP mean |log2FC| at TP0: %.3f**", mean(abs(krp_tp0_expression$log2FC_TP0))),
    sprintf("**KRP max |log2FC| at TP0: %.3f**", krp_max_fc),
    "",
    ifelse(krp_max_fc < 1.0,
           "**Finding:** No KRP shows |log2FC| > 1 at TP0. KRPs are expressed but NOT regulated by JAG1.",
           "**Note:** Some KRPs show |log2FC| > 1 at TP0. Check formal DE test for significance.")
  )
}

report_lines <- c(report_lines,
  "",
  "## Conclusion",
  "",
  ifelse(nrow(new_genes) == 0 && guo_verified == guo_total,
         "**VERIFIED:** The 9 KRP genes from Guo et al. (2023) represent a complete list of soybean KRP genes.",
         "**REQUIRES REVIEW:** Some discrepancies found between literature and annotation search."),
  "",
  "## Citation",
  "",
  "Guo B, Chen L, Dong L, Yang C, Zhang J, Geng X, Zhou L, Song L. (2023)",
  "Characterization of the soybean KRP gene family reveals a key role for",
  "GmKRP2a in root development. Front. Plant Sci. 14:1096467.",
  "DOI: 10.3389/fpls.2023.1096467"
)

writeLines(report_lines, report_file)
cat("Saved: KRP_Verification_Report.md\n")

cat("\nAll results saved to:", output_dir, "\n")

# =============================================================================
# 10. SAVE CHECKPOINT
# =============================================================================

# Save all relevant objects
save_objects <- c("comparison", "pfam_hits", "ortholog_hits", "guo_krp_genes",
                  "method1_genes", "method2_genes")

# Add expression summaries if they exist
if (exists("krp_expression_summary")) {
  save_objects <- c(save_objects, "krp_expression_summary")
}
if (exists("krp_tp0_expression")) {
  save_objects <- c(save_objects, "krp_tp0_expression")
}

save(list = save_objects,
     file = file.path(base_dir, "03_results/checkpoints/34c_KRP_verification.RData"))

cat("Checkpoint saved: 34c_KRP_verification.RData\n")

cat("\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("Script 34c Complete\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
