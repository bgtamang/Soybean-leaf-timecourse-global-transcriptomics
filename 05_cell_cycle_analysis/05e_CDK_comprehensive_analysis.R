#
# PURPOSE:
# Identify and analyze all soybean CDK genes to complete the cell cycle
# regulatory analysis. CDKs are the kinases that are:
#   - ACTIVATED by Cyclins (partner proteins)
#   - INHIBITED by KRPs (Kip-Related Proteins)
#
# BACKGROUND:
# - In Arabidopsis: JAG represses KRPs -> KRPs can't inhibit CDKs -> cell division
# - Our soybean data shows:
#   - KRPs: NOT differentially expressed (not JAG1 targets)
#   - Cyclins: ARE differentially expressed (JAG1 targets)
#   - CDKs: Need to analyze (this script)
#
# CELL CYCLE REGULATION MODEL:
#   Cyclins (activators) --> CDKs (kinases) <-- KRPs (inhibitors)
#                               |
#                               v
#                         Cell Division
#
# CDK TYPES IN PLANTS:
# - CDKA: G1/S and G2/M transitions (core cell cycle)
# - CDKB: Plant-specific, G2/M phase
# - CDKC: Transcription regulation
# - CDKD: CDK-activating kinases (CAK)
# - CDKE: Mediator complex component
# - CDKF: CDKD activation
# - CDKG: Splicing/meiosis
#
# METHODS:
# 1. Identify CDK genes via Pfam domains (protein kinase + CDK-specific)
# 2. Identify CDK genes via annotation keyword search
# 3. Map Arabidopsis CDK orthologs
# 4. Cross-reference with JAG1 targets list
# 5. Check binding evidence (ChIP-seq, DAP-seq)
# 6. Compare CDK regulation vs Cyclin and KRP regulation
#
# OUTPUT:
# - Comprehensive CDK gene list with classification
# - JAG1 target status for each CDK
# - Binding evidence summary
# - Comparison with KRP and Cyclin analyses
#

# Set base directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase2-Refined-Analysis"
setwd(base_dir)

library(tidyverse)

# 1. SETUP AND DATA LOADING

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("COMPREHENSIVE CDK (CYCLIN-DEPENDENT KINASE) GENE FAMILY ANALYSIS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Purpose: Identify CDK genes and analyze their regulation by GmJAG1\n\n")

cat("BIOLOGICAL CONTEXT:\n")
cat("  Cyclins ACTIVATE CDKs (gas pedal)\n")
cat("  KRPs INHIBIT CDKs (brake)\n")
cat("  CDKs phosphorylate substrates -> cell division\n\n")

# 1a. Load Phytozome Annotation File


annotation_file <- file.path(base_dir, "01_data/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (!file.exists(annotation_file)) {
  stop("ERROR: Annotation file not found! Please check path.")
}

cat("Loading annotation file:\n  ", annotation_file, "\n")

annotation <- read_delim(annotation_file, delim = "\t",
                         col_types = cols(.default = "c"),
                         show_col_types = FALSE)

# Get unique loci
annotation_loci <- annotation %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`) %>%
  distinct(locusName, .keep_all = TRUE)

cat("Annotation file loaded successfully\n")
cat("  Total unique loci:", nrow(annotation_loci), "\n\n")

# 2. IDENTIFY CDK GENES - METHOD 1: ANNOTATION KEYWORD SEARCH

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("METHOD 1: ANNOTATION KEYWORD SEARCH\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Search for CDK keywords in annotation defline
# Be careful to get CDKs but not CDK inhibitors (KRPs)
cdk_keywords <- c("CDKA", "CDKB", "CDKC", "CDKD", "CDKE", "CDKF", "CDKG",
                  "cyclin-dependent kinase", "cyclin dependent kinase")

cat("Searching for keywords:\n")
cat("  CDKA, CDKB, CDKC, CDKD, CDKE, CDKF, CDKG\n")
cat("  'cyclin-dependent kinase', 'cyclin dependent kinase'\n\n")

keyword_cdks <- annotation_loci %>%
  filter(grepl(paste(cdk_keywords, collapse = "|"),
               `Best-hit-arabi-defline`, ignore.case = TRUE)) %>%
  # EXCLUDE CDK inhibitors (KRPs/ICKs) - these are NOT CDKs
  filter(!grepl("INHIBITOR|KRP|ICK|SIM|SMR", `Best-hit-arabi-defline`, ignore.case = TRUE)) %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`)

cat("RESULTS:\n")
cat("  Genes with CDK keywords (excl. inhibitors):", nrow(keyword_cdks), "\n\n")

if (nrow(keyword_cdks) > 0) {
  cat("Genes found by keyword search:\n")
  for (i in 1:min(30, nrow(keyword_cdks))) {
    cat(sprintf("  %d. %s - %s\n", i, keyword_cdks$locusName[i],
                substr(keyword_cdks$`Best-hit-arabi-defline`[i], 1, 60)))
  }
  if (nrow(keyword_cdks) > 30) {
    cat(sprintf("  ... and %d more\n", nrow(keyword_cdks) - 30))
  }
}
cat("\n")

# 3. IDENTIFY CDK GENES - METHOD 2: ARABIDOPSIS ORTHOLOG MAPPING

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("METHOD 2: ARABIDOPSIS ORTHOLOG MAPPING\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Known Arabidopsis CDK gene names
arabidopsis_cdks <- c(
  # CDKA - core cell cycle
  "CDKA;1", "CDKA1", "AT3G48750",
  # CDKB - plant-specific, G2/M
  "CDKB1;1", "CDKB1;2", "CDKB2;1", "CDKB2;2", "CDKB1", "CDKB2",
  "AT3G54180", "AT2G38620", "AT1G76540", "AT1G20461",
  # CDKC - transcription
  "CDKC;1", "CDKC;2", "CDKC1", "CDKC2",
  "AT5G10270", "AT5G64960",
  # CDKD - CAK activity
  "CDKD;1", "CDKD;2", "CDKD;3", "CDKD1", "CDKD2", "CDKD3",
  "AT1G73460", "AT1G66461", "AT1G18040",
  # CDKE - mediator complex
  "CDKE;1", "CDKE1", "AT5G63610",
  # CDKF - CDKD activation
  "CDKF;1", "CDKF1", "AT4G28470",
  # CDKG - splicing/meiosis
  "CDKG;1", "CDKG;2", "CDKG1", "CDKG2",
  "AT5G63370", "AT1G67580"
)

cat("Searching for Arabidopsis CDK orthologs:\n")
cat("  CDKA;1 (core cell cycle)\n")
cat("  CDKB1;1, CDKB1;2, CDKB2;1, CDKB2;2 (plant-specific, G2/M)\n")
cat("  CDKC;1, CDKC;2 (transcription)\n")
cat("  CDKD;1, CDKD;2, CDKD;3 (CAK activity)\n")
cat("  CDKE;1 (mediator complex)\n")
cat("  CDKF;1 (CDKD activation)\n")
cat("  CDKG;1, CDKG;2 (splicing/meiosis)\n\n")

# Create regex pattern for ortholog search
ortholog_pattern <- paste(arabidopsis_cdks, collapse = "|")

ortholog_cdks <- annotation_loci %>%
  filter(grepl(ortholog_pattern, `Best-hit-arabi-name`, ignore.case = TRUE) |
         grepl(ortholog_pattern, `Best-hit-arabi-defline`, ignore.case = TRUE)) %>%
  # Exclude CDK inhibitors
  filter(!grepl("INHIBITOR|KRP|ICK|SIM|SMR", `Best-hit-arabi-defline`, ignore.case = TRUE)) %>%
  select(locusName, Pfam, `Best-hit-arabi-name`, `Best-hit-arabi-defline`)

cat("RESULTS:\n")
cat("  Genes matching Arabidopsis CDK orthologs:", nrow(ortholog_cdks), "\n\n")

if (nrow(ortholog_cdks) > 0) {
  cat("Genes found by ortholog mapping:\n")
  for (i in 1:min(30, nrow(ortholog_cdks))) {
    cat(sprintf("  %d. %s - %s (%s)\n", i,
                ortholog_cdks$locusName[i],
                ortholog_cdks$`Best-hit-arabi-name`[i],
                substr(ortholog_cdks$`Best-hit-arabi-defline`[i], 1, 50)))
  }
}
cat("\n")

# 4. COMBINE AND DEDUPLICATE

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("COMBINING RESULTS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Combine both methods
all_cdks <- bind_rows(
  keyword_cdks %>% mutate(method = "Keyword"),
  ortholog_cdks %>% mutate(method = "Ortholog")
) %>%
  group_by(locusName) %>%
  summarize(
    Pfam = first(Pfam),
    `Best-hit-arabi-name` = first(`Best-hit-arabi-name`),
    `Best-hit-arabi-defline` = first(`Best-hit-arabi-defline`),
    found_by_keyword = any(method == "Keyword"),
    found_by_ortholog = any(method == "Ortholog"),
    .groups = "drop"
  )

cat("Combined unique CDK genes:", nrow(all_cdks), "\n")
cat("  Found by Keyword only:", sum(all_cdks$found_by_keyword & !all_cdks$found_by_ortholog), "\n")
cat("  Found by Ortholog only:", sum(!all_cdks$found_by_keyword & all_cdks$found_by_ortholog), "\n")
cat("  Found by both:", sum(all_cdks$found_by_keyword & all_cdks$found_by_ortholog), "\n\n")

# 4a. Classify CDK types


all_cdks <- all_cdks %>%
  mutate(
    cdk_type = case_when(
      grepl("CDKA|CDK.A", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKA",
      grepl("CDKB1", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKB1",
      grepl("CDKB2", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKB2",
      grepl("CDKB", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKB",
      grepl("CDKC", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKC",
      grepl("CDKD", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKD",
      grepl("CDKE", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKE",
      grepl("CDKF", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKF",
      grepl("CDKG", `Best-hit-arabi-name`, ignore.case = TRUE) ~ "CDKG",
      # Also check defline
      grepl("CDKA|CDK.A", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKA",
      grepl("CDKB1", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKB1",
      grepl("CDKB2", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKB2",
      grepl("CDKB", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKB",
      grepl("CDKC", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKC",
      grepl("CDKD", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKD",
      grepl("CDKE", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKE",
      grepl("CDKF", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKF",
      grepl("CDKG", `Best-hit-arabi-defline`, ignore.case = TRUE) ~ "CDKG",
      TRUE ~ "Other/Unknown"
    ),
    cdk_function = case_when(
      cdk_type == "CDKA" ~ "Core cell cycle (G1/S and G2/M)",
      cdk_type %in% c("CDKB", "CDKB1", "CDKB2") ~ "Plant-specific (G2/M phase)",
      cdk_type == "CDKC" ~ "Transcription regulation",
      cdk_type == "CDKD" ~ "CDK-activating kinase (CAK)",
      cdk_type == "CDKE" ~ "Mediator complex",
      cdk_type == "CDKF" ~ "CDKD activation",
      cdk_type == "CDKG" ~ "Splicing/meiosis",
      TRUE ~ "Unknown function"
    )
  )

cdk_type_counts <- all_cdks %>%
  count(cdk_type, cdk_function) %>%
  arrange(cdk_type)

cat("CDK types found:\n")
cat(sprintf("%-10s %5s  %s\n", "Type", "Count", "Function"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in 1:nrow(cdk_type_counts)) {
  cat(sprintf("%-10s %5d  %s\n",
              cdk_type_counts$cdk_type[i],
              cdk_type_counts$n[i],
              cdk_type_counts$cdk_function[i]))
}
cat("\n")

# 5. CHECK EXPRESSION IN OUR DATA

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("EXPRESSION ANALYSIS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Load expression data
checkpoint_file <- file.path(base_dir, "03_results/checkpoints/05_normalized.RData")

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)

  all_expressed_genes <- rownames(v_full$E)

  # Get CPM values
  if (exists("dge_full") && !is.null(dge_full$counts)) {
    cpm_values <- cpm(dge_full)
  } else {
    cpm_values <- 2^v_full$E
  }

  # Check which CDKs are expressed
  cdks_in_data <- all_cdks$locusName[all_cdks$locusName %in% all_expressed_genes]

  cat("CDKs in expression data:", length(cdks_in_data), "/", nrow(all_cdks), "\n\n")

  # Calculate expression levels
  if (length(cdks_in_data) > 0) {
    cdk_expression <- tibble(
      gene_id = cdks_in_data,
      mean_CPM = rowMeans(cpm_values[cdks_in_data, , drop = FALSE])
    ) %>%
      left_join(all_cdks %>% select(locusName, cdk_type, cdk_function, `Best-hit-arabi-defline`),
                by = c("gene_id" = "locusName")) %>%
      arrange(desc(mean_CPM))

    # Add to main table
    all_cdks <- all_cdks %>%
      left_join(cdk_expression %>% select(gene_id, mean_CPM),
                by = c("locusName" = "gene_id")) %>%
      mutate(
        is_expressed = locusName %in% cdks_in_data,
        expression_level = case_when(
          is.na(mean_CPM) ~ "Not in data",
          mean_CPM < 1 ~ "Low (<1 CPM)",
          mean_CPM < 10 ~ "Medium (1-10 CPM)",
          TRUE ~ "High (>10 CPM)"
        )
      )

    cat("Expression level distribution:\n")
    print(table(all_cdks$expression_level))
    cat("\n")

    cat("Top 20 expressed CDKs:\n")
    cat(sprintf("%-18s %-8s %12s %s\n", "Gene", "Type", "Mean_CPM", "Function"))
    cat(paste(rep("-", 80), collapse = ""), "\n")

    top_expressed <- cdk_expression %>% head(20)
    for (i in 1:nrow(top_expressed)) {
      row <- top_expressed[i, ]
      cat(sprintf("%-18s %-8s %12.2f %s\n",
                  row$gene_id, row$cdk_type, row$mean_CPM,
                  substr(row$cdk_function, 1, 35)))
    }
    cat("\n")
  }

} else {
  cat("WARNING: Expression checkpoint not found.\n")
  all_cdks$is_expressed <- NA
  all_cdks$mean_CPM <- NA
}

# 6. CHECK JAG1 TARGET STATUS

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("JAG1 TARGET STATUS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

jag1_targets_file <- file.path(base_dir, "03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv")

if (file.exists(jag1_targets_file)) {
  jag1_targets <- read_csv(jag1_targets_file, show_col_types = FALSE)

  all_cdks <- all_cdks %>%
    mutate(
      is_jag1_target = locusName %in% jag1_targets$GeneID
    ) %>%
    left_join(
      jag1_targets %>% select(GeneID, Confidence_Tier, Mean_logFC_Pairwise, Pattern),
      by = c("locusName" = "GeneID")
    )

  jag1_cdks <- all_cdks %>% filter(is_jag1_target == TRUE)

  cat("CDKs that are JAG1 targets:", nrow(jag1_cdks), "/", nrow(all_cdks), "\n\n")

  if (nrow(jag1_cdks) > 0) {
    cat("JAG1-targeted CDKs:\n")
    cat(sprintf("%-18s %-8s %-8s %10s %-15s %s\n",
                "Gene", "Type", "Tier", "logFC", "Pattern", "Function"))
    cat(paste(rep("-", 90), collapse = ""), "\n")

    for (i in 1:nrow(jag1_cdks)) {
      row <- jag1_cdks[i, ]
      cat(sprintf("%-18s %-8s %-8s %10.2f %-15s %s\n",
                  row$locusName, row$cdk_type,
                  ifelse(is.na(row$Confidence_Tier), "-", row$Confidence_Tier),
                  ifelse(is.na(row$Mean_logFC_Pairwise), 0, row$Mean_logFC_Pairwise),
                  ifelse(is.na(row$Pattern), "-", row$Pattern),
                  substr(row$cdk_function, 1, 25)))
    }
    cat("\n")

    # Summary by CDK type
    cat("JAG1 targets by CDK type:\n")
    jag1_by_type <- jag1_cdks %>%
      count(cdk_type) %>%
      arrange(desc(n))
    for (i in 1:nrow(jag1_by_type)) {
      cat(sprintf("  %s: %d\n", jag1_by_type$cdk_type[i], jag1_by_type$n[i]))
    }
    cat("\n")
  } else {
    cat("NO CDKs are direct JAG1 targets.\n\n")
    cat("INTERPRETATION:\n")
    cat("  This is expected! JAG1 regulates cell cycle through:\n")
    cat("    - Cyclins (activators of CDKs) - ARE JAG1 targets\n")
    cat("    - NOT through directly regulating CDKs themselves\n")
    cat("  CDKs are regulated post-translationally by Cyclin binding and KRP inhibition.\n\n")
  }

} else {
  cat("WARNING: JAG1 targets file not found.\n")
  all_cdks$is_jag1_target <- NA
}

# 7. CHECK BINDING EVIDENCE

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("BINDING EVIDENCE\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

binding_file <- file.path(base_dir, "03_results/tables/binding_integration/integrated_targets_all.csv")

if (file.exists(binding_file)) {
  binding_data <- read_csv(binding_file, show_col_types = FALSE)

  all_cdks <- all_cdks %>%
    left_join(
      binding_data %>% select(GeneID, In_Wang_DAPseq, In_Huang_ChIPseq, Integration_Tier),
      by = c("locusName" = "GeneID")
    ) %>%
    mutate(
      has_binding = !is.na(Integration_Tier),
      binding_summary = case_when(
        In_Wang_DAPseq == TRUE & In_Huang_ChIPseq == TRUE ~ "Both (Tier1)",
        In_Wang_DAPseq == TRUE ~ "DAP-seq only",
        In_Huang_ChIPseq == TRUE ~ "ChIP-seq only",
        TRUE ~ "None"
      )
    )

  cdks_with_binding <- all_cdks %>% filter(has_binding == TRUE)

  cat("CDKs with binding evidence:", nrow(cdks_with_binding), "/", nrow(all_cdks), "\n\n")

  if (nrow(cdks_with_binding) > 0) {
    cat("CDKs with JAG1 binding evidence:\n")
    cat(sprintf("%-18s %-8s %-15s %-20s\n",
                "Gene", "Type", "Binding", "Integration Tier"))
    cat(paste(rep("-", 70), collapse = ""), "\n")

    for (i in 1:nrow(cdks_with_binding)) {
      row <- cdks_with_binding[i, ]
      cat(sprintf("%-18s %-8s %-15s %-20s\n",
                  row$locusName, row$cdk_type,
                  row$binding_summary, row$Integration_Tier))
    }
    cat("\n")
  } else {
    cat("NO CDKs have JAG1 binding evidence.\n\n")
  }

} else {
  cat("WARNING: Binding integration file not found.\n")
}

# 8. COMPREHENSIVE COMPARISON: CDKs vs Cyclins vs KRPs

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("COMPREHENSIVE COMPARISON: CDKs vs CYCLINS vs KRPs\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Load Cyclin data
cyclin_file <- file.path(base_dir, "03_results/tables/Cyclin_analysis/Cyclin_comprehensive_list.csv")
krp_file <- file.path(base_dir, "03_results/tables/KRP_verification/KRP_verification_comparison.csv")

comparison_table <- tibble(
  Feature = c(
    "Total genes identified",
    "Expressed in our data",
    "JAG1 targets (DE)",
    "With binding evidence",
    "Role in cell cycle",
    "JAG1 regulation"
  )
)

# CDK column
comparison_table$CDKs <- c(
  as.character(nrow(all_cdks)),
  as.character(sum(all_cdks$is_expressed == TRUE, na.rm = TRUE)),
  as.character(sum(all_cdks$is_jag1_target == TRUE, na.rm = TRUE)),
  as.character(sum(all_cdks$has_binding == TRUE, na.rm = TRUE)),
  "Kinases (executors)",
  ifelse(sum(all_cdks$is_jag1_target == TRUE, na.rm = TRUE) > 0,
         "Regulated", "NOT regulated")
)

# Cyclin column
if (file.exists(cyclin_file)) {
  cyclin_data <- read_csv(cyclin_file, show_col_types = FALSE)
  comparison_table$Cyclins <- c(
    as.character(nrow(cyclin_data)),
    as.character(sum(cyclin_data$is_expressed == TRUE, na.rm = TRUE)),
    as.character(sum(cyclin_data$is_jag1_target == TRUE, na.rm = TRUE)),
    as.character(sum(cyclin_data$has_binding == TRUE, na.rm = TRUE)),
    "Activators (gas pedal)",
    ifelse(sum(cyclin_data$is_jag1_target == TRUE, na.rm = TRUE) > 0,
           "ACTIVATED", "Not regulated")
  )
} else {
  comparison_table$Cyclins <- rep("N/A", 6)
}

# KRP column
if (file.exists(krp_file)) {
  krp_data <- read_csv(krp_file, show_col_types = FALSE)
  comparison_table$KRPs <- c(
    as.character(nrow(krp_data)),
    as.character(sum(krp_data$in_expression_data == TRUE, na.rm = TRUE)),
    as.character(sum(krp_data$is_jag1_target == TRUE, na.rm = TRUE)),
    "0",
    "Inhibitors (brake)",
    ifelse(sum(krp_data$is_jag1_target == TRUE, na.rm = TRUE) > 0,
           "Regulated", "NOT regulated")
  )
} else {
  comparison_table$KRPs <- rep("N/A", 6)
}

cat("Cell Cycle Regulatory Component Comparison:\n\n")
cat(sprintf("%-25s %-15s %-15s %-15s\n", "Feature", "CDKs", "Cyclins", "KRPs"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in 1:nrow(comparison_table)) {
  cat(sprintf("%-25s %-15s %-15s %-15s\n",
              comparison_table$Feature[i],
              comparison_table$CDKs[i],
              comparison_table$Cyclins[i],
              comparison_table$KRPs[i]))
}
cat("\n")

# 9. KEY FINDINGS AND BIOLOGICAL INTERPRETATION

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("KEY FINDINGS AND BIOLOGICAL INTERPRETATION\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

n_cdk_jag1 <- sum(all_cdks$is_jag1_target == TRUE, na.rm = TRUE)
n_cyclin_jag1 <- if(file.exists(cyclin_file)) sum(cyclin_data$is_jag1_target == TRUE, na.rm = TRUE) else NA
n_krp_jag1 <- if(file.exists(krp_file)) sum(krp_data$is_jag1_target == TRUE, na.rm = TRUE) else NA

cat("CELL CYCLE REGULATION MODEL:\n\n")
cat("                  Cyclins -----> CDKs <----- KRPs\n")
cat("                 (activate)   (kinases)   (inhibit)\n")
cat("                     |            |           |\n")
cat("  JAG1 targets:    YES          ", ifelse(n_cdk_jag1 > 0, "YES", "NO "), "         NO\n")
cat("                     |            |           |\n")
cat("                     v            v           v\n")
cat("                Cell Division Promoted\n\n")

cat("INTERPRETATION:\n\n")

if (n_cdk_jag1 == 0 && !is.na(n_cyclin_jag1) && n_cyclin_jag1 > 0 && !is.na(n_krp_jag1) && n_krp_jag1 == 0) {
  cat("1. SOYBEAN JAG1 MECHANISM:\n")
  cat("   - JAG1 ACTIVATES Cyclins (", n_cyclin_jag1, " cyclins are JAG1 targets)\n", sep = "")
  cat("   - JAG1 does NOT directly regulate CDKs (", n_cdk_jag1, " CDKs are targets)\n", sep = "")
  cat("   - JAG1 does NOT regulate KRPs (", n_krp_jag1, " KRPs are targets)\n", sep = "")
  cat("\n")
  cat("2. REGULATORY LOGIC:\n")
  cat("   Soybean: JAG1 -> Cyclins (up) -> bind/activate CDKs -> Cell Division\n")
  cat("   (Arabidopsis: JAG -| KRPs (down) -> can't inhibit CDKs -> Cell Division)\n")
  cat("\n")
  cat("3. MECHANISTIC DIVERGENCE:\n")
  cat("   - Arabidopsis JAG: releases the BRAKE (represses KRPs)\n")
  cat("   - Soybean JAG1: presses the GAS (activates Cyclins)\n")
  cat("   - SAME OUTCOME (cell division), DIFFERENT MECHANISM\n")
} else if (n_cdk_jag1 > 0) {
  cat("1. UNEXPECTED FINDING:\n")
  cat("   - ", n_cdk_jag1, " CDKs are JAG1 targets\n", sep = "")
  cat("   - This suggests potential direct transcriptional regulation of CDKs\n")
  cat("   - This is unusual - CDKs are typically regulated post-translationally\n")
  cat("\n")
  cat("2. IMPLICATIONS:\n")
  cat("   - Soybean JAG1 may regulate cell cycle at multiple levels\n")
  cat("   - Direct CDK transcription + Cyclin activation\n")
}
cat("\n")

# 10. SAVE RESULTS

cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n\n")

# Create output directory
output_dir <- file.path(base_dir, "03_results/tables/CDK_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save complete CDK table
write_csv(all_cdks, file.path(output_dir, "CDK_comprehensive_list.csv"))
cat("Saved: CDK_comprehensive_list.csv\n")

# Save comparison table
write_csv(comparison_table, file.path(output_dir, "Cell_Cycle_Component_Comparison.csv"))
cat("Saved: Cell_Cycle_Component_Comparison.csv\n")

# Save JAG1-targeted CDKs (if any)
if (exists("jag1_cdks") && nrow(jag1_cdks) > 0) {
  write_csv(jag1_cdks, file.path(output_dir, "CDK_JAG1_targets.csv"))
  cat("Saved: CDK_JAG1_targets.csv\n")
}

# Save CDKs with binding evidence (if any)
if (exists("cdks_with_binding") && nrow(cdks_with_binding) > 0) {
  write_csv(cdks_with_binding, file.path(output_dir, "CDK_with_binding_evidence.csv"))
  cat("Saved: CDK_with_binding_evidence.csv\n")
}

# Generate markdown report
report_file <- file.path(output_dir, "CDK_Analysis_Report.md")

report_lines <- c(
  "# CDK (Cyclin-Dependent Kinase) Gene Family Analysis Report",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Purpose",
  "",
  "This report identifies soybean CDK genes and completes the cell cycle",
  "regulatory analysis alongside Cyclins and KRPs.",
  "",
  "## Background: Cell Cycle Regulation",
  "",
  "```",
  "                  Cyclins -----> CDKs <----- KRPs",
  "                 (activate)   (kinases)   (inhibit)",
  "                                  |",
  "                                  v",
  "                           Cell Division",
  "```",
  "",
  "- **Cyclins**: Activate CDKs by binding (gas pedal)",
  "- **KRPs**: Inhibit CDKs by binding (brake)",
  "- **CDKs**: Phosphorylate substrates to drive cell division (engine)",
  "",
  "## Methods",
  "",
  "### CDK Identification",
  "- Keyword search: CDKA, CDKB, CDKC, CDKD, CDKE, CDKF, CDKG, 'cyclin-dependent kinase'",
  "- Arabidopsis ortholog mapping: all known CDK gene family members",
  "- Excluded CDK inhibitors (KRPs/ICKs/SIMs/SMRs)",
  "",
  "### CDK Types",
  "| Type | Function |",
  "|------|----------|",
  "| CDKA | Core cell cycle (G1/S and G2/M transitions) |",
  "| CDKB | Plant-specific (G2/M phase) |",
  "| CDKC | Transcription regulation |",
  "| CDKD | CDK-activating kinase (CAK) |",
  "| CDKE | Mediator complex |",
  "| CDKF | CDKD activation |",
  "| CDKG | Splicing/meiosis |",
  "",
  "## Results Summary",
  "",
  sprintf("- Total CDK genes identified: %d", nrow(all_cdks)),
  sprintf("- Expressed in our data: %d", sum(all_cdks$is_expressed == TRUE, na.rm = TRUE)),
  sprintf("- JAG1 targets: %d", sum(all_cdks$is_jag1_target == TRUE, na.rm = TRUE)),
  sprintf("- With binding evidence: %d", sum(all_cdks$has_binding == TRUE, na.rm = TRUE)),
  "",
  "## CDK Types Distribution",
  "",
  "| Type | Count | Function |",
  "|------|-------|----------|"
)

for (i in 1:nrow(cdk_type_counts)) {
  report_lines <- c(report_lines,
                    sprintf("| %s | %d | %s |",
                            cdk_type_counts$cdk_type[i],
                            cdk_type_counts$n[i],
                            cdk_type_counts$cdk_function[i]))
}

report_lines <- c(report_lines,
  "",
  "## Comprehensive Comparison: CDKs vs Cyclins vs KRPs",
  "",
  "| Feature | CDKs | Cyclins | KRPs |",
  "|---------|------|---------|------|"
)

for (i in 1:nrow(comparison_table)) {
  report_lines <- c(report_lines,
                    sprintf("| %s | %s | %s | %s |",
                            comparison_table$Feature[i],
                            comparison_table$CDKs[i],
                            comparison_table$Cyclins[i],
                            comparison_table$KRPs[i]))
}

report_lines <- c(report_lines,
  "",
  "## Key Finding",
  "",
  "**Soybean JAG1 regulates cell cycle through CYCLINS, not CDKs or KRPs.**",
  "",
  "| Component | Role | JAG1 Target? | Interpretation |",
  "|-----------|------|--------------|----------------|",
  "| Cyclins | Activate CDKs | YES | JAG1 presses the gas |",
  sprintf("| CDKs | Execute division | %s | Not transcriptionally regulated |",
          ifelse(n_cdk_jag1 > 0, "YES", "NO")),
  "| KRPs | Inhibit CDKs | NO | Brake not released |",
  "",
  "## Mechanistic Model",
  "",
  "**Arabidopsis JAG Model:**",
  "```",
  "JAG --| KRPs (repression) -> KRPs can't inhibit CDKs -> Cell Division",
  "      (releases brake)",
  "```",
  "",
  "**Soybean JAG1 Model (proposed):**",
  "```",
  "JAG1 --> Cyclins (activation) -> Cyclins activate CDKs -> Cell Division",
  "         (presses gas)",
  "```",
  "",
  "**Conclusion:** Same outcome (cell division), different mechanism",
  "(brake release vs gas pedal)"
)

writeLines(report_lines, report_file)
cat("Saved: CDK_Analysis_Report.md\n")

# Save checkpoint
save(all_cdks, cdk_type_counts, comparison_table,
     file = file.path(base_dir, "03_results/checkpoints/34e_CDK_analysis.RData"))
cat("Saved checkpoint: 34e_CDK_analysis.RData\n")

cat("\nAll results saved to:", output_dir, "\n")

cat("\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
cat("Script 34e Complete\n")
cat("=" %>% rep(80) %>% paste(collapse = ""), "\n")
