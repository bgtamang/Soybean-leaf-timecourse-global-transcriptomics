# Script 42: Gene ID Conversion (Wm82 → ZH13) for Single-Cell Integration
# Converts gene IDs using SoyBase pangene correspondence table (Glycine.pan5.MKRS)
# Required before running the Python analysis scripts on Fan et al. 2025 snRNA-seq data

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 42: GENE ID CONVERSION (Wm82 <-> ZH13)\n")
cat("  For single-cell transcriptomics integration\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

sci_dir  <- "single_cell_integration"
data_dir <- file.path(sci_dir, "data")
res_dir  <- file.path(sci_dir, "results")
dir.create(file.path(data_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

# ===== LOAD PACKAGES =====

library(readxl)

# ===== PARAMETERS =====

# SoyBase pangene table URL
PANGENE_URL <- "https://data.soybase.org/Glycine/GENUS/pangenes/Glycine.pan5.MKRS/Glycine.pan5.MKRS.table_ref_lines.tsv.gz"

# Fan et al. supplementary (shoot apex cluster markers)
fan_supp <- file.path(base_dir, "Literature/Fan_et_al_2025/mmc6.xlsx")

cat("Pangene source: Glycine.pan5.MKRS\n")
cat("Fan et al. supplementary:", fan_supp, "\n\n")


# ===== SECTION 1: DOWNLOAD PANGENE TABLE =====

cat("========================================\n")
cat("SECTION 1: DOWNLOAD PANGENE TABLE\n")
cat("========================================\n\n")

pangene_file <- file.path(data_dir, "pangene_ref_lines.tsv.gz")

if (!file.exists(pangene_file)) {
  cat("Downloading SoyBase pangene table...\n")
  download.file(PANGENE_URL, pangene_file, mode = "wb")
  cat("Downloaded.\n\n")
} else {
  cat("Pangene table already exists, skipping download.\n\n")
}


# ===== SECTION 2: PARSE PANGENE TABLE =====

cat("========================================\n")
cat("SECTION 2: PARSE PANGENE TABLE\n")
cat("========================================\n\n")

pangene <- read.delim(
  gzfile(pangene_file),
  header = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cat("Pangene table:", nrow(pangene), "rows x", ncol(pangene), "columns\n")

# Find the Wm82 a2.v1 and ZH13 columns
wm82_cols <- grep("Wm82", names(pangene), value = TRUE)
zh13_cols <- grep("Zh13|ZH13", names(pangene), value = TRUE)
wm82_col <- wm82_cols[grep("gnm2", wm82_cols)][1]  # Wm82.gnm2.ann1 = Glyma.XXG format
zh13_col <- zh13_cols[1]
cat("Wm82 column:", wm82_col, "\n")
cat("ZH13 column:", zh13_col, "\n")

# Build clean lookup
pangene_id_col <- names(pangene)[1]
lookup <- pangene[, c(pangene_id_col, wm82_col, zh13_col)]
names(lookup) <- c("pangene", "Wm82", "Zh13")
lookup <- lookup[lookup$Wm82 != "" & lookup$Zh13 != "" &
                 !is.na(lookup$Wm82) & !is.na(lookup$Zh13), ]
lookup <- lookup[lookup$Wm82 != "NONE" & lookup$Zh13 != "NONE", ]

# Flag 1:1 vs multi-gene rows
lookup$is_1to1 <- !grepl(",", lookup$Wm82) & !grepl(",", lookup$Zh13)

cat("\nRows with both Wm82 + ZH13:", nrow(lookup), "\n")
cat("Clean 1:1 mappings:", sum(lookup$is_1to1), "\n")
cat("Multi-gene rows:", sum(!lookup$is_1to1), "\n\n")


# ===== SECTION 3: CONVERSION FUNCTIONS =====

cat("========================================\n")
cat("SECTION 3: DEFINE CONVERSION FUNCTIONS\n")
cat("========================================\n\n")

convert_wm82_to_zh13 <- function(glyma_ids, lookup_table) {
  lookup_1to1 <- lookup_table[lookup_table$is_1to1, ]
  result <- data.frame(
    Wm82 = glyma_ids,
    Zh13 = NA_character_,
    pangene = NA_character_,
    mapping_type = NA_character_,
    stringsAsFactors = FALSE
  )

  # Pass 1: exact 1:1
  idx <- match(glyma_ids, lookup_1to1$Wm82)
  found <- !is.na(idx)
  result$Zh13[found] <- lookup_1to1$Zh13[idx[found]]
  result$pangene[found] <- lookup_1to1$pangene[idx[found]]
  result$mapping_type[found] <- "1:1"

  # Pass 2: grep multi-gene rows for remainder
  remaining <- which(!found)
  if (length(remaining) > 0) {
    lookup_multi <- lookup_table[!lookup_table$is_1to1, ]
    for (j in remaining) {
      gene <- glyma_ids[j]
      hit <- lookup_multi[grepl(gene, lookup_multi$Wm82, fixed = TRUE), ]
      if (nrow(hit) > 0) {
        result$Zh13[j] <- hit$Zh13[1]
        result$pangene[j] <- hit$pangene[1]
        result$mapping_type[j] <- "multi"
      }
    }
  }

  return(result)
}

convert_zh13_to_wm82 <- function(zh13_ids, lookup_table) {
  lookup_1to1 <- lookup_table[lookup_table$is_1to1, ]
  result <- data.frame(
    Zh13 = zh13_ids,
    Wm82 = NA_character_,
    pangene = NA_character_,
    mapping_type = NA_character_,
    stringsAsFactors = FALSE
  )

  idx <- match(zh13_ids, lookup_1to1$Zh13)
  found <- !is.na(idx)
  result$Wm82[found] <- lookup_1to1$Wm82[idx[found]]
  result$pangene[found] <- lookup_1to1$pangene[idx[found]]
  result$mapping_type[found] <- "1:1"

  remaining <- which(!found)
  if (length(remaining) > 0) {
    lookup_multi <- lookup_table[!lookup_table$is_1to1, ]
    for (j in remaining) {
      gene <- zh13_ids[j]
      hit <- lookup_multi[grepl(gene, lookup_multi$Zh13, fixed = TRUE), ]
      if (nrow(hit) > 0) {
        result$Wm82[j] <- hit$Wm82[1]
        result$pangene[j] <- hit$pangene[1]
        result$mapping_type[j] <- "multi"
      }
    }
  }

  return(result)
}

cat("Defined: convert_wm82_to_zh13(), convert_zh13_to_wm82()\n\n")


# ===== SECTION 4: VERIFY WITH KNOWN GENES =====

cat("========================================\n")
cat("SECTION 4: VERIFICATION\n")
cat("========================================\n\n")

test_genes <- c("Glyma.20G116200", "Glyma.10G273800")
test_names <- c("JAG1", "JAG2")

for (i in seq_along(test_genes)) {
  gene <- test_genes[i]
  hit <- lookup[grepl(gene, lookup$Wm82), ]
  if (nrow(hit) > 0) {
    cat(sprintf("  %s (%s) -> %s  [%s]\n",
                test_names[i], gene, hit$Zh13[1],
                ifelse(hit$is_1to1[1], "1:1", "multi")))
  } else {
    cat(sprintf("  %s (%s) -> NOT FOUND\n", test_names[i], gene))
  }
}
cat("\n")


# ===== SECTION 5: CONVERT ALL GENE LISTS =====

cat("========================================\n")
cat("SECTION 5: BATCH CONVERSION\n")
cat("========================================\n\n")

# --- 5A: JAG1 targets (1,567 genes) ---
cat("--- JAG1 targets ---\n")
jag1_targets <- read.csv("03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv")
cat("  Loaded:", nrow(jag1_targets), "genes\n")

jag1_conv <- convert_wm82_to_zh13(jag1_targets$GeneID, lookup)
jag1_conv$Confidence_Tier <- jag1_targets$Confidence_Tier
jag1_conv$Pattern <- jag1_targets$Pattern
cat(sprintf("  Converted: %d / %d (%.1f%%)\n",
            sum(!is.na(jag1_conv$Zh13)), nrow(jag1_conv),
            100 * sum(!is.na(jag1_conv$Zh13)) / nrow(jag1_conv)))

write.csv(jag1_conv, file.path(data_dir, "gene_lists/jag1_targets_1567_converted.csv"),
          row.names = FALSE)

# --- 5B: High-confidence targets (79 genes) ---
cat("\n--- High-confidence targets ---\n")
hc_targets <- read.csv("03_results/tables/38_binding_comprehensive/high_confidence_targets.csv")
cat("  Loaded:", nrow(hc_targets), "genes\n")

hc_conv <- convert_wm82_to_zh13(hc_targets$GeneID, lookup)
hc_conv$DE_Tier <- hc_targets$DE_Tier
hc_conv$Evidence_Class <- hc_targets$Evidence_Class
cat(sprintf("  Converted: %d / %d (%.1f%%)\n",
            sum(!is.na(hc_conv$Zh13)), nrow(hc_conv),
            100 * sum(!is.na(hc_conv$Zh13)) / nrow(hc_conv)))

write.csv(hc_conv, file.path(data_dir, "gene_lists/high_confidence_79_converted.csv"),
          row.names = FALSE)

# --- 5C: KRPs (9 genes) ---
cat("\n--- KRP cell cycle inhibitors ---\n")
krps <- read.csv("03_results/tables/KRP_verification/KRP_expression_summary.csv")
cat("  Loaded:", nrow(krps), "genes\n")

krp_conv <- convert_wm82_to_zh13(krps$gene_id, lookup)
krp_conv$name <- krps$name
cat(sprintf("  Converted: %d / %d\n", sum(!is.na(krp_conv$Zh13)), nrow(krp_conv)))

write.csv(krp_conv, file.path(data_dir, "gene_lists/krps_9_converted.csv"),
          row.names = FALSE)

# --- 5D: Cyclins ---
cat("\n--- Cyclins (comprehensive) ---\n")
cyclins <- read.csv("03_results/tables/Cyclin_analysis/Cyclin_comprehensive_list.csv")
cat("  Loaded:", nrow(cyclins), "genes\n")

cyclin_conv <- convert_wm82_to_zh13(cyclins$locusName, lookup)
cyclin_conv$cyclin_type <- cyclins$cyclin_type
cyclin_conv$is_jag1_target <- cyclins$is_jag1_target
cyclin_conv$is_expressed <- cyclins$is_expressed
cat(sprintf("  Converted: %d / %d (%.1f%%)\n",
            sum(!is.na(cyclin_conv$Zh13)), nrow(cyclin_conv),
            100 * sum(!is.na(cyclin_conv$Zh13)) / nrow(cyclin_conv)))

write.csv(cyclin_conv, file.path(data_dir, "gene_lists/cyclins_converted.csv"),
          row.names = FALSE)

# --- 5E: Cyclin JAG1 targets (7 genes) ---
cat("\n--- Cyclin JAG1 targets ---\n")
cyclin_jag <- read.csv("03_results/tables/Cyclin_analysis/Cyclin_JAG1_targets.csv")
cat("  Loaded:", nrow(cyclin_jag), "genes\n")

cyclin_jag_conv <- convert_wm82_to_zh13(cyclin_jag$locusName, lookup)
cyclin_jag_conv$cyclin_type <- cyclin_jag$cyclin_type
cat(sprintf("  Converted: %d / %d\n", sum(!is.na(cyclin_jag_conv$Zh13)), nrow(cyclin_jag_conv)))

write.csv(cyclin_jag_conv, file.path(data_dir, "gene_lists/cyclin_jag1_targets_converted.csv"),
          row.names = FALSE)


# ===== SECTION 6: CONVERT FAN ET AL. MARKERS (ZH13 → Wm82) =====

cat("\n========================================\n")
cat("SECTION 6: FAN ET AL. MARKER CONVERSION\n")
cat("========================================\n\n")

if (file.exists(fan_supp)) {
  fan_markers <- read_excel(fan_supp, sheet = "Supplemental Table 5B", skip = 1)
  cat("Fan markers loaded:", nrow(fan_markers), "\n")

  fan_gene_col <- names(fan_markers)[1]
  fan_conv <- convert_zh13_to_wm82(fan_markers[[fan_gene_col]], lookup)
  fan_conv$Cell_cluster <- fan_markers$`Cell cluster`
  fan_conv$Cell_type <- fan_markers$`Cell type`
  fan_conv$Arabidopsis_Gene <- fan_markers$`Arabidopsis Gene`
  fan_conv$Description <- fan_markers$Description
  fan_conv$avg_log2FC <- fan_markers$`Average log fold change`
  fan_conv$p_val_adj <- fan_markers$`Adjusted p value`

  cat(sprintf("  Converted: %d / %d (%.1f%%)\n",
              sum(!is.na(fan_conv$Wm82)), nrow(fan_conv),
              100 * sum(!is.na(fan_conv$Wm82)) / nrow(fan_conv)))

  write.csv(fan_conv, file.path(data_dir, "fan2025_shoot_apex_markers_converted.csv"),
            row.names = FALSE)
} else {
  cat("Fan et al. supplementary file not found, skipping.\n")
  cat("  Expected at:", fan_supp, "\n")
}


# ===== SECTION 7: SAVE LOOKUP TABLE =====

cat("\n========================================\n")
cat("SECTION 7: SAVE LOOKUP TABLE\n")
cat("========================================\n\n")

write.csv(lookup, file.path(data_dir, "gene_id_conversion_full_lookup.csv"),
          row.names = FALSE)
cat("Saved: gene_id_conversion_full_lookup.csv\n")


# ===== SUMMARY =====

cat("\n")
cat("================================================================\n")
cat("  GENE ID CONVERSION COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Pangene rows: %d total, %d with both Wm82+ZH13\n",
            nrow(pangene), nrow(lookup)))
cat(sprintf("  1:1 mappings: %d\n", sum(lookup$is_1to1)))
cat("\n")
cat(sprintf("  JAG1 targets:     %d / %d converted (%.1f%%)\n",
            sum(!is.na(jag1_conv$Zh13)), nrow(jag1_conv),
            100 * sum(!is.na(jag1_conv$Zh13)) / nrow(jag1_conv)))
cat(sprintf("  High-confidence:  %d / %d converted\n",
            sum(!is.na(hc_conv$Zh13)), nrow(hc_conv)))
cat(sprintf("  KRPs:             %d / %d converted\n",
            sum(!is.na(krp_conv$Zh13)), nrow(krp_conv)))
cat(sprintf("  Cyclins:          %d / %d converted\n",
            sum(!is.na(cyclin_conv$Zh13)), nrow(cyclin_conv)))
cat("\n")
cat("  Output:", data_dir, "\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
