# Updated: 2026-01-25
# Now extracts FULL peak coordinates, not just gene lists

rm(list = ls())
gc()

cat("\n")
cat("  EXTRACT PUBLISHED BINDING DATA - WITH PEAK COORDINATES\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("01_data/published", recursive = TRUE, showWarnings = FALSE)
dir.create("01_data/published/peaks", recursive = TRUE, showWarnings = FALSE)


if (!requireNamespace("readxl", quietly = TRUE)) {
  cat("Installing readxl package...\n")
  install.packages("readxl")
}
library(readxl)


cat("  Extracting FULL peak coordinates\n")

huang_file <- file.path(base_dir, "Phase2-Refined-Analysis/Manuscript/Literature",
                        "Huang et al 2021 Chip Seq/1-s2.0-S0888754320320723-mmc2.xlsx")

# Alternative path if not found
if (!file.exists(huang_file)) {
  huang_file <- file.path(base_dir, "Huang et al 2021 Chip Seq",
                          "1-s2.0-S0888754320320723-mmc2.xlsx")
}

if (file.exists(huang_file)) {
  cat("Reading Huang et al. 2021 ChIP-seq data...\n\n")

  # Read Table S3: GmJAG1 binding sites with FULL data
  cat("Reading Table S3 (ChIP-seq binding sites with coordinates)...\n")
  huang_s3 <- read_excel(huang_file, sheet = "Table S3", skip = 2)

  # Set proper column names
  colnames(huang_s3) <- c("Chromosome", "Start", "End", "GmJBS_ID", "Peak_Width",
                          "GeneID", "Gene_Direction", "Genomic_Region",
                          "Distance_to_TSS", "FPKM_Leaf", "FPKM_SD", "Chromatin_State")

  cat("  Total peak-gene associations:", nrow(huang_s3), "\n")
  cat("  Unique peaks (GmJBS IDs):", length(unique(huang_s3$GmJBS_ID)), "\n")
  cat("  Unique genes:", length(unique(huang_s3$GeneID)), "\n\n")

  cat("Genomic Region Distribution:\n")
  region_table <- table(huang_s3$Genomic_Region)
  for (region in names(region_table)) {
    pct <- round(region_table[region] / nrow(huang_s3) * 100, 1)
    cat(sprintf("  %s: %d (%.1f%%)\n", region, region_table[region], pct))
  }

  cat("\nBinding Location Analysis:\n")

  # Count genes by binding region
  promoter_genes <- unique(huang_s3$GeneID[huang_s3$Genomic_Region == "Promoter"])
  genic_genes <- unique(huang_s3$GeneID[huang_s3$Genomic_Region == "Genic"])
  downstream_genes <- unique(huang_s3$GeneID[huang_s3$Genomic_Region == "Downstream"])
  intergenic_genes <- unique(huang_s3$GeneID[huang_s3$Genomic_Region == "InterGenicRegion"])

  cat(sprintf("  Genes with PROMOTER binding: %d\n", length(promoter_genes)))
  cat(sprintf("  Genes with GENIC (gene body) binding: %d\n", length(genic_genes)))
  cat(sprintf("  Genes with DOWNSTREAM binding: %d\n", length(downstream_genes)))
  cat(sprintf("  Genes with INTERGENIC binding: %d\n", length(intergenic_genes)))

  # Identify genes with ONLY gene body binding (no promoter)
  only_genic <- setdiff(genic_genes, promoter_genes)
  cat(sprintf("\n  >>> CRITICAL: %d genes have binding ONLY in gene body (not promoter)\n",
              length(only_genic)))
  cat("  >>> These are MISSED by promoter-only motif searches!\n\n")


  # 1. Save complete peak-gene association table
  write.csv(huang_s3, "01_data/published/peaks/Huang_2021_ChIPseq_peaks_full.csv",
            row.names = FALSE)
  cat("Saved: 01_data/published/peaks/Huang_2021_ChIPseq_peaks_full.csv\n")

  # 2. Create gene-level summary with binding region info
  huang_gene_summary <- data.frame(
    GeneID = unique(huang_s3$GeneID),
    stringsAsFactors = FALSE
  )

  # Add binding region flags
  huang_gene_summary$Has_Promoter_Binding <- huang_gene_summary$GeneID %in% promoter_genes
  huang_gene_summary$Has_Genic_Binding <- huang_gene_summary$GeneID %in% genic_genes
  huang_gene_summary$Has_Downstream_Binding <- huang_gene_summary$GeneID %in% downstream_genes
  huang_gene_summary$Has_Intergenic_Binding <- huang_gene_summary$GeneID %in% intergenic_genes

  # Categorize primary binding location
  huang_gene_summary$Primary_Binding_Region <- "Multiple"
  huang_gene_summary$Primary_Binding_Region[huang_gene_summary$Has_Promoter_Binding &
                                             !huang_gene_summary$Has_Genic_Binding] <- "Promoter_Only"
  huang_gene_summary$Primary_Binding_Region[!huang_gene_summary$Has_Promoter_Binding &
                                             huang_gene_summary$Has_Genic_Binding] <- "Genic_Only"
  huang_gene_summary$Primary_Binding_Region[huang_gene_summary$Has_Promoter_Binding &
                                             huang_gene_summary$Has_Genic_Binding] <- "Promoter_and_Genic"

  # Count peaks per gene
  peaks_per_gene <- table(huang_s3$GeneID)
  huang_gene_summary$N_Peaks <- peaks_per_gene[huang_gene_summary$GeneID]

  # Get closest peak distance to TSS for each gene
  closest_dist <- aggregate(abs(as.numeric(huang_s3$Distance_to_TSS)),
                            by = list(GeneID = huang_s3$GeneID),
                            FUN = min, na.rm = TRUE)
  colnames(closest_dist)[2] <- "Closest_Peak_Distance"
  huang_gene_summary <- merge(huang_gene_summary, closest_dist, by = "GeneID", all.x = TRUE)

  # Save gene-level summary
  write.csv(huang_gene_summary, "01_data/published/Huang_2021_ChIPseq_gene_summary.csv",
            row.names = FALSE)
  cat("Saved: 01_data/published/Huang_2021_ChIPseq_gene_summary.csv\n")

  # 3. Also save legacy format for backward compatibility
  huang_output <- data.frame(
    GeneID = huang_gene_summary$GeneID,
    ChIPseq_binding = 1,
    stringsAsFactors = FALSE
  )
  write.csv(huang_output, "01_data/published/Huang_2021_ChIPseq_targets.csv", row.names = FALSE)
  cat("Saved: 01_data/published/Huang_2021_ChIPseq_targets.csv (legacy format)\n")

  # Summary statistics
  cat("\nGene-level binding summary:\n")
  cat(sprintf("  Total genes with ChIP-seq binding: %d\n", nrow(huang_gene_summary)))
  cat(sprintf("  Promoter only: %d\n", sum(huang_gene_summary$Primary_Binding_Region == "Promoter_Only")))
  cat(sprintf("  Genic (gene body) only: %d\n", sum(huang_gene_summary$Primary_Binding_Region == "Genic_Only")))
  cat(sprintf("  Both promoter and genic: %d\n", sum(huang_gene_summary$Primary_Binding_Region == "Promoter_and_Genic")))
  cat(sprintf("  Multiple regions: %d\n", sum(huang_gene_summary$Primary_Binding_Region == "Multiple")))

} else {
  cat("ERROR: Huang supplementary file not found at:\n")
  cat("  ", huang_file, "\n")
}


cat("\n========================================\n")
cat("  Extracting FULL peak coordinates\n")

wang_file <- file.path(base_dir, "Phase2-Refined-Analysis/Manuscript/Literature",
                       "Wang et al DAP-Seq/plants-3024718-supplementary.xlsx")

# Alternative path
if (!file.exists(wang_file)) {
  wang_file <- file.path(base_dir, "Wang et al DAP-Seq", "plants-3024718-supplementary.xlsx")
}

if (file.exists(wang_file)) {
  cat("Reading Wang et al. 2024 DAP-seq data...\n\n")

  # Read Table S2: DAP-seq enriched peaks
  cat("Reading Table S2 (DAP-seq peaks with coordinates)...\n")
  wang_s2 <- read_excel(wang_file, sheet = "Table S2", skip = 1)
  colnames(wang_s2) <- c("Chromosome", "Start", "End", "Peak_ID", "GeneID",
                         "Strand", "Genomic_Region", "Distance_to_TSS", "GmJAG1_Peak")

  cat("  Total peak-gene associations:", nrow(wang_s2), "\n")
  cat("  Unique peaks:", length(unique(wang_s2$Peak_ID)), "\n")

  # Filter to valid gene IDs
  wang_s2 <- wang_s2[!is.na(wang_s2$GeneID) & grepl("Glyma", wang_s2$GeneID), ]
  cat("  Unique genes:", length(unique(wang_s2$GeneID)), "\n\n")

  # Genomic region distribution
  cat("Genomic Region Distribution:\n")
  if ("Genomic_Region" %in% colnames(wang_s2)) {
    region_table <- table(wang_s2$Genomic_Region)
    for (region in names(region_table)) {
      pct <- round(region_table[region] / nrow(wang_s2) * 100, 1)
      cat(sprintf("  %s: %d (%.1f%%)\n", region, region_table[region], pct))
    }
  }

  # Save full peak data
  write.csv(wang_s2, "01_data/published/peaks/Wang_2024_DAPseq_peaks_full.csv",
            row.names = FALSE)
  cat("\nSaved: 01_data/published/peaks/Wang_2024_DAPseq_peaks_full.csv\n")

  # Read Table S4: Motifs position at gene promoter
  cat("\nReading Table S4 (Motif positions)...\n")
  wang_s4 <- read_excel(wang_file, sheet = "Table S4", skip = 1)
  colnames(wang_s4) <- c("Motif_ID", "GeneID", "Start", "End", "Strand", "Pvalue", "Match_Seq")

  wang_s4 <- wang_s4[!is.na(wang_s4$GeneID) & grepl("Glyma", wang_s4$GeneID), ]
  cat("  Motif occurrences:", nrow(wang_s4), "\n")
  cat("  Unique genes with motifs:", length(unique(wang_s4$GeneID)), "\n")

  # Save motif data
  write.csv(wang_s4, "01_data/published/peaks/Wang_2024_DAPseq_motifs.csv",
            row.names = FALSE)
  cat("Saved: 01_data/published/peaks/Wang_2024_DAPseq_motifs.csv\n")

  # Create gene-level summary
  wang_peak_genes <- unique(wang_s2$GeneID)
  wang_motif_genes <- unique(wang_s4$GeneID)
  wang_all_genes <- unique(c(wang_peak_genes, wang_motif_genes))

  wang_gene_summary <- data.frame(
    GeneID = wang_all_genes,
    Has_DAPseq_Peak = wang_all_genes %in% wang_peak_genes,
    Has_Motif = wang_all_genes %in% wang_motif_genes,
    stringsAsFactors = FALSE
  )

  # Count peaks per gene
  if (length(wang_peak_genes) > 0) {
    peaks_per_gene <- table(wang_s2$GeneID)
    wang_gene_summary$N_Peaks <- peaks_per_gene[wang_gene_summary$GeneID]
    wang_gene_summary$N_Peaks[is.na(wang_gene_summary$N_Peaks)] <- 0
  }

  write.csv(wang_gene_summary, "01_data/published/Wang_2024_DAPseq_gene_summary.csv",
            row.names = FALSE)
  cat("Saved: 01_data/published/Wang_2024_DAPseq_gene_summary.csv\n")

  # Legacy format
  wang_output <- data.frame(
    GeneID = wang_all_genes,
    DAPseq_peak = as.integer(wang_all_genes %in% wang_peak_genes),
    Motif_in_promoter = as.integer(wang_all_genes %in% wang_motif_genes),
    stringsAsFactors = FALSE
  )
  write.csv(wang_output, "01_data/published/Wang_2024_DAPseq_targets.csv", row.names = FALSE)
  cat("Saved: 01_data/published/Wang_2024_DAPseq_targets.csv (legacy format)\n")

  cat("\nGene-level summary:\n")
  cat(sprintf("  Total genes: %d\n", nrow(wang_gene_summary)))
  cat(sprintf("  With DAP-seq peak: %d\n", sum(wang_gene_summary$Has_DAPseq_Peak)))
  cat(sprintf("  With motif: %d\n", sum(wang_gene_summary$Has_Motif)))

} else {
  cat("ERROR: Wang supplementary file not found at:\n")
  cat("  ", wang_file, "\n")
}


cat("\n========================================\n")
cat("COMPARISON: ChIP-seq vs DAP-seq\n")

if (exists("huang_gene_summary") && exists("wang_gene_summary")) {

  huang_genes <- huang_gene_summary$GeneID
  wang_genes <- wang_gene_summary$GeneID

  overlap <- intersect(huang_genes, wang_genes)
  only_chipseq <- setdiff(huang_genes, wang_genes)
  only_dapseq <- setdiff(wang_genes, huang_genes)

  cat("Gene overlap:\n")
  cat(sprintf("  ChIP-seq (Huang, in vivo): %d genes\n", length(huang_genes)))
  cat(sprintf("  DAP-seq (Wang, in vitro): %d genes\n", length(wang_genes)))
  cat(sprintf("  Overlap (both methods): %d genes (%.1f%% of DAP-seq)\n",
              length(overlap), length(overlap)/length(wang_genes)*100))
  cat(sprintf("  Only ChIP-seq: %d genes\n", length(only_chipseq)))
  cat(sprintf("  Only DAP-seq: %d genes\n", length(only_dapseq)))

  cat("\nInterpretation:\n")
  cat("  - ChIP-seq captures in vivo binding (includes indirect/co-factor mediated)\n")
  cat("  - DAP-seq captures in vitro binding (direct DNA-protein interaction)\n")
  cat("  - Low overlap suggests many ChIP-seq targets are indirect\n")
}


cat("\n========================================\n")
cat("OUTPUT FILES CREATED\n")

cat("Full peak coordinate data:\n")
cat("  1. 01_data/published/peaks/Huang_2021_ChIPseq_peaks_full.csv\n")
cat("  2. 01_data/published/peaks/Wang_2024_DAPseq_peaks_full.csv\n")
cat("  3. 01_data/published/peaks/Wang_2024_DAPseq_motifs.csv\n")

cat("\nGene-level summaries:\n")
cat("  4. 01_data/published/Huang_2021_ChIPseq_gene_summary.csv\n")
cat("  5. 01_data/published/Wang_2024_DAPseq_gene_summary.csv\n")

cat("\nLegacy format (for backward compatibility):\n")
cat("  6. 01_data/published/Huang_2021_ChIPseq_targets.csv\n")
cat("  7. 01_data/published/Wang_2024_DAPseq_targets.csv\n")

