# Integrates Huang et al. 2021 ChIP-seq and Wang et al. 2024 DAP-seq with DE results
# Uses actual peak coordinates, not just gene lists

rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 38: COMPREHENSIVE BINDING DATA INTEGRATION\n")
cat("  Integrating ChIP-seq (Huang 2021) + DAP-seq (Wang 2024)\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/38_binding_comprehensive", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/38_binding_comprehensive", recursive = TRUE, showWarnings = FALSE)
required_packages <- c("dplyr", "tidyr", "ggplot2", "readxl", "GenomicRanges", "rtracklayer")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GenomicRanges", "rtracklayer")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  cat("  Loaded:", pkg, "\n")
}
cat("\n")


# Promoter window for peak-to-gene association
PROMOTER_UPSTREAM <- 2000   # Huang et al. used 2kb for JTG definition
PROMOTER_DOWNSTREAM <- 0    # TSS only

# GmJAG1 binding motifs (all three validated by EMSA)
MOTIFS <- list(
  M0 = list(seq = "GTTGGA", source = "Huang 2021 ChIP-seq", method = "in vivo"),
  M1 = list(seq = "ACGCCACT", source = "Wang 2024 DAP-seq", method = "in vitro"),
  M2 = list(seq = "ACTGGCAG", source = "Wang 2024 DAP-seq", method = "in vitro")
)

cat("Analysis parameters:\n")
cat("  Promoter window:", PROMOTER_UPSTREAM, "bp upstream of TSS\n")
cat("  M0 motif:", MOTIFS$M0$seq, "(", MOTIFS$M0$source, ")\n")
cat("  M1 motif:", MOTIFS$M1$seq, "(", MOTIFS$M1$source, ")\n")
cat("  M2 motif:", MOTIFS$M2$seq, "(", MOTIFS$M2$source, ")\n\n")
gff_file <- "01_data/genome/Gmax_880_Wm82.a6.v1.gene.gff3"

if (file.exists(gff_file)) {
  cat("Loading gene annotation...\n")
  gff <- import(gff_file)
  genes <- gff[gff$type == "gene"]
  cat("  Total genes:", length(genes), "\n\n")

  # Create gene TSS table
  gene_tss <- data.frame(
    GeneID = genes$ID,
    Chr = as.character(seqnames(genes)),
    Strand = as.character(strand(genes)),
    TSS = ifelse(as.character(strand(genes)) == "+", start(genes), end(genes)),
    stringsAsFactors = FALSE
  )
  # Clean gene IDs (remove version numbers if present)
  gene_tss$GeneID <- gsub("\\.v.*", "", gene_tss$GeneID)
  gene_tss$GeneID <- gsub("\\.Wm.*", "", gene_tss$GeneID)

} else {
  stop("GFF file not found: ", gff_file)
}
load("03_results/checkpoints/14_JAG1_targets.RData")

jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("JAG1 targets from RNA-seq:", nrow(jag1_targets), "\n")
cat("  Gold tier:", sum(jag1_targets$Confidence_Tier == "Gold"), "\n")
cat("  Silver tier:", sum(jag1_targets$Confidence_Tier == "Silver"), "\n")
cat("  Bronze tier:", sum(jag1_targets$Confidence_Tier == "Bronze"), "\n\n")

target_genes <- unique(jag1_targets$GeneID)
huang_file <- file.path(base_dir, "Huang et al 2021 Chip Seq", "1-s2.0-S0888754320320723-mmc2.xlsx")

if (file.exists(huang_file)) {
  cat("Reading Huang et al. 2021 ChIP-seq data...\n")

  # Read Table S3: GmJAG1 binding sites (peaks)
  huang_peaks <- read_excel(huang_file, sheet = "Table S3", skip = 2)

  # Standardize column names
  if (ncol(huang_peaks) >= 10) {
    colnames(huang_peaks)[1:10] <- c("Chromosome", "Start", "End", "GmJBS_ID", "GmJBS_Width",
                                      "Associated_Gene", "Strand", "Location", "Distance", "Score")
  }

  cat("  Total ChIP-seq peaks:", nrow(huang_peaks), "\n")

  # Filter to valid peaks
  huang_peaks <- huang_peaks[!is.na(huang_peaks$Chromosome) & !is.na(huang_peaks$Start), ]
  cat("  Valid peaks:", nrow(huang_peaks), "\n")

  # Analyze peak locations
  if ("Location" %in% colnames(huang_peaks)) {
    cat("\n  Peak location distribution:\n")
    print(table(huang_peaks$Location))
  }

  # Filter to promoter peaks (within 2kb of TSS)
  if ("Distance" %in% colnames(huang_peaks)) {
    huang_peaks$Distance <- as.numeric(huang_peaks$Distance)
    huang_promoter_peaks <- huang_peaks[abs(huang_peaks$Distance) <= PROMOTER_UPSTREAM, ]
    cat("\n  Peaks within", PROMOTER_UPSTREAM, "bp of TSS:", nrow(huang_promoter_peaks), "\n")
  } else {
    huang_promoter_peaks <- huang_peaks
  }

  # Extract genes with ChIP-seq peaks in promoter
  huang_genes <- unique(huang_promoter_peaks$Associated_Gene)
  huang_genes <- huang_genes[!is.na(huang_genes) & grepl("Glyma", huang_genes)]
  cat("  Genes with promoter ChIP-seq peaks:", length(huang_genes), "\n")

  has_huang <- TRUE

} else {
  cat("WARNING: Huang et al. supplementary file not found.\n")
  cat("  Expected:", huang_file, "\n")
  has_huang <- FALSE
  huang_genes <- character(0)
}


cat("\n========================================\n")

wang_file <- file.path(base_dir, "Wang et al DAP-Seq", "plants-3024718-supplementary.xlsx")

if (file.exists(wang_file)) {
  cat("Reading Wang et al. 2024 DAP-seq data...\n")

  # Read Table S2: DAP-seq enriched peaks
  wang_peaks <- read_excel(wang_file, sheet = "Table S2", skip = 1)

  if (ncol(wang_peaks) >= 9) {
    colnames(wang_peaks)[1:9] <- c("Chr_ID", "Start", "End", "Peak_ID", "Associated_Gene",
                                    "Strand", "Location", "Distance", "GmJAG1_Peak")
  }

  cat("  Total DAP-seq peaks:", nrow(wang_peaks), "\n")

  # Analyze peak locations
  if ("Location" %in% colnames(wang_peaks)) {
    cat("\n  Peak location distribution:\n")
    print(table(wang_peaks$Location))
  }

  # Key insight from Wang et al.: Only 2.4% overlap with ChIP-seq
  # Most DAP-seq peaks are in distal intergenic regions (direct binding)

  # Filter to promoter peaks
  wang_promoter_peaks <- wang_peaks[grepl("Promoter|promoter", wang_peaks$Location, ignore.case = TRUE), ]
  cat("\n  Promoter DAP-seq peaks:", nrow(wang_promoter_peaks), "\n")

  # Extract genes with DAP-seq peaks
  wang_genes_all <- unique(wang_peaks$Associated_Gene)
  wang_genes_all <- wang_genes_all[!is.na(wang_genes_all) & grepl("Glyma", wang_genes_all)]

  wang_genes_promoter <- unique(wang_promoter_peaks$Associated_Gene)
  wang_genes_promoter <- wang_genes_promoter[!is.na(wang_genes_promoter) & grepl("Glyma", wang_genes_promoter)]

  cat("  Genes with ANY DAP-seq peak:", length(wang_genes_all), "\n")
  cat("  Genes with PROMOTER DAP-seq peak:", length(wang_genes_promoter), "\n")

  # Read Table S4: Motif positions in promoters (M0, M1, M2)
  cat("\nReading Table S4 (Motif positions)...\n")
  wang_motifs <- read_excel(wang_file, sheet = "Table S4", skip = 1)

  if (ncol(wang_motifs) >= 7) {
    colnames(wang_motifs) <- c("Motif_ID", "Gene_ID", "Start", "End", "Strand", "Pvalue", "Match_Seq")
  }

  # Genes with each motif in promoter
  wang_m0_genes <- unique(wang_motifs$Gene_ID[wang_motifs$Motif_ID == "M0" | grepl("M0|GTTGGA", wang_motifs$Motif_ID)])
  wang_m1_genes <- unique(wang_motifs$Gene_ID[wang_motifs$Motif_ID == "M1" | grepl("M1|ACGCCACT", wang_motifs$Motif_ID)])
  wang_m2_genes <- unique(wang_motifs$Gene_ID[wang_motifs$Motif_ID == "M2" | grepl("M2|ACTGGCAG", wang_motifs$Motif_ID)])

  cat("  Genes with M0 motif:", length(wang_m0_genes), "\n")
  cat("  Genes with M1 motif:", length(wang_m1_genes), "\n")
  cat("  Genes with M2 motif:", length(wang_m2_genes), "\n")

  has_wang <- TRUE

} else {
  cat("WARNING: Wang et al. supplementary file not found.\n")
  cat("  Expected:", wang_file, "\n")
  has_wang <- FALSE
  wang_genes_all <- character(0)
  wang_genes_promoter <- character(0)
}


cat("\n========================================\n")

if (has_huang && has_wang) {
  # Gene-level overlap
  overlap_chipseq_dapseq <- intersect(huang_genes, wang_genes_all)
  cat("Overlap between binding datasets (gene level):\n")
  cat("  ChIP-seq genes (Huang):", length(huang_genes), "\n")
  cat("  DAP-seq genes (Wang):", length(wang_genes_all), "\n")
  cat("  Overlap:", length(overlap_chipseq_dapseq), "\n")
  cat("  % of DAP-seq in ChIP-seq:", round(100 * length(overlap_chipseq_dapseq) / length(wang_genes_all), 1), "%\n")

  cat("\n  Note: Wang et al. reported only 2.4% peak overlap\n")
  cat("  This suggests ChIP-seq captures indirect binding via co-factors\n")
  cat("  while DAP-seq captures direct GmJAG1-DNA interactions\n\n")
}
# Create comprehensive integration table
integrated <- data.frame(
  GeneID = target_genes,
  DE_Tier = jag1_targets$Confidence_Tier[match(target_genes, jag1_targets$GeneID)],
  stringsAsFactors = FALSE
)

# Add logFC if available
if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  integrated$logFC <- jag1_targets$Mean_logFC_Pairwise[match(target_genes, jag1_targets$GeneID)]
}

# Add binding evidence
integrated$ChIPseq_Huang <- integrated$GeneID %in% huang_genes
integrated$DAPseq_Wang <- integrated$GeneID %in% wang_genes_all
integrated$DAPseq_Promoter <- integrated$GeneID %in% wang_genes_promoter

# Add motif evidence (from Wang Table S4 if available)
if (has_wang && exists("wang_m0_genes")) {
  integrated$Has_M0_Wang <- integrated$GeneID %in% wang_m0_genes
  integrated$Has_M1_Wang <- integrated$GeneID %in% wang_m1_genes
  integrated$Has_M2_Wang <- integrated$GeneID %in% wang_m2_genes
} else {
  integrated$Has_M0_Wang <- FALSE
  integrated$Has_M1_Wang <- FALSE
  integrated$Has_M2_Wang <- FALSE
}

# Calculate evidence scores
integrated$N_Binding_Types <- as.integer(integrated$ChIPseq_Huang) +
                               as.integer(integrated$DAPseq_Wang)

integrated$N_Motifs <- as.integer(integrated$Has_M0_Wang) +
                        as.integer(integrated$Has_M1_Wang) +
                        as.integer(integrated$Has_M2_Wang)
# Multi-evidence tiers (combining in vivo and in vitro)
integrated$Evidence_Class <- "DE_only"

# Tier 1: Both ChIP-seq AND DAP-seq (rarest, strongest)
tier1_idx <- integrated$ChIPseq_Huang & integrated$DAPseq_Wang
integrated$Evidence_Class[tier1_idx] <- "ChIPseq_AND_DAPseq"

# Tier 2a: ChIP-seq only (in vivo binding)
tier2a_idx <- integrated$ChIPseq_Huang & !integrated$DAPseq_Wang
integrated$Evidence_Class[tier2a_idx] <- "ChIPseq_only"

# Tier 2b: DAP-seq only (in vitro direct binding)
tier2b_idx <- !integrated$ChIPseq_Huang & integrated$DAPseq_Wang
integrated$Evidence_Class[tier2b_idx] <- "DAPseq_only"

cat("Evidence classification:\n")
print(table(integrated$Evidence_Class))
cat("\n")

# Combined priority score
# DE tier: Gold=3, Silver=2, Bronze=1
# Binding: ChIP+DAP=3, ChIP=2, DAP=2, motif_only=1, none=0

de_score <- c("Gold" = 3, "Silver" = 2, "Bronze" = 1)
integrated$DE_Score <- de_score[integrated$DE_Tier]

binding_score <- case_when(
  integrated$ChIPseq_Huang & integrated$DAPseq_Wang ~ 3,
  integrated$ChIPseq_Huang ~ 2,
  integrated$DAPseq_Wang ~ 2,
  integrated$N_Motifs > 0 ~ 1,
  TRUE ~ 0
)
integrated$Binding_Score <- binding_score

integrated$Total_Score <- integrated$DE_Score + integrated$Binding_Score

# Sort by total score
integrated <- integrated[order(-integrated$Total_Score, -integrated$DE_Score), ]

cat("Priority score distribution:\n")
print(table(integrated$Total_Score))
cat("\n")
# Highest confidence: Gold/Silver DE + both binding datasets
top_tier <- integrated[integrated$DE_Tier %in% c("Gold", "Silver") &
                        integrated$Evidence_Class == "ChIPseq_AND_DAPseq", ]
cat("HIGHEST CONFIDENCE (Gold/Silver DE + ChIP + DAP):", nrow(top_tier), "genes\n")
if (nrow(top_tier) > 0) {
  print(top_tier[, c("GeneID", "DE_Tier", "logFC", "Evidence_Class", "Total_Score")])
}

# High confidence: Gold/Silver DE + any binding
high_tier <- integrated[integrated$DE_Tier %in% c("Gold", "Silver") &
                         integrated$N_Binding_Types >= 1, ]
cat("\nHIGH CONFIDENCE (Gold/Silver DE + any binding):", nrow(high_tier), "genes\n")

# Medium confidence: Any DE tier + both binding
medium_tier <- integrated[integrated$Evidence_Class == "ChIPseq_AND_DAPseq", ]
cat("\nMEDIUM-HIGH (Any DE + both binding):", nrow(medium_tier), "genes\n")


cat("\n========================================\n")

motif_by_class <- integrated %>%
  group_by(Evidence_Class) %>%
  summarise(
    N = n(),
    Pct_M0 = round(mean(Has_M0_Wang) * 100, 1),
    Pct_M1 = round(mean(Has_M1_Wang) * 100, 1),
    Pct_M2 = round(mean(Has_M2_Wang) * 100, 1),
    Pct_Any_Motif = round(mean(N_Motifs > 0) * 100, 1),
    Mean_logFC = round(mean(logFC, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Motif presence by evidence class:\n\n")
print(as.data.frame(motif_by_class))


cat("\n========================================\n")

# Plot 1: Evidence class distribution
png("03_results/figures/38_binding_comprehensive/evidence_class_distribution.png",
    width = 900, height = 600, res = 120)

evidence_counts <- as.data.frame(table(integrated$Evidence_Class))
colnames(evidence_counts) <- c("Class", "Count")
evidence_counts$Class <- factor(evidence_counts$Class,
                                 levels = c("ChIPseq_AND_DAPseq", "ChIPseq_only",
                                           "DAPseq_only", "DE_only"))

ggplot(evidence_counts, aes(x = Class, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  scale_fill_manual(values = c("ChIPseq_AND_DAPseq" = "#FFD700",
                                "ChIPseq_only" = "#3498DB",
                                "DAPseq_only" = "#E74C3C",
                                "DE_only" = "#95A5A6")) +
  labs(title = "GmJAG1 Target Classification by Binding Evidence",
       subtitle = "ChIP-seq (Huang 2021, in vivo) vs DAP-seq (Wang 2024, in vitro)",
       x = "Evidence Class",
       y = "Number of DE Genes") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
cat("Saved: evidence_class_distribution.png\n")

# Plot 2: DE tier vs binding evidence heatmap
png("03_results/figures/38_binding_comprehensive/tier_binding_heatmap.png",
    width = 800, height = 600, res = 120)

tier_binding <- integrated %>%
  group_by(DE_Tier, Evidence_Class) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Evidence_Class, values_from = N, values_fill = 0)

tier_binding_long <- integrated %>%
  group_by(DE_Tier, Evidence_Class) %>%
  summarise(N = n(), .groups = "drop")

tier_binding_long$DE_Tier <- factor(tier_binding_long$DE_Tier, levels = c("Gold", "Silver", "Bronze"))
tier_binding_long$Evidence_Class <- factor(tier_binding_long$Evidence_Class,
                                            levels = c("ChIPseq_AND_DAPseq", "ChIPseq_only",
                                                      "DAPseq_only", "DE_only"))

ggplot(tier_binding_long, aes(x = Evidence_Class, y = DE_Tier, fill = N)) +
  geom_tile(color = "white") +
  geom_text(aes(label = N), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "#E74C3C") +
  labs(title = "DE Tier vs Binding Evidence",
       x = "Binding Evidence",
       y = "DE Confidence Tier",
       fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
cat("Saved: tier_binding_heatmap.png\n")

# Plot 3: Priority score vs logFC
png("03_results/figures/38_binding_comprehensive/priority_vs_logFC.png",
    width = 800, height = 600, res = 120)

ggplot(integrated, aes(x = factor(Total_Score), y = logFC, fill = DE_Tier)) +
  geom_boxplot(outlier.shape = 21) +
  scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Fold Change by Priority Score",
       subtitle = "Higher score = more evidence",
       x = "Priority Score (DE + Binding)",
       y = "log2 Fold Change (Narrow/Broad)") +
  theme_minimal()

dev.off()
cat("Saved: priority_vs_logFC.png\n")

# Plot 4: Venn-style summary
png("03_results/figures/38_binding_comprehensive/evidence_summary.png",
    width = 1000, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Panel A: Binding overlap
n_chipseq <- sum(integrated$ChIPseq_Huang)
n_dapseq <- sum(integrated$DAPseq_Wang)
n_both <- sum(integrated$ChIPseq_Huang & integrated$DAPseq_Wang)
n_neither <- sum(!integrated$ChIPseq_Huang & !integrated$DAPseq_Wang)

barplot(c(n_chipseq, n_dapseq, n_both, n_neither),
        names.arg = c("ChIP-seq\nonly", "DAP-seq\nonly", "Both", "Neither"),
        col = c("#3498DB", "#E74C3C", "#FFD700", "#95A5A6"),
        main = "A. Binding Evidence Among DE Genes",
        ylab = "Number of Genes")

# Panel B: Total score distribution
score_counts <- table(integrated$Total_Score)
barplot(score_counts,
        col = colorRampPalette(c("#95A5A6", "#FFD700"))(length(score_counts)),
        main = "B. Priority Score Distribution",
        xlab = "Priority Score",
        ylab = "Number of Genes")

dev.off()
cat("Saved: evidence_summary.png\n")


cat("\n========================================\n")

# Save full integrated table
write.csv(integrated,
          "03_results/tables/38_binding_comprehensive/integrated_all_evidence.csv",
          row.names = FALSE)
cat("Saved: integrated_all_evidence.csv (", nrow(integrated), " genes)\n")

# Save high-confidence targets
if (nrow(high_tier) > 0) {
  write.csv(high_tier,
            "03_results/tables/38_binding_comprehensive/high_confidence_targets.csv",
            row.names = FALSE)
  cat("Saved: high_confidence_targets.csv (", nrow(high_tier), " genes)\n")
}

# Save top candidates (both binding types)
if (nrow(medium_tier) > 0) {
  write.csv(medium_tier,
            "03_results/tables/38_binding_comprehensive/dual_binding_targets.csv",
            row.names = FALSE)
  cat("Saved: dual_binding_targets.csv (", nrow(medium_tier), " genes)\n")
}

# Save summary statistics
summary_stats <- list(
  n_de_targets = nrow(integrated),
  n_chipseq = sum(integrated$ChIPseq_Huang),
  n_dapseq = sum(integrated$DAPseq_Wang),
  n_both_binding = sum(integrated$ChIPseq_Huang & integrated$DAPseq_Wang),
  n_high_confidence = nrow(high_tier),
  motif_by_class = motif_by_class,
  evidence_counts = table(integrated$Evidence_Class)
)

save(integrated, summary_stats, high_tier, medium_tier,
     file = "03_results/checkpoints/38_binding_comprehensive.RData")
cat("Saved: 38_binding_comprehensive.RData\n")


cat("\n================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("INTEGRATION SUMMARY:\n")

cat("RNA-seq DE targets:", nrow(integrated), "\n\n")

cat("BINDING EVIDENCE:\n")
cat("  ChIP-seq (Huang 2021, in vivo):", sum(integrated$ChIPseq_Huang),
    "(", round(100*mean(integrated$ChIPseq_Huang), 1), "%)\n")
cat("  DAP-seq (Wang 2024, in vitro):", sum(integrated$DAPseq_Wang),
    "(", round(100*mean(integrated$DAPseq_Wang), 1), "%)\n")
cat("  BOTH binding types:", sum(integrated$ChIPseq_Huang & integrated$DAPseq_Wang),
    "(", round(100*mean(integrated$ChIPseq_Huang & integrated$DAPseq_Wang), 1), "%)\n\n")

cat("KEY INSIGHT:\n")
cat("  ChIP-seq (in vivo) captures indirect binding via co-factors\n")
cat("  DAP-seq (in vitro) captures direct GmJAG1-DNA binding\n")
cat("  Genes with BOTH have strongest evidence for direct regulation\n\n")

cat("EVIDENCE CLASSIFICATION:\n")
print(table(integrated$Evidence_Class))

cat("\nOUTPUT FILES:\n")
cat("  - 03_results/tables/38_binding_comprehensive/*.csv\n")
cat("  - 03_results/figures/38_binding_comprehensive/*.png\n")
cat("  - 03_results/checkpoints/38_binding_comprehensive.RData\n\n")

cat("FOR MANUSCRIPT:\n")
cat("  Use 'high_confidence_targets.csv' for validated target discussion\n")
cat("  Use 'dual_binding_targets.csv' for strongest mechanistic evidence\n\n")
