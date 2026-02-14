# Script 25a: Promoter Window Validation
# Phase 1: Validate motif distribution using published gene lists from Huang et al. and Wang et al.
# Using OUR genome reference (Wm82.a6.v1) to confirm consistency

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 25a: PROMOTER WINDOW VALIDATION\n")
cat("  Phase 1: Validate using published gene lists\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/promoter/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/25a_promoter_validation", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD PACKAGES =====
cat("Loading packages...\n")
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# ===== PARAMETERS =====
PROMOTER_LENGTH <- 3000  # Using 3000 bp to capture 100% of promoter binding sites

# ===== DEFINE MOTIFS =====
motifs <- list(
  M0 = list(seq = "GTTGGA", source = "ChIP-seq (Huang 2021)", core = TRUE),
  M1 = list(seq = "ACGCCACT", source = "DAP-seq (Wang 2024)", core = FALSE),
  M2 = list(seq = "ACTGGCAG", source = "DAP-seq (Wang 2024)", core = FALSE)
)

cat("\nMotifs to scan:\n")
for (m in names(motifs)) {
  cat(sprintf("  %s: %s (%s)\n", m, motifs[[m]]$seq, motifs[[m]]$source))
}

# ===== LOAD GENOME AND ANNOTATION =====
cat("\n========================================\n")
cat("LOADING GENOME AND ANNOTATION\n")
cat("========================================\n\n")

genome_file <- "01_data/genome/Gmax_880_v6.0.fa"
gff_file <- "01_data/genome/Gmax_880_Wm82.a6.v1.gene.gff3"

cat("Loading genome:", genome_file, "\n")
genome <- readDNAStringSet(genome_file)
cat("  Chromosomes:", length(genome), "\n")

cat("Loading annotation:", gff_file, "\n")
gff <- import(gff_file)
genes <- gff[gff$type == "gene"]
cat("  Genes:", length(genes), "\n\n")

# ===== LOAD PUBLISHED GENE LISTS =====
cat("========================================\n")
cat("LOADING PUBLISHED GENE LISTS\n")
cat("========================================\n\n")

# --- Huang et al. 2021 ChIP-seq ---
huang_file <- file.path(base_dir, "Huang et al 2021 Chip Seq/1-s2.0-S0888754320320723-mmc2.xlsx")
cat("Loading Huang et al. 2021 data...\n")
huang_data <- read_excel(huang_file, sheet = "Table S3", skip = 2)
colnames(huang_data) <- c("Chromosome", "Start", "End", "GmJBS_ID", "GmJBS_width",
                          "Gene_ID", "Gene_Direction", "Genomic_Region",
                          "Distance_to_TSS", "FPKM_Leaf", "FPKM_SD", "Chromatin_State")

# Get unique genes with promoter binding
huang_promoter_genes <- huang_data %>%
  filter(Genomic_Region == "Promoter") %>%
  pull(Gene_ID) %>%
  unique()
cat("  Huang et al. promoter-bound genes:", length(huang_promoter_genes), "\n")

# Get all genes with binding
huang_all_genes <- huang_data %>%
  pull(Gene_ID) %>%
  unique()
cat("  Huang et al. all bound genes:", length(huang_all_genes), "\n")

# --- Wang et al. 2024 DAP-seq ---
wang_file <- file.path(base_dir, "Wang et al DAP-Seq/plants-3024718-supplementary.xlsx")
cat("\nLoading Wang et al. 2024 data...\n")
wang_data <- read_excel(wang_file, sheet = "Table S4", skip = 1)
colnames(wang_data) <- c("Motif_ID", "Promoter_Location", "Motif_Start",
                         "Motif_End", "Orientation", "P_value", "Match_Sequence")

# Extract gene IDs from promoter location
wang_data$Gene_ID <- gsub("_Promoter.*", "", wang_data$Promoter_Location)
wang_genes <- unique(wang_data$Gene_ID)
cat("  Wang et al. genes with motifs:", length(wang_genes), "\n\n")

# ===== FUNCTION TO EXTRACT PROMOTER =====
extract_promoter <- function(gene_id, genes_gr, genome, upstream = 2000) {
  # Find gene - try multiple matching strategies
  gene_match <- genes_gr[grepl(paste0("^", gene_id), genes_gr$ID) |
                         grepl(paste0("^", gene_id), genes_gr$Name)]

  if (length(gene_match) == 0) {
    # Try without version number
    gene_id_base <- gsub("\\.v.*", "", gene_id)
    gene_match <- genes_gr[grepl(paste0("^", gene_id_base), genes_gr$ID)]
  }

  if (length(gene_match) == 0) return(NULL)

  gene <- gene_match[1]
  chr <- as.character(seqnames(gene))
  strand_val <- as.character(strand(gene))

  if (!chr %in% names(genome)) return(NULL)

  chr_seq <- genome[[chr]]
  chr_len <- length(chr_seq)

  if (strand_val == "+") {
    prom_end <- start(gene) - 1
    prom_start <- max(1, prom_end - upstream + 1)
  } else {
    prom_start <- end(gene) + 1
    prom_end <- min(chr_len, prom_start + upstream - 1)
  }

  if (prom_start > prom_end || prom_start < 1 || prom_end > chr_len) return(NULL)

  prom_seq <- subseq(chr_seq, start = prom_start, end = prom_end)
  if (strand_val == "-") prom_seq <- reverseComplement(prom_seq)

  return(list(
    sequence = as.character(prom_seq),
    length = nchar(as.character(prom_seq)),
    strand = strand_val
  ))
}

# ===== FUNCTION TO FIND MOTIF POSITIONS =====
find_motif_positions <- function(sequence, motif) {
  seq_obj <- DNAString(sequence)
  motif_obj <- DNAString(motif)
  rc_motif <- reverseComplement(motif_obj)

  # Forward matches
  fwd_matches <- matchPattern(motif_obj, seq_obj, fixed = FALSE)
  fwd_pos <- start(fwd_matches)

  # Reverse complement matches
  rc_matches <- matchPattern(rc_motif, seq_obj, fixed = FALSE)
  rc_pos <- start(rc_matches)

  all_positions <- c(fwd_pos, rc_pos)

  if (length(all_positions) == 0) return(NULL)

  # Convert to distance from TSS
  # Position 1 in sequence = -PROMOTER_LENGTH from TSS
  # Position PROMOTER_LENGTH = -1 from TSS (immediately upstream)
  seq_len <- nchar(sequence)
  distance_from_tss <- all_positions - seq_len - 1  # Negative values

  return(distance_from_tss)
}

# ===== PHASE 1A: VALIDATE HUANG ET AL. GENES =====
cat("========================================\n")
cat("PHASE 1A: VALIDATING HUANG ET AL. GENES\n")
cat("========================================\n\n")

cat("Extracting promoters for Huang et al. promoter-bound genes...\n")
cat("Using", PROMOTER_LENGTH, "bp upstream window\n\n")

huang_results <- data.frame()
huang_promoters <- list()

pb <- txtProgressBar(min = 0, max = length(huang_promoter_genes), style = 3)
for (i in seq_along(huang_promoter_genes)) {
  gene_id <- huang_promoter_genes[i]
  prom <- extract_promoter(gene_id, genes, genome, PROMOTER_LENGTH)

  if (!is.null(prom)) {
    huang_promoters[[gene_id]] <- prom$sequence

    # Scan for each motif
    for (m in names(motifs)) {
      positions <- find_motif_positions(prom$sequence, motifs[[m]]$seq)
      if (!is.null(positions)) {
        for (pos in positions) {
          huang_results <- rbind(huang_results, data.frame(
            Source = "Huang_2021",
            Gene_ID = gene_id,
            Motif = m,
            Position = pos,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nHuang et al. validation results:\n")
cat("  Promoters extracted:", length(huang_promoters), "/", length(huang_promoter_genes), "\n")
cat("  Total motif occurrences found:", nrow(huang_results), "\n")

# ===== PHASE 1B: VALIDATE WANG ET AL. GENES =====
cat("\n========================================\n")
cat("PHASE 1B: VALIDATING WANG ET AL. GENES\n")
cat("========================================\n\n")

cat("Extracting promoters for Wang et al. genes...\n")

wang_results <- data.frame()
wang_promoters <- list()

pb <- txtProgressBar(min = 0, max = length(wang_genes), style = 3)
for (i in seq_along(wang_genes)) {
  gene_id <- wang_genes[i]
  prom <- extract_promoter(gene_id, genes, genome, PROMOTER_LENGTH)

  if (!is.null(prom)) {
    wang_promoters[[gene_id]] <- prom$sequence

    # Scan for each motif
    for (m in names(motifs)) {
      positions <- find_motif_positions(prom$sequence, motifs[[m]]$seq)
      if (!is.null(positions)) {
        for (pos in positions) {
          wang_results <- rbind(wang_results, data.frame(
            Source = "Wang_2024",
            Gene_ID = gene_id,
            Motif = m,
            Position = pos,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nWang et al. validation results:\n")
cat("  Promoters extracted:", length(wang_promoters), "/", length(wang_genes), "\n")
cat("  Total motif occurrences found:", nrow(wang_results), "\n")

# ===== COMBINE RESULTS =====
all_results <- rbind(huang_results, wang_results)

# ===== ANALYZE DISTRIBUTION =====
cat("\n========================================\n")
cat("POSITION DISTRIBUTION ANALYSIS\n")
cat("========================================\n\n")

# Summary by source and motif
distribution_summary <- all_results %>%
  mutate(abs_position = abs(Position)) %>%
  group_by(Source, Motif) %>%
  summarise(
    N_occurrences = n(),
    N_genes = n_distinct(Gene_ID),
    Mean_dist = round(mean(abs_position), 0),
    Median_dist = round(median(abs_position), 0),
    Pct_within_500 = round(sum(abs_position <= 500) / n() * 100, 1),
    Pct_within_1000 = round(sum(abs_position <= 1000) / n() * 100, 1),
    Pct_within_1500 = round(sum(abs_position <= 1500) / n() * 100, 1),
    Pct_within_2000 = round(sum(abs_position <= 2000) / n() * 100, 1),
    Pct_within_3000 = round(sum(abs_position <= 3000) / n() * 100, 1),
    .groups = "drop"
  )

cat("Distribution summary:\n\n")
print(as.data.frame(distribution_summary))

# Overall by source
overall_summary <- all_results %>%
  mutate(abs_position = abs(Position)) %>%
  group_by(Source) %>%
  summarise(
    N_occurrences = n(),
    N_genes = n_distinct(Gene_ID),
    Pct_within_500 = round(sum(abs_position <= 500) / n() * 100, 1),
    Pct_within_1000 = round(sum(abs_position <= 1000) / n() * 100, 1),
    Pct_within_1500 = round(sum(abs_position <= 1500) / n() * 100, 1),
    Pct_within_2000 = round(sum(abs_position <= 2000) / n() * 100, 1),
    Pct_within_3000 = round(sum(abs_position <= 3000) / n() * 100, 1),
    .groups = "drop"
  )

cat("\n\nOverall distribution by source:\n\n")
print(as.data.frame(overall_summary))

# ===== VISUALIZATIONS =====
cat("\n========================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("========================================\n\n")

# Plot 1: Distribution histogram by source
png("03_results/figures/25a_promoter_validation/position_distribution_by_source.png",
    width = 1000, height = 600, res = 120)

ggplot(all_results, aes(x = Position, fill = Source)) +
  geom_histogram(binwidth = 100, alpha = 0.7, position = "identity") +
  facet_wrap(~Source, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = c(-1000, -2000, -3000), linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Huang_2021" = "#3498DB", "Wang_2024" = "#E74C3C")) +
  labs(
    title = "Motif Position Distribution - Validation with Our Genome (Wm82.a6.v1)",
    subtitle = paste("Using", PROMOTER_LENGTH, "bp promoter window"),
    x = "Position relative to TSS (bp)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()
cat("Saved: position_distribution_by_source.png\n")

# Plot 2: Density comparison
png("03_results/figures/25a_promoter_validation/position_density_comparison.png",
    width = 900, height = 500, res = 120)

ggplot(all_results, aes(x = Position, fill = Source, color = Source)) +
  geom_density(alpha = 0.3, size = 1) +
  geom_vline(xintercept = c(-1000, -2000, -3000), linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Huang_2021" = "#3498DB", "Wang_2024" = "#E74C3C")) +
  scale_color_manual(values = c("Huang_2021" = "#3498DB", "Wang_2024" = "#E74C3C")) +
  labs(
    title = "Motif Position Density - Comparing Huang (ChIP-seq) vs Wang (DAP-seq)",
    x = "Position relative to TSS (bp)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

dev.off()
cat("Saved: position_density_comparison.png\n")

# Plot 3: By motif type
png("03_results/figures/25a_promoter_validation/position_by_motif.png",
    width = 1000, height = 700, res = 120)

ggplot(all_results, aes(x = Position, fill = Motif)) +
  geom_histogram(binwidth = 100, alpha = 0.7) +
  facet_grid(Source ~ Motif) +
  geom_vline(xintercept = -1000, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71")) +
  labs(
    title = "Motif Position Distribution by Type and Source",
    subtitle = "Red dashed line = 1000 bp upstream",
    x = "Position relative to TSS (bp)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()
cat("Saved: position_by_motif.png\n")

# Plot 4: Cumulative distribution
png("03_results/figures/25a_promoter_validation/cumulative_distribution.png",
    width = 800, height = 500, res = 120)

all_results %>%
  mutate(abs_position = abs(Position)) %>%
  ggplot(aes(x = abs_position, color = Source)) +
  stat_ecdf(size = 1.2) +
  geom_vline(xintercept = c(1000, 2000, 3000), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Huang_2021" = "#3498DB", "Wang_2024" = "#E74C3C")) +
  labs(
    title = "Cumulative Distribution of Motif Distance from TSS",
    subtitle = "Validated using our genome reference (Wm82.a6.v1)",
    x = "Distance from TSS (bp)",
    y = "Cumulative proportion"
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  annotate("text", x = 1000, y = 0.3, label = "1000 bp", angle = 90, vjust = -0.5)

dev.off()
cat("Saved: cumulative_distribution.png\n")

# ===== SAVE RESULTS =====
cat("\n========================================\n")
cat("SAVING RESULTS\n")
cat("========================================\n\n")

write.csv(all_results,
          "03_results/tables/promoter/validation/motif_positions_validation.csv",
          row.names = FALSE)
cat("Saved: motif_positions_validation.csv\n")

write.csv(distribution_summary,
          "03_results/tables/promoter/validation/distribution_summary.csv",
          row.names = FALSE)
cat("Saved: distribution_summary.csv\n")

write.csv(overall_summary,
          "03_results/tables/promoter/validation/overall_summary.csv",
          row.names = FALSE)
cat("Saved: overall_summary.csv\n")

# Save checkpoint
save(
  huang_results, wang_results, all_results,
  distribution_summary, overall_summary,
  huang_promoters, wang_promoters,
  file = "03_results/checkpoints/25a_promoter_validation.RData"
)
cat("Saved: 25a_promoter_validation.RData\n")

# ===== FINAL SUMMARY =====
cat("\n================================================================\n")
cat("  SCRIPT 25a COMPLETE: PROMOTER WINDOW VALIDATION\n")
cat("================================================================\n\n")

cat("VALIDATION SUMMARY:\n")
cat("-------------------\n")
cat("Genome reference: Wm82.a6.v1\n")
cat("Promoter window tested:", PROMOTER_LENGTH, "bp\n\n")

cat("HUANG ET AL. 2021 (ChIP-seq):\n")
huang_sum <- overall_summary[overall_summary$Source == "Huang_2021", ]
if (nrow(huang_sum) > 0) {
  cat(sprintf("  Motif occurrences: %d in %d genes\n",
              huang_sum$N_occurrences, huang_sum$N_genes))
  cat(sprintf("  Within 1000 bp: %.1f%%\n", huang_sum$Pct_within_1000))
  cat(sprintf("  Within 2000 bp: %.1f%%\n", huang_sum$Pct_within_2000))
  cat(sprintf("  Within 3000 bp: %.1f%%\n", huang_sum$Pct_within_3000))
}

cat("\nWANG ET AL. 2024 (DAP-seq):\n")
wang_sum <- overall_summary[overall_summary$Source == "Wang_2024", ]
if (nrow(wang_sum) > 0) {
  cat(sprintf("  Motif occurrences: %d in %d genes\n",
              wang_sum$N_occurrences, wang_sum$N_genes))
  cat(sprintf("  Within 1000 bp: %.1f%%\n", wang_sum$Pct_within_1000))
  cat(sprintf("  Within 2000 bp: %.1f%%\n", wang_sum$Pct_within_2000))
  cat(sprintf("  Within 3000 bp: %.1f%%\n", wang_sum$Pct_within_3000))
}

cat("\n\nCOMPARISON TO PUBLISHED VALUES:\n")
cat("Published (Huang promoter sites): 68.7% within 1000 bp, 93.5% within 2000 bp, ~100% within 3000 bp\n")
cat("Published (Wang): 47.1% within 1000 bp, 100% within 2000 bp (used 2000 bp window)\n")
cat("\nCheck if our validation matches these distributions.\n")
cat("If consistent, proceed to Phase 2 with 3000 bp window for complete coverage.\n\n")

cat("OUTPUT FILES:\n")
cat("  - 03_results/tables/promoter/validation/*.csv\n")
cat("  - 03_results/figures/25a_promoter_validation/*.png\n")
cat("  - 03_results/checkpoints/25a_promoter_validation.RData\n\n")

cat("NEXT: If validation successful, run Phase 2 on JAG1 targets\n")
