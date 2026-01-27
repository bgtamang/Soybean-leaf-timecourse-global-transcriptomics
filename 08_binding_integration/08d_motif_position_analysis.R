# Analyze WHERE in the promoter GmJAG1 binding motifs are located

rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 25c: MOTIF POSITION DISTRIBUTION ANALYSIS\n")

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/promoter", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/25c_motif_positions", recursive = TRUE, showWarnings = FALSE)

library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyr)

library(GenomicRanges)
library(rtracklayer)

cat("Loading checkpoint data...\n")
load("03_results/checkpoints/25b_promoter_motifs.RData")

# Get target gene IDs from the checkpoint
target_genes <- target_motif_counts$GeneID
cat("  Target genes:", length(target_genes), "\n")

# Load genome and annotation
cat("Loading genome and annotation...\n")
genome_file <- "01_data/genome/Gmax_880_v6.0.fa"
gff_file <- "01_data/genome/Gmax_880_Wm82.a6.v1.gene.gff3"

genome <- readDNAStringSet(genome_file)
gff <- import(gff_file)
genes <- gff[gff$type == "gene"]
cat("  Genome loaded:", length(genome), "chromosomes\n")
cat("  Genes loaded:", length(genes), "\n\n")

PROMOTER_LENGTH <- 3000  # Updated to capture 100% of promoter binding sites

extract_promoter <- function(gene_id, genes_gr, genome, upstream = 3000) {
  gene_match <- genes_gr[grepl(gene_id, genes_gr$ID) | grepl(gene_id, genes_gr$Name)]
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

  return(as.character(prom_seq))
}

cat("Extracting promoter sequences...\n")
target_promoters <- list()
pb <- txtProgressBar(min = 0, max = length(target_genes), style = 3)
for (i in seq_along(target_genes)) {
  prom <- extract_promoter(target_genes[i], genes, genome, PROMOTER_LENGTH)
  if (!is.null(prom)) target_promoters[[target_genes[i]]] <- prom
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\n  Promoters extracted:", length(target_promoters), "\n\n")

motifs <- list(
  M0 = list(seq = "GTTGGA", source = "ChIP-seq (Huang 2021)"),
  M1 = list(seq = "ACGCCACT", source = "DAP-seq (Wang 2024)"),
  M2 = list(seq = "ACTGGCAG", source = "DAP-seq (Wang 2024)")
)

find_motif_positions <- function(sequence, motif, promoter_length = 3000) {
  # Returns positions relative to TSS (negative values)
  # Position -1 is immediately upstream of TSS
  # Position -3000 is 3000 bp upstream

  seq_obj <- DNAString(sequence)
  motif_obj <- DNAString(motif)
  rc_motif <- reverseComplement(motif_obj)

  # Find forward matches
  fwd_matches <- matchPattern(motif_obj, seq_obj, fixed = FALSE)
  fwd_positions <- start(fwd_matches)

 # Find reverse complement matches
  rc_matches <- matchPattern(rc_motif, seq_obj, fixed = FALSE)
  rc_positions <- start(rc_matches)

  # Combine all positions
  all_positions <- c(fwd_positions, rc_positions)

  if (length(all_positions) == 0) {
    return(NULL)
  }

  # Convert to distance from TSS
  # Position 1 in sequence = -promoter_length from TSS
  # Position promoter_length in sequence = -1 from TSS
  seq_len <- nchar(sequence)
  distance_from_tss <- all_positions - seq_len - 1

  return(distance_from_tss)
}

cat("SCANNING FOR MOTIF POSITIONS\n")

# Store all motif positions
all_positions <- data.frame()

for (gene_id in names(target_promoters)) {
  seq <- target_promoters[[gene_id]]

  for (motif_name in names(motifs)) {
    positions <- find_motif_positions(seq, motifs[[motif_name]]$seq, PROMOTER_LENGTH)

    if (!is.null(positions) && length(positions) > 0) {
      df <- data.frame(
        GeneID = gene_id,
        Motif = motif_name,
        Position = positions,
        Source = motifs[[motif_name]]$source,
        stringsAsFactors = FALSE
      )
      all_positions <- rbind(all_positions, df)
    }
  }
}

cat("Total motif occurrences found:", nrow(all_positions), "\n")
cat("\nMotif occurrences by type:\n")
print(table(all_positions$Motif))

cat("\n========================================\n")
cat("POSITION STATISTICS\n")

position_summary <- all_positions %>%
  group_by(Motif, Source) %>%
  summarise(
    N_occurrences = n(),
    N_genes = n_distinct(GeneID),
    Mean_position = round(mean(Position), 1),
    Median_position = round(median(Position), 1),
    SD_position = round(sd(Position), 1),
    Min_position = min(Position),
    Max_position = max(Position),
    Pct_within_500bp = round(sum(Position >= -500) / n() * 100, 1),
    Pct_within_1000bp = round(sum(Position >= -1000) / n() * 100, 1),
    Pct_within_2000bp = round(sum(Position >= -2000) / n() * 100, 1),
    Pct_within_3000bp = round(sum(Position >= -3000) / n() * 100, 1),
    .groups = "drop"
  )

cat("Position summary by motif:\n\n")
print(as.data.frame(position_summary))

cat("\n========================================\n")
cat("GENERATING VISUALIZATIONS\n")

# Plot 1: Distribution histogram by motif
png("03_results/figures/25c_motif_positions/motif_position_histogram.png",
    width = 1000, height = 600, res = 120)

ggplot(all_positions, aes(x = Position, fill = Motif)) +
  geom_histogram(binwidth = 100, alpha = 0.7, position = "identity") +
  facet_wrap(~Motif, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71")) +
  geom_vline(xintercept = c(-1000, -2000), linetype = "dashed", color = "gray40") +
  labs(
    title = "Distribution of GmJAG1 Binding Motifs in Target Promoters",
    subtitle = paste("Position relative to TSS (", PROMOTER_LENGTH, " bp promoter window)", sep = ""),
    x = "Position relative to TSS (bp)",
    y = "Count",
    fill = "Motif"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

dev.off()
cat("Saved: motif_position_histogram.png\n")

# Plot 2: Density plot comparing all motifs
png("03_results/figures/25c_motif_positions/motif_position_density.png",
    width = 900, height = 500, res = 120)

ggplot(all_positions, aes(x = Position, fill = Motif, color = Motif)) +
  geom_density(alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71"),
                    labels = c("M0 (ChIP-seq)", "M1 (DAP-seq)", "M2 (DAP-seq)")) +
  scale_color_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71"),
                     labels = c("M0 (ChIP-seq)", "M1 (DAP-seq)", "M2 (DAP-seq)")) +
  geom_vline(xintercept = c(-1000, -2000, -3000), linetype = "dashed", color = "gray50") +
  annotate("text", x = -1000, y = Inf, label = "-1000 bp", vjust = 2, size = 3) +
  annotate("text", x = -2000, y = Inf, label = "-2000 bp", vjust = 2, size = 3) +
  labs(
    title = "Positional Distribution of GmJAG1 Binding Motifs",
    subtitle = "Comparing ChIP-seq (in vivo) vs DAP-seq (in vitro) derived motifs",
    x = "Position relative to TSS (bp)",
    y = "Density",
    fill = "Motif",
    color = "Motif"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

dev.off()
cat("Saved: motif_position_density.png\n")

# Plot 3: Box plot of positions
png("03_results/figures/25c_motif_positions/motif_position_boxplot.png",
    width = 700, height = 500, res = 120)

ggplot(all_positions, aes(x = Motif, y = Position, fill = Motif)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = c(-1000, -2000), linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71")) +
  labs(
    title = "Motif Position Distribution by Type",
    subtitle = paste(PROMOTER_LENGTH, "bp promoter window"),
    x = "Motif",
    y = "Position relative to TSS (bp)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_flip()

dev.off()
cat("Saved: motif_position_boxplot.png\n")

# Plot 4: Cumulative distribution
png("03_results/figures/25c_motif_positions/motif_position_cumulative.png",
    width = 800, height = 500, res = 120)

ggplot(all_positions, aes(x = Position, color = Motif)) +
  stat_ecdf(size = 1.2) +
  scale_color_manual(values = c("M0" = "#3498DB", "M1" = "#E74C3C", "M2" = "#2ECC71"),
                     labels = c("M0 (ChIP-seq)", "M1 (DAP-seq)", "M2 (DAP-seq)")) +
  geom_vline(xintercept = c(-1000, -2000, -3000), linetype = "dashed", color = "gray50") +
  labs(
    title = "Cumulative Distribution of Motif Positions",
    subtitle = paste("Proportion of motifs within given distance of TSS (", PROMOTER_LENGTH, " bp window)", sep = ""),
    x = "Position relative to TSS (bp)",
    y = "Cumulative proportion",
    color = "Motif"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

dev.off()
cat("Saved: motif_position_cumulative.png\n")

cat("\n========================================\n")
cat("BINNED POSITION ANALYSIS\n")

# Create bins (500 bp bins for 3000 bp window)
all_positions$Bin <- cut(all_positions$Position,
                         breaks = c(-3000, -2500, -2000, -1500, -1000, -500, 0),
                         labels = c("-3000 to -2500", "-2500 to -2000", "-2000 to -1500",
                                    "-1500 to -1000", "-1000 to -500", "-500 to TSS"),
                         include.lowest = TRUE)

bin_summary <- all_positions %>%
  group_by(Motif, Bin) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Motif) %>%
  mutate(Percent = round(Count / sum(Count) * 100, 1)) %>%
  ungroup()

cat("Motif distribution by promoter region:\n\n")
print(as.data.frame(bin_summary %>%
                      select(Motif, Bin, Count, Percent) %>%
                      pivot_wider(names_from = Motif, values_from = c(Count, Percent))))

# Plot 5: Stacked bar by bin
png("03_results/figures/25c_motif_positions/motif_position_bins.png",
    width = 1000, height = 500, res = 120)

ggplot(bin_summary, aes(x = Motif, y = Percent, fill = Bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  geom_text(aes(label = ifelse(Percent >= 5, paste0(round(Percent), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2.5) +
  labs(
    title = "Distribution of Motifs Across Promoter Regions",
    subtitle = paste("Percentage of motifs in each 500 bp bin (", PROMOTER_LENGTH, " bp window)", sep = ""),
    x = "Motif",
    y = "Percentage (%)",
    fill = "Region"
  ) +
  theme_minimal() +
  coord_flip()

dev.off()
cat("Saved: motif_position_bins.png\n")

cat("\n========================================\n")
cat("SAVING RESULTS\n")

write.csv(all_positions,
          "03_results/tables/promoter/motif_positions_all.csv",
          row.names = FALSE)
cat("Saved: motif_positions_all.csv\n")

write.csv(position_summary,
          "03_results/tables/promoter/motif_position_summary.csv",
          row.names = FALSE)
cat("Saved: motif_position_summary.csv\n")

write.csv(bin_summary,
          "03_results/tables/promoter/motif_position_bins.csv",
          row.names = FALSE)
cat("Saved: motif_position_bins.csv\n")

cat("\n================================================================\n")

cat("KEY FINDINGS:\n")
cat("Promoter window:", PROMOTER_LENGTH, "bp\n\n")
for (m in unique(position_summary$Motif)) {
  row <- position_summary[position_summary$Motif == m, ]
  cat(sprintf("%s (%s):\n", m, row$Source))
  cat(sprintf("  - %d occurrences in %d genes\n", row$N_occurrences, row$N_genes))
  cat(sprintf("  - Mean position: %.0f bp from TSS\n", row$Mean_position))
  cat(sprintf("  - %.1f%% within 1000 bp\n", row$Pct_within_1000bp))
  cat(sprintf("  - %.1f%% within 2000 bp\n", row$Pct_within_2000bp))
  cat(sprintf("  - %.1f%% within 3000 bp\n\n", row$Pct_within_3000bp))
}

cat("OUTPUT FILES:\n")
cat("  - 03_results/tables/promoter/motif_positions_all.csv\n")
cat("  - 03_results/tables/promoter/motif_position_summary.csv\n")
cat("  - 03_results/tables/promoter/motif_position_bins.csv\n")
cat("  - 03_results/figures/25c_motif_positions/*.png\n\n")
