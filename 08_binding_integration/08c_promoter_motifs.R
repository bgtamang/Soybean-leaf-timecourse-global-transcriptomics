# Script 25b: Promoter Motif Analysis
# Main promoter motif scanning and enrichment analysis

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 25b: PROMOTER MOTIF ANALYSIS\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/promoter", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/25b_promoter", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "Biostrings",
  "GenomicRanges",
  "rtracklayer"
)

# Check and install if needed
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("Biostrings", "GenomicRanges", "rtracklayer")) {
      cat("Installing", pkg, "from Bioconductor...\n")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  All packages loaded\n\n")

# ===== DEFINE PARAMETERS =====

PROMOTER_LENGTH <- 3000  # bp upstream of TSS (validated from Huang et al. 2021 ChIP-seq)
MIN_MOTIF_SCORE <- 0.8   # For PWM matching (if used)

# ===== DEFINE GmJAG1 BINDING MOTIFS =====

cat("========================================\n")
cat("SECTION 1: DEFINE GmJAG1 BINDING MOTIFS\n")
cat("========================================\n\n")

# These are the experimentally validated GmJAG1 binding motifs
gmjag1_motifs <- list(
  M0 = list(
    name = "GmJAG1_M0",
    sequence = "GTTGGA",
    source = "Huang et al. 2021 ChIP-seq",
    method = "In vivo binding (ChIPmentation)"
  ),
  M1 = list(
    name = "GmJAG1_M1",
    sequence = "ACGCCACT",
    source = "Wang et al. 2024 DAP-seq",
    method = "In vitro binding"
  ),
  M2 = list(
    name = "GmJAG1_M2",
    sequence = "ACTGGCAG",
    source = "Wang et al. 2024 DAP-seq",
    method = "In vitro binding"
  )
)

cat("GmJAG1 binding motifs:\n")
for (m in names(gmjag1_motifs)) {
  cat(sprintf("  %s: %s (%s)\n",
              gmjag1_motifs[[m]]$name,
              gmjag1_motifs[[m]]$sequence,
              gmjag1_motifs[[m]]$source))
}
cat("\n")

# Create motif summary table
motif_table <- data.frame(
  Motif_ID = c("M0", "M1", "M2"),
  Motif_Name = sapply(gmjag1_motifs, function(x) x$name),
  Sequence = sapply(gmjag1_motifs, function(x) x$sequence),
  Length = sapply(gmjag1_motifs, function(x) nchar(x$sequence)),
  Source = sapply(gmjag1_motifs, function(x) x$source),
  stringsAsFactors = FALSE
)
write.csv(motif_table, "03_results/tables/promoter/GmJAG1_motifs.csv", row.names = FALSE)

# ===== LOAD GENOME AND ANNOTATION =====

cat("========================================\n")
cat("SECTION 2: LOAD GENOME AND ANNOTATION\n")
cat("========================================\n\n")

# File paths
genome_file <- "01_data/genome/Gmax_880_v6.0.fa"
gff_file <- "01_data/genome/Gmax_880_Wm82.a6.v1.gene.gff3"

# Check files exist
if (!file.exists(genome_file)) {
  stop("Genome file not found: ", genome_file)
}
if (!file.exists(gff_file)) {
  stop("GFF file not found: ", gff_file)
}

cat("Loading genome sequence...\n")
cat("  File:", genome_file, "\n")
genome <- readDNAStringSet(genome_file)
cat("  Chromosomes loaded:", length(genome), "\n")
cat("  Total genome size:", sum(width(genome)), "bp\n\n")

cat("Loading gene annotation...\n")
cat("  File:", gff_file, "\n")
gff <- import(gff_file)
cat("  Total features:", length(gff), "\n")

# Filter to genes only
genes <- gff[gff$type == "gene"]
cat("  Genes found:", length(genes), "\n\n")

# ===== LOAD JAG1 TARGETS =====

cat("========================================\n")
cat("SECTION 3: LOAD JAG1 TARGETS\n")
cat("========================================\n\n")

load("03_results/checkpoints/14_JAG1_targets.RData")

jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("JAG1 targets loaded:", nrow(jag1_targets), "\n")

# Count by tier
tier_counts <- table(jag1_targets$Confidence_Tier)
cat("\nTargets by tier:\n")
print(tier_counts)
cat("\n")

# Get target gene IDs
target_genes <- unique(jag1_targets$GeneID)

# ===== EXTRACT PROMOTER SEQUENCES =====

cat("========================================\n")
cat("SECTION 4: EXTRACT PROMOTER SEQUENCES\n")
cat("========================================\n\n")

cat("Extracting", PROMOTER_LENGTH, "bp upstream of TSS...\n\n")

# Function to extract promoter for a gene
extract_promoter <- function(gene_id, genes_gr, genome, upstream = 3000) {
  # Find the gene
  gene_match <- genes_gr[grepl(gene_id, genes_gr$ID) | grepl(gene_id, genes_gr$Name)]

  if (length(gene_match) == 0) {
    return(NULL)
  }

  gene <- gene_match[1]
  chr <- as.character(seqnames(gene))
  strand_val <- as.character(strand(gene))

  # Check if chromosome exists in genome
  if (!chr %in% names(genome)) {
    return(NULL)
  }

  chr_seq <- genome[[chr]]
  chr_len <- length(chr_seq)

  if (strand_val == "+") {
    # Promoter is upstream of start
    prom_end <- start(gene) - 1
    prom_start <- max(1, prom_end - upstream + 1)
  } else {
    # Promoter is downstream of end (on reverse strand)
    prom_start <- end(gene) + 1
    prom_end <- min(chr_len, prom_start + upstream - 1)
  }

  if (prom_start > prom_end || prom_start < 1 || prom_end > chr_len) {
    return(NULL)
  }

  # Extract sequence
  prom_seq <- subseq(chr_seq, start = prom_start, end = prom_end)

  # Reverse complement if on minus strand
  if (strand_val == "-") {
    prom_seq <- reverseComplement(prom_seq)
  }

  return(as.character(prom_seq))
}

# Extract promoters for all targets
cat("Extracting promoters for JAG1 targets...\n")
target_promoters <- list()
successful <- 0
failed <- 0

pb <- txtProgressBar(min = 0, max = length(target_genes), style = 3)
for (i in seq_along(target_genes)) {
  gene_id <- target_genes[i]
  prom <- extract_promoter(gene_id, genes, genome, PROMOTER_LENGTH)
  if (!is.null(prom)) {
    target_promoters[[gene_id]] <- prom
    successful <- successful + 1
  } else {
    failed <- failed + 1
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nPromoter extraction results:\n")
cat("  Successful:", successful, "\n")
cat("  Failed:", failed, "\n")
cat("  Success rate:", round(successful / length(target_genes) * 100, 1), "%\n\n")

# Extract promoters for background (all genes)
cat("Extracting promoters for background genes...\n")
all_gene_ids <- unique(genes$ID)
# Remove version numbers if present
all_gene_ids <- gsub("\\..*", "", all_gene_ids)

# Sample background (use all or subset)
set.seed(42)
if (length(all_gene_ids) > 10000) {
  bg_gene_ids <- sample(all_gene_ids, 10000)
} else {
  bg_gene_ids <- all_gene_ids
}

# Remove targets from background
bg_gene_ids <- setdiff(bg_gene_ids, target_genes)

bg_promoters <- list()
pb <- txtProgressBar(min = 0, max = length(bg_gene_ids), style = 3)
for (i in seq_along(bg_gene_ids)) {
  gene_id <- bg_gene_ids[i]
  prom <- extract_promoter(gene_id, genes, genome, PROMOTER_LENGTH)
  if (!is.null(prom)) {
    bg_promoters[[gene_id]] <- prom
  }
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nBackground promoters extracted:", length(bg_promoters), "\n\n")

# ===== SCAN FOR MOTIFS =====

cat("========================================\n")
cat("SECTION 5: SCAN FOR GmJAG1 MOTIFS\n")
cat("========================================\n\n")

# Function to count motif occurrences
count_motif <- function(sequence, motif, allow_rc = TRUE) {
  seq_obj <- DNAString(sequence)
  motif_obj <- DNAString(motif)

  # Forward strand
  fwd_matches <- countPattern(motif_obj, seq_obj, fixed = FALSE)

  # Reverse complement
  if (allow_rc) {
    rc_motif <- reverseComplement(motif_obj)
    rc_matches <- countPattern(rc_motif, seq_obj, fixed = FALSE)
    return(fwd_matches + rc_matches)
  }

  return(fwd_matches)
}

# Scan all target promoters
cat("Scanning target promoters for GmJAG1 motifs...\n")
target_motif_counts <- data.frame(
  GeneID = names(target_promoters),
  Promoter_Length = sapply(target_promoters, nchar),
  M0_count = NA,
  M1_count = NA,
  M2_count = NA,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(target_motif_counts))) {
  seq <- target_promoters[[target_motif_counts$GeneID[i]]]
  target_motif_counts$M0_count[i] <- count_motif(seq, gmjag1_motifs$M0$sequence)
  target_motif_counts$M1_count[i] <- count_motif(seq, gmjag1_motifs$M1$sequence)
  target_motif_counts$M2_count[i] <- count_motif(seq, gmjag1_motifs$M2$sequence)
}

# Calculate derived metrics
target_motif_counts$Total_motifs <- target_motif_counts$M0_count +
                                     target_motif_counts$M1_count +
                                     target_motif_counts$M2_count
target_motif_counts$Has_M0 <- target_motif_counts$M0_count > 0
target_motif_counts$Has_M1 <- target_motif_counts$M1_count > 0
target_motif_counts$Has_M2 <- target_motif_counts$M2_count > 0
target_motif_counts$Has_any_motif <- target_motif_counts$Total_motifs > 0
target_motif_counts$N_motif_types <- (target_motif_counts$Has_M0 +
                                       target_motif_counts$Has_M1 +
                                       target_motif_counts$Has_M2)

# Add tier information
target_motif_counts$Tier <- jag1_targets$Confidence_Tier[match(
  target_motif_counts$GeneID, jag1_targets$GeneID)]

cat("Target motif scanning complete\n\n")

# Scan background promoters
cat("Scanning background promoters for comparison...\n")
bg_motif_counts <- data.frame(
  GeneID = names(bg_promoters),
  Promoter_Length = sapply(bg_promoters, nchar),
  M0_count = NA,
  M1_count = NA,
  M2_count = NA,
  stringsAsFactors = FALSE
)

pb <- txtProgressBar(min = 0, max = nrow(bg_motif_counts), style = 3)
for (i in seq_len(nrow(bg_motif_counts))) {
  seq <- bg_promoters[[bg_motif_counts$GeneID[i]]]
  bg_motif_counts$M0_count[i] <- count_motif(seq, gmjag1_motifs$M0$sequence)
  bg_motif_counts$M1_count[i] <- count_motif(seq, gmjag1_motifs$M1$sequence)
  bg_motif_counts$M2_count[i] <- count_motif(seq, gmjag1_motifs$M2$sequence)
  setTxtProgressBar(pb, i)
}
close(pb)

bg_motif_counts$Total_motifs <- bg_motif_counts$M0_count +
                                 bg_motif_counts$M1_count +
                                 bg_motif_counts$M2_count
bg_motif_counts$Has_any_motif <- bg_motif_counts$Total_motifs > 0

cat("\nBackground motif scanning complete\n\n")

# ===== ENRICHMENT ANALYSIS =====

cat("========================================\n")
cat("SECTION 6: MOTIF ENRICHMENT ANALYSIS\n")
cat("========================================\n\n")

# Calculate enrichment for each motif
enrichment_results <- data.frame(
  Motif = c("M0", "M1", "M2", "Any_motif"),
  Motif_Sequence = c(gmjag1_motifs$M0$sequence,
                     gmjag1_motifs$M1$sequence,
                     gmjag1_motifs$M2$sequence,
                     "Any"),
  Targets_with_motif = NA,
  Targets_total = nrow(target_motif_counts),
  Target_percent = NA,
  Background_with_motif = NA,
  Background_total = nrow(bg_motif_counts),
  Background_percent = NA,
  Fold_enrichment = NA,
  P_value = NA,
  stringsAsFactors = FALSE
)

# M0 enrichment
enrichment_results$Targets_with_motif[1] <- sum(target_motif_counts$Has_M0)
enrichment_results$Background_with_motif[1] <- sum(bg_motif_counts$M0_count > 0)

# M1 enrichment
enrichment_results$Targets_with_motif[2] <- sum(target_motif_counts$Has_M1)
enrichment_results$Background_with_motif[2] <- sum(bg_motif_counts$M1_count > 0)

# M2 enrichment
enrichment_results$Targets_with_motif[3] <- sum(target_motif_counts$Has_M2)
enrichment_results$Background_with_motif[3] <- sum(bg_motif_counts$M2_count > 0)

# Any motif enrichment
enrichment_results$Targets_with_motif[4] <- sum(target_motif_counts$Has_any_motif)
enrichment_results$Background_with_motif[4] <- sum(bg_motif_counts$Has_any_motif)

# Calculate percentages and enrichment
for (i in 1:4) {
  enrichment_results$Target_percent[i] <- round(
    enrichment_results$Targets_with_motif[i] / enrichment_results$Targets_total[i] * 100, 1)
  enrichment_results$Background_percent[i] <- round(
    enrichment_results$Background_with_motif[i] / enrichment_results$Background_total[i] * 100, 1)

  # Fold enrichment
  target_rate <- enrichment_results$Target_percent[i]
  bg_rate <- enrichment_results$Background_percent[i]
  enrichment_results$Fold_enrichment[i] <- round(target_rate / max(bg_rate, 0.1), 2)

  # Fisher's exact test
  cont_table <- matrix(c(
    enrichment_results$Targets_with_motif[i],
    enrichment_results$Targets_total[i] - enrichment_results$Targets_with_motif[i],
    enrichment_results$Background_with_motif[i],
    enrichment_results$Background_total[i] - enrichment_results$Background_with_motif[i]
  ), nrow = 2)

  fisher_result <- fisher.test(cont_table)
  enrichment_results$P_value[i] <- fisher_result$p.value
}

enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
enrichment_results$Significant <- enrichment_results$FDR < 0.05

cat("MOTIF ENRICHMENT RESULTS:\n")
cat("========================\n\n")
print(enrichment_results[, c("Motif", "Motif_Sequence", "Target_percent",
                              "Background_percent", "Fold_enrichment", "P_value")])
cat("\n")

write.csv(enrichment_results, "03_results/tables/promoter/motif_enrichment.csv", row.names = FALSE)

# ===== TIER-SPECIFIC ANALYSIS =====

cat("========================================\n")
cat("SECTION 7: TIER-SPECIFIC MOTIF ANALYSIS\n")
cat("========================================\n\n")

tier_motif_summary <- target_motif_counts %>%
  group_by(Tier) %>%
  summarise(
    N_genes = n(),
    Pct_with_M0 = round(mean(Has_M0) * 100, 1),
    Pct_with_M1 = round(mean(Has_M1) * 100, 1),
    Pct_with_M2 = round(mean(Has_M2) * 100, 1),
    Pct_with_any = round(mean(Has_any_motif) * 100, 1),
    Mean_total_motifs = round(mean(Total_motifs), 2),
    Mean_motif_types = round(mean(N_motif_types), 2),
    .groups = "drop"
  )

cat("Motif presence by confidence tier:\n\n")
print(as.data.frame(tier_motif_summary))
cat("\n")

write.csv(tier_motif_summary, "03_results/tables/promoter/tier_motif_summary.csv", row.names = FALSE)

# ===== HIGH PRIORITY TARGETS =====

cat("========================================\n")
cat("SECTION 8: HIGH PRIORITY TARGETS\n")
cat("========================================\n\n")

# Targets with multiple motif types
high_motif_targets <- target_motif_counts[target_motif_counts$N_motif_types >= 2, ]
cat("Targets with 2+ motif types:", nrow(high_motif_targets), "\n")

# Targets with all 3 motifs
all_motif_targets <- target_motif_counts[target_motif_counts$N_motif_types == 3, ]
cat("Targets with ALL 3 motif types:", nrow(all_motif_targets), "\n\n")

if (nrow(all_motif_targets) > 0) {
  cat("Genes with all 3 GmJAG1 motifs:\n")
  print(all_motif_targets[order(-all_motif_targets$Total_motifs),
                          c("GeneID", "Tier", "M0_count", "M1_count", "M2_count", "Total_motifs")])
}

write.csv(target_motif_counts, "03_results/tables/promoter/target_motif_counts.csv", row.names = FALSE)

# ===== VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 9: VISUALIZATIONS\n")
cat("========================================\n\n")

# Plot 1: Motif enrichment bar plot
png("03_results/figures/25b_promoter/motif_enrichment.png",
    width = 800, height = 600, res = 120)

plot_data <- enrichment_results[1:3, ]  # Exclude "Any"
plot_data$Motif <- factor(plot_data$Motif, levels = c("M0", "M1", "M2"))

barplot_data <- rbind(
  data.frame(Motif = plot_data$Motif, Percent = plot_data$Target_percent, Group = "JAG1 Targets"),
  data.frame(Motif = plot_data$Motif, Percent = plot_data$Background_percent, Group = "Background")
)

ggplot(barplot_data, aes(x = Motif, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0(Percent, "%")),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("JAG1 Targets" = "#E74C3C", "Background" = "#95A5A6")) +
  labs(title = "GmJAG1 Motif Enrichment in Target Promoters",
       subtitle = paste("Promoter region:", PROMOTER_LENGTH, "bp upstream of TSS"),
       x = "Motif",
       y = "Genes with Motif (%)",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "top")

dev.off()
cat("Saved: motif_enrichment.png\n")

# Plot 2: Tier-specific motif presence
png("03_results/figures/25b_promoter/tier_motif_presence.png",
    width = 900, height = 600, res = 120)

tier_plot_data <- tier_motif_summary %>%
  pivot_longer(cols = starts_with("Pct_with_M"),
               names_to = "Motif", values_to = "Percent") %>%
  mutate(Motif = gsub("Pct_with_", "", Motif))

tier_plot_data$Tier <- factor(tier_plot_data$Tier, levels = c("Gold", "Silver", "Bronze"))

ggplot(tier_plot_data, aes(x = Tier, y = Percent, fill = Motif)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0(Percent, "%")),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "GmJAG1 Motif Presence by Target Confidence Tier",
       x = "Confidence Tier",
       y = "Genes with Motif (%)",
       fill = "Motif") +
  theme_minimal() +
  theme(legend.position = "top")

dev.off()
cat("Saved: tier_motif_presence.png\n")

# Plot 3: Total motif count distribution
png("03_results/figures/25b_promoter/motif_count_distribution.png",
    width = 800, height = 600, res = 120)

ggplot(target_motif_counts, aes(x = Total_motifs, fill = Tier)) +
  geom_histogram(binwidth = 1, position = "stack", color = "white") +
  scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")) +
  labs(title = "Distribution of GmJAG1 Motif Counts in Target Promoters",
       x = "Total Motif Count (M0 + M1 + M2)",
       y = "Number of Genes",
       fill = "Tier") +
  theme_minimal()

dev.off()
cat("Saved: motif_count_distribution.png\n")

# Plot 4: Venn-style summary
png("03_results/figures/25b_promoter/motif_presence_summary.png",
    width = 800, height = 600, res = 120)

par(mar = c(5, 6, 4, 2))
motif_summary <- data.frame(
  Motif = c("M0\n(GTTGGA)", "M1\n(ACGCCACT)", "M2\n(ACTGGCAG)", "Any Motif"),
  Percent = c(
    mean(target_motif_counts$Has_M0) * 100,
    mean(target_motif_counts$Has_M1) * 100,
    mean(target_motif_counts$Has_M2) * 100,
    mean(target_motif_counts$Has_any_motif) * 100
  )
)

barplot(motif_summary$Percent,
        names.arg = motif_summary$Motif,
        col = c("#3498DB", "#E74C3C", "#2ECC71", "#9B59B6"),
        main = "GmJAG1 Binding Motifs in Target Promoters",
        ylab = "Targets with Motif (%)",
        ylim = c(0, max(motif_summary$Percent) * 1.2),
        las = 1)

# Add percentage labels
text(x = 1:4 * 1.2 - 0.5, y = motif_summary$Percent + 2,
     labels = paste0(round(motif_summary$Percent, 1), "%"), cex = 0.9)

# Add enrichment fold
abline(h = mean(bg_motif_counts$Has_any_motif) * 100, lty = 2, col = "gray")
text(4.5, mean(bg_motif_counts$Has_any_motif) * 100 + 1, "Background", cex = 0.8, col = "gray40")

dev.off()
cat("Saved: motif_presence_summary.png\n")

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 10: SAVE CHECKPOINT\n")
cat("========================================\n\n")

promoter_analysis <- list(
  motifs = gmjag1_motifs,
  target_motif_counts = target_motif_counts,
  bg_motif_counts = bg_motif_counts,
  enrichment_results = enrichment_results,
  tier_summary = tier_motif_summary,
  parameters = list(
    promoter_length = PROMOTER_LENGTH,
    n_targets = nrow(target_motif_counts),
    n_background = nrow(bg_motif_counts)
  )
)

save(
  promoter_analysis,
  target_motif_counts,
  bg_motif_counts,
  enrichment_results,
  tier_motif_summary,
  gmjag1_motifs,
  file = "03_results/checkpoints/25b_promoter_motifs.RData"
)

cat("Checkpoint saved: 25b_promoter_motifs.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 25b COMPLETE: PROMOTER MOTIF ANALYSIS\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("ANALYSIS SUMMARY:\n")
cat("  Promoter length analyzed:", PROMOTER_LENGTH, "bp\n")
cat("  Target promoters scanned:", nrow(target_motif_counts), "\n")
cat("  Background promoters scanned:", nrow(bg_motif_counts), "\n\n")

cat("MOTIF ENRICHMENT RESULTS:\n")
for (i in 1:3) {
  cat(sprintf("  %s (%s): %.1f%% targets vs %.1f%% background (%.1fx enrichment, p=%.2e)\n",
              enrichment_results$Motif[i],
              enrichment_results$Motif_Sequence[i],
              enrichment_results$Target_percent[i],
              enrichment_results$Background_percent[i],
              enrichment_results$Fold_enrichment[i],
              enrichment_results$P_value[i]))
}

cat("\nKEY FINDINGS:\n")
cat("  Targets with any GmJAG1 motif:", sum(target_motif_counts$Has_any_motif),
    "(", round(mean(target_motif_counts$Has_any_motif) * 100, 1), "%)\n")
cat("  Targets with 2+ motif types:", sum(target_motif_counts$N_motif_types >= 2), "\n")
cat("  Targets with ALL 3 motifs:", sum(target_motif_counts$N_motif_types == 3), "\n\n")

cat("OUTPUT FILES:\n")
cat("  - 03_results/tables/promoter/GmJAG1_motifs.csv\n")
cat("  - 03_results/tables/promoter/motif_enrichment.csv\n")
cat("  - 03_results/tables/promoter/tier_motif_summary.csv\n")
cat("  - 03_results/tables/promoter/target_motif_counts.csv\n")
cat("  - 03_results/figures/25b_promoter/*.png\n\n")

cat("NEXT: Run Script 25c for motif position distribution analysis\n\n")
