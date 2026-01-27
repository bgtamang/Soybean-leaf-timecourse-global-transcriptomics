
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 11: DIFFERENTIAL EXPRESSION ANALYSIS\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/11_DE", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "dplyr",
  "tidyr"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")
load("03_results/checkpoints/10_DE_setup.RData")
cat("Loaded checkpoint from script 10\n")
cat("  DGEList:", nrow(dge), "genes x", ncol(dge), "samples\n")
cat("  Design matrix:", nrow(design), "x", ncol(design), "\n")
cat("  Contrast matrix:", ncol(contrast_matrix), "contrasts\n\n")
cat("Estimating dispersions with robust=TRUE...\n")
cat("  (This accounts for outlier genes)\n\n")

# Estimate dispersions
dge <- estimateDisp(dge, design, robust = TRUE)

# Report dispersion values
cat("Dispersion estimates:\n")
cat("  Common dispersion:", round(dge$common.dispersion, 4), "\n")
cat("  BCV (sqrt of common dispersion):", round(sqrt(dge$common.dispersion), 4), "\n")
cat("  Trend method:", dge$trend.method, "\n\n")

# --- Plot 1: BCV plot ---
cat("Creating BCV plot...\n")

png("03_results/figures/11_DE/BCV_plot.png",
    width = 8, height = 6, units = "in", res = 300)

plotBCV(dge, main = "Biological Coefficient of Variation")

dev.off()
cat("  Saved: 03_results/figures/11_DE/BCV_plot.png\n\n")
cat("Fitting quasi-likelihood GLM with robust=TRUE...\n")

# Fit the model
fit <- glmQLFit(dge, design, robust = TRUE)

cat("  Model fit complete\n")

# Prior df can be a vector with robust=TRUE, so report median
if (is.numeric(fit$prior.df)) {
  if (length(fit$prior.df) == 1) {
    cat("  Prior df:", round(fit$prior.df, 2), "\n\n")
  } else {
    cat("  Prior df (median):", round(median(fit$prior.df), 2), "\n\n")
  }
} else {
  cat("  Prior df: (robust estimation used)\n\n")
}

# --- Plot 2: QL dispersion plot ---
cat("Creating QL dispersion plot...\n")

png("03_results/figures/11_DE/QL_dispersion.png",
    width = 8, height = 6, units = "in", res = 300)

plotQLDisp(fit, main = "Quasi-Likelihood Dispersion")

dev.off()
cat("  Saved: 03_results/figures/11_DE/QL_dispersion.png\n\n")
cat("Testing", ncol(contrast_matrix), "contrasts...\n")
cat("  Running both glmQLFTest (standard) and glmTreat (FC threshold)...\n\n")

# Define fold-change threshold for glmTreat (log2(1.2) ≈ 0.263)
LFC_THRESHOLD <- log2(1.2)
cat("  glmTreat LFC threshold:", round(LFC_THRESHOLD, 3), "(FC = 1.2)\n\n")

# Store all results
de_results <- list()        # Standard QL F-test results
de_results_treat <- list()  # glmTreat results (with FC threshold)

de_summary <- data.frame(
  Contrast = character(),
  Total_DE = integer(),
  Up = integer(),
  Down = integer(),
  Total_DE_treat = integer(),
  Up_treat = integer(),
  Down_treat = integer(),
  stringsAsFactors = FALSE
)

# Test each contrast
for (i in 1:ncol(contrast_matrix)) {
  contrast_name <- colnames(contrast_matrix)[i]

  # Run standard QL F-test
  qlf <- glmQLFTest(fit, contrast = contrast_matrix[, i])

  # Run glmTreat with FC threshold (more stringent)
  treat <- glmTreat(fit, contrast = contrast_matrix[, i], lfc = LFC_THRESHOLD)

  # Get all results (not just top)
  results <- topTags(qlf, n = Inf, sort.by = "none")$table
  results_treat <- topTags(treat, n = Inf, sort.by = "none")$table

  # Add gene names as column
  results$Gene <- rownames(results)
  results_treat$Gene <- rownames(results_treat)

  # Store results
  de_results[[contrast_name]] <- results
  de_results_treat[[contrast_name]] <- results_treat

  # Count significant genes - Standard (FDR < 0.05, |logFC| > 1)
  sig <- results$FDR < 0.05 & abs(results$logFC) > 1
  up <- sig & results$logFC > 0
  down <- sig & results$logFC < 0

  # Count significant genes - Treat (FDR < 0.05, tests FC > threshold)
  sig_treat <- results_treat$FDR < 0.05 & abs(results_treat$logFC) > LFC_THRESHOLD
  up_treat <- sig_treat & results_treat$logFC > 0
  down_treat <- sig_treat & results_treat$logFC < 0

  de_summary <- rbind(de_summary, data.frame(
    Contrast = contrast_name,
    Total_DE = sum(sig),
    Up = sum(up),
    Down = sum(down),
    Total_DE_treat = sum(sig_treat),
    Up_treat = sum(up_treat),
    Down_treat = sum(down_treat),
    stringsAsFactors = FALSE
  ))

  # Progress update every 10 contrasts
  if (i %% 10 == 0 || i == ncol(contrast_matrix)) {
    cat("  Completed", i, "of", ncol(contrast_matrix), "contrasts\n")
  }
}

cat("\nDE testing complete!\n\n")
cat("Creating TestResults objects for Venn diagrams and summaries...\n")

# Create decision matrices using decideTests
# Standard tests (FDR < 0.05)
edgeR_coded <- do.call(cbind, lapply(de_results, function(x) {
  # Manual decision: -1, 0, 1
  decision <- rep(0, nrow(x))
  decision[x$FDR < 0.05 & x$logFC > 1] <- 1
  decision[x$FDR < 0.05 & x$logFC < -1] <- -1
  return(decision)
}))

# Treat tests (more stringent)
edgeR_tr_coded <- do.call(cbind, lapply(names(de_results_treat), function(contrast_name) {
  x <- de_results_treat[[contrast_name]]
  decision <- rep(0, nrow(x))
  decision[x$FDR < 0.05 & x$logFC > LFC_THRESHOLD] <- 1
  decision[x$FDR < 0.05 & x$logFC < -LFC_THRESHOLD] <- -1
  return(decision)
}))

# Set row and column names
rownames(edgeR_coded) <- rownames(de_results[[1]])
colnames(edgeR_coded) <- names(de_results)
rownames(edgeR_tr_coded) <- rownames(de_results_treat[[1]])
colnames(edgeR_tr_coded) <- names(de_results_treat)

# Convert to TestResults objects (for vennDiagram compatibility)
edgeR_coded <- new("TestResults", edgeR_coded)
edgeR_tr_coded <- new("TestResults", edgeR_tr_coded)

# Set labels
attr(edgeR_coded, "labels") <- c("Down", "NotSig", "Up")
attr(edgeR_tr_coded, "labels") <- c("Down", "NotSig", "Up")

cat("  Standard TestResults:", nrow(edgeR_coded), "genes x", ncol(edgeR_coded), "contrasts\n")
cat("  Treat TestResults:", nrow(edgeR_tr_coded), "genes x", ncol(edgeR_tr_coded), "contrasts\n\n")

# Summary of TestResults
cat("Summary of DE genes (glmTreat, FC > 1.2):\n")
print(summary(edgeR_tr_coded)[, 1:min(10, ncol(edgeR_tr_coded))])
cat("  ... (showing first 10 contrasts)\n\n")
# Add contrast type from description
de_summary <- merge(de_summary, contrast_desc[, c("Contrast", "Type")],
                    by = "Contrast", all.x = TRUE)

# Summary by contrast type
cat("DE genes by contrast type (FDR < 0.05, |logFC| > 1):\n\n")

summary_by_type <- de_summary %>%
  group_by(Type) %>%
  summarize(
    N_contrasts = n(),
    Mean_DE = round(mean(Total_DE)),
    Total_Up = sum(Up),
    Total_Down = sum(Down),
    .groups = "drop"
  )

print(as.data.frame(summary_by_type))

# Top contrasts by number of DE genes
cat("\n\nTop 10 contrasts by number of DE genes:\n")
top_contrasts <- de_summary %>%
  arrange(desc(Total_DE)) %>%
  head(10)
print(top_contrasts[, c("Contrast", "Total_DE", "Up", "Down", "Type")])

# Key contrasts (Narrow vs Broad at TP0)
cat("\n\nKey contrasts for JAG1 analysis (Narrow vs Broad at TP0):\n")
key_results <- de_summary[de_summary$Contrast %in% key_contrasts, ]
if (nrow(key_results) > 0) {
  print(key_results[, c("Contrast", "Total_DE", "Up", "Down")])
}

# Pooled comparison
cat("\n\nPooled Narrow vs Broad at TP0:\n")
pooled_result <- de_summary[de_summary$Contrast == "NarrowvsBroad_TP0", ]
if (nrow(pooled_result) > 0) {
  print(pooled_result[, c("Contrast", "Total_DE", "Up", "Down")])
}


cat("\n========================================\n")

# --- Plot 3: DE gene counts barplot ---
cat("Creating DE summary barplot...\n")

# Prepare data for plotting (temporal vs TP0 contrasts)
temporal_de <- de_summary[de_summary$Type == "Temporal_vsTP0", ]

if (nrow(temporal_de) > 0) {
  # Extract line and timepoint from contrast name
  temporal_de <- temporal_de %>%
    mutate(
      Line = gsub("_TP[0-9]vsTP0", "", Contrast),
      Timepoint = gsub(".*_(TP[0-9])vsTP0", "\\1", Contrast)
    )

  # Reshape for plotting
  temporal_long <- temporal_de %>%
    select(Line, Timepoint, Up, Down) %>%
    pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
    mutate(Count = ifelse(Direction == "Down", -Count, Count))

  p1 <- ggplot(temporal_long, aes(x = Timepoint, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "identity") +
    facet_wrap(~Line, ncol = 2) +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
    geom_hline(yintercept = 0, color = "black") +
    theme_bw(base_size = 12) +
    labs(title = "Differentially Expressed Genes (vs TP0)",
         subtitle = "FDR < 0.05, |logFC| > 1",
         x = "Timepoint",
         y = "Number of DE Genes") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"))

  ggsave("03_results/figures/11_DE/DE_barplot_temporal.png", p1,
         width = 10, height = 8, dpi = 300)
  cat("  Saved: 03_results/figures/11_DE/DE_barplot_temporal.png\n")
}

# --- Plot 4: Between-line DE at TP0 ---
cat("Creating between-line DE barplot...\n")

between_de <- de_summary[de_summary$Type == "Between_line_TP0", ]

if (nrow(between_de) > 0) {
  between_de$Contrast <- factor(between_de$Contrast, levels = between_de$Contrast)

  between_long <- between_de %>%
    select(Contrast, Up, Down) %>%
    pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
    mutate(Count = ifelse(Direction == "Down", -Count, Count))

  p2 <- ggplot(between_long, aes(x = Contrast, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
    geom_hline(yintercept = 0, color = "black") +
    theme_bw(base_size = 12) +
    labs(title = "DE Genes: Narrow vs Broad Lines at TP0",
         subtitle = "Positive = Higher in Narrow (potential JAG1 targets)",
         x = "Comparison",
         y = "Number of DE Genes") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("03_results/figures/11_DE/DE_barplot_between_lines.png", p2,
         width = 10, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/11_DE/DE_barplot_between_lines.png\n")
}


cat("\n========================================\n")

# Save DE summary
write.csv(de_summary, "03_results/tables/DE_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/DE_summary.csv\n")

# Save checkpoint (includes both standard and treat results)
save(dge, design, contrast_matrix, contrast_desc, fit,
     de_results, de_results_treat,  # Both result types
     de_summary,
     edgeR_coded, edgeR_tr_coded,   # TestResults objects
     LFC_THRESHOLD,
     key_contrasts, broad_lines, narrow_lines,
     file = "03_results/checkpoints/11_DE_results.RData")
cat("Saved: 03_results/checkpoints/11_DE_results.RData\n")
cat("  Includes: de_results (standard), de_results_treat (FC threshold)\n")
cat("  Includes: edgeR_coded, edgeR_tr_coded (TestResults objects)\n")


cat("\n========================================\n")
cat("SESSION INFO\n")

print(sessionInfo())


cat("\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Contrasts tested:", ncol(contrast_matrix), "\n")
cat("    - Total DE genes across all contrasts:", sum(de_summary$Total_DE), "\n")
cat("    - Mean DE genes per contrast:", round(mean(de_summary$Total_DE)), "\n")
cat("\n")
cat("  Key findings (Narrow vs Broad at TP0):\n")
for (kc in key_contrasts) {
  res <- de_summary[de_summary$Contrast == kc, ]
  if (nrow(res) > 0) {
    cat("    -", kc, ":", res$Total_DE, "DE (", res$Up, "up,", res$Down, "down)\n")
  }
}
cat("\n")
