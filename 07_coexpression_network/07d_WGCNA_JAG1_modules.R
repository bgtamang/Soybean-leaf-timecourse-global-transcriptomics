
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 21: JAG1 MODULE ANALYSIS & HUB GENES\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/21_WGCNA_JAG1", recursive = TRUE, showWarnings = FALSE)
required_packages <- c(
  "WGCNA",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "reshape2"
)

invisible(lapply(required_packages, library, character.only = TRUE))
enableWGCNAThreads()
cat("  Packages loaded\n\n")
load("03_results/checkpoints/20_WGCNA_traits.RData")
cat("Loaded WGCNA trait data\n")

# Load JAG1 targets from checkpoint 14
load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n")

# target_table contains all genes with Confidence_Tier column
# Filter to get actual targets
jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("JAG1 targets:", nrow(jag1_targets), "\n\n")

# Key gene IDs
JAG1_ID <- "Glyma.20G116200"
JAG2_ID <- "Glyma.10G273800"

cat("Key genes:\n")
cat("  JAG1:", JAG1_ID, "\n")
cat("  JAG2:", JAG2_ID, "\n\n")
# Calculate module membership (correlation of gene expression with module eigengene)
# kME = correlation between gene expression and module eigengene

MEs_mat <- as.matrix(MEs)

# Calculate kME for all genes
# datExpr is samples x genes, MEs_mat is samples x modules
# cor() will correlate each gene (column) with each module eigengene (column)
kMEs <- cor(datExpr, MEs_mat, use = "pairwise.complete.obs")
colnames(kMEs) <- paste0("kME_", gsub("ME", "", colnames(MEs)))

cat("Module membership (kME) calculated\n")
cat("  Genes:", nrow(kMEs), "\n")
cat("  Modules:", ncol(kMEs), "\n\n")
# Find JAG1's module
genes_in_wgcna <- colnames(datExpr)

if (JAG1_ID %in% genes_in_wgcna) {
  JAG1_module <- net$colors[genes_in_wgcna == JAG1_ID]
  JAG1_module_color <- labels2colors(JAG1_module)

  cat("JAG1 module assignment:\n")
  cat("  Module number:", JAG1_module, "\n")
  cat("  Module color:", JAG1_module_color, "\n")

  # Get JAG1's kME values
  JAG1_kME <- kMEs[JAG1_ID, ]
  cat("\n  JAG1 module membership (top 5):\n")
  top_kME <- sort(JAG1_kME, decreasing = TRUE)[1:5]
  for (i in 1:5) {
    cat(sprintf("    %s: %.3f\n", names(top_kME)[i], top_kME[i]))
  }

  # JAG1's membership in its own module
  own_module_kme <- paste0("kME_", JAG1_module)
  cat(sprintf("\n  JAG1's kME in own module (%s): %.3f\n",
              JAG1_module_color, JAG1_kME[own_module_kme]))

} else {
  cat("JAG1 not found in WGCNA gene set\n")
  cat("This may occur if JAG1 didn't pass variance/FDR filters\n")
  JAG1_module <- NA
  JAG1_module_color <- "unknown"
}

# Check JAG2
if (JAG2_ID %in% genes_in_wgcna) {
  JAG2_module <- net$colors[genes_in_wgcna == JAG2_ID]
  JAG2_module_color <- labels2colors(JAG2_module)
  cat("\nJAG2 module assignment:\n")
  cat("  Module number:", JAG2_module, "\n")
  cat("  Module color:", JAG2_module_color, "\n")
} else {
  cat("\nJAG2 not found in WGCNA gene set\n")
  JAG2_module <- NA
}


cat("\n========================================\n")

# Get JAG1 targets
jag1_target_genes <- unique(jag1_targets$GeneID)
cat("Total JAG1 targets:", length(jag1_target_genes), "\n")

# Find which targets are in WGCNA
targets_in_wgcna <- intersect(jag1_target_genes, genes_in_wgcna)
cat("JAG1 targets in WGCNA:", length(targets_in_wgcna), "\n\n")

# Get module assignments for targets
target_modules <- net$colors[match(targets_in_wgcna, genes_in_wgcna)]
names(target_modules) <- targets_in_wgcna

# Count targets per module
target_counts <- table(target_modules)

# Calculate enrichment using Fisher's exact test
all_module_counts <- table(net$colors)
n_wgcna_genes <- length(net$colors)
n_targets <- length(target_modules)

enrichment_results <- data.frame()

for (mod in names(all_module_counts)) {
  # 2x2 contingency table
  in_module_target <- sum(target_modules == mod)
  in_module_not_target <- all_module_counts[mod] - in_module_target
  not_module_target <- n_targets - in_module_target
  not_module_not_target <- n_wgcna_genes - n_targets - in_module_not_target

  cont_table <- matrix(c(in_module_target, in_module_not_target,
                         not_module_target, not_module_not_target),
                       nrow = 2)

  # Fisher's exact test
  fisher_result <- fisher.test(cont_table)

  enrichment_results <- rbind(enrichment_results, data.frame(
    Module = as.numeric(mod),
    Module_Color = labels2colors(as.numeric(mod)),
    Module_Size = as.numeric(all_module_counts[mod]),
    Targets_in_Module = in_module_target,
    Expected = round(n_targets * all_module_counts[mod] / n_wgcna_genes, 1),
    Fold_Enrichment = round(in_module_target / (n_targets * all_module_counts[mod] / n_wgcna_genes), 2),
    P_value = fisher_result$p.value,
    Odds_Ratio = fisher_result$estimate
  ))
}

# Adjust p-values
enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
enrichment_results <- enrichment_results[order(enrichment_results$P_value), ]

# Identify significantly enriched modules
sig_enriched <- enrichment_results[enrichment_results$FDR < 0.05 & enrichment_results$Fold_Enrichment > 1, ]

cat("Modules significantly enriched for JAG1 targets (FDR < 0.05):\n")
if (nrow(sig_enriched) > 0) {
  print(sig_enriched[, c("Module", "Module_Color", "Module_Size", "Targets_in_Module",
                         "Fold_Enrichment", "FDR")])
} else {
  cat("  No significantly enriched modules\n")
}

write.csv(enrichment_results, "03_results/tables/WGCNA/JAG1_module_enrichment.csv",
          row.names = FALSE)
cat("\nSaved: JAG1_module_enrichment.csv\n")


cat("\n========================================\n")

# Hub genes: high kME in their module (kME > 0.8)
hub_threshold <- 0.8

hub_genes <- data.frame()

for (mod in unique(net$colors)) {
  if (mod == 0) next  # Skip grey module

  mod_genes <- genes_in_wgcna[net$colors == mod]
  mod_color <- labels2colors(mod)
  kme_col <- paste0("kME_", mod)

  # Get kME for module genes
  mod_kME <- kMEs[mod_genes, kme_col]

  # Identify hubs
  hubs <- mod_kME[mod_kME > hub_threshold]

  if (length(hubs) > 0) {
    hub_df <- data.frame(
      GeneID = names(hubs),
      Module = mod,
      Module_Color = mod_color,
      kME = as.numeric(hubs),
      stringsAsFactors = FALSE
    )
    hub_genes <- rbind(hub_genes, hub_df)
  }
}

hub_genes <- hub_genes[order(-hub_genes$kME), ]

cat("Hub gene summary (kME > ", hub_threshold, "):\n", sep = "")
cat("  Total hub genes:", nrow(hub_genes), "\n")

# Summarize by module
hub_summary <- aggregate(GeneID ~ Module + Module_Color, data = hub_genes, FUN = length)
colnames(hub_summary)[3] <- "N_Hubs"
hub_summary <- hub_summary[order(-hub_summary$N_Hubs), ]

cat("\nHub genes per module (top 10):\n")
print(head(hub_summary, 10))

# Check if any JAG1 targets are hubs
target_hubs <- hub_genes[hub_genes$GeneID %in% jag1_target_genes, ]
cat("\nJAG1 targets that are hub genes:", nrow(target_hubs), "\n")

if (nrow(target_hubs) > 0) {
  cat("\nTop 10 JAG1 target hubs:\n")
  print(head(target_hubs, 10))
}

write.csv(hub_genes, "03_results/tables/WGCNA/hub_genes.csv", row.names = FALSE)
cat("\nSaved: hub_genes.csv\n")

# UPDATED: Use Line_LW_Ratio instead of Narrow_TP0 for consistency with Figure 5

cat("\n========================================\n")

# Use Line_LW_Ratio (actual measured L:W ratios per line) if available
# This is consistent with Figure 5 module-trait correlations
if ("Line_LW_Ratio" %in% colnames(moduleTraitCor)) {
  pheno_trait <- "Line_LW_Ratio"
  cat("Using Line_LW_Ratio (actual measured phenotype) for module selection\n\n")
} else if ("V1_Ratio" %in% colnames(moduleTraitCor)) {
  pheno_trait <- "V1_Ratio"
  cat("Using V1_Ratio for module selection\n\n")
} else {
  pheno_trait <- "Narrow"
  cat("WARNING: Line_LW_Ratio not found, falling back to Narrow (binary)\n\n")
}

pheno_cor <- moduleTraitCor[, pheno_trait]
pheno_pval <- moduleTraitPvalue[, pheno_trait]

narrow_modules <- data.frame(
  Module = gsub("ME", "", names(pheno_cor)),
  Correlation = pheno_cor,
  P_value = pheno_pval,
  stringsAsFactors = FALSE
)
narrow_modules$FDR <- p.adjust(narrow_modules$P_value, method = "BH")
narrow_modules$Significant <- narrow_modules$FDR < 0.05

# Add module info
narrow_modules$Module_Size <- module_info$Size[match(narrow_modules$Module, as.character(module_info$Module))]
narrow_modules$Module_Color <- labels2colors(as.numeric(narrow_modules$Module))

# Add JAG1 target enrichment
narrow_modules$Target_Enrichment <- enrichment_results$Fold_Enrichment[match(
  as.numeric(narrow_modules$Module), enrichment_results$Module)]
narrow_modules$Target_Enrichment_FDR <- enrichment_results$FDR[match(
  as.numeric(narrow_modules$Module), enrichment_results$Module)]

narrow_modules <- narrow_modules[order(narrow_modules$P_value), ]

cat("Modules significantly correlated with", pheno_trait, "(FDR < 0.05):\n")
sig_narrow <- narrow_modules[narrow_modules$Significant, ]
if (nrow(sig_narrow) > 0) {
  print(sig_narrow[, c("Module", "Module_Color", "Correlation", "P_value", "Target_Enrichment")])
} else {
  cat("  No significantly correlated modules\n")
}

# Also keep Narrow_TP0 correlation for reference (legacy)
if ("Narrow_TP0" %in% colnames(moduleTraitCor)) {
  narrow_modules$Narrow_TP0_Cor <- moduleTraitCor[match(paste0("ME", narrow_modules$Module),
                                                         rownames(moduleTraitCor)), "Narrow_TP0"]
}


# Check if phenotype traits are available
if (exists("phenotype_available") && phenotype_available && "V1_Ratio" %in% colnames(traitData)) {
  cat("\n========================================\n")

  v1_ratio_cor <- moduleTraitCor[, "V1_Ratio"]
  v1_ratio_pval <- moduleTraitPvalue[, "V1_Ratio"]

  pheno_modules <- data.frame(
    Module = gsub("ME", "", names(v1_ratio_cor)),
    V1_Ratio_Cor = v1_ratio_cor,
    V1_Ratio_Pval = v1_ratio_pval,
    stringsAsFactors = FALSE
  )

  # Add Narrow_TP0 correlation for comparison
  pheno_modules$Narrow_TP0_Cor <- narrow_modules$Correlation[match(pheno_modules$Module, narrow_modules$Module)]

  pheno_modules$V1_FDR <- p.adjust(pheno_modules$V1_Ratio_Pval, method = "BH")
  pheno_modules$Significant <- pheno_modules$V1_FDR < 0.05

  # Add target enrichment
  pheno_modules$Target_Enrichment <- enrichment_results$Fold_Enrichment[match(
    as.numeric(pheno_modules$Module), enrichment_results$Module)]

  pheno_modules <- pheno_modules[order(pheno_modules$V1_Ratio_Pval), ]

  cat("Modules significantly correlated with V1 Leaf Ratio (p < 0.05):\n")
  sig_pheno <- pheno_modules[pheno_modules$Significant, ]
  if (nrow(sig_pheno) > 0) {
    print(sig_pheno[, c("Module", "V1_Ratio_Cor", "V1_Ratio_Pval", "Narrow_TP0_Cor", "Target_Enrichment")])
  } else {
    cat("  No significantly correlated modules\n")
  }

  # Identify modules enriched for JAG1 targets AND correlated with V1_Ratio
  jag1_pheno_modules <- pheno_modules[pheno_modules$Significant &
                                        !is.na(pheno_modules$Target_Enrichment) &
                                        pheno_modules$Target_Enrichment > 1.5, ]

  if (nrow(jag1_pheno_modules) > 0) {
    cat("\nModules enriched for JAG1 targets AND correlated with V1_Ratio:\n")
    print(jag1_pheno_modules[, c("Module", "V1_Ratio_Cor", "Target_Enrichment")])
  }
}


cat("\n========================================\n")

# Plot 1: JAG1 target enrichment across modules
png("03_results/figures/21_WGCNA_JAG1/target_enrichment.png",
    width = 1000, height = 600, res = 120)

# Prepare data - exclude grey module
plot_data <- enrichment_results[enrichment_results$Module != 0, ]
plot_data <- plot_data[order(-plot_data$Fold_Enrichment), ]
plot_data <- head(plot_data, 20)  # Top 20 modules

par(mar = c(7, 4, 4, 2))
barplot(
  plot_data$Fold_Enrichment,
  names.arg = plot_data$Module_Color,
  col = ifelse(plot_data$FDR < 0.05, "red", "gray"),
  las = 2,
  main = "JAG1 Target Enrichment by Module",
  ylab = "Fold Enrichment",
  border = "white"
)
abline(h = 1, lty = 2, col = "blue")
legend("topright", legend = c("FDR < 0.05", "Not significant"),
       fill = c("red", "gray"))

dev.off()
cat("Saved: target_enrichment.png\n")

# Plot 2: Module-trait correlation vs target enrichment
png("03_results/figures/21_WGCNA_JAG1/correlation_vs_enrichment.png",
    width = 800, height = 800, res = 120)

plot_data <- merge(narrow_modules, enrichment_results,
                   by.x = "Module", by.y = "Module", all.x = TRUE)

plot(plot_data$Correlation, log2(plot_data$Fold_Enrichment + 0.01),
     pch = 16,
     col = ifelse(plot_data$Significant, "red",
                  ifelse(plot_data$Target_Enrichment_FDR < 0.05, "blue", "gray")),
     cex = log10(plot_data$Module_Size + 1),
     xlab = "Correlation with Narrow_TP0",
     ylab = "log2(JAG1 Target Fold Enrichment)",
     main = "Module-Trait Correlation vs JAG1 Target Enrichment")
abline(h = 0, v = 0, lty = 2, col = "gray")

# Label significant points
sig_points <- plot_data[plot_data$Significant | plot_data$Target_Enrichment_FDR < 0.05, ]
if (nrow(sig_points) > 0) {
  text(sig_points$Correlation, log2(sig_points$Fold_Enrichment + 0.01),
       labels = sig_points$Module_Color.x, pos = 3, cex = 0.7)
}

legend("topleft",
       legend = c("Sig. Narrow_TP0 corr", "Sig. target enrichment", "Not significant"),
       col = c("red", "blue", "gray"), pch = 16)

dev.off()
cat("Saved: correlation_vs_enrichment.png\n")

# Plot 3: Hub gene kME distribution in top modules
png("03_results/figures/21_WGCNA_JAG1/hub_gene_distribution.png",
    width = 1000, height = 600, res = 120)

top_modules_hubs <- head(hub_summary$Module, 9)

par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

for (mod in top_modules_hubs) {
  mod_color <- labels2colors(mod)
  mod_genes <- genes_in_wgcna[net$colors == mod]
  kme_col <- paste0("kME_", mod)

  kme_vals <- kMEs[mod_genes, kme_col]

  hist(kme_vals, breaks = 20, col = mod_color, border = "white",
       main = paste("Module", mod, "(", mod_color, ")"),
       xlab = "kME", xlim = c(-0.5, 1))
  abline(v = hub_threshold, col = "red", lty = 2, lwd = 2)

  # Mark JAG1 targets
  targets_in_mod <- intersect(mod_genes, jag1_target_genes)
  if (length(targets_in_mod) > 0) {
    target_kme <- kMEs[targets_in_mod, kme_col]
    rug(target_kme, col = "blue", lwd = 2)
  }
}

dev.off()
cat("Saved: hub_gene_distribution.png\n")

# Plot 4: Top hub genes expression heatmap
if (nrow(hub_genes) > 0) {
  png("03_results/figures/21_WGCNA_JAG1/top_hubs_heatmap.png",
      width = 1000, height = 800, res = 120)

  # Get top 50 hub genes
  top_hubs <- head(hub_genes$GeneID, min(50, nrow(hub_genes)))
  hub_expr <- datExpr[, top_hubs, drop = FALSE]

  # Sample annotation
  sample_anno <- data.frame(
    Leaf_type = targets$Leaf_type,
    Timepoint = targets$Timepoint,
    row.names = rownames(datExpr)
  )

  # Gene annotation (module color)
  gene_anno <- data.frame(
    Module = hub_genes$Module_Color[1:length(top_hubs)],
    row.names = top_hubs
  )

  # Create color mapping for modules
  unique_colors <- unique(gene_anno$Module)
  module_color_map <- setNames(unique_colors, unique_colors)

  anno_colors <- list(
    Leaf_type = c("Broad" = "#E69F00", "Narrow" = "#56B4E9"),
    Timepoint = c("TP0" = "#1b9e77", "TP1" = "#d95f02", "TP2" = "#7570b3",
                  "TP3" = "#e7298a", "TP4" = "#66a61e"),
    Module = module_color_map
  )

  pheatmap(
    t(hub_expr),
    annotation_col = sample_anno,
    annotation_row = gene_anno,
    annotation_colors = anno_colors,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    main = paste("Top", length(top_hubs), "Hub Genes Expression"),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "row"
  )

  dev.off()
  cat("Saved: top_hubs_heatmap.png\n")
}


cat("\n========================================\n")

if (!is.na(JAG1_module) && JAG1_module != 0) {
  # Get all genes in JAG1's module
  jag1_mod_genes <- genes_in_wgcna[net$colors == JAG1_module]

  cat("JAG1's module (", JAG1_module_color, ") analysis:\n")
  cat("  Total genes:", length(jag1_mod_genes), "\n")

  # Module eigengene pattern
  me_col <- paste0("ME", JAG1_module)
  if (me_col %in% colnames(MEs)) {
    cat("\n  Module eigengene correlation with traits:\n")
    cat(sprintf("    Narrow: r = %.3f, p = %.4f\n",
                moduleTraitCor[me_col, "Narrow"],
                moduleTraitPvalue[me_col, "Narrow"]))
    cat(sprintf("    TP0: r = %.3f, p = %.4f\n",
                moduleTraitCor[me_col, "TP0"],
                moduleTraitPvalue[me_col, "TP0"]))
    cat(sprintf("    Narrow_TP0: r = %.3f, p = %.4f\n",
                moduleTraitCor[me_col, "Narrow_TP0"],
                moduleTraitPvalue[me_col, "Narrow_TP0"]))
  }

  # JAG1 targets in module
  targets_in_jag1_mod <- intersect(jag1_mod_genes, jag1_target_genes)
  cat("\n  JAG1 targets in this module:", length(targets_in_jag1_mod), "\n")

  # Hub genes in module
  hubs_in_jag1_mod <- hub_genes[hub_genes$Module == JAG1_module, ]
  cat("  Hub genes in this module:", nrow(hubs_in_jag1_mod), "\n")

  if (nrow(hubs_in_jag1_mod) > 0) {
    cat("\n  Top hubs in JAG1's module:\n")
    print(head(hubs_in_jag1_mod[, c("GeneID", "kME")], 10))
  }

  # Plot JAG1 module expression
  png("03_results/figures/21_WGCNA_JAG1/JAG1_module_expression.png",
      width = 1000, height = 600, res = 120)

  par(mfrow = c(1, 2))

  # Module eigengene
  me_data <- data.frame(
    ME = MEs[, me_col],
    Timepoint = factor(targets$Timepoint, levels = paste0("TP", 0:4)),
    Leaf_type = targets$Leaf_type
  )

  means <- aggregate(ME ~ Timepoint + Leaf_type, data = me_data, FUN = mean)
  broad_means <- means[means$Leaf_type == "Broad", ]
  narrow_means <- means[means$Leaf_type == "Narrow", ]

  ylim_range <- range(c(broad_means$ME, narrow_means$ME))
  ylim_range <- ylim_range + c(-0.1, 0.1) * diff(ylim_range)

  plot(1:5, broad_means$ME, type = "b", pch = 16, col = "#E69F00",
       ylim = ylim_range, xlab = "Timepoint", ylab = "Module Eigengene",
       main = paste0("JAG1's Module (", JAG1_module_color, ") Eigengene"),
       xaxt = "n", lwd = 2)
  axis(1, at = 1:5, labels = paste0("TP", 0:4))
  lines(1:5, narrow_means$ME, type = "b", pch = 17, col = "#56B4E9", lwd = 2)
  abline(h = 0, lty = 2, col = "gray")
  legend("topright", legend = c("Broad", "Narrow"),
         col = c("#E69F00", "#56B4E9"), pch = c(16, 17), lwd = 2)

  # JAG1 expression for comparison
  if (JAG1_ID %in% rownames(logCPM)) {
    jag1_data <- data.frame(
      Expression = logCPM[JAG1_ID, ],
      Timepoint = factor(targets$Timepoint, levels = paste0("TP", 0:4)),
      Leaf_type = targets$Leaf_type
    )

    jag1_means <- aggregate(Expression ~ Timepoint + Leaf_type, data = jag1_data, FUN = mean)
    broad_jag1 <- jag1_means[jag1_means$Leaf_type == "Broad", ]
    narrow_jag1 <- jag1_means[jag1_means$Leaf_type == "Narrow", ]

    ylim_range <- range(c(broad_jag1$Expression, narrow_jag1$Expression))
    ylim_range <- ylim_range + c(-0.5, 0.5)

    plot(1:5, broad_jag1$Expression, type = "b", pch = 16, col = "#E69F00",
         ylim = ylim_range, xlab = "Timepoint", ylab = "log2 CPM",
         main = "JAG1 (Glyma.20G116200) Expression",
         xaxt = "n", lwd = 2)
    axis(1, at = 1:5, labels = paste0("TP", 0:4))
    lines(1:5, narrow_jag1$Expression, type = "b", pch = 17, col = "#56B4E9", lwd = 2)
    legend("topright", legend = c("Broad", "Narrow"),
           col = c("#E69F00", "#56B4E9"), pch = c(16, 17), lwd = 2)
  }

  dev.off()
  cat("\nSaved: JAG1_module_expression.png\n")

} else {
  cat("JAG1 is not assigned to a non-grey module\n")
}


cat("\n========================================\n")

# Create master table of module-JAG1 relationships
module_jag1_summary <- data.frame(
  Module = as.character(module_info$Module),
  Module_Color = module_info$Color,
  Module_Size = module_info$Size,
  stringsAsFactors = FALSE
)

# Add trait correlations
module_jag1_summary$Cor_Narrow <- moduleTraitCor[paste0("ME", module_jag1_summary$Module), "Narrow"]
module_jag1_summary$Cor_TP0 <- moduleTraitCor[paste0("ME", module_jag1_summary$Module), "TP0"]
module_jag1_summary$Cor_Narrow_TP0 <- moduleTraitCor[paste0("ME", module_jag1_summary$Module), "Narrow_TP0"]

# Add target enrichment
module_jag1_summary$JAG1_Targets <- enrichment_results$Targets_in_Module[match(
  as.numeric(module_jag1_summary$Module), enrichment_results$Module)]
module_jag1_summary$Target_Enrichment <- enrichment_results$Fold_Enrichment[match(
  as.numeric(module_jag1_summary$Module), enrichment_results$Module)]
module_jag1_summary$Enrichment_FDR <- enrichment_results$FDR[match(
  as.numeric(module_jag1_summary$Module), enrichment_results$Module)]

# Add hub gene count
hub_counts <- aggregate(GeneID ~ Module, data = hub_genes, FUN = length)
colnames(hub_counts)[2] <- "N_Hubs"
module_jag1_summary$N_Hubs <- hub_counts$N_Hubs[match(
  as.numeric(module_jag1_summary$Module), hub_counts$Module)]
module_jag1_summary$N_Hubs[is.na(module_jag1_summary$N_Hubs)] <- 0

# Mark JAG1's module
module_jag1_summary$Contains_JAG1 <- as.numeric(module_jag1_summary$Module) == JAG1_module

write.csv(module_jag1_summary, "03_results/tables/WGCNA/module_JAG1_summary.csv",
          row.names = FALSE)
cat("Saved: module_JAG1_summary.csv\n")


cat("\n========================================\n")

# Prepare objects to save
objects_to_save <- c(
  "kMEs",
  "hub_genes",
  "enrichment_results",
  "narrow_modules",
  "module_jag1_summary",
  "JAG1_module",
  "JAG1_module_color",
  "target_modules",
  "jag1_targets",
  "jag1_target_genes",
  # Keep previous objects
  "traitData",
  "moduleTraitCor",
  "moduleTraitPvalue",
  "MEs",
  "module_info",
  "datExpr",
  "net",
  "wgcna_prep",
  "module_colors",
  "targets",
  "logCPM"
)

# Add phenotype objects if available
if (exists("phenotype_available") && phenotype_available) {
  objects_to_save <- c(objects_to_save, "phenotype_available")
  if (exists("pheno_modules")) {
    objects_to_save <- c(objects_to_save, "pheno_modules")
  }
  if (exists("phenotype_summary")) {
    objects_to_save <- c(objects_to_save, "phenotype_summary")
  }
}

save(
  list = objects_to_save,
  file = "03_results/checkpoints/21_WGCNA_JAG1.RData"
)

cat("Checkpoint saved: 21_WGCNA_JAG1.RData\n")


cat("\n================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("KEY RESULTS:\n")
cat("  - JAG1 module:", ifelse(is.na(JAG1_module), "Not in WGCNA", paste(JAG1_module, "(", JAG1_module_color, ")")), "\n")
cat("  - Total hub genes identified:", nrow(hub_genes), "\n")
cat("  - Modules enriched for JAG1 targets (FDR < 0.05):", sum(enrichment_results$FDR < 0.05 & enrichment_results$Fold_Enrichment > 1), "\n")
cat("  - Modules correlated with L:W phenotype (FDR < 0.05):", sum(narrow_modules$Significant), "\n")
cat("  - JAG1 targets that are hub genes:", nrow(target_hubs), "\n")

# Add phenotype correlation results if available
if (exists("pheno_modules") && exists("sig_pheno")) {
  cat("  - Modules correlated with V1_Ratio (p < 0.05):", nrow(sig_pheno), "\n")
}
cat("\n")

cat("OUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/21_WGCNA_JAG1.RData\n")
cat("  - Tables: 03_results/tables/WGCNA/\n")
cat("    - hub_genes.csv\n")
cat("    - JAG1_module_enrichment.csv\n")
cat("    - module_JAG1_summary.csv\n")
cat("  - Figures: 03_results/figures/21_WGCNA_JAG1/\n")
cat("    - target_enrichment.png\n")
cat("    - correlation_vs_enrichment.png\n")
cat("    - hub_gene_distribution.png\n")
cat("    - top_hubs_heatmap.png\n")
cat("    - JAG1_module_expression.png\n\n")

