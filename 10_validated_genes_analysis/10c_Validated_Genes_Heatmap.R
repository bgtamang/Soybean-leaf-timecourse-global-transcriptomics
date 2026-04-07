# =============================================================================
# Figure 7: Validated GmJAG1 Targets — Heatmap + Bar/Dot
# Self-contained: reads from checkpoints + TableS30 + phenotype correlations
# =============================================================================

rm(list = ls())
gc()

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase3-Refined-Analysis"
setwd(base_dir)

fig7_dir <- "Manuscript/Figure/Figure7"
dir.create(fig7_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD VALIDATED GENES
# =============================================================================

targets <- read_csv("03_results/tables/supplementary_v6/TableS30_HighConfidence_79_Targets.csv",
                    show_col_types = FALSE)
cat("Loaded", nrow(targets), "validated targets\n")

# r(LW) already in TableS30
targets <- targets %>% rename(r_LW = Cor_LW_Ratio)

# =============================================================================
# 2. PER-TIMEPOINT FOLD CHANGES (Narrow vs Broad at each timepoint)
# =============================================================================

# Use edgeR GLM fit coefficients to compute NarrowVsBroad at each TP
load("03_results/checkpoints/11_DE_results.RData")

gene_ids <- targets$GeneID

# fit$coefficients has columns for each group: GenotypeLine_TPx
# Narrow lines: PI612713B, PI547745
# Broad lines: PI532462A, LD112170
coef_mat <- fit$coefficients  # genes x groups
coef_cols <- colnames(coef_mat)

cat("Design groups:", paste(coef_cols, collapse = ", "), "\n")

# Compute Narrow vs Broad logFC at each timepoint
tp_labels <- paste0("TP", 0:4)
narrow_lines <- c("PI612713B", "PI547745")
broad_lines <- c("PI532462A", "LD112170")

fc_matrix <- matrix(NA, nrow = length(gene_ids), ncol = 5,
                    dimnames = list(gene_ids, paste0("FC_", tp_labels)))

all_genes <- rownames(coef_mat)

for (tp_idx in seq_along(tp_labels)) {
  tp <- tp_labels[tp_idx]

  # Find narrow columns for this TP
  narrow_cols <- paste0(narrow_lines, "_", tp)
  narrow_cols <- narrow_cols[narrow_cols %in% coef_cols]

  # Find broad columns for this TP
  broad_cols <- paste0(broad_lines, "_", tp)
  broad_cols <- broad_cols[broad_cols %in% coef_cols]

  if (length(narrow_cols) > 0 && length(broad_cols) > 0) {
    narrow_mean <- rowMeans(coef_mat[, narrow_cols, drop = FALSE])
    broad_mean <- rowMeans(coef_mat[, broad_cols, drop = FALSE])
    logfc <- narrow_mean - broad_mean

    matched <- match(gene_ids, all_genes)
    fc_matrix[, tp_idx] <- logfc[matched]
  }
}

tp_fc <- as.data.frame(fc_matrix)
tp_fc$GeneID <- gene_ids

cat("\nNarrow vs Broad logFC at each timepoint:\n")
for (col in colnames(fc_matrix)) {
  vals <- fc_matrix[, col]
  cat(sprintf("  %s: %d non-NA, range [%.2f, %.2f]\n",
              col, sum(!is.na(vals)), min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

# Clean up
rm(fit, dge, de_results, de_results_treat, contrast_matrix, design,
   edgeR_coded, edgeR_tr_coded, de_summary, contrast_desc)
gc()

# =============================================================================
# 3. CURATED FUNCTIONAL CATEGORIES (from Phase2 + new genes)
# =============================================================================

category_map <- tribble(
  ~GeneID,              ~Functional_Category,       ~Short_Function,
  # --- Hormone signaling ---
  "Glyma.02G061400",    "Hormone (Auxin)",          "NPH3",
  "Glyma.11G168428",    "Hormone (Auxin)",          "CYP83B1",
  "Glyma.04G013800",    "Hormone (JA)",             "JAZ1",
  "Glyma.11G225500",    "Hormone (Cytokinin)",      "cZOGT",
  "Glyma.15G263332",    "Hormone (GA)",             "KS",
  "Glyma.15G196500",    "Hormone (Light/GA)",       "PhyB",
  "Glyma.04G180400",    "Hormone (ABA)",            "RD22",
  "Glyma.05G204800",    "Hormone (ABA)",            "OSM34",
  "Glyma.19G133900",    "Hormone (BR)",             "BPG2",
  # --- Defense/NBS-LRR ---
  "Glyma.03G052800",    "Defense (NBS-LRR)",        "NBS-LRR",
  "Glyma.06G310000",    "Defense (NBS-LRR)",        "NBS-LRR",
  "Glyma.16G213800",    "Defense (NBS-LRR)",        "NBS-LRR",
  "Glyma.16G215000",    "Defense (NBS-LRR)",        "NBS-LRR",
  "Glyma.14G204600",    "Defense (NB-ARC)",         "NB-ARC",
  "Glyma.14G205000",    "Defense (NB-ARC)",         "NB-ARC",
  "Glyma.14G205300",    "Defense (NB-ARC)",         "NB-ARC",
  # --- Receptor kinase / Signaling ---
  "Glyma.19G167300",    "Receptor kinase",          "MIK2",
  "Glyma.18G273100",    "Receptor kinase",          "ANXUR1",
  "Glyma.12G198000",    "Receptor kinase",          "RLK",
  "Glyma.06G184400",    "Receptor kinase",          "LRR-RLK",
  "Glyma.16G202200",    "Receptor kinase",          "Ser/Thr kinase",
  "Glyma.16G018000",    "Circadian clock",          "PRR5",
  # --- Transcription / Chromatin ---
  "Glyma.02G226700",    "Transcription/Chromatin",  "ELP5",
  "Glyma.11G251900",    "Transcription/Chromatin",  "CAMTA4",
  "Glyma.12G022600",    "Transcription/Chromatin",  "SCL23",
  "Glyma.09G276900",    "Transcription/Chromatin",  "TBP2",
  "Glyma.20G062300",    "Transcription/Chromatin",  "ADA2a",
  "Glyma.02G189700",    "Epigenetics",              "AGO4",
  # --- Transport ---
  "Glyma.03G171600",    "Transport",                "KEA1",
  "Glyma.05G102600",    "Transport",                "CAX2",
  "Glyma.11G097000",    "Transport",                "LHT8",
  "Glyma.17G192000",    "Transport",                "AAP6",
  # --- Protein folding / ER ---
  "Glyma.02G093200",    "Protein folding/ER",       "HSP70",
  "Glyma.13G357700",    "Protein folding/ER",       "PDI",
  "Glyma.14G178800",    "Protein folding/ER",       "DnaJ1",
  "Glyma.16G098700",    "Protein folding/ER",       "ERDJ3A",
  # --- Vesicle trafficking ---
  "Glyma.14G199600",    "Vesicle/Secretory",        "EXO70B1",
  "Glyma.01G115100",    "Vesicle/Secretory",        "CLC",
  # --- RNA processing / Translation ---
  "Glyma.13G136500",    "RNA processing",           "RBP",
  "Glyma.17G178001",    "Translation",              "RPL25",
  "Glyma.09G000400",    "RNA processing",           "PPR",
  "Glyma.10G043900",    "RNA processing",           "PPR",
  "Glyma.03G000300",    "RNA processing",           "PPR",
  # --- Metabolism ---
  "Glyma.03G073702",    "Metabolism",               "GSTU8",
  "Glyma.03G119400",    "Metabolism",               "ECH1",
  "Glyma.14G058602",    "Metabolism",               "2OG-FeII",
  "Glyma.06G137000",    "Metabolism",               "Fe2OG",
  "Glyma.19G029700",    "Metabolism",               "UGT",
  "Glyma.03G250300",    "Metabolism",               "FUT",
  "Glyma.14G118550",    "Metabolism",               "OFUT",
  "Glyma.03G031200",    "Metabolism",               "CYP",
  # --- Organellar ---
  "Glyma.15G250600",    "Chloroplast",              "THF1",
  "Glyma.05G181700",    "Mitochondria",             "Mitofilin",
  "Glyma.06G059200",    "Mitochondria",             "MIOREX",
  "Glyma.06G059402",    "Mitochondria",             "MIOREX",
  "Glyma.13G099400",    "Mitochondria",             "Metaxin",
  # --- Cell cycle ---
  "Glyma.06G062500",    "Cell cycle",               "CDT1",
  # --- Stress/Signaling/Other ---
  "Glyma.03G041000",    "Stress response",          "LURP1",
  "Glyma.05G127400",    "Cell death",               "PDCD5",
  "Glyma.06G267300",    "Signaling",                "ADPRC",
  "Glyma.06G304100",    "Unknown (DUF247)",         "DUF247",
  "Glyma.06G304600",    "Unknown (DUF247)",         "DUF247",
  "Glyma.12G230967",    "Signaling",                "NHL",
  "Glyma.15G051900",    "Metal homeostasis",        "HMA",
  "Glyma.20G220200",    "Protein turnover",         "F-box",
  "Glyma.07G022602",    "Metabolism",               "XIP",
  "Glyma.03G173300",    "Metabolism",               "Misc",
  "Glyma.03G255900",    "Metabolism",               "Misc",
  "Glyma.18G195200",    "Metabolism",               "Misc",
  "Glyma.06G044000",    "Metabolism",               "Misc",
  "Glyma.17G079900",    "Metabolism",               "Misc",
  "Glyma.19G020700",    "Metabolism",               "Misc",
  # --- No annotation ---
  "Glyma.02G181350",    "No annotation",            "unknown",
  "Glyma.15G240132",    "No annotation",            "unknown",
  "Glyma.18G245400",    "No annotation",            "unknown",
  "Glyma.13G000600",    "No annotation",            "unknown",
  # --- New genes in Phase3 (not in Phase2 72) ---
  "Glyma.18G273600",    "Receptor kinase",          "RLK",
  "Glyma.13G053700",    "Signaling",                "Kinase",
  "Glyma.06G261400",    "Signaling",                "ADPRC"
)

# Merge categories
targets <- targets %>%
  left_join(category_map, by = "GeneID")

# For any genes still missing a category, use Arabi_Defline or mark unknown
targets <- targets %>%
  mutate(
    Functional_Category = ifelse(is.na(Functional_Category), "No annotation", Functional_Category),
    Short_Function = ifelse(is.na(Short_Function), "", Short_Function)
  )

# Create gene labels
targets <- targets %>%
  mutate(Gene_Label = ifelse(Short_Function != "",
                              paste0(GeneID, " (", Short_Function, ")"),
                              GeneID))

# Handle duplicate labels
dup_check <- targets$Gene_Label[duplicated(targets$Gene_Label)]
if (length(dup_check) > 0) {
  for (lab in unique(dup_check)) {
    idx <- which(targets$Gene_Label == lab)
    for (i in seq_along(idx)) {
      targets$Gene_Label[idx[i]] <- paste0(targets$GeneID[idx[i]], " (", targets$Short_Function[idx[i]], "-", i, ")")
    }
  }
}

# Broad categories for row splitting
targets <- targets %>%
  mutate(Broad_Category = case_when(
    grepl("Hormone", Functional_Category) ~ "Hormone signaling",
    grepl("Defense|NB", Functional_Category) ~ "Defense/NBS-LRR",
    grepl("Receptor|Circadian|Signaling", Functional_Category) & !grepl("Hormone", Functional_Category) ~ "Receptor kinase/Signaling",
    grepl("Transcription|Epigenetics|Chromatin", Functional_Category) ~ "Transcription/Chromatin",
    grepl("Transport", Functional_Category) ~ "Transport",
    grepl("Protein folding|Protein turnover", Functional_Category) ~ "Protein folding/ER",
    grepl("Vesicle|Secretory", Functional_Category) ~ "Vesicle trafficking",
    grepl("RNA|Translation", Functional_Category) ~ "RNA processing",
    grepl("Metabolism", Functional_Category) ~ "Metabolism",
    grepl("Chloroplast|Mitochondria", Functional_Category) ~ "Organellar",
    grepl("Cell cycle", Functional_Category) ~ "Cell cycle",
    grepl("Stress|Cell death|Metal", Functional_Category) ~ "Stress/Defense",
    grepl("No annotation|Unknown|DUF", Functional_Category) ~ "Unannotated",
    TRUE ~ "Other"
  ))

cat("\n=== Category Summary ===\n")
cat_summary <- targets %>% count(Broad_Category) %>% arrange(desc(n))
for (i in seq_len(nrow(cat_summary))) {
  cat(sprintf("  %-28s %d\n", cat_summary$Broad_Category[i], cat_summary$n[i]))
}

# Warn about uncategorized genes
uncategorized <- targets %>% filter(Functional_Category == "No annotation" & Short_Function == "")
if (nrow(uncategorized) > 0) {
  cat("\nWARNING:", nrow(uncategorized), "genes have no category:\n")
  cat(uncategorized$GeneID, sep = "\n")
}

cat_colors <- c(
  "Hormone signaling"        = "#E41A1C",
  "Defense/NBS-LRR"          = "#377EB8",
  "Receptor kinase/Signaling"= "#4DAF4A",
  "Transcription/Chromatin"  = "#984EA3",
  "Transport"                = "#FF7F00",
  "Protein folding/ER"       = "#A65628",
  "Vesicle trafficking"      = "#F781BF",
  "RNA processing"           = "#999999",
  "Metabolism"               = "#66C2A5",
  "Organellar"               = "#8DA0CB",
  "Cell cycle"               = "#FC8D62",
  "Stress/Defense"           = "#E78AC3",
  "Unannotated"              = "#BDBDBD",
  "Other"                    = "#E5C494"
)

# =============================================================================
# 4. MERGE TIMEPOINT FCs
# =============================================================================

dat <- targets %>%
  left_join(tp_fc, by = "GeneID")

cat("\nFinal data:", nrow(dat), "genes with timepoint FCs\n")

# =============================================================================
# 5. FIGURE A: Bar + Dot Plot
# =============================================================================

cat("\n=== Building Figure A: Bar + Dot Plot ===\n")

plot_dat <- dat %>%
  arrange(FC_TP1) %>%
  mutate(Gene_Label = factor(Gene_Label, levels = Gene_Label))

p_bar <- ggplot(plot_dat, aes(x = FC_TP1, y = Gene_Label, fill = Broad_Category)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = cat_colors, name = "Functional\nCategory") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = expression(log[2]*"FC at TP1 (mean narrow genotypes vs TP0)"), y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 5.5, color = "black"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 2, 5, 5)
  )

p_dot <- ggplot(plot_dat, aes(x = 1, y = Gene_Label)) +
  geom_point(aes(size = abs(r_LW),
                 color = ifelse(r_LW >= 0, "Positive", "Negative")),
             alpha = 0.7) +
  scale_size_continuous(range = c(0.5, 4), name = "|r(LW)|",
                        breaks = c(0.5, 0.7, 0.9)) +
  scale_color_manual(values = c("Positive" = "#D73027", "Negative" = "#4575B4"),
                     name = "Direction") +
  labs(x = "r(LW)", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 0)
  )

p_combined <- p_bar + p_dot +
  plot_layout(widths = c(4, 1)) +
  plot_annotation(
    title = paste0(nrow(dat), " Validated GmJAG1 Targets: Expression Change at TP1 and Leaf Shape Correlation"),
    subtitle = expression("Bars: "*log[2]*"FC at TP1 (mean narrow genotypes vs TP0) | Dots: |r| with leaf length:width ratio"),
    theme = theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  )

ggsave(file.path(fig7_dir, "Fig7A_BarDot_TP1.png"),
       p_combined, width = 12, height = 18, dpi = 300)
cat("Saved: Fig7A_BarDot_TP1.png\n")

# =============================================================================
# 6. FIGURE B: ComplexHeatmap
# =============================================================================

cat("\n=== Building Figure B: ComplexHeatmap ===\n")

dat_sorted <- dat %>%
  arrange(Broad_Category, desc(FC_TP1))

mat <- dat_sorted %>%
  select(FC_TP0, FC_TP1, FC_TP2, FC_TP3, FC_TP4) %>%
  as.matrix()
rownames(mat) <- dat_sorted$Gene_Label
colnames(mat) <- c("TP0", "TP1", "TP2", "TP3", "TP4")

# Non-linear color ramp: more resolution in 0-5 range, full scale to 15
col_fun <- colorRamp2(c(-5, -2, 0, 2, 5, 10, 15),
                       c("#2166AC", "#92C5DE", "white", "#FDDBC7", "#D6604D", "#B2182B", "#67001F"))

de_tier <- dat_sorted$DE_Tier
tier_colors <- c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")

bind_class <- dat_sorted$Binding_Class
bind_colors <- c("ChIP_Only" = "#1B9E77", "ChIP_AND_DAP" = "#D95F02", "DAP_Only" = "#7570B3")

r_vals <- dat_sorted$r_LW

n_peaks <- as.numeric(dat_sorted$ChIP_N_Peaks)
n_peaks[is.na(n_peaks)] <- 0

ha_left <- rowAnnotation(
  DE_Tier = de_tier,
  Binding = bind_class,
  col = list(
    DE_Tier = tier_colors,
    Binding = bind_colors
  ),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 7),
  annotation_width = unit(c(0.3, 0.3), "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    DE_Tier = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)),
    Binding = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))
  )
)

ha_right <- rowAnnotation(
  `r(LW)` = anno_barplot(r_vals,
                          gp = gpar(fill = ifelse(r_vals >= 0, "#D73027", "#4575B4")),
                          width = unit(4, "cm"),
                          border = FALSE,
                          axis_param = list(gp = gpar(fontsize = 8))),
  `ChIP\nPeaks` = anno_barplot(n_peaks,
                                gp = gpar(fill = "#636363"),
                                width = unit(3, "cm"),
                                border = FALSE,
                                axis_param = list(gp = gpar(fontsize = 6))),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8)
)

ht <- Heatmap(mat,
              name = "log2FC",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 0,
              column_names_centered = TRUE,
              width = unit(10, "cm"),
              left_annotation = ha_left,
              right_annotation = ha_right,
              row_split = factor(dat_sorted$Broad_Category,
                                levels = names(sort(table(dat_sorted$Broad_Category),
                                                     decreasing = TRUE))),
              row_title_gp = gpar(fontsize = 7),
              row_title_rot = 0,
              row_gap = unit(1, "mm"),
              border = TRUE,
              heatmap_legend_param = list(
                title = "log2FC",
                title_gp = gpar(fontsize = 8),
                labels_gp = gpar(fontsize = 8),
                legend_height = unit(3, "cm")
              ),
              column_title = paste0(nrow(dat_sorted), " Validated GmJAG1 Targets: Temporal Expression and Leaf Shape Correlation"),
              column_title_gp = gpar(fontsize = 11, fontface = "bold")
)

png(file.path(fig7_dir, "Fig7B_ComplexHeatmap.png"),
    width = 14, height = 20, units = "in", res = 300)
draw(ht, padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

cat("Saved: Fig7B_ComplexHeatmap.png\n")
cat("\n=== FIGURE 7 COMPLETE ===\n")
