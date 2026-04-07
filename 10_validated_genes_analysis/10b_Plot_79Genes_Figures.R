# =============================================================================
# Two Figures for Validated GmJAG1 Targets
# Figure A: Bar + Dot plot (logFC at TP1 + r(LW) circle)
# Figure B: ComplexHeatmap (logFC across timepoints + annotations)
# =============================================================================

rm(list = ls())
gc()

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))

out_dir <- "Manuscript/Figure/Figure6/Validated_Genes_Analysis"

# =============================================================================
# LOAD DATA
# =============================================================================

annot <- read_csv(file.path(out_dir, "TableS_Validated_Genes_AllTargets.csv"),
                  show_col_types = FALSE)

tp_fc <- read_csv(file.path(out_dir, "validated_genes_TP1_FC.csv"),
                  show_col_types = FALSE)

# Merge
dat <- annot %>%
  left_join(tp_fc, by = "GeneID")

cat("Merged data:", nrow(dat), "genes\n")

# =============================================================================
# DEFINE BROAD CATEGORIES AND COLORS
# =============================================================================

dat <- dat %>%
  mutate(Broad_Category = case_when(
    grepl("Hormone", Functional_Category) ~ "Hormone signaling",
    grepl("Defense|NB", Functional_Category) ~ "Defense/NBS-LRR",
    grepl("Receptor|Circadian", Functional_Category) ~ "Receptor kinase/Signaling",
    grepl("Transcription|Epigenetics", Functional_Category) ~ "Transcription/Chromatin",
    grepl("Transport", Functional_Category) ~ "Transport",
    grepl("Protein folding", Functional_Category) ~ "Protein folding/ER",
    grepl("Vesicle|Secretory", Functional_Category) ~ "Vesicle trafficking",
    grepl("RNA|Translation", Functional_Category) ~ "RNA processing",
    grepl("Metabolism", Functional_Category) ~ "Metabolism",
    grepl("Chloroplast|Mitochondria", Functional_Category) ~ "Organellar",
    grepl("No annotation|Unknown", Functional_Category) ~ "Unannotated",
    TRUE ~ "Other"
  ))

# Category colors
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
  "Unannotated"              = "#BDBDBD",
  "Other"                    = "#E5C494"
)

# Create gene labels: ID + short function name
dat <- dat %>%
  mutate(
    Func_Short = case_when(
      grepl("RD22", Short_Function) ~ "RD22",
      grepl("NPH3", Short_Function) ~ "NPH3",
      grepl("OSM34|Osmotin", Short_Function) ~ "OSM34",
      grepl("AXR1", Short_Function) ~ "AXR1",
      grepl("JAZ1|TIFY", Short_Function) ~ "JAZ1",
      grepl("ELP5|Elongator", Short_Function) ~ "ELP5",
      grepl("AGO4|Argonaute", Short_Function) ~ "AGO4",
      grepl("MIK2", Short_Function) ~ "MIK2",
      grepl("SCL23|SCARECROW", Short_Function) ~ "SCL23",
      grepl("CAMTA", Short_Function) ~ "CAMTA4",
      grepl("PhyB|Phytochrome B", Short_Function) ~ "PhyB",
      grepl("CYP83B1", Short_Function) ~ "CYP83B1",
      grepl("ent-Kaurene|kaurene", Short_Function) ~ "KS",
      grepl("EXO70", Short_Function) ~ "EXO70B1",
      grepl("PDI|disulfide", Short_Function) ~ "PDI",
      grepl("zeatin|glucosyltransferase.*cytokinin", Short_Function, ignore.case = TRUE) ~ "cZOGT",
      grepl("Zeatin", Short_Function) ~ "cZOGT",
      grepl("HSP70|Heat shock 70", Short_Function) ~ "HSP70",
      grepl("PRR5|Pseudo-response", Short_Function) ~ "PRR5",
      grepl("TBP2|TATA", Short_Function) ~ "TBP2",
      grepl("ADA2|Transcriptional adapter", Short_Function) ~ "ADA2a",
      grepl("CDT1", Short_Function) ~ "CDT1",
      grepl("ANXUR", Short_Function) ~ "ANXUR",
      grepl("Clathrin", Short_Function) ~ "CLC",
      grepl("BPG2|brassinazole", Short_Function) ~ "BPG2",
      grepl("TIR-NBS|NB-ARC|DISEASE", Short_Function, ignore.case = TRUE) ~ "NBS-LRR",
      grepl("DnaJ.*ERDJ3", Short_Function) ~ "ERDJ3A",
      grepl("DnaJ.*mitoch", Short_Function) ~ "DnaJ1",
      grepl("AAP6|Amino acid permease 6", Short_Function) ~ "AAP6",
      grepl("LHT8|Lysine histidine", Short_Function) ~ "LHT8",
      grepl("CAX2|cation.*exchanger", Short_Function, ignore.case = TRUE) ~ "CAX2",
      grepl("KEA1|antiporter", Short_Function) ~ "KEA1",
      grepl("RPL25|Ribosomal", Short_Function) ~ "RPL25",
      grepl("PPR|Pentatrico", Short_Function) ~ "PPR",
      grepl("RRM|RNA-binding", Short_Function) ~ "RBP",
      grepl("MIOREX", Short_Function) ~ "MIOREX",
      grepl("Mitofilin", Short_Function) ~ "Mitofilin",
      grepl("Metaxin", Short_Function) ~ "Metaxin",
      grepl("THF1|Thylakoid", Short_Function) ~ "THF1",
      grepl("F-box", Short_Function) ~ "F-box",
      grepl("NHL", Short_Function) ~ "NHL",
      grepl("LURP", Short_Function) ~ "LURP1",
      grepl("PDCD5|cell death", Short_Function) ~ "PDCD5",
      grepl("GST|Glutathione", Short_Function) ~ "GSTU8",
      grepl("ECH1|Dienoyl", Short_Function) ~ "ECH1",
      grepl("HMA", Short_Function) ~ "HMA",
      grepl("DUF247", Short_Function) ~ "DUF247",
      grepl("ADP-ribosyl", Short_Function) ~ "ADPRC",
      grepl("O-fucosyl", Short_Function) ~ "OFUT",
      grepl("Fucosyl.*glyco", Short_Function) ~ "FUT",
      grepl("linamarin", Short_Function) ~ "UGT",
      grepl("2-oxo.*oxygenase", Short_Function) ~ "2OG-FeII",
      grepl("Fe2OG", Short_Function) ~ "Fe2OG",
      grepl("Cytochrome P450.*oxidative", Short_Function) ~ "CYP",
      grepl("glycolate", Short_Function) ~ "PLGG1",
      grepl("No Arabidopsis", Short_Function) ~ "unknown",
      grepl("Xylanase", Short_Function) ~ "XIP",
      TRUE ~ ""
    ),
    Gene_Label = ifelse(Func_Short != "",
                        paste0(GeneID, " (", Func_Short, ")"),
                        GeneID)
  )

# For NBS-LRR duplicates, add a suffix
nbs_idx <- which(dat$Func_Short == "NBS-LRR")
if (length(nbs_idx) > 1) {
  for (i in seq_along(nbs_idx)) {
    dat$Gene_Label[nbs_idx[i]] <- paste0(dat$GeneID[nbs_idx[i]], " (NBS-LRR-", i, ")")
  }
}

# Same for PPR duplicates
ppr_idx <- which(dat$Func_Short == "PPR")
if (length(ppr_idx) > 1) {
  for (i in seq_along(ppr_idx)) {
    dat$Gene_Label[ppr_idx[i]] <- paste0(dat$GeneID[ppr_idx[i]], " (PPR-", i, ")")
  }
}

# Same for MIOREX duplicates
mio_idx <- which(dat$Func_Short == "MIOREX")
if (length(mio_idx) > 1) {
  for (i in seq_along(mio_idx)) {
    dat$Gene_Label[mio_idx[i]] <- paste0(dat$GeneID[mio_idx[i]], " (MIOREX-", i, ")")
  }
}

# Same for DUF247 duplicates
duf_idx <- which(dat$Func_Short == "DUF247")
if (length(duf_idx) > 1) {
  for (i in seq_along(duf_idx)) {
    dat$Gene_Label[duf_idx[i]] <- paste0(dat$GeneID[duf_idx[i]], " (DUF247-", i, ")")
  }
}

# Same for unknown duplicates
unk_idx <- which(dat$Func_Short == "unknown")
if (length(unk_idx) > 1) {
  for (i in seq_along(unk_idx)) {
    dat$Gene_Label[unk_idx[i]] <- paste0(dat$GeneID[unk_idx[i]], " (unknown-", i, ")")
  }
}

# =============================================================================
# FIGURE A: Bar + Dot Plot
# =============================================================================

cat("\n=== Building Figure A: Bar + Dot Plot ===\n")

# Sort by FC_TP1
plot_dat <- dat %>%
  arrange(FC_TP1) %>%
  mutate(Gene_Label = factor(Gene_Label, levels = Gene_Label))

# Panel 1: Horizontal bar plot of FC_TP1
p_bar <- ggplot(plot_dat, aes(x = FC_TP1, y = Gene_Label, fill = Broad_Category)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = cat_colors, name = "Functional\nCategory") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = expression(log[2]*"FC at TP1"), y = NULL) +
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

# Panel 2: Dot plot of r(LW)
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

# Combine
library(patchwork)

p_combined <- p_bar + p_dot +
  plot_layout(widths = c(4, 1)) +
  plot_annotation(
    title = paste0(nrow(dat), " Validated GmJAG1 Targets: Expression Change at TP1 and Leaf Shape Correlation"),
    subtitle = expression("Bars: "*log[2]*"FC at TP1 (peak GmJAG1 activity) | Dots: |r| with leaf length:width ratio"),
    theme = theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  )

ggsave(file.path(out_dir, "Fig_Validated_Genes_BarDot_TP1.png"),
       p_combined, width = 12, height = 16, dpi = 300)

cat("Saved: Fig_Validated_Genes_BarDot_TP1.png\n")

# =============================================================================
# FIGURE B: ComplexHeatmap
# =============================================================================

cat("\n=== Building Figure B: ComplexHeatmap ===\n")

# Sort by functional category then by FC_TP1 within category
dat_sorted <- dat %>%
  arrange(Broad_Category, desc(FC_TP1))

# Heatmap matrix: logFC across timepoints
mat <- dat_sorted %>%
  select(FC_TP0, FC_TP1, FC_TP2, FC_TP3, FC_TP4) %>%
  as.matrix()
rownames(mat) <- dat_sorted$Gene_Label
colnames(mat) <- c("TP0", "TP1", "TP2", "TP3", "TP4")

# Color function for heatmap - diverging
max_val <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-5, 0, 5, 10, 15),
                       c("#2166AC", "white", "#FDDBC7", "#D6604D", "#67001F"))

# Row annotations
# 1. Functional category
row_cat <- dat_sorted$Broad_Category
cat_cols_ha <- cat_colors[unique(row_cat)]

# 2. r(LW) as bar
r_vals <- dat_sorted$r_LW

# 3. DE Tier
tier_colors <- c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")
de_tier <- dat_sorted$DE_Tier

# 4. Binding class
bind_colors <- c("ChIP_Only" = "#1B9E77", "ChIP_AND_DAP" = "#D95F02", "DAP_Only" = "#7570B3")
bind_class <- dat_sorted$Binding_Class

# 5. Number of ChIP peaks
n_peaks <- as.numeric(dat_sorted$ChIP_Peaks)
n_peaks[is.na(n_peaks)] <- 0

# Left annotation: DE tier + binding (no category color since rows are already split by category)
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

# Right annotation: r(LW) barplot + ChIP peaks
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

# Draw heatmap
ht <- Heatmap(mat,
              name = "log2FC",
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
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

# Save PNG only
png(file.path(out_dir, "Fig_Validated_Genes_ComplexHeatmap.png"),
    width = 14, height = 18, units = "in", res = 300)
draw(ht, padding = unit(c(5, 5, 5, 5), "mm"),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()

cat("Saved: Fig_Validated_Genes_ComplexHeatmap.png\n")

cat("\n=== COMPLETE ===\n")
