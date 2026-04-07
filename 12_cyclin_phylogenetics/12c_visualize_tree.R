#!/usr/bin/env Rscript
# ============================================================
# 03_visualize_tree.R
# Publication-quality D-type cyclin phylogenetic tree (ggtree)
# GmJAG1 Manuscript - Task 02
#
# Shows: JAG1 targets are in CYCD1/CYCD6, not CYCD3
# ============================================================

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(treeio)
    library(ggplot2)
    library(dplyr)
    library(readr)
})

# ============================================================
# PATHS
# ============================================================

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
    base_dir <- dirname(dirname(normalizePath(script_path)))
} else {
    base_dir <- file.path(
        "C:/Users/bgtam/OneDrive - University of Illinois - Urbana",
        "USDA-UIUC WORKS/Bishal/Projects/Data Collections/MANUSCRIPTS",
        "5. Soybean RNASeq leaf timeseries/Soybean-RNASEQ",
        "Phase3-Refined-Analysis/Manuscript/TASK_02_Cyclin_Phylogenetics"
    )
}

cat("Base directory:", base_dir, "\n")

# Windows MAX_PATH workaround
temp_data <- "C:/Users/bgtam/cycd_temp"
if (dir.exists(temp_data)) {
    tree_file    <- file.path(temp_data, "CYCD_tree.treefile")
    metadata_csv <- file.path(temp_data, "CYCD_sequence_metadata.csv")
} else {
    tree_file    <- file.path(base_dir, "03_tree", "CYCD_tree.treefile")
    metadata_csv <- file.path(base_dir, "05_reports", "CYCD_sequence_metadata.csv")
}

out_dir <- file.path(base_dir, "04_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (nchar(file.path(out_dir, "CYCD_phylogeny_annotated.pdf")) > 255) {
    out_dir <- temp_data
}

out_pdf <- file.path(out_dir, "CYCD_phylogeny_annotated.pdf")
out_png <- file.path(out_dir, "CYCD_phylogeny_annotated.png")

# ============================================================
# LOAD DATA
# ============================================================

cat("Loading data...\n")
tree <- read.tree(tree_file)
meta <- read_csv(metadata_csv, show_col_types = FALSE)

# Fix label mismatch: IQ-TREE replaced ; with _
meta_join <- meta %>% mutate(tree_label = gsub(";", "_", label))

# ============================================================
# PRUNE: drop poplar, trim rice to 1 per subclade
# Result: Soybean (39) + Arabidopsis (10) + Medicago (15) +
#         Rice reps (6) = 70 tips
# ============================================================

cat("Pruning tree...\n")

poplar_tips <- grep("^Ptr_", tree$tip.label, value = TRUE)

rice_keep <- c(
    "Osa_LOC_Os08g32540_CYCD1",
    "Osa_LOC_Os07g42860_CYCD2",
    "Osa_LOC_Os09g02360_CYCD3",
    "Osa_LOC_Os09g29100_CYCD4",
    "Osa_LOC_Os03g42070_CYCD5",
    "Osa_LOC_Os07g37010_CYCD6"
)
rice_drop <- setdiff(grep("^Osa_", tree$tip.label, value = TRUE), rice_keep)

tree <- drop.tip(tree, c(poplar_tips, rice_drop))
meta_join <- meta_join %>% filter(tree_label %in% tree$tip.label)
cat("  Pruned to", Ntip(tree), "tips\n")

# ============================================================
# ROOT on rice CYCD3 (longest branch, clear outgroup)
# ============================================================

tree <- root(tree, outgroup = "Osa_LOC_Os09g02360_CYCD3", resolve.root = TRUE)

# ============================================================
# TIP ANNOTATIONS
# ============================================================

tip_df <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
    left_join(
        meta_join %>% select(tree_label, gene_id, species, subclade,
                              JAG1_target, JAG_reference),
        by = c("label" = "tree_label")
    ) %>%
    mutate(
        display_label = case_when(
            species == "Arabidopsis" ~ subclade,
            TRUE ~ gene_id
        ),
        subclade_major = case_when(
            grepl("CYCD1", subclade) ~ "CYCD1",
            grepl("CYCD[24]", subclade) ~ "CYCD2/4",
            grepl("CYCD3", subclade) ~ "CYCD3",
            grepl("CYCD5", subclade) ~ "CYCD5",
            grepl("CYCD6", subclade) ~ "CYCD6",
            grepl("CYCD7", subclade) ~ "CYCD7",
            TRUE ~ NA_character_
        ),
        species_name = case_when(
            species == "Arabidopsis" ~ "A. thaliana",
            species == "Soybean"     ~ "G. max",
            species == "Medicago"    ~ "M. truncatula",
            species == "Rice"        ~ "O. sativa",
            TRUE ~ species
        ),
        highlight = case_when(
            JAG1_target == "YES"   ~ "JAG1_target",
            JAG_reference == "YES" ~ "JAG_ref",
            TRUE                   ~ "none"
        )
    )

# ============================================================
# FIND SUBCLADE MRCA NODES (with outlier removal)
# ============================================================

cat("Finding subclade clades...\n")

get_clade_tips <- function(tree, node) {
    extract.clade(tree, node)$tip.label
}

find_pure_clade <- function(tree, tips, min_purity = 0.65, min_keep = 0.5) {
    if (length(tips) < 2) return(NULL)
    all_dists <- cophenetic.phylo(tree)
    current_tips <- tips
    max_remove <- ceiling(length(tips) * (1 - min_keep))
    for (i in seq_len(max_remove + 1)) {
        if (length(current_tips) < 2) break
        node <- getMRCA(tree, current_tips)
        if (is.null(node)) break
        desc <- get_clade_tips(tree, node)
        purity <- sum(desc %in% current_tips) / length(desc)
        if (purity >= min_purity) {
            return(list(node = node, purity = purity,
                        n_match = sum(desc %in% current_tips),
                        n_total = length(desc)))
        }
        sub_d <- all_dists[current_tips, current_tips]
        outlier <- names(which.max(rowMeans(sub_d)))
        current_tips <- setdiff(current_tips, outlier)
    }
    return(NULL)
}

subclade_nodes <- list()
for (sc in unique(na.omit(tip_df$subclade_major))) {
    sc_tips <- tip_df$label[tip_df$subclade_major == sc]
    if (length(sc_tips) < 2) next

    node <- getMRCA(tree, sc_tips)
    if (is.null(node)) next

    desc_tips <- get_clade_tips(tree, node)
    purity <- sum(desc_tips %in% sc_tips) / length(desc_tips)

    if (purity >= 0.55) {
        subclade_nodes[[sc]] <- list(node = node, purity = purity)
        cat(sprintf("  %s: node=%d (%.0f%% pure)\n", sc, node, purity * 100))
    } else {
        result <- find_pure_clade(tree, sc_tips)
        if (!is.null(result)) {
            subclade_nodes[[sc]] <- result
            cat(sprintf("  %s: node=%d (%.0f%% pure, after outlier removal)\n",
                        sc, result$node, result$purity * 100))
        }
    }
}

# ============================================================
# COLOR PALETTES
# ============================================================

subclade_fills <- c(
    "CYCD1"   = "#FFF3E0",
    "CYCD2/4" = "#E8F5E9",
    "CYCD3"   = "#E3F2FD",
    "CYCD5"   = "#F3E5F5",
    "CYCD6"   = "#FCE4EC",
    "CYCD7"   = "#FFFDE7"
)

subclade_bar_colors <- c(
    "CYCD1"   = "#E65100",
    "CYCD2/4" = "#2E7D32",
    "CYCD3"   = "#1565C0",
    "CYCD5"   = "#6A1B9A",
    "CYCD6"   = "#C62828",
    "CYCD7"   = "#F9A825"
)

species_colors <- c(
    "G. max"        = "#D32F2F",
    "A. thaliana"   = "#1976D2",
    "M. truncatula" = "#388E3C",
    "O. sativa"     = "#F57C00"
)

# ============================================================
# BUILD FIGURE
# ============================================================

cat("Building figure...\n")

# --- Base tree: thick branches ---
p <- ggtree(tree, layout = "rectangular", linewidth = 0.8,
            color = "grey20", ladderize = TRUE) %<+% tip_df

# --- Subclade highlights (tight around tips) ---
for (sc in names(subclade_nodes)) {
    p <- p + geom_hilight(
        node = subclade_nodes[[sc]]$node,
        fill = subclade_fills[sc],
        alpha = 0.5,
        extend = 0.35
    )
}

# --- Subclade clade labels (bars placed well past highlights) ---
for (sc in names(subclade_nodes)) {
    p <- p + geom_cladelab(
        node = subclade_nodes[[sc]]$node,
        label = sc,
        fontsize = 4.5,
        fontface = "bold",
        color = subclade_bar_colors[sc],
        barcolor = subclade_bar_colors[sc],
        barsize = 2.5,
        offset = 0.55,
        offset.text = 0.03,
        align = TRUE,
        angle = 0,
        hjust = -0.05
    )
}

# --- Tip labels: large, species-colored ---
p <- p +
    geom_tiplab(
        aes(color = species_name, label = display_label),
        size = 2.8,
        align = FALSE,
        offset = 0.05
    ) +
    scale_color_manual(
        values = species_colors,
        name = "Species",
        guide = guide_legend(order = 1, override.aes = list(size = 4))
    )

# --- JAG1 target markers (large red stars) ---
p <- p +
    geom_tippoint(
        data = function(d) d[d$isTip & !is.na(d$highlight) &
                               d$highlight == "JAG1_target", ],
        shape = 8, size = 2.25, color = "#D32F2F", stroke = 1
    )

# --- AtCYCD3;3 reference (large blue diamond) ---
p <- p +
    geom_tippoint(
        data = function(d) d[d$isTip & !is.na(d$highlight) &
                               d$highlight == "JAG_ref", ],
        shape = 18, size = 5, color = "#1565C0"
    )

# --- Bootstrap at key nodes (SH-aLRT >= 90 only) ---
p <- p +
    geom_text2(
        data = function(d) {
            d$bs_val <- suppressWarnings(
                as.numeric(sub("/.*", "", d$label))
            )
            d[!d$isTip & !is.na(d$bs_val) & d$bs_val >= 90, ]
        },
        aes(label = round(
            suppressWarnings(as.numeric(sub("/.*", "", label)))
        )),
        size = 2.2, hjust = 1.2, vjust = -0.4, color = "grey35"
    )

# --- X-axis space for labels ---
tree_fort <- fortify(tree)
max_x <- max(tree_fort$x, na.rm = TRUE)

# --- Theme ---
p <- p +
    theme_tree2() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.82, 0.55),
        legend.background = element_rect(fill = "white", color = "grey80",
                                          linewidth = 0.3),
        legend.text = element_text(size = 10, face = "italic"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.margin = margin(5, 8, 5, 8),
        plot.title = element_text(size = 14, face = "bold", hjust = 0,
                                  margin = margin(b = 4)),
        plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0,
                                     margin = margin(b = 8)),
        plot.caption = element_text(size = 8, color = "grey40", hjust = 0,
                                    margin = margin(t = 10)),
        plot.margin = margin(15, 15, 15, 15),
        axis.text.x = element_text(size = 7, color = "grey40")
    ) +
    labs(
        title = "D-type Cyclin Phylogeny",
        subtitle = paste0(
            "Q.plant+I+G4 | 1,000 UFBoot + SH-aLRT | ",
            Ntip(tree), " sequences (4 species) | ",
            "Rooted on O. sativa"
        ),
        caption = paste0(
            "*  GmJAG1-regulated D-type cyclins     ",
            "<>  AtCYCD3;3 (Arabidopsis JAG target; Schiessl et al. 2014)     ",
            "SH-aLRT >= 90 shown at nodes"
        )
    ) +
    xlim(NA, max_x * 2.2)

# ============================================================
# SAVE
# ============================================================

cat("Saving figures...\n")

fig_height <- max(9, Ntip(tree) * 0.17)
fig_width <- 12

ggsave(out_pdf, p, width = fig_width, height = fig_height,
       units = "in", limitsize = FALSE)
cat("  PDF:", out_pdf, "\n")

ggsave(out_png, p, width = fig_width, height = fig_height,
       units = "in", dpi = 300, limitsize = FALSE)
cat("  PNG:", out_png, "\n")

cat(sprintf("\nFigure: %.1f x %.1f inches\n", fig_width, fig_height))
cat("Done!\n")
