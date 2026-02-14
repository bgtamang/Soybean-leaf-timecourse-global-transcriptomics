# Script 10: Differential Expression Setup

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 10: DIFFERENTIAL EXPRESSION SETUP\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",
  "limma",
  "dplyr"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded checkpoint from script 06\n\n")

# Use primary dataset for main analysis
# dge_corrected contains all samples with batch correction
# targets_primary contains the sample metadata
dge <- dge_corrected
targets <- targets_primary

# Ensure Sample column exists
if (!"Sample" %in% colnames(targets)) {
  targets$Sample <- rownames(targets)
}

# Create group column if it doesn't exist (Line_Timepoint format)
if (!"group" %in% colnames(targets)) {
  targets$group <- paste(targets$Line, targets$Timepoint, sep = "_")
}

cat("Using PRIMARY dataset:\n")
cat("  Samples:", ncol(dge), "\n")
cat("  Genes:", nrow(dge), "\n\n")

# ===== SECTION 2: VERIFY EXPERIMENTAL GROUPS =====

cat("========================================\n")
cat("SECTION 2: EXPERIMENTAL GROUPS\n")
cat("========================================\n\n")

# Check group structure
cat("Experimental groups:\n")
print(table(targets$group))

cat("\nGroup breakdown:\n")
cat("  Lines:", paste(unique(targets$Line), collapse = ", "), "\n")
cat("  Timepoints:", paste(sort(unique(targets$Timepoint)), collapse = ", "), "\n")
cat("  Leaf types:", paste(unique(targets$Leaf_type), collapse = ", "), "\n\n")

# Ensure group is a factor with proper levels
targets$group <- factor(targets$group)
group_levels <- levels(targets$group)

cat("Group levels (", length(group_levels), " total):\n", sep = "")
print(group_levels)

# ===== SECTION 3: CREATE DESIGN MATRIX =====

cat("\n========================================\n")
cat("SECTION 3: DESIGN MATRIX\n")
cat("========================================\n\n")

# Create design matrix with one coefficient per group (no intercept)
design <- model.matrix(~ 0 + group, data = targets)

# Clean column names (remove "group" prefix)
colnames(design) <- gsub("^group", "", colnames(design))

# Set rownames to sample names
rownames(design) <- targets$Sample

cat("Design matrix created:\n")
cat("  Dimensions:", nrow(design), "samples x", ncol(design), "groups\n")
cat("  Columns:", paste(colnames(design), collapse = ", "), "\n\n")

# Save design matrix
write.csv(design, "03_results/tables/design_matrix.csv")
cat("Saved: 03_results/tables/design_matrix.csv\n\n")

# ===== SECTION 4: CREATE CONTRAST MATRIX =====

cat("========================================\n")
cat("SECTION 4: CONTRAST MATRIX\n")
cat("========================================\n\n")

cat("Creating contrasts for all relevant comparisons...\n\n")

# Get unique lines and timepoints
lines <- unique(targets$Line)
timepoints <- sort(unique(targets$Timepoint))

# Build contrast list dynamically
contrast_list <- list()
contrast_desc <- data.frame(
  Contrast = character(),
  Description = character(),
  Type = character(),
  stringsAsFactors = FALSE
)

# --- Type 1: Temporal contrasts (TPx vs TP0) within each line ---
cat("Type 1: Temporal contrasts (vs TP0)...\n")

for (line in lines) {
  for (tp in timepoints[timepoints != "TP0"]) {
    contrast_name <- paste0(line, "_", tp, "vsTP0")
    group_tp <- paste0(line, "_", tp)
    group_tp0 <- paste0(line, "_TP0")

    if (group_tp %in% colnames(design) && group_tp0 %in% colnames(design)) {
      contrast_list[[contrast_name]] <- paste0(group_tp, " - ", group_tp0)
      contrast_desc <- rbind(contrast_desc, data.frame(
        Contrast = contrast_name,
        Description = paste0(line, ": ", tp, " vs TP0"),
        Type = "Temporal_vsTP0"
      ))
    }
  }
}

# --- Type 2: Temporal contrasts (consecutive timepoints) ---
cat("Type 2: Temporal contrasts (consecutive)...\n")

tp_pairs <- list(
  c("TP2", "TP1"), c("TP3", "TP1"), c("TP4", "TP1"),
  c("TP3", "TP2"), c("TP4", "TP2"),
  c("TP4", "TP3")
)

for (line in lines) {
  for (pair in tp_pairs) {
    tp_a <- pair[1]
    tp_b <- pair[2]
    contrast_name <- paste0(line, "_", tp_a, "vs", tp_b)
    group_a <- paste0(line, "_", tp_a)
    group_b <- paste0(line, "_", tp_b)

    if (group_a %in% colnames(design) && group_b %in% colnames(design)) {
      contrast_list[[contrast_name]] <- paste0(group_a, " - ", group_b)
      contrast_desc <- rbind(contrast_desc, data.frame(
        Contrast = contrast_name,
        Description = paste0(line, ": ", tp_a, " vs ", tp_b),
        Type = "Temporal_consecutive"
      ))
    }
  }
}

# --- Type 3: Between-line comparisons at TP0 ---
cat("Type 3: Between-line comparisons at TP0...\n")

# Broad vs Narrow at TP0
broad_lines <- c("PI532462A", "LD112170")
narrow_lines <- c("PI612713B", "PI547745")

for (broad in broad_lines) {
  for (narrow in narrow_lines) {
    contrast_name <- paste0(narrow, "vs", broad, "_TP0")
    group_narrow <- paste0(narrow, "_TP0")
    group_broad <- paste0(broad, "_TP0")

    if (group_narrow %in% colnames(design) && group_broad %in% colnames(design)) {
      contrast_list[[contrast_name]] <- paste0(group_narrow, " - ", group_broad)
      contrast_desc <- rbind(contrast_desc, data.frame(
        Contrast = contrast_name,
        Description = paste0("TP0: ", narrow, " (Narrow) vs ", broad, " (Broad)"),
        Type = "Between_line_TP0"
      ))
    }
  }
}

# Within leaf-type comparisons at TP0
# Broad vs Broad
contrast_list[["PI532462Avs LD112170_TP0"]] <- "PI532462A_TP0 - LD112170_TP0"
contrast_desc <- rbind(contrast_desc, data.frame(
  Contrast = "PI532462AvsLD112170_TP0",
  Description = "TP0: PI532462A vs LD112170 (both Broad)",
  Type = "Within_leaftype_TP0"
))

# Narrow vs Narrow
contrast_list[["PI547745vsPI612713B_TP0"]] <- "PI547745_TP0 - PI612713B_TP0"
contrast_desc <- rbind(contrast_desc, data.frame(
  Contrast = "PI547745vsPI612713B_TP0",
  Description = "TP0: PI547745 vs PI612713B (both Narrow)",
  Type = "Within_leaftype_TP0"
))

# --- Type 4: Pooled Leaf type comparison at TP0 ---
cat("Type 4: Pooled leaf type comparison...\n")

# Average of narrow lines vs average of broad lines at TP0
contrast_list[["NarrowvsBroad_TP0"]] <-
  "(PI612713B_TP0 + PI547745_TP0)/2 - (PI532462A_TP0 + LD112170_TP0)/2"
contrast_desc <- rbind(contrast_desc, data.frame(
  Contrast = "NarrowvsBroad_TP0",
  Description = "TP0: Narrow (pooled) vs Broad (pooled)",
  Type = "Pooled_leaftype"
))

# Create contrast matrix using makeContrasts
cat("\nBuilding contrast matrix...\n")

# Convert list to makeContrasts call
contrast_strings <- unlist(contrast_list)
contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)

# Set proper column names
colnames(contrast_matrix) <- names(contrast_list)

cat("Contrast matrix created:\n")
cat("  Dimensions:", nrow(contrast_matrix), "groups x", ncol(contrast_matrix), "contrasts\n")
cat("  Total contrasts:", ncol(contrast_matrix), "\n\n")

# Summarize by type
cat("Contrasts by type:\n")
print(table(contrast_desc$Type))

# ===== SECTION 5: KEY CONTRASTS FOR JAG1 ANALYSIS =====

cat("\n========================================\n")
cat("SECTION 5: KEY CONTRASTS FOR JAG1\n")
cat("========================================\n\n")

# Identify the key contrasts for JAG1 target identification
# These are the Narrow vs Broad comparisons at TP0

key_contrasts <- contrast_desc$Contrast[contrast_desc$Type == "Between_line_TP0"]

cat("Key contrasts for JAG1 target identification:\n")
cat("  (Narrow vs Broad at TP0 - genes UP in narrow are potential JAG1 targets)\n\n")

for (kc in key_contrasts) {
  desc <- contrast_desc$Description[contrast_desc$Contrast == kc]
  cat("  -", kc, ":", desc, "\n")
}

cat("\n  Plus pooled comparison:\n")
cat("  - NarrowvsBroad_TP0 : TP0: Narrow (pooled) vs Broad (pooled)\n")

# ===== SECTION 6: SAVE OUTPUTS =====

cat("\n========================================\n")
cat("SECTION 6: SAVE OUTPUTS\n")
cat("========================================\n\n")

# Save contrast matrix
write.csv(contrast_matrix, "03_results/tables/contrast_matrix.csv")
cat("Saved: 03_results/tables/contrast_matrix.csv\n")

# Save contrast descriptions
write.csv(contrast_desc, "03_results/tables/contrast_descriptions.csv", row.names = FALSE)
cat("Saved: 03_results/tables/contrast_descriptions.csv\n")

# Save checkpoint
save(dge, targets, design, contrast_matrix, contrast_list, contrast_desc,
     key_contrasts, broad_lines, narrow_lines,
     file = "03_results/checkpoints/10_DE_setup.RData")
cat("Saved: 03_results/checkpoints/10_DE_setup.RData\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 10: DE SETUP - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Design matrix:", nrow(design), "samples x", ncol(design), "groups\n")
cat("    - Contrast matrix:", ncol(contrast_matrix), "contrasts\n")
cat("    - Temporal (vs TP0):", sum(contrast_desc$Type == "Temporal_vsTP0"), "\n")
cat("    - Temporal (consecutive):", sum(contrast_desc$Type == "Temporal_consecutive"), "\n")
cat("    - Between-line (TP0):", sum(contrast_desc$Type == "Between_line_TP0"), "\n")
cat("    - Within leaf-type:", sum(contrast_desc$Type == "Within_leaftype_TP0"), "\n")
cat("    - Pooled leaf-type:", sum(contrast_desc$Type == "Pooled_leaftype"), "\n")
cat("\n")
cat("  Next: Run 11_DE_analysis.R\n")
cat("================================================================\n")
