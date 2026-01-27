# Build Two Tables for 72 Validated GmJAG1 Targets
# 1. Supplementary Table: All 72 genes with corrected functional categories
# 2. Main Text Table: ~16 literature-supported genes with citations

rm(list = ls())
gc()

library(tidyverse)

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))

out_dir <- "Manuscript/Figure/Figure6/72_Genes_Analysis"

# 1. LOAD DATA

targets <- read_csv("03_results/tables/supplementary_v6/TableS30_HighConfidence_72_Targets.csv",
                    show_col_types = FALSE)

pheno <- read_csv("Manuscript/Figure/Figure6/Fig6D_All_Target_Phenotype_Correlations.csv",
                  show_col_types = FALSE)

hub <- read_csv("03_results/tables/WGCNA/hub_genes.csv", show_col_types = FALSE)

cat("Loaded", nrow(targets), "targets,", nrow(pheno), "phenotype correlations\n")

# 2. MERGE

merged <- targets %>%
  left_join(pheno %>% select(GeneID, Cor_LW_Ratio, Pval_LW_Ratio),
            by = "GeneID") %>%
  left_join(hub %>% select(GeneID, Module_Color, kME),
            by = "GeneID")

# 3. CORRECTED FUNCTIONAL CATEGORIZATION

# Manual categorization based on actual Arabidopsis ortholog knowledge
# This replaces the flawed keyword-based approach

category_map <- tribble(
  ~GeneID,              ~Functional_Category,       ~Short_Function,
  # --- Hormone signaling ---
  "Glyma.06G305800",    "Hormone (Auxin)",          "NEDD8-activating enzyme AXR1; auxin perception via SCF-TIR1",
  "Glyma.02G061400",    "Hormone (Auxin)",          "NPH3 family; phototropin-mediated auxin redistribution",
  "Glyma.11G168428",    "Hormone (Auxin)",          "CYP83B1; cytochrome P450 in auxin homeostasis (IAOx pathway)",
  "Glyma.04G013800",    "Hormone (JA)",             "JAZ1/TIFY10A; jasmonate signaling repressor",
  "Glyma.11G225500",    "Hormone (Cytokinin)",      "Zeatin O-glucosyltransferase; cytokinin inactivation",
  "Glyma.15G263332",    "Hormone (GA)",             "ent-Kaurene synthase; gibberellin biosynthesis",
  "Glyma.15G196500",    "Hormone (Light/GA)",       "Phytochrome B; light perception, GA signaling suppression",
  "Glyma.04G180400",    "Hormone (ABA)",            "Dehydration-responsive protein RD22; ABA signaling effector",
  "Glyma.05G204800",    "Hormone (ABA)",            "Osmotin-like protein OSM34; ABA receptor interactor (PR5 family)",
  "Glyma.19G133900",    "Hormone (BR)",             "BPG2; brassinazole-insensitive pale green 2, chloroplastic",
  # --- Defense/NBS-LRR ---
  "Glyma.03G052800",    "Defense (NBS-LRR)",        "TIR-NBS-LRR disease resistance protein",
  "Glyma.06G310000",    "Defense (NBS-LRR)",        "TIR-NBS-LRR disease resistance protein",
  "Glyma.16G087100",    "Defense (NBS-LRR)",        "TIR-NBS-LRR disease resistance protein",
  "Glyma.16G213800",    "Defense (NBS-LRR)",        "TIR-NBS-LRR disease resistance protein (with LRR domain)",
  "Glyma.16G215000",    "Defense (NBS-LRR)",        "TIR-NBS-LRR disease resistance protein",
  "Glyma.14G204600",    "Defense (NB-ARC)",         "NB-ARC domain protein; innate immunity",
  "Glyma.14G205000",    "Defense (NB-ARC)",         "NB-ARC domain protein; innate immunity",
  "Glyma.14G205300",    "Defense (NB-ARC)",         "NB-ARC domain protein; innate immunity",
  # --- Receptor kinase / Signaling ---
  "Glyma.19G167300",    "Receptor kinase",          "MIK2/LRR-KISS; cell wall integrity sensing RLK",
  "Glyma.18G273100",    "Receptor kinase",          "ANXUR1-related RLK; CrRLK1L family, cell wall integrity",
  "Glyma.12G198000",    "Receptor kinase",          "Receptor-like Ser/Thr protein kinase",
  "Glyma.06G184400",    "Receptor kinase",          "LRR receptor-like kinase",
  "Glyma.16G202200",    "Receptor kinase",          "Ser/Thr protein kinase with retrotransposon domain",
  "Glyma.16G018000",    "Circadian clock",          "Pseudo-response regulator 5 (PRR5); circadian clock component",
  # --- Transcription / Chromatin ---
  "Glyma.02G226700",    "Transcription/Chromatin",  "Elongator complex protein 5 (ELP5); RNAPII transcription elongation",
  "Glyma.11G251900",    "Transcription/Chromatin",  "Calmodulin-binding transcription activator 4 (CAMTA4)",
  "Glyma.12G022600",    "Transcription/Chromatin",  "SCARECROW-LIKE 23 (SCL23); GRAS family TF, vascular patterning",
  "Glyma.09G276900",    "Transcription/Chromatin",  "TATA-box binding protein 2 (TBP2); basal transcription",
  "Glyma.20G062300",    "Transcription/Chromatin",  "Transcriptional adapter 2-alpha (ADA2a); histone acetylation",
  "Glyma.02G189700",    "Epigenetics",              "ARGONAUTE 4 (AGO4); siRNA-directed DNA methylation (RdDM)",
  # --- Transport ---
  "Glyma.03G171600",    "Transport",                "K+/H+ antiporter 1 (KEA1); potassium/proton exchange",
  "Glyma.05G102600",    "Transport",                "Vacuolar cation/proton exchanger 2 (CAX2)",
  "Glyma.11G097000",    "Transport",                "Lysine histidine transporter-like 8 (LHT8); amino acid transport",
  "Glyma.17G192000",    "Transport",                "Amino acid permease 6 (AAP6)",
  "Glyma.09G211200",    "Transport",                "Plastidal glycolate/glycerate translocator 1; photorespiration",
  # --- Protein folding / ER ---
  "Glyma.02G093200",    "Protein folding/ER",       "Heat shock 70 kDa protein 4 (HSP70-4); chaperone",
  "Glyma.13G357700",    "Protein folding/ER",       "Protein disulfide isomerase-like 1-5 (PDI4); ER protein folding",
  "Glyma.14G178800",    "Protein folding/ER",       "DnaJ chaperone 1, mitochondrial; co-chaperone",
  "Glyma.16G098700",    "Protein folding/ER",       "DnaJ protein ERDJ3A; ER-localized co-chaperone",
  # --- Vesicle trafficking / Secretory ---
  "Glyma.14G199600",    "Vesicle/Secretory",        "Exocyst complex component EXO70B1; polarized secretion",
  "Glyma.01G115100",    "Vesicle/Secretory",        "Clathrin light chain; endocytosis/vesicle trafficking",
  # --- RNA processing / Translation ---
  "Glyma.13G136500",    "RNA processing",           "RNA-binding protein (RRM domain); post-transcriptional regulation",
  "Glyma.17G178001",    "Translation",              "Ribosomal protein L25; large subunit",
  "Glyma.09G000400",    "RNA processing",           "Pentatricopeptide repeat (PPR-DYW); organellar RNA editing",
  "Glyma.10G043900",    "RNA processing",           "Pentatricopeptide repeat (PPR); organellar RNA processing",
  # --- Metabolism / Enzyme ---
  "Glyma.03G073702",    "Metabolism",               "Glutathione S-transferase U8 (GSTU8); detoxification",
  "Glyma.03G119400",    "Metabolism",               "Dienoyl-CoA isomerase (ECH1); fatty acid beta-oxidation",
  "Glyma.14G058602",    "Metabolism",               "2-oxoglutarate/Fe(II)-dependent oxygenase; oxidative modification",
  "Glyma.06G137000",    "Metabolism",               "Fe2OG dioxygenase domain protein",
  "Glyma.19G029700",    "Metabolism",               "UDP-glucosyltransferase; linamarin synthase",
  "Glyma.03G250300",    "Metabolism",               "Fucosyltransferase; glycoprotein modification",
  "Glyma.14G118550",    "Metabolism",               "O-fucosyltransferase family protein",
  # --- Chloroplast / Photosynthesis ---
  "Glyma.15G250600",    "Chloroplast",              "Thylakoid formation 1 (THF1); chloroplast biogenesis",
  # --- Cell cycle ---
  "Glyma.06G062500",    "Cell cycle",               "DNA replication factor CDT1; replication licensing",
  # --- Other with annotation ---
  "Glyma.03G041000",    "Stress response",          "LURP-ONE-related; defense/pathogen responsive",
  "Glyma.05G127400",    "Cell death",               "Programmed cell death protein 5 (PDCD5)",
  "Glyma.05G181700",    "Mitochondria",             "Mitofilin (IMMT); mitochondrial inner membrane organization",
  "Glyma.06G059200",    "Mitochondria",             "MIOREX complex component 2; mitochondrial ribosome",
  "Glyma.06G059402",    "Mitochondria",             "MIOREX complex component 2; mitochondrial ribosome",
  "Glyma.06G267300",    "Signaling",                "ADP-ribosyl cyclase; cyclic ADP-ribose hydrolase, Ca2+ signaling",
  "Glyma.06G304100",    "Unknown (DUF247)",         "Plant protein of unknown function (DUF247)",
  "Glyma.06G304600",    "Unknown (DUF247)",         "DUF247 with RNA-binding; unknown function",
  "Glyma.12G230967",    "Signaling",                "NHL domain-containing protein; innate immunity",
  "Glyma.13G099400",    "Mitochondria",             "Metaxin (MTX); mitochondrial outer membrane protein import",
  "Glyma.15G051900",    "Metal homeostasis",        "HMA domain protein; heavy metal associated",
  "Glyma.20G220200",    "Protein turnover",         "F-box protein PP2-B10-related; SCF ubiquitin ligase",
  "Glyma.03G000300",    "RNA processing",           "PPR repeat protein; organellar RNA processing",
  "Glyma.03G031200",    "Metabolism",               "Cytochrome P450; oxidative metabolism",
  # --- No annotation ---
  "Glyma.02G181350",    "No annotation",            "No Arabidopsis ortholog identified",
  "Glyma.15G240132",    "No annotation",            "No Arabidopsis ortholog identified",
  "Glyma.18G245400",    "No annotation",            "No Arabidopsis ortholog identified"
)

# Merge categories
merged <- merged %>%
  left_join(category_map, by = "GeneID")

# Check for missing
missing <- merged %>% filter(is.na(Functional_Category))
if (nrow(missing) > 0) {
  cat("WARNING: Missing categories for:\n")
  print(missing$GeneID)
}

# 4. SUPPLEMENTARY TABLE: ALL 72 GENES

# Clean up Arabidopsis defline
parse_defline <- function(defline) {
  if (is.na(defline) || defline == "") return(NA_character_)
  desc <- gsub("^\\(\\d+ of \\d+\\)\\s*", "", defline)
  desc <- gsub("PTHR\\d+(:SF\\d+)?\\s*-\\s*", "", desc)
  desc <- gsub("^[\\d\\.]+(/[\\d\\\\.]+)*\\s*-\\s*", "", desc)
  desc <- trimws(desc)
  if (nchar(desc) > 100) desc <- paste0(substr(desc, 1, 97), "...")
  return(desc)
}

supp_table <- merged %>%
  mutate(
    Arabi_Short = sapply(Best_hit_arabi_defline, parse_defline)
  ) %>%
  select(
    GeneID,
    Arabidopsis_Ortholog = Best_hit_arabi_name,
    Functional_Category,
    Short_Function,
    DE_Tier,
    logFC = Mean_logFC,
    r_LW = Cor_LW_Ratio,
    Binding_Class,
    ChIP_Peaks = ChIP_N_Peaks,
    ChIP_Region = ChIP_Primary_Region,
    Final_Tier,
    WGCNA_Module = Module_Color,
    Pfam
  ) %>%
  mutate(logFC = round(as.numeric(logFC), 2),
         r_LW = round(r_LW, 2)) %>%
  arrange(Functional_Category, desc(abs(logFC)))

write_csv(supp_table, file.path(out_dir, "TableS_72Genes_AllTargets.csv"))
cat("\nSaved supplementary table:", nrow(supp_table), "genes\n")

# Print category summary
cat("\n=== FUNCTIONAL CATEGORY SUMMARY ===\n")
cat_summary <- supp_table %>%
  mutate(Broad_Category = case_when(
    grepl("Hormone", Functional_Category) ~ "Hormone signaling",
    grepl("Defense|NB", Functional_Category) ~ "Defense/NBS-LRR",
    grepl("Receptor|Signaling|Circadian", Functional_Category) ~ "Receptor kinase/Signaling",
    grepl("Transcription|Epigenetics", Functional_Category) ~ "Transcription/Chromatin",
    grepl("Transport", Functional_Category) ~ "Transport",
    grepl("Protein folding", Functional_Category) ~ "Protein folding/ER",
    grepl("Vesicle|Secretory", Functional_Category) ~ "Vesicle trafficking",
    grepl("RNA|Translation", Functional_Category) ~ "RNA processing/Translation",
    grepl("Metabolism", Functional_Category) ~ "Metabolism/Enzyme",
    grepl("Chloroplast", Functional_Category) ~ "Chloroplast",
    grepl("Mitochondria", Functional_Category) ~ "Mitochondria",
    grepl("Cell cycle", Functional_Category) ~ "Cell cycle",
    grepl("No annotation|Unknown", Functional_Category) ~ "Unknown/Unannotated",
    TRUE ~ "Other"
  )) %>%
  group_by(Broad_Category) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(desc(N))

for (i in seq_len(nrow(cat_summary))) {
  cat(sprintf("  %-30s %d (%.0f%%)\n",
              cat_summary$Broad_Category[i],
              cat_summary$N[i],
              100 * cat_summary$N[i] / nrow(supp_table)))
}

# 5. MAIN TEXT TABLE: LITERATURE-SUPPORTED GENES

# Only genes with direct literature support on Arabidopsis ortholog
lit_genes <- tribble(
  ~GeneID,            ~Key_Citation,
  "Glyma.04G180400",  "Abe et al. 2003 Plant Cell",
  "Glyma.19G167300",  "Van der Does et al. 2017 PLoS Genet",
  "Glyma.02G061400",  "Kozuka et al. 2013 Plant Cell Physiol",
  "Glyma.11G225500",  "Kudo et al. 2012 Plant Physiol",
  "Glyma.05G204800",  "Park & Kim 2021 Int J Mol Sci",
  "Glyma.06G305800",  "Stirnberg et al. 1999 Plant Physiol",
  "Glyma.04G013800",  "Gonzalez et al. 2015 Plant Cell",
  "Glyma.15G263332",  "Gao et al. 2020 Hortic Res",
  "Glyma.15G196500",  "Wu et al. 2011 PLoS One",
  "Glyma.11G168428",  "Bak & Feyereisen 2001 Plant Physiol",
  "Glyma.12G022600",  "Cui et al. 2014 Plant J",
  "Glyma.02G226700",  "Nelissen et al. 2005 PNAS",
  "Glyma.11G251900",  "Galon et al. 2010 Planta",
  "Glyma.02G189700",  "Zilberman et al. 2003 Science",
  "Glyma.14G199600",  "Hala et al. 2008 Plant Cell",
  "Glyma.13G357700",  "Lu & Christopher 2008 Mol Genet Genomics"
)

main_table <- merged %>%
  inner_join(lit_genes, by = "GeneID") %>%
  mutate(logFC = round(as.numeric(Mean_logFC), 2),
         r_LW = round(Cor_LW_Ratio, 2)) %>%
  select(
    GeneID,
    Arabidopsis_Ortholog = Best_hit_arabi_name,
    Functional_Category,
    Short_Function,
    DE_Tier,
    logFC,
    r_LW,
    Binding_Class,
    ChIP_Peaks = ChIP_N_Peaks,
    Final_Tier,
    Key_Citation
  ) %>%
  arrange(Functional_Category, desc(abs(logFC)))

write_csv(main_table, file.path(out_dir, "Table_MainText_LiteratureSupported.csv"))
cat("\nSaved main text table:", nrow(main_table), "literature-supported genes\n")

cat("\n=== MAIN TEXT TABLE PREVIEW ===\n")
for (i in seq_len(nrow(main_table))) {
  r <- main_table[i, ]
  cat(sprintf("%-18s | %s | %-22s | logFC=%6.2f | r=%.2f | %s | %s\n",
              r$GeneID, r$Arabidopsis_Ortholog,
              r$Functional_Category,
              r$logFC, r$r_LW,
              r$Final_Tier,
              r$Key_Citation))
}

# 6. CATEGORY SUMMARY TABLE

broad_summary <- supp_table %>%
  mutate(Broad_Category = case_when(
    grepl("Hormone", Functional_Category) ~ "Hormone signaling",
    grepl("Defense|NB", Functional_Category) ~ "Defense/NBS-LRR",
    grepl("Receptor|Signaling|Circadian", Functional_Category) ~ "Receptor kinase/Signaling",
    grepl("Transcription|Epigenetics", Functional_Category) ~ "Transcription/Chromatin",
    grepl("Transport", Functional_Category) ~ "Transport",
    grepl("Protein folding", Functional_Category) ~ "Protein folding/ER",
    grepl("Vesicle|Secretory", Functional_Category) ~ "Vesicle trafficking",
    grepl("RNA|Translation", Functional_Category) ~ "RNA processing/Translation",
    grepl("Metabolism", Functional_Category) ~ "Metabolism/Enzyme",
    grepl("Chloroplast", Functional_Category) ~ "Chloroplast",
    grepl("Mitochondria", Functional_Category) ~ "Mitochondria",
    grepl("Cell cycle", Functional_Category) ~ "Cell cycle",
    grepl("No annotation|Unknown", Functional_Category) ~ "Unknown/Unannotated",
    TRUE ~ "Other"
  )) %>%
  group_by(Broad_Category) %>%
  summarise(
    N = n(),
    Pct = round(100 * n() / nrow(supp_table), 1),
    Example_Genes = paste(head(GeneID, 3), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(N))

write_csv(broad_summary, file.path(out_dir, "Table_Functional_Categories_Corrected.csv"))
cat("\nSaved corrected category summary\n")

cat("\n=== COMPLETE ===\n")
