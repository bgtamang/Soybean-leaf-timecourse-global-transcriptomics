#          to provide specific binding evidence for manuscript

rm(list = ls())
gc()

cat("\n")
cat("  DETAILED PEAK EXTRACTION FOR TOP TARGETS\n")
cat("  Extracting specific binding coordinates for key genes\n")

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))

library(dplyr)
library(tidyr)

dir.create("03_results/tables/binding_integration/detailed_peaks",
           recursive = TRUE, showWarnings = FALSE)


cat("Loading data...\n")

# Full ChIP-seq peaks with coordinates
chipseq_peaks <- read.csv("01_data/published/peaks/Huang_2021_ChIPseq_peaks_full.csv",
                          stringsAsFactors = FALSE)
cat("  ChIP-seq peaks:", nrow(chipseq_peaks), "peak-gene associations\n")

# DAP-seq peaks
dapseq_peaks <- read.csv("01_data/published/peaks/Wang_2024_DAPseq_peaks_full.csv",
                         stringsAsFactors = FALSE)
cat("  DAP-seq peaks:", nrow(dapseq_peaks), "peak-gene associations\n")

# Integrated binding results
integrated <- read.csv("03_results/tables/binding_integration/integrated_binding_analysis.csv",
                       stringsAsFactors = FALSE)
cat("  Integrated targets:", nrow(integrated), "genes\n\n")
# Use ALL DE genes - no special filtering
genes_of_interest <- integrated$GeneID

cat("Total DE genes to analyze:", length(genes_of_interest), "\n")
cat("  - With ChIP-seq binding:", sum(integrated$Has_ChIPseq), "\n")
cat("  - With DAP-seq binding:", sum(integrated$Has_DAPseq), "\n")
cat("  - With both:", sum(integrated$Binding_Class == "ChIP_AND_DAP"), "\n")


cat("\n========================================\n")

# Function to extract peak details for a gene
extract_peak_details <- function(gene_id, peaks_df, source_name) {
  gene_peaks <- peaks_df %>%
    filter(GeneID == gene_id)

  if (nrow(gene_peaks) == 0) {
    # Return consistent columns even when no peaks
    return(data.frame(
      GeneID = gene_id,
      Source = source_name,
      N_Peaks = 0,
      Regions = NA_character_,
      Closest_Distance = NA_real_,
      Peak_Details = "No peaks",
      stringsAsFactors = FALSE
    ))
  }

  # Summarize peaks
  peak_summary <- gene_peaks %>%
    mutate(
      Peak_Info = paste0(Genomic_Region, " (", Distance_to_TSS, "bp from TSS)")
    ) %>%
    summarize(
      N_Peaks = n(),
      Regions = paste(unique(Genomic_Region), collapse = ", "),
      Closest_Distance = min(abs(as.numeric(Distance_to_TSS)), na.rm = TRUE),
      Peak_Details = paste(Peak_Info, collapse = "; ")
    )

  data.frame(
    GeneID = gene_id,
    Source = source_name,
    N_Peaks = peak_summary$N_Peaks,
    Regions = peak_summary$Regions,
    Closest_Distance = peak_summary$Closest_Distance,
    Peak_Details = peak_summary$Peak_Details,
    stringsAsFactors = FALSE
  )
}

# Extract ChIP-seq details for genes of interest
cat("Extracting ChIP-seq peak details...\n")
chipseq_details <- do.call(rbind, lapply(genes_of_interest, function(g) {
  extract_peak_details(g, chipseq_peaks, "ChIP-seq_Huang2021")
}))

# Extract DAP-seq details
cat("Extracting DAP-seq peak details...\n")
dapseq_details <- do.call(rbind, lapply(genes_of_interest, function(g) {
  extract_peak_details(g, dapseq_peaks, "DAP-seq_Wang2024")
}))


cat("\n========================================\n")

# Use all integrated results (all DE genes)
comprehensive <- integrated %>%
  select(GeneID, DE_Tier, Mean_logFC, Binding_Class, Final_Tier, Priority_Score,
         ChIP_Promoter, ChIP_Genic, ChIP_Primary_Region, ChIP_Closest_Distance)

# Add ChIP-seq peak details
chipseq_summary <- chipseq_details %>%
  filter(N_Peaks > 0) %>%
  select(GeneID, ChIP_N_Peaks = N_Peaks, ChIP_Regions = Regions,
         ChIP_Peak_Details = Peak_Details)

comprehensive <- left_join(comprehensive, chipseq_summary, by = "GeneID")

# Add DAP-seq peak details
dapseq_summary <- dapseq_details %>%
  filter(N_Peaks > 0) %>%
  select(GeneID, DAP_N_Peaks = N_Peaks, DAP_Peak_Details = Peak_Details)

comprehensive <- left_join(comprehensive, dapseq_summary, by = "GeneID")

# Sort by priority
comprehensive <- comprehensive %>%
  arrange(desc(Priority_Score), desc(abs(Mean_logFC)))

cat("Comprehensive table created:", nrow(comprehensive), "genes\n")


cat("\n========================================\n")

# Create a clean summary for top targets
manuscript_summary <- comprehensive %>%
  filter(Priority_Score >= 5) %>%
  mutate(
    Binding_Evidence = case_when(
      Binding_Class == "ChIP_AND_DAP" ~ paste0("ChIP-seq (", ChIP_Regions, ") + DAP-seq"),
      Binding_Class == "ChIP_Only" ~ paste0("ChIP-seq (", ChIP_Regions, ")"),
      Binding_Class == "DAP_Only" ~ "DAP-seq only",
      TRUE ~ "DE only"
    ),
    Distance_Summary = ifelse(!is.na(ChIP_Closest_Distance),
                              paste0(ChIP_Closest_Distance, " bp from TSS"),
                              "N/A")
  ) %>%
  select(GeneID, DE_Tier, logFC = Mean_logFC, Binding_Evidence,
         Binding_Location = ChIP_Primary_Region, Distance_Summary,
         N_ChIP_Peaks = ChIP_N_Peaks, Priority_Score)

cat("Top targets with binding details:\n")
print(head(manuscript_summary, 20))


cat("\n========================================\n")

# Create supplementary table with ALL DE genes and their binding details
supp_table <- integrated %>%
  left_join(
    chipseq_peaks %>%
      group_by(GeneID) %>%
      summarize(
        ChIP_Peak_Count = n(),
        ChIP_Binding_Regions = paste(unique(Genomic_Region), collapse = "; "),
        ChIP_Distances_to_TSS = paste(Distance_to_TSS, collapse = "; "),
        ChIP_Peak_IDs = paste(unique(GmJBS_ID), collapse = "; "),
        ChIP_Chromosomes = paste(unique(Chromosome), collapse = "; "),
        .groups = "drop"
      ),
    by = "GeneID"
  ) %>%
  left_join(
    dapseq_peaks %>%
      group_by(GeneID) %>%
      summarize(
        DAP_Peak_Count = n(),
        DAP_Binding_Regions = paste(unique(Genomic_Region), collapse = "; "),
        DAP_Distances_to_TSS = paste(Distance_to_TSS, collapse = "; "),
        .groups = "drop"
      ),
    by = "GeneID"
  ) %>%
  select(
    GeneID,
    DE_Tier,
    Mean_logFC,
    Binding_Class,
    Final_Tier,
    Priority_Score,
    # ChIP-seq details
    Has_ChIPseq,
    ChIP_Peak_Count,
    ChIP_Primary_Region,
    ChIP_Closest_Distance,
    ChIP_Binding_Regions,
    ChIP_Distances_to_TSS,
    ChIP_Peak_IDs,
    # DAP-seq details
    Has_DAPseq,
    DAP_Peak_Count,
    DAP_Binding_Regions,
    DAP_Distances_to_TSS
  ) %>%
  arrange(desc(Priority_Score), desc(abs(Mean_logFC)))

# Fill NAs
supp_table$ChIP_Peak_Count[is.na(supp_table$ChIP_Peak_Count)] <- 0
supp_table$DAP_Peak_Count[is.na(supp_table$DAP_Peak_Count)] <- 0

cat("Supplementary table created for ALL", nrow(supp_table), "DE genes\n")
cat("  - With ChIP-seq binding:", sum(supp_table$Has_ChIPseq), "\n")
cat("  - With DAP-seq binding:", sum(supp_table$Has_DAPseq), "\n")
cat("  - With both:", sum(supp_table$Binding_Class == "ChIP_AND_DAP"), "\n\n")

# Summary of binding locations for genes WITH ChIP-seq
cat("ChIP-seq binding location summary (for DE genes with binding):\n")
chipseq_de_genes <- supp_table %>% filter(Has_ChIPseq == TRUE)
cat("  Promoter binding:", sum(integrated$ChIP_Promoter[integrated$Has_ChIPseq], na.rm = TRUE), "\n")
cat("  Gene body binding:", sum(integrated$ChIP_Genic[integrated$Has_ChIPseq], na.rm = TRUE), "\n")
cat("  Downstream binding:", sum(integrated$ChIP_Downstream[integrated$Has_ChIPseq], na.rm = TRUE), "\n")


cat("\n========================================\n")

# Save FULL supplementary table (this is the key output!)
write.csv(supp_table,
          "03_results/tables/binding_integration/detailed_peaks/Supplementary_Table_Binding_Evidence.csv",
          row.names = FALSE)
cat("Saved: Supplementary_Table_Binding_Evidence.csv (ALL DE genes)\n")

# Save comprehensive table
write.csv(comprehensive,
          "03_results/tables/binding_integration/detailed_peaks/comprehensive_binding_details.csv",
          row.names = FALSE)
cat("Saved: comprehensive_binding_details.csv\n")

# Save manuscript summary
write.csv(manuscript_summary,
          "03_results/tables/binding_integration/detailed_peaks/manuscript_binding_summary.csv",
          row.names = FALSE)
cat("Saved: manuscript_binding_summary.csv\n")


# Save full peak details for all genes of interest
all_peak_details <- rbind(chipseq_details, dapseq_details)
write.csv(all_peak_details,
          "03_results/tables/binding_integration/detailed_peaks/all_peak_details.csv",
          row.names = FALSE)
cat("Saved: all_peak_details.csv\n")


cat("\n========================================\n")

cat("Use these as templates for describing binding evidence:\n\n")

# For top Tier1 targets
tier1_genes <- comprehensive %>% filter(Binding_Class == "ChIP_AND_DAP")

if (nrow(tier1_genes) > 0) {
  cat("TIER 1 (ChIP-seq + DAP-seq) examples:\n")
  for (i in 1:min(3, nrow(tier1_genes))) {
    g <- tier1_genes[i, ]
    cat(sprintf('  "%s showed significant differential expression (logFC = %.2f, %s tier) ',
                g$GeneID, g$Mean_logFC, g$DE_Tier))
    cat(sprintf('with GmJAG1 binding evidence from both ChIP-seq (%s, %d bp from TSS) ',
                g$ChIP_Primary_Region, g$ChIP_Closest_Distance))
    cat('and DAP-seq, indicating direct transcriptional regulation."\n\n')
  }
}

# For ChIP-seq only targets with promoter binding
promoter_targets <- comprehensive %>%
  filter(Binding_Class == "ChIP_Only", ChIP_Promoter == TRUE) %>%
  head(3)

if (nrow(promoter_targets) > 0) {
  cat("PROMOTER BINDING examples:\n")
  for (i in 1:nrow(promoter_targets)) {
    g <- promoter_targets[i, ]
    cat(sprintf('  "%s (logFC = %.2f) has a GmJAG1 ChIP-seq peak %d bp upstream of the TSS, ',
                g$GeneID, g$Mean_logFC, g$ChIP_Closest_Distance))
    cat('consistent with promoter-mediated regulation."\n\n')
  }
}

# For gene body binding
genic_targets <- comprehensive %>%
  filter(ChIP_Primary_Region == "Genic_Only") %>%
  head(2)

if (nrow(genic_targets) > 0) {
  cat("GENE BODY BINDING examples:\n")
  for (i in 1:nrow(genic_targets)) {
    g <- genic_targets[i, ]
    cat(sprintf('  "%s shows GmJAG1 binding within the gene body (%d bp from TSS), ',
                g$GeneID, g$ChIP_Closest_Distance))
    cat('suggesting potential regulation through intragenic binding."\n\n')
  }
}


cat("\n========================================\n")

cat("Output files in: 03_results/tables/binding_integration/detailed_peaks/\n")
cat("  1. Supplementary_Table_Binding_Evidence.csv - ALL DE genes with binding details\n")
cat("  2. comprehensive_binding_details.csv - Full details for all DE genes\n")
cat("  3. manuscript_binding_summary.csv - Clean summary for manuscript\n")
cat("  4. all_peak_details.csv - Raw peak data for all genes\n\n")
