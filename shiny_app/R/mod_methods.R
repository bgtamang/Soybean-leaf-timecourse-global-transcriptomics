# Methods module - Documentation and about

methods_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(8,
        card(
          card_header("Study Overview"),
          card_body(
            h4("RNA-seq Analysis of GmJAG1 Transcriptional Targets in Soybean Leaf Development"),
            hr(),
            h5("Background"),
            p("This study investigates the transcriptional targets of the GmJAG1 transcription ",
              "factor in soybean (Glycine max) leaf development. GmJAG1 is a C2H2 zinc finger ",
              "protein that regulates leaf shape through control of cell proliferation."),
            p("Two narrow-leaf genotypes (PI 547745 and PI 612713B) carry the D9H mutation in ",
              "the EAR motif of GmJAG1, while two broad-leaf genotypes (LD11-2170 and PI 532462A) ",
              "have functional GmJAG1."),
            hr(),
            h5("Experimental Design"),
            tags$ul(
              tags$li(tags$b("Genotypes:"), " 4 (2 broad-leaf, 2 narrow-leaf)"),
              tags$li(tags$b("Timepoints:"), " 5 (TP1-TP5, shoot apex through leaf expansion)"),
              tags$li(tags$b("Replicates:"), " 3 biological replicates per condition"),
              tags$li(tags$b("Total Samples:"), " 60"),
              tags$li(tags$b("Sequencing:"), " Illumina NovaSeq 6000, 100bp single-end")
            )
          )
        ),
        card(
          card_header("Analysis Pipeline"),
          card_body(
            h5("1. Read Quantification"),
            p("Raw reads were quantified using Salmon v1.10.0 with selective alignment against ",
              "the soybean reference transcriptome (Wm82.a6.v1)."),
            h5("2. Normalization & Batch Correction"),
            p("Transcript-level abundance estimates were summarized to gene-level counts using ",
              "tximport. Batch effects between 2021 and 2022 sequencing runs were corrected ",
              "using ComBat-seq. TMM normalization was applied using edgeR."),
            h5("3. Differential Expression"),
            p("Differential expression analysis was performed using the edgeR quasi-likelihood (QL) ",
              "framework. Genes with FDR < 0.05 and |log2FC| > 1 were considered differentially expressed."),
            h5("4. JAG1 Target Identification"),
            p("JAG1 targets were identified through a tiered approach based on consistency ",
              "across pairwise genotype comparisons at TP1:"),
            tags$ul(
              tags$li(tags$b("Gold:"), " DE in all 4 pairwise comparisons (76 genes)"),
              tags$li(tags$b("Silver:"), " DE in 3 of 4 comparisons (69 genes)"),
              tags$li(tags$b("Bronze:"), " DE in 2+ comparisons or pooled contrast (1,422 genes)")
            ),
            h5("5. Co-expression Analysis"),
            p("Weighted gene co-expression network analysis (WGCNA) identified 22 co-expression ",
              "modules (21 named + grey/unassigned) from 35,941 expressed genes. Module-trait correlation analysis linked ",
              "modules to leaf shape (broad vs. narrow) and developmental stage. JAG1 target ",
              "enrichment was assessed per module using Fisher's exact test."),
            h5("6. Functional Enrichment"),
            p("GO enrichment analysis was performed using clusterProfiler to identify ",
              "biological processes and molecular functions associated with JAG1 targets ",
              "at each confidence tier (Gold, Silver, Bronze)."),
            h5("7. Cell Cycle Machinery"),
            p("Cell cycle gene families were curated with Pfam domain validation: ",
              "101 cyclins (PF00134/PF02984), 46 CDKs (PF00069/PF07714), 9 KRP inhibitors, ",
              "12 TPR co-repressor complex genes, and 15 HDAC complex genes (Yang et al. 2018). ",
              "Expression patterns were analyzed at TP1 (shoot apex) to characterize ",
              "JAG1-mediated cell proliferation control."),
            h5("8. Hormone Pathway Analysis"),
            p("Hormone-related genes (459 total) were mapped to 8 sub-pathways using KEGG ",
              "Orthology (KO) assignments: auxin, cytokinin, gibberellin, ABA, ethylene, ",
              "brassinosteroid, jasmonic acid, and salicylic acid. Enrichment of JAG1 targets ",
              "within each pathway was tested using Fisher's exact test."),
            h5("9. Multi-Evidence Integration"),
            p("A 79-gene high-confidence set was identified by integrating DE evidence, ",
              "ChIP-seq/DAP-seq binding data (O'Malley et al. 2016), phenotype correlations, ",
              "temporal expression patterns, and WGCNA module membership. Genes required ",
              "convergent support from multiple independent evidence lines."),
            h5("10. Single-Cell Spatial Validation"),
            p("Published soybean snRNA-seq data (Fan et al. 2025, Molecular Plant) from ",
              "shoot apical meristem (9,869 cells) and leaf (10,084 cells) were used to ",
              "validate spatial co-localization of JAG1 and its targets in the Leaf ",
              "Primordium 1 (LP1) cluster.")
          )
        )
      ),
      column(4,
        card(
          card_header("Software Versions"),
          card_body(
            tags$table(class = "table table-sm",
              tags$tr(tags$td("R"), tags$td("4.5.2")),
              tags$tr(tags$td("Salmon"), tags$td("1.10.0")),
              tags$tr(tags$td("tximport"), tags$td("1.30.0")),
              tags$tr(tags$td("edgeR"), tags$td("4.6.3")),
              tags$tr(tags$td("limma"), tags$td("3.58.1")),
              tags$tr(tags$td("WGCNA"), tags$td("1.73")),
              tags$tr(tags$td("clusterProfiler"), tags$td("4.16.0")),
              tags$tr(tags$td("ComBat-seq"), tags$td("sva 3.50.0"))
            )
          )
        ),
        card(
          card_header("Reference Genome"),
          card_body(
            p(tags$b("Species:"), " Glycine max"),
            p(tags$b("Assembly:"), " Wm82.a6.v1"),
            p(tags$b("Source:"), " Phytozome v13"),
            p(tags$b("Genes:"), " 48,359 total"),
            p(tags$b("After Filtering:"), " 35,941 (CPM > 1 in 3+ samples)")
          )
        ),
        card(
          card_header("Key Thresholds"),
          card_body(
            tags$table(class = "table table-sm",
              tags$tr(tags$td("DE FDR"), tags$td("< 0.05")),
              tags$tr(tags$td("DE |log2FC|"), tags$td("> 1")),
              tags$tr(tags$td("Expression filter"), tags$td("CPM > 1 in 3+ samples")),
              tags$tr(tags$td("WGCNA minModuleSize"), tags$td("30")),
              tags$tr(tags$td("WGCNA mergeCutHeight"), tags$td("0.20"))
            )
          )
        ),
        card(
          card_header("Contact"),
          card_body(
            p(tags$b("Author:"), " Bishal Tamang"),
            p(tags$b("Institution:"), " University of Illinois Urbana-Champaign"),
            p(tags$b("Department:"), " Crop Sciences"),
            hr(),
            p(class = "small text-muted", "For questions about this dashboard or data, please contact the author.")
          )
        ),
        card(
          card_header("Data Availability"),
          card_body(
            p("Raw sequencing data and processed files are available at:"),
            p(tags$b("NCBI GEO:"), " ",
              tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317596",
                     "GSE317596", target = "_blank")),
            p(tags$b("GitHub:"), " ",
              tags$a(href = "https://github.com/bgtamang/Soybean-leaf-timecourse-global-transcriptomics",
                     "bgtamang/Soybean-leaf-timecourse-global-transcriptomics", target = "_blank"))
          )
        )
      )
    )
  )
}

methods_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # No server logic needed for static content
  })
}
