# GmJAG1 Transcriptional Target Analysis Pipeline

Analysis scripts for identifying and characterizing transcriptional targets of GmJAG1 in soybean using comparative transcriptomics across genotypes with contrasting leaf morphology.

## Citation

Tamang BG, Kramer C, Ainsworth EA. Natural leaf shape variation reveals diverse transcriptional targets of GmJAG1 during soybean leaf development. *XX* (2026).

## Data Availability

### RNA-seq data
Raw FASTQ files and processed count matrices: [NCBI GEO GSE317596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317596)

### Reference genome
Glycine max Williams 82 genome assembly Gmax_880_v6.0 and annotation Wm82.a6.v1 from [Phytozome](https://phytozome-next.jgi.doe.gov/) (Espina et al., 2024). Required files:
- `Gmax_880_v6.0.fa.gz` (genome assembly)
- `Gmax_880_Wm82.a6.v1.transcript.fa.gz` (transcript sequences, for Salmon index)
- `Gmax_880_Wm82.a6.v1.P14.annotation_info.txt` (gene annotations, Arabidopsis orthologs, Pfam domains)

### Published binding data (for 08_binding_integration)
- **Huang et al. 2021** ChIP-seq supplementary Table S2: `1-s2.0-S0888754320320723-mmc2.xlsx` from [Huang et al. (2021) *Genomics* 113:1304-1316](https://doi.org/10.1016/j.ygeno.2020.12.033)
- **Wang et al. 2024** DAP-seq supplementary data: `plants-3024718-supplementary.xlsx` from [Wang et al. (2024) *Plants* 13:1024](https://doi.org/10.3390/plants13071024)

### Published snRNA-seq data (for 11_single_cell_integration)
- **Fan et al. 2025** soybean shoot apex and leaf single-nucleus RNA-seq: Pre-processed h5ad files (`SAM.scRNA.h5ad`, `Leaf.scRNA.h5ad`) from [SoyOmics](https://ngdc.cncb.ac.cn/soyomics). Reference: [Fan et al. (2025) *Molecular Plant* 18:669-689](https://doi.org/10.1016/j.molp.2025.02.003)
- **SoyBase pangene table**: `Glycine.pan5.MKRS.table_ref_lines.tsv.gz` from [SoyBase](https://www.soybase.org/) (downloaded automatically by 11a script)

## Setup

All R scripts define a `base_dir` variable at the top. Before running, update this path in each script to point to your local project directory:

```r
base_dir <- "/path/to/your/project"
```

Shell scripts in `01_preprocessing/` contain SLURM directives for HPC execution. Update paths and SLURM parameters for your computing environment.

Python scripts in `11_single_cell_integration/` define a `BASE_DIR` variable similarly. Update this and ensure the conda environment has the required packages (see below).

Place downloaded data as follows:
```
<base_dir>/
├── data/genome/                       # Phytozome genome and transcript FASTA files
├── data/fastq/                        # Raw FASTQ files from GEO
├── Manuscript/Literature/
│   ├── Huang et al 2021 Chip Seq/     # Huang ChIP-seq supplementary Excel
│   └── Wang et al DAP-Seq/            # Wang DAP-seq supplementary Excel
├── single_cell_integration/data/      # snRNA-seq h5ad files and pangene table
```

## Software Requirements

- **R** v4.5.2
- **Salmon** v1.10.0 (read quantification)
- **Python** 3.10 (single-cell integration step 11, phylogenetics step 12)
- **MAFFT** v7.520 (multiple sequence alignment, step 12)
- **trimAl** v1.4.1 (alignment trimming, step 12)
- **IQ-TREE2** v2.2.2.6 (phylogenetic inference, step 12)

### Python Packages (for 11_single_cell_integration and 12_cyclin_phylogenetics)

| Package | Purpose |
|---------|---------|
| scanpy | Single-cell data I/O and analysis |
| pandas | Data manipulation |
| numpy | Array operations |
| scipy | Statistical tests (Mann-Whitney, Spearman, Pearson) |

### R Packages

| Package | Version | Purpose |
|---------|---------|---------|
| limma | 3.64.3 | Differential expression |
| edgeR | 4.6.3 | TMM normalization, filtering |
| tximport | 1.36.1 | Salmon output import |
| sva (ComBat-seq) | 3.56.0 | Batch correction |
| WGCNA | 1.73 | Co-expression network |
| clusterProfiler | 4.16.0 | GO enrichment |
| enrichplot | — | GO enrichment visualization |
| AnnotationDbi | — | Annotation database interface |
| GO.db | — | Gene Ontology database |
| vegan | 2.7-2 | PERMANOVA |
| ComplexHeatmap | 2.24.1 | Heatmap visualization |
| pheatmap | 1.0.13 | Heatmap visualization |
| ggplot2 | 4.0.1 | Plotting |
| patchwork | — | Plot composition |
| UpSetR | — | UpSet plot visualization |
| GenomicRanges | — | Genomic interval operations |
| Biostrings | — | Sequence handling |
| rtracklayer | — | Genome annotation I/O |
| readxl | — | Excel file import |
| tidyverse | — | Data manipulation |
| scales | — | Axis scaling |
| RColorBrewer | — | Color palettes |
| jsonlite | — | JSON I/O |
| ggtree | — | Phylogenetic tree visualization |

## Directory Structure

Scripts are organized in execution order:

```
01_preprocessing/          Salmon indexing, read mapping, quantification
02_data_processing/        Data import, metadata, QC, batch correction, normalization
03_exploratory_analysis/   PCA, variance partitioning, expression overview
04_differential_expression/ DE analysis (limma-voom), tiered target classification,
                           temporal dynamics
05_cell_cycle_analysis/    KRP, cyclin, and CDK gene family analyses
06_functional_enrichment/  GO enrichment, KEGG pathways, hormone pathway enrichment,
                           functional integration across evidence types
07_coexpression_network/   WGCNA network construction, module-trait correlation,
                           JAG1 target module enrichment
08_binding_integration/    ChIP-seq/DAP-seq peak extraction, promoter motif analysis,
                           binding evidence integration
09_multi_evidence_integration/ Multi-layer evidence filtering, phenotype correlation,
                           cross-validation
10_validated_genes_analysis/ Functional categorization, visualization, and heatmap
                           of 79 high-confidence targets (Tier 1: 76, Tier 2: 69,
                           Tier 3: 1,422; total: 1,567 tiered targets)
11_single_cell_integration/ Gene ID conversion (Wm82↔ZH13), hypothesis-driven
                           pathway validation and data-driven target mapping
                           using published snRNA-seq (Fan et al. 2025)
12_cyclin_phylogenetics/   D-type cyclin sequence extraction, MAFFT alignment,
                           IQ-TREE2 ML phylogeny, ggtree visualization
shiny_app/                 Interactive Shiny dashboard for exploring results
```

## Execution Order

Folders are numbered in the order they should be run. Within each folder, scripts are numbered sequentially (e.g., 02a, 02b, 02c...). Key dependencies:

1. **01_preprocessing** requires raw FASTQ files and Salmon installation
2. **02-03** process and explore quantified data
3. **04** produces DE results used by all downstream analyses
4. **05-08** can be run independently after step 4
5. **06f** (functional integration) optionally incorporates WGCNA results from step 7
6. **09** integrates results from steps 4, 6f, 7, and 8
7. **10** uses the 79 high-confidence targets from step 9
8. **11** uses target gene lists from step 4 and high-confidence targets from step 9; requires Python 3.10+ with scanpy
9. **12** requires MAFFT, trimAl, and IQ-TREE2; can be run independently after step 5

## Experimental Design

Four soybean genotypes (2 narrow-leaf with D9H mutation in GmJAG1, 2 broad-leaf with functional GmJAG1) sampled at 5 developmental timepoints (TP1-TP5, meristem to expanding leaf), 3 biological replicates each, 60 samples total.

## Interactive Dashboard

The `shiny_app/` folder contains an interactive Shiny application for exploring the RNA-seq results, including gene expression profiles, differential expression, GO enrichment, WGCNA modules, cell cycle analysis, and phenotype correlations.

**Live dashboard:** [https://bgtamang.shinyapps.io/GmJAG1-Dashboard/](https://bgtamang.shinyapps.io/GmJAG1-Dashboard/)

To run locally:

```r
shiny::runApp("shiny_app")
```

## License

MIT License
