# Proteomics Differential Expression & Network Analysis Pipeline

An R pipeline for processing, analyzing, and visualizing label-free quantitative proteomics data. Designed for mice tear fluid proteomics (Nondiabetic vs. Diabetic). This pipeline was created for Differential Expression Analysis and STRING network visualization for the manuscript titled "..." . This pipeline was based on https://github.com/Iaguilaror/low_BMD_in_PMWMX 

# Pipeline Overview

The workflow is divided into 7 sequential modules, each handling a specific step of the proteomics analysis:

*   Module 01: Differential Expression Analysis (DE) - Calculates statistics, Fold Changes, p-values (ANOVA), and generates Volcano Plots.
*   Module 02: Heatmaps- Calculates Z-scores and generates heatmaps from DE results (UP/DOWN).
*   Module 03: Venn Diagrams & Signature Extraction- Identifies intersecting and exclusive proteins between comparisons to define specific biomarker signatures.
*   Module 04: Lollipop Plots- Visualizes exclusively altered proteins (Log2FC).
*   Module 05: Protein-Protein Interaction (PPI) Networks- Maps Fold Change data onto STRING networks using `igraph` and `ggraph`.
*   Module 06: Composite Network & Enrichment plots- Assembles PPI networks and GO/InterPro enrichment dot plots.
*   Module 07: Supplementary Tables Generator- Formats STRING enrichment outputs into  Excel workbooks.

# Prerequisites & Installation

This pipeline requires **R (>= 4.1.0)**. 
Please install the following CRAN and Bioconductor packages before running the scripts:


# CRAN Packages

```R

install.packages(c("tidyverse", "janitor", "matrixStats", "ggplot2", "cowplot", 
                   "ggsci", "ggrepel", "openxlsx", "readxl", "ggvenn", 
                   "patchwork", "igraph", "ggraph", "tidygraph"))
```
# Bioconductor Packages

```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

# Data Requirements
To run this pipeline, you will need:

-Raw Quantification Data: CSV files exported from your proteomics software (e.g., Progenesis). Ensure it contains columns for accession, peptides, fold change, and ANOVA p-values.

-Sample Sheet: A .csv containing sample metadata (e.g., column sample and column condition).

-STRING Database Exports (For Modules 5-7):

-Network edges: Exported as string_interactions.tsv.

-Enrichment files: Exported as TSV (e.g., enrichment.Process.tsv, enrichment.InterPro.tsv, enrichment.all.tsv).

# Usage
The scripts are designed to be run interactively. Execute the script in RStudio. You do not need to hardcode file paths; the pipeline uses ``choose.files()`` and ``choose.dir()`` to open native file dialogs.

Note: The file selection functions are optimized for Windows. Mac/Linux users may need to replace these functions with file.choose() or provide absolute paths.

# Important Notes & Customization
If you are adapting this pipeline for your own dataset, please review the following parameters:

Module 01: The significance thresholds are set to p-value < 0.05 and Fold Change = 1.5 (ratio_threshold <- log2(1.5)). Change these in the Define Statistical Parameters section if needed.

Module 02: Column extraction currently restricts to 3 replicates per condition (cols_cond1[1:3]). Modify this line if your N > 3.

Module 04: The Lollipop plot height is optimized (height = 9) for ~40 genes. If your signature has significantly more or fewer proteins, adjust the ggsave() height parameter to prevent overlapping text.

Module 06: Contains manual gene name mapping (e.g., Tf -> Trf) to correct specific identifier mismatches between the raw dataset and the STRING database. Users must adjust or remove this case_when mapping to fit their own datasets.

# License
--- 

# Citation
If you use this pipeline in your research, please cite the associated article:
---
