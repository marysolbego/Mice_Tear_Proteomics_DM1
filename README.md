# Proteomics Differential Expression & Network Analysis Pipeline

An R pipeline for processing, analyzing, and visualizing label-free quantitative proteomics data. Designed for mice tear fluid proteomics (Nondiabetic vs. Diabetic). This pipeline was created for Differential Expression Analysis and STRING network visualization for the manuscript titled "Tear Fluid as an Early Molecular Diagnostic Window in a Type 1 Diabetic Mouse Model" . This pipeline was based on https://github.com/Iaguilaror/low_BMD_in_PMWMX 

# Pipeline Overview

The workflow is divided into 3 main R scripts, each for an specific phase of the proteomic analysis: 

1. Differential Expression Analysis (DEA): Calculates descriptive statistics, Log2 Fold Changes, and ANOVA p-values. Generates annotated Volcano Plots and outputs results in Excel/CSV formats. Heatmaps: Plots heatmaps for significant proteins (UP/DOWN).
Venn Diagrams: Identifies intersecting and altered proteins across longitudinal comparisons to define specific biomarker signatures.
Lollipop Plots:Visualizes the altered proteins (Log2FC) of the specific diabetic signature.
Multivariate Statistics: Performs condition-specific Principal Component Analysis (PCA) and a Combined Longitudinal PCA with PERMANOVA to track disease progression trajectories.

2. PPI Network Analysis for Diabetic Signature: Integrates the diabetic biomarker signature (from the previous step) with the STRING interaction database.

3. PPI Network & Enrichment Analysis : Builds the network 1-month Nondiabetic vs. 1-month Diabetic.
Enrichment Dot Plots: Visualizes Gene Ontology (GO - Biological Process) and InterPro Functional Domains enrichment.
Supplementary Tables: Extracts the STRING enrichment TSV outputs in a table. 




# Prerequisites & Installation

This pipeline requires **R (>= 4.1.0)**. 
Please install the following CRAN and Bioconductor packages before running the scripts:


## CRAN Packages

```R

install.packages(c("tidyverse", "janitor", "matrixStats", "ggplot2", "cowplot", 
                   "ggsci", "ggrepel", "openxlsx", "readxl", "ggvenn", 
                   "patchwork", "igraph", "ggraph", "tidygraph", "vegan", "ggpubr"))
```
## Bioconductor Packages

```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

# Data Requirements
To run this pipeline, you will need:

1. Raw_data:
CSV files exported from your proteomics software (e.g., Progenesis). Ensure it contains columns for accession, peptides, fold change, and ANOVA p-values. File names must follow the convention: `ConditionA_vs_ConditionB_total_proteins.csv`.
2.Sample Sheet (Raw_data/sample_sheet.csv):
A metadata CSV containing at least two columns: muestra (sample ID matching raw data headers) and condition.
3. STRING Database Exports:
Network edges: Exported from STRING as string_interactions.tsv.
Enrichment files:Exported from STRING as TSV files (e.g., enrichment.Process.tsv, enrichment.InterPro.tsv, enrichment.all.tsv).


---

## License

---

## Citation

## If you use this pipeline in your research, please cite the following article:

