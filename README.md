## Functional genomics and tumor microenvironment analysis reveal prognostic biological subtypes in Mantle cell lymphoma

This repository contains the code for all manuscript figures, organized by directory (Figure1/, Figure2/, …). Figures 1–4, 6, and 7 were generated using a combination of Linux tools and R, while Figure 5 (tumor microenvironment analysis) was produced using R and Python.

Within each figure directory:  
- `Data/` contains all input files required to generate the figure.  
- `*.R` or `*.ipynb` scripts implement the analysis.  
- Below is the mapping of figures to code and data sources.

| Figure | Location in Repo | Main Script(s) | Input Data Type | Output Data Type |
|--------|-----------------|----------------|------------|--------|
| 1 | `MCL_2024/Figure1/` | Figure_1.R | Filtered Maf files for Coh-1 and Coh-2, Clinical Metadata, GRCh38 Aligned bam files available at **dbGap study accession: phs003849.v1.p1 (n=267)**, **Note:** Please use BC_23_cases.txt and BC_23_clinical.tsv for 23 external cases as sequencing data is not deposited for these cases under this accession| Fig_1a.pdf, Fig_1b.pdf, Fig_1c.pdf
| 2 | `MCL_2024/Figure2/` | Figure_2.R | Copy number segmentation files for Coh-1 and Coh-2, Clinical Metadata, RNA-Seq data (n=47) Coh-1 | Fig_2a.pdf, Fig_2b.pdf, Fig_2c.pdf, Fig_2d.pdf, Fig_2e.pdf, Fig_2f.pdf, Fig_2g.pdf|
| 3 | `MCL_2024/Figure3/` | Figure_3.R | Mutation and Copy number aberrations across 153 cases from Coh-1 (n=153), Clinical Metadata, RNA-Seq (n=47), Cibersortx Immune Deconvolution profiles (n=47) | Fig_3a.pdf, Fig_3b.pdf, Fig_3c.pdf, Fig_3d.pdf |
| 4 | `MCL_2024/Figure4/` | Figure_4.R | Mutation and Copy number aberrations across 153 cases from Coh-1 and Coh-2 (n=290), Clinical Metadata | Fig_4a.pdf, Fig_4b.pdf, Fig_4c.pdf|
| 5 | `MCL_2024/Figure5/` | TME_analysis.ipynb | Tumor microenvironment (IMC processed data) | Fig_5a.pdf, Fig_5b.pdf, Fig_5c.pdf, Fig_5d.pdf, Fig_5e.pdf, Fig_5f.pdf, Fig_5g.pdf|
| 6 | `MCL_2024/Figure6/` | Figure_6.R | CUT&RUN analysis, Differential gene and binding analysis using DESeq2, Functional annotation of Differential genes and promoters | Fig_6b.pdf, Fig_6c.pdf, Fig_6d.pdf, Fig_6e.pdf, Fig_6f.pdf, Fig_6g.pdf |
| 7 | `MCL_2024/Figure7/` | Figure_7.R | CUT&RUN bigWigs – normalized coverage tracks (RPKM/CPM) representing histone modification or transcription factor binding signals across the genome | Fig_7b.pdf|



## Data Paths for users
To run these scripts on your own system, you must edit the data_dir variable at the top of each script to point to the location where you have downloaded the corresponding input data.
- User: edit this path to where you downloaded the data
- data_dir <- "Figure1/Data"   # default (relative path inside this repo)
- Example (if data are stored elsewhere):
- data_dir <- "/path/to/my/local/data"
- Load file
data <- read.csv(file.path(data_dir, "input.csv"))
