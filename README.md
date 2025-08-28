## Functional genomics and tumor microenvironment analysis reveal prognostic biological subtypes in Mantle cell lymphoma

This repository contains the code used to generate figures in the manuscript. Each figure in this repository has its own directory (e.g., `Figure1/`, `Figure2/`, …).  
Within each figure directory:  
- `Data/` contains all input files required to generate the figure.  
- `*.R` or `*.py` scripts implement the analysis.  
- `output/` contains the resulting plots or tables (as shown in the paper).
- Below is the mapping of figures to code and data sources.

| Figure | Location in Repo | Main Script(s) | Input Data Type | Output Data Type |
|--------|-----------------|----------------|------------|--------|
| 1 | `MCL_2024/Figure1/` | Figure_1.R | Filtered Maf files for Coh-1 and Coh-2, Clinical Metadata, GRCh38 Aligned bam files available at **dbGap study accession: phs003849.v1.p1 (n=267)**, **Note:** Please use BC_23_cases.txt and BC_23_clinical.tsv for 23 external cases as sequencing data is not deposited for these cases under this accession| Fig_1a.pdf, Fig_1b.pdf, Fig_1c.pdf
| 2 | `MCL_2024/Figure2/` | Figure_2.R | Copy number segmentation files for Coh-1 and Coh-2, Clinical Metadata, RNA-Seq data (n=47) Coh-1 | Fig_2a.pdf, Fig_2b.pdf, Fig_2c.pdf, Fig_2d.pdf, Fig_2e.pdf, Fig_2f.pdf, Fig_2g.pdf|
| 3 | `MCL_2024/Figure3/` | Figure_3.R | Mutation and Copy number aberrations across 153 cases from Coh-1 (n=153), Clinical Metadata, RNA-Seq (n=47), Cibersortx Immune Deconvolution profiles (n=47) | Fig_3a.pdf, Fig_3b.pdf, Fig_3c.pdf, Fig_3d.pdf |
| 4 | `MCL_2024/Figure4/` | Figure_4.R | Mutation and Copy number aberrations across 153 cases from Coh-1 and Coh-2 (n=290), Clinical Metadata | Fig_4a.pdf, Fig_4b.pdf, Fig_4c.pdf|
| 5 | `MCL_2024/Figure5/` | Figure_5.R | Tumor microenvironment (IMC processed data) | Supplementary Tables / processed IMC | `Figure5/output/Figure5.pdf` |
| 6 | `MCL_2024/Figure6/` | Figure_6.R | Wet-lab validation (qPCR / Western blot plots) | Supplementary Data | `Figure6/output/Figure6.pdf` |
| 7 | `MCL_2024/Figure7/` | Figure_7.R | Validation (e.g. CUT&RUN in MCL lines) | GEO **GSE271594**, **GSE271503** | `Figure7/output/Figure7.pdf` |
