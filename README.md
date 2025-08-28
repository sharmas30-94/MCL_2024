## Figure-to-Script Mapping: Functional genomics and tumor microenvironment analysis reveal prognostic biological subtypes in Mantle cell lymphoma

This repository contains the code used to generate figures in the manuscript.  
Each figure has its own directory with scripts, input files, and outputs.  
Below is the mapping of figures to code and data sources.

| Figure | Location in Repo | Main Script(s) | Input Data Type | Output Data Type |
|--------|-----------------|----------------|------------|--------|
| 1 | `MCL_2024/Figure1/Figure_1.R` | Figure_1.R | Filtered Maf files for Coh-1 and Coh-2, Clinical Metadata, **dbGap study accession: phs003849.v1.p1 (n=267) **, Note: Please use BC_23_cases.txt and BC_23_clinical.tsv for 23 external cases as sequencing data is not deposited under this accession| Fig_1a.pdf, Fig_1b.pdf, Fig_1c.pdf
| 2 | `MCL_2024/Figure2/Figure_2.R` | Figure_2.R | Copy number segmentation files for Coh-1 and Coh-2, Clinical Metadata | Fig_2a.pdf, Fig_2b.pdf, Fig_2c.pdf, Fig_2d.pdf, Fig_2e.pdf, Fig_2f.pdf, Fig_2g.pdf|
| 3 | `MCL_2024/Figure3/Figure_3.R` | Figure_3.R | Tabulated Mutation and Copy numner aberrations across 152 cases | ... |
| 4 | `MCL_2024/Figure4/Figure_4.R` | ... | ... | ... |
| 5 | `MCL_2024/Figure_5.R` | Tumor microenvironment (IMC processed data) | Supplementary Tables / processed IMC | `Figure5/output/Figure5.pdf` |
| 6 | `MCL_2024Figure6/Figure_6.R` | Wet-lab validation (qPCR / Western blot plots) | Supplementary Data | `Figure6/output/Figure6.pdf` |
| 7 | `MCL_2024Figure7/Figure_7.R` | Validation (e.g. CUT&RUN in MCL lines) | GEO **GSE271594**, **GSE271503** | `Figure7/output/Figure7.pdf` |
