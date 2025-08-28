## Figure-to-Script Mapping: Functional genomics and tumor microenvironment analysis reveal prognostic biological subtypes in Mantle cell lymphoma

This repository contains the code used to generate figures in the manuscript.  
Each figure has its own directory with scripts, input files, and outputs.  
Below is the mapping of figures to code and data sources.

| Figure | Location in Repo | Main Script(s) | Input Data Type | Output Data Type |
|--------|-----------------|----------------|------------|--------|
| 1 | `MCL_2024/Figure1/Figure_1.R` | Figure_1.R | Filtered Maf files for Coh-1 and Coh-2| Figure_1A.pdf, Figure_1B.pdf, Figure_1C.pdf
| 2 | `MCL_2024/Figure2/Figure_2.R` | ... | ... | ... |
| 3 | `MCL_2024/Figure3/Figure_3.R` | ... | ... | ... |
| 4 | `MCL_2024/Figure4/Figure_4.R` | ... | ... | ... |
| 5 | `MCL_2024/Figure_5.R` | Tumor microenvironment (IMC processed data) | Supplementary Tables / processed IMC | `Figure5/output/Figure5.pdf` |
| 6 | `MCL_2024Figure6/Figure_6.R` | Wet-lab validation (qPCR / Western blot plots) | Supplementary Data | `Figure6/output/Figure6.pdf` |
| 7 | `MCL_2024Figure7/Figure_7.R` | Validation (e.g. CUT&RUN in MCL lines) | GEO **GSE271594**, **GSE271503** | `Figure7/output/Figure7.pdf` |
