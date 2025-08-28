## Figure-to-Script Mapping: Functional genomics and tumor microenvironment analysis reveal prognostic biological subtypes in Mantle cell lymphoma

This repository contains the code used to generate figures in the manuscript.  
Each figure has its own directory with scripts, input files, and outputs.  
Below is the mapping of figures to code and data sources.

| Figure | Location in Repo | Main Script(s) | Input Data | Output |
|--------|-----------------|----------------|------------|--------|
| 1 | `Figure1/Figure_1.R` | R script performing [genomics analysis step: e.g. differential expression, mutation frequency, CNA plot] | dbGaP **phs003849.v1.p1** BAMs (Coh-1, Coh-2); processed RNA-seq counts (**GSE271664**) | `Figure1/output/Figure1.pdf` (reproduces panel from the paper) |
| 2 | `Figure2/Figure_2.R` | ... | ... | ... |
| 3 | `Figure3/Figure_3.R` | ... | ... | ... |
| 4 | `Figure4/Figure_4.R` | ... | ... | ... |
| 5 | `Figure5/Figure_5.R` | Tumor microenvironment (IMC processed data) | Supplementary Tables / processed IMC | `Figure5/output/Figure5.pdf` |
| 6 | `Figure6/Figure_6.R` | Wet-lab validation (qPCR / Western blot plots) | Supplementary Data | `Figure6/output/Figure6.pdf` |
| 7 | `Figure7/Figure_7.R` | Validation (e.g. CUT&RUN in MCL lines) | GEO **GSE271594**, **GSE271503** | `Figure7/output/Figure7.pdf` |
