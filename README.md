# BIOCLOCK

This repository contains the code used for the **BIOCLOCK** pilot study. The manuscript is currently being prepared for submission and will be available soon.

## Scripts

This repository contains three main scripts:
- 1_preprocessing.R: This script performs quality control on the DNA methylation data, applies the 101 preprocessing pipelines, and estimates the epigenetic ages and blood leukocyte composition.
- 2_analysis.R: This script contains the statistical analyses perfomed on the epigenetic ages and other phenotypic variables.
- 3_figures.R: This script contains the code used to generate the figures in the manuscript.

It also contains three addtitional scripts:
- function_runPipelines: This script is used to run the 101 preprocessing pipelines and is adapted from this repository: https://github.com/anilpsori/_pipelines_and_biomarkers/tree/main.
- Biolearn_EpiAge.py: This script is used to calculate GrimAge epigenetic age using the Biolearn python package.
- runcalcPCClock.R: This script is used to calculate epigenetic age according to the principal component clocks and is adapted from this repository: https://github.com/MorganLevineLab/PC-Clocks.

## Packages

The following package and software versions were used:

## Directory structure

The scripts rely on the following directory structure:

```plaintext
.
├── Data
│   ├── Infinium
│   ├── Objects
│   └── Pheno
├── Images
│   ├── Paper
│   ├── PrePro
│   └── QC
├── PC_Clocks
├── Results
│   ├── Analysis
│   └── PrePro
└── Scripts
    └── Additional
```
