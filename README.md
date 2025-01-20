# BIOCLOCK

This repository contains the code used for the **BIOCLOCK** pilot study. The manuscript is currently being prepared for submission and will be available soon.

## Scripts

This repository contains three main scripts:
- [1_preprocessing.R](./Scripts/1_preprocessing.R): This script performs quality control on the DNA methylation data, applies the 101 preprocessing pipelines, and estimates the epigenetic ages and blood leukocyte composition.
- [2_analysis.R](./Scripts/2_analysis.R): This script contains the statistical analyses perfomed on the epigenetic ages and other phenotypic variables.
- [3_figures.R](./Scripts/3_figures.R): This script contains the code used to generate the figures in the manuscript.

It also contains four addtitional scripts:
- [function_runPipelines.R](./Scripts/Additional/function_runPipelines.R): This script is used to run the 101 preprocessing pipelines and is adapted from this repository: https://github.com/anilpsori/_pipelines_and_biomarkers/tree/main.
- [Biolearn_EpiAge.py](./Scripts/Additional/Biolearn_EpiAge.py): This script is used to calculate the Horvath1, Horvath2, GrimAge1 and GrimAge2 epigenetic ages using the Biolearn python package (https://github.com/bio-learn/biolearn).
- [runcalcPCClock.R](./Scripts/Additional/runcalcPCClock.R): This script is used to calculate epigenetic age according to the principal component clocks (PCHorvath1, PCHorvath2, PCGrimAge) and is adapted from this repository: https://github.com/MorganLevineLab/PC-Clocks.
- [plottingFunctions.R](./Scripts/Additional/plottingFunctions.R): This script contains additional plotting fuctions used in [3_figures.R](./Scripts/3_figures.R).

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
├── Results
│   ├── Analysis
│   └── PrePro
└── Scripts
    └── Additional
        └── PC_Clocks
```

## Environment Information

Platform: Linux-5.15.0-76-generic-x86_64-with-glibc2.35 

OS: Linux (Ubuntu 22.04.5 LTS)

Machine: x86_64 x86_64-conda-linux-gnu (64-bit)

For reproducibility, the mamba/conda environments (one for R and one for Python) are available in the folder [Environments](./Environments). There is both a completely resolved environment and a minimal environment for each. In case of system incompatibilities when using the complete environment, use the minimal ones instead.

The following GitHub packages are installed in the R session using Bioconductor as they are not available in the conda package repository:
 - illuminaHumanMethylationEPICv2manifest (GitHub release 0.99.1)
 - illuminaHumanMethylationEPICv2anno.20a1.hg38 (GitHub release 0.99.0)

The PC-CLocks repository (https://github.com/MorganLevineLab/PC-Clocks) will need to be cloned into the Scripts/Additional folder and named "PC_Clocks" in order to calculate the PC clocks.

## DOI

Zenodo DOI for this repository:

DOI: 10.5281/zenodo.14674514
