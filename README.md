# DNA_meth_illumina_pipeline

## Introduction


Welcome to the DNA Methylation Analysis Pipeline for Infinium Data repository! This pipeline is designed to facilitate the analysis of DNA methylation data generated from the Illumina HumanMethylation450 BeadChip platform. The Infinium assay is widely used for studying DNA methylation patterns at a genome-wide scale, providing valuable insights into various biological processes and disease mechanisms.

This repository aims to provide a comprehensive and user-friendly pipeline for processing, analyzing, and interpreting DNA methylation data obtained from the Infinium platform.

### Workflow

* Quality control: Perform data quality assessment to identify any potential issues or biases in the raw data.
* Preprocessing: Apply normalization techniques to minimize technical variations and batch effects within and across samples.
* Differential methylation analysis: Identify differentially methylated regions (DMRs) or individual CpG sites associated with various conditions or phenotypes of interest.
Functional interpretation: Perform functional enrichment analysis to gain insights into the biological processes and pathways related to the identified DMRs.
Visualization: Generate informative plots and figures to visualize the DNA methylation patterns and differential methylation results.

## Installation and Dependencies
To use this pipeline, you will need to have the following software and packages installed:

```r


# BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

# minfi
BiocManager::install("minfi")

# Illumina manifest
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# factoextra
install.packages("factoextra")

# qqman
install.packages("qqman")

# gplots
install.packages("gplots")

# future.apply (optional to parallelize)
install.packages("future.apply")

```
