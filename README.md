# Single Cell Analyses

This repository showcases several scripts for single-cell RNA-seq analyses using the Seurat package and other R libraries.

## Table of Contents

1. [Introduction](#introduction)
2. [Scripts](#scripts)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Contact](#contact)

## Introduction

This project demonstrates a typical workflow for single-cell RNA-seq and spatial transcriptomics data analysis, including data loading, quality control, data integration, clustering, visualization, gene ontology analysis, trajectory analysis, and reference mapping for cell type annotation.

## Scripts
- `ReferenceMappingAnlaysis.R` : Script for using Seurat to annotate cell labels by mapping query datasets to reference datasets.
- `SpatialAnalysis.R` : Script for using Seurat to analyze spatial transcriptomic datasets.
- `GeneOntologyAnalysis.R` : Script to analyze gene ontology and identify enriched biological processes.
- `TrajectoryAnalaysis.R` : Script to  infer cell lineage relationships and developmental pathways.

## Installation

To run the scripts in this repository, you need to have R and the required packages installed.

### Prerequisites

Make sure you have R installed. You can download it from [CRAN](https://cran.r-project.org/).

### Install Required Packages

Install the required R packages by running:

```r
install.packages(c(
  "Seurat",
  "data.table",
  "sctransform",
  "ggplot2",
  "dplyr",
  "tidyr",
  "dittoSeq",
  "presto",
  "tictoc",
  "ComplexHeatmap",
  "stringr",
  "tidyverse",
  "RColorBrewer",
  "monocle3",
  "clusterProfiler",
  "AnnotationDbi",
  "DOSE",
  "cowplot"
))

# Bioconductor packages need to be installed separately
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Mm.eg.db"))

