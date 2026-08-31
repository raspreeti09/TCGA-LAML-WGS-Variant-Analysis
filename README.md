# TCGA-LAML WGS Variant Analysis

## Overview

This project presents a comprehensive whole-genome sequencing (WGS) variant analysis workflow for **acute myeloid leukemia (AML)** using publicly available data from **The Cancer Genome Atlas (TCGA)**.

The pipeline is designed to identify and characterize somatic genomic alterations in **TCGA-LAML** samples and investigate their potential biological and clinical relevance.

## Objectives

The major objectives of this project are to:

* Characterize the somatic mutation landscape of TCGA-LAML samples
* Identify frequently mutated genes
* Identify potential driver genes
* Analyze mutational signatures
* Investigate gene-level co-mutation patterns
* Characterize copy number variations (CNVs)
* Correlate genomic alterations with clinical characteristics
* Perform survival analysis
* Generate summary visualizations of the genomic landscape

## Analysis Workflow

The analysis is organized into a sequential R-based workflow:

```text
TCGA-LAML WGS Data
        ↓
Data Acquisition
        ↓
Data Preprocessing & Quality Control
        ↓
Mutation Landscape Analysis
        ↓
Driver Gene Identification
        ↓
Mutational Signature Analysis
        ↓
Co-mutation Analysis
        ↓
Copy Number Variation Analysis
        ↓
Clinical Correlation
        ↓
Survival Analysis
        ↓
Summary Visualization
```

## Pipeline Scripts

| Script                       | Description                                                      |
| ---------------------------- | ---------------------------------------------------------------- |
| `00_run_pipeline.R`          | Main script for running the complete analysis workflow           |
| `01_setup.R`                 | Environment setup and package configuration                      |
| `02_data_acquisition.R`      | Acquisition of TCGA-LAML genomic and clinical data               |
| `03_data_preprocessing.R`    | Data cleaning, preprocessing and preparation                     |
| `04_mutation_landscape.R`    | Analysis and visualization of the mutation landscape             |
| `05_driver_genes.R`          | Identification and characterization of potential driver genes    |
| `06_mutational_signatures.R` | Analysis of mutational signatures                                |
| `07_comutation_analysis.R`   | Analysis of co-occurring and mutually exclusive mutations        |
| `08_cnv_analysis.R`          | Copy number variation analysis                                   |
| `09_clinical_correlation.R`  | Correlation of genomic alterations with clinical characteristics |
| `10_survival_analysis.R`     | Survival analysis based on genomic and clinical variables        |
| `11_summary_visualization.R` | Generation of summary plots and visualizations                   |

## Dataset

The project focuses on **TCGA-LAML (Acute Myeloid Leukemia)** genomic data generated through The Cancer Genome Atlas project.

The analysis uses publicly available genomic and clinical information associated with TCGA-LAML samples.

## Tools and Technologies

* **R**
* **Bioconductor**
* **TCGA data resources**
* Statistical analysis and visualization packages in R
* Genomic variant analysis tools and packages

## Reproducibility

The analysis has been divided into modular R scripts to make the workflow easier to understand, reproduce, and modify.

The recommended execution order is:

```text
00 → 01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11
```

Alternatively, the complete workflow can be executed using:

```text
00_run_pipeline.R
```

## Project Structure

```text
TCGA-LAML-WGS-Variant-Analysis/
│
├── README.md
│
├── 00_run_pipeline.R
├── 01_setup.R
├── 02_data_acquisition.R
├── 03_data_preprocessing.R
├── 04_mutation_landscape.R
├── 05_driver_genes.R
├── 06_mutational_signatures.R
├── 07_comutation_analysis.R
├── 08_cnv_analysis.R
├── 09_clinical_correlation.R
├── 10_survival_analysis.R
└── 11_summary_visualization.R
```

## Expected Outputs

The workflow is intended to generate:

* Mutation frequency summaries
* Mutation landscape visualizations
* Driver gene results
* Mutational signature profiles
* Co-mutation analysis results
* CNV summaries and visualizations
* Clinical association results
* Kaplan-Meier survival plots
* Summary figures describing the genomic landscape of TCGA-LAML

## Applications

This workflow can be used as a reproducible framework for exploring the genomic characteristics of AML and investigating relationships between genomic alterations and clinical outcomes.

## Author

**Ras Preeti Sharma**

M.Sc. Systems Biology and Bioinformatics
Panjab University, Chandigarh

---

*This project is intended for research and educational purposes.*
