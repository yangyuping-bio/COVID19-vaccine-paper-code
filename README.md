# COVID19-vaccine-paper-code

This repository contains the full analysis code used in our study on immune responses to the **inactivated SARS-CoV-2 vaccine (BBIBP-CorV)**, including multi-omics data analysis such as single-cell transcriptomics, antibody titers, metabolomics, and immune receptor repertoire.

## 📁 Modules Overview

1. **1_scRNA_analysis/**  
   Single-cell RNA-seq data preprocessing, clustering, annotation, and downstream analysis.

2. **2_antibody_analysis/**  
   Antibody titer dynamics analysis and visualization across multiple time points.

3. **3_metabolomics_analysis/**  
   Preprocessing and statistical analysis of metabolomics data across vaccination stages.

4. **4_TCR_analysis/**  
   T-cell receptor (TCR) repertoire profiling and clonal expansion analysis.

5. **5_BCR_analysis/**  
   B-cell receptor (BCR) repertoire profiling, including VDJ usage and clonal tracking.

6. **6_vaccine_vs_infection/**  
   Comparative analysis between vaccine-induced responses and natural infection.

## 🧪 Data

The code assumes input data structured according to the folder-level README in each module. If you're interested in accessing the raw data, please contact the corresponding author or refer to the related publication.

## 🛠️ Environment

This project was primarily implemented in **R**, with dependencies including `Seurat`, `tidyverse`, `GSEABase`, `ggplot2`, and others as listed in the scripts. Python-based scripts are also used where necessary.

## 📄 License

MIT License.

---

If you use this codebase, please cite our paper (coming soon).
