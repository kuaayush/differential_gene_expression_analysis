# Comparative Analysis of Differential Gene Expression in HNSC and CESC

This repository contains code, plots, and documentation from a self-initiated project exploring differential gene expression between **Head and Neck Squamous Cell Carcinoma (HNSC)** and **Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC)** using bulk RNA-seq data from **The Cancer Genome Atlas (TCGA)**.

## 📌 Project Overview

This project aims to identify **shared and distinct gene expression signatures** between two HPV-driven cancers — HNSC and CESC — and uncover **biological pathways and processes** that distinguish them using differential gene expression analysis (DGEA) and functional enrichment analysis.

## 🧬 Objective

* Compare bulk RNA-seq profiles of HNSC and CESC from TCGA.
* Identify differentially expressed genes (DEGs).
* Perform **Principal Component Analysis**, generate **MA plots**, **Volcano plots**, **Gene Ontology (GO)**, and **Pathway Enrichment Analysis** to interpret biological relevance.

## 🧪 Dataset

* **Source**: TCGA (downloaded via `TCGAbiolinks` in R)
* **Type**: Bulk RNA-seq
* **Samples**: Tumor tissues from HNSC and CESC patients
* **Notes**: Bulk RNA-seq provides an average expression profile across all cells in a sample, making it suitable for identifying dominant transcriptional trends.

## 🧭 Pipeline Summary

```
Data Download → Preprocessing → Normalization → PCA → Heatmap →
Differential Gene Expression Analysis (DESeq2) →
Visualization (MA & Volcano Plots) →
GO and Pathway Enrichment Analysis
```

## 🧠 Biological Insights

Even though both cancers are associated with **HPV**, they manifest in distinct anatomical regions and contexts. By comparing their transcriptomes:

* PCA revealed initial clustering by sequencing depth (PC1), corrected by normalization.
* Post-normalization PCA revealed clear separation by cancer type, indicating distinct expression profiles.
* DEGs were identified using DESeq2 and analyzed for biological significance.
* Enrichment analyses highlighted **cancer-specific pathways** and **biological processes**.

All generated plots including PCA, MA plots, volcano plots, heatmaps, GO and pathway enrichment results are provided in the [**Figures**](./figures) folder.

## 📚 References

* *Comprehensive Analysis of Differential Gene Expression to Identify Common Gene Signatures in Multiple Cancers* \[PMC]
* *Systems-level differential gene expression analysis reveals new candidate biomarkers for breast cancer* \[Nature]

## 📎 Notes

* Batch effect correction was performed to remove technical variation.
* This project was developed for learning purposes and is open to feedback and collaboration.
