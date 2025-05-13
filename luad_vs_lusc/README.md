# Differential Gene Expression Analysis: LUAD vs LUSC (TCGA)

This is my redo of the tutorial at [Diving into Genetics and Genomics](https://divingintogeneticsandgenomics.com/post/pca-tcga/)

This repository contains an independent reimplementation of the PCA and differential gene expression analysis workflow originally described in the blog post [**"PCA on TCGA data"**](https://divingintogeneticsandgenomics.com/post/pca-tcga/), adapted to compare **Lung Adenocarcinoma (LUAD)** and **Lung Squamous Cell Carcinoma (LUSC)** using **TCGA bulk RNA-seq** data.

---

## 🧬 Overview

The goal of this analysis is to explore transcriptomic differences between LUAD and LUSC, two major subtypes of non-small cell lung cancer (NSCLC), and to identify genes that are differentially expressed between them.

---

## 📌 Main Steps

1. **Data Retrieval**

   * Bulk RNA-seq HTSeq-counts and clinical metadata downloaded using the `TCGAbiolinks` package.

2. **Data Preprocessing**

   * Filtering low-count genes.
   * Annotating gene symbols using `org.Hs.eg.db`.

3. **Normalization and Transformation**

   * Normalized count data using DESeq2.
   * Variance stabilizing transformation (VST) for PCA.

4. **Principal Component Analysis (PCA)**

   * Initial PCA revealed batch effects related to sequencing depth.
   * Normalization corrected technical biases.
   * Final PCA successfully separated LUAD and LUSC along principal components.

5. **Differential Gene Expression Analysis (DGEA)**

   * DESeq2 used to identify genes significantly upregulated or downregulated in LUAD vs. LUSC.
   * Volcano and MA plots generated to visualize significant DEGs.

---

## 📊 Results Snapshot

* PCA shows clear separation between LUAD and LUSC samples.
* Top DEGs identified based on log2 fold change and adjusted p-value.
* MA and Volcano plots provide a global view of differential expression.

---

## 🛠️ Tools Used

* R (v4.x)
* [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* `ggplot2`, `EnhancedVolcano`, `org.Hs.eg.db`, `dplyr`

---

## 📌 Notes

This repository is intended for educational and exploratory purposes. It replicates the original PCA tutorial and extends it to differential gene expression analysis using best practices.

---

## 📬 Contact

Feel free to raise issues or suggestions [here](https://github.com/kuaayush/differential_gene_expression_analysis/issues) or reach out directly.
