 

# After following the tutorial by Dr. Tommy Tang to compare transcriptome profiles of LUAD and LUSC, 
# now I am applying the codes to compare the transcriptome 
# profiles of Head and Neck Cancer vs Cervical Cancer. 

# To compare Head and Neck Cancer (HNSC) vs. Cervical Cancer (CESC), I will be retrieving gene expression data from 
# TCGA for both cancer types.


library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)


# Step 1: Query and Download Raw Counts for HNSC & CESC

# Load The Cancer Genome Atlas (TCGA) data for Head and Neck Squamous Cell Carcinoma (HNSC)
HNSC_query <- GDCquery(project = "TCGA-HNSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(HNSC_query)
HNSC_data <- GDCprepare(HNSC_query)
# This returns a summarizedExperiment object
saveRDS(HNSC_data, "/Users/aayushojha/Documents/DGE/tcga/blog_data2/HNSC_SummarizedExperiment.rds")
# Read the saved RDS file
HNSC_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/blog_data2/HNSC_SummarizedExperiment.rds")


# Load The Cancer Genome Atlas (TCGA) data for Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC)
CESC_query <- GDCquery(project = "TCGA-CESC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(CESC_query)
CESC_data <- GDCprepare(CESC_query)
#This returns a summarizedExperiment object
saveRDS(CESC_data, "/Users/aayushojha/Documents/DGE/tcga/blog_data2/CESC_SummarizedExperiment.rds")
# Read the saved RDS file
CESC_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/blog_data2/CESC_SummarizedExperiment.rds")



# #There are many metadata for each sample 
# colData(HNSC_data) %>%
#   colnames() %>%
#   tail()
# 
# # CIMP methylation subtypes
# colData(HNSC_data) %>%
#   as.data.frame() %>%
#   janitor::tabyl('paper_Methylation')





# Step 2: Extract and Preprocess Raw Counts

# Extract raw count matrices
HNSC_mat <- assay(HNSC_data)
HNSC_mat[1:5, 1:5]
dim(HNSC_mat)

CESC_mat <- assay(CESC_data)
CESC_mat[1:5, 1:5]
dim(CESC_mat)

# # Add labels
# HNSC_counts$CancerType <- "HNSC"
# CESC_counts$CancerType <- "CESC"
# # Combine datasets
# combined_counts <- rbind(HNSC_counts, CESC_counts)


# Convert ensembl id to gene symbols 
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# remove the version number in the end of the ENSEMBL id
HNSC_genes<- rownames(HNSC_mat) %>%
  tibble::enframe() %>%
  mutate(ENSEMBL =stringr::str_replace(value, "\\.[0-9]+", ""))
head(HNSC_genes)


# There are duplicated gene symbols for different ENSEMBL id
# Check to see if there are any duplicated ENSEMBL
clusterProfiler::bitr(HNSC_genes$ENSEMBL,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
  janitor::get_dupes(SYMBOL) %>%
  head()


# keep only one of the duplicated Ids 
HNSC_gene_map<- clusterProfiler::bitr(HNSC_genes$ENSEMBL,
                                           fromType = "ENSEMBL",
                                           toType = "SYMBOL",
                                           OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL,.keep_all = TRUE)

HNSC_gene_map<- HNSC_gene_map %>%
  left_join(HNSC_genes)

head(HNSC_gene_map)



# Subset the gene expression matrices and replace the rownames to gene symbol.

CESC_mat <- CESC_mat[HNSC_gene_map$value, ]
row.names(CESC_mat) <- HNSC_gene_map$SYMBOL
CESC_mat[1:5, 1:5]

HNSC_mat<- HNSC_mat[HNSC_gene_map$value, ]
row.names(HNSC_mat)<- HNSC_gene_map$SYMBOL
HNSC_mat[1:5, 1:5]

dim(CESC_mat)
dim(HNSC_mat)



# Combine HNSC and CESC carcinoma. 
# PCA plot use top variable genes.

# double check the genes are the same
all.equal(rownames(HNSC_mat), rownames(CESC_mat))


hpv_mat <- cbind(CESC_mat, HNSC_mat)
hpv_meta <- data.frame(cancer_type = c(rep( "CESC", ncol(CESC_mat)), 
                                             rep("HNSC", ncol(HNSC_mat))))

dim(hpv_mat)


# PCA plot using ggplot2

library(ggplot2)
library(ggfortify)

# select the top 1000 most variable genes 
hpv_gene_idx<- order(rowVars(hpv_mat), decreasing = TRUE)[1:1000]

hpv_mat_sub <- hpv_mat[hpv_gene_idx, ]

hpv_pca_res <- prcomp(t(hpv_mat_sub), scale. = TRUE)

autoplot(hpv_pca_res, data = hpv_meta , color ="cancer_type") +
  scale_color_manual(values = c("blue", "red")) +
  ggtitle("HNSC vs CESC")


# HNSC and CESC samples are not separated in the first PC.
# Let’s see how sequencing depth is correlated with first PC
# read https://divingintogeneticsandgenomics.com/post/pca-in-action/ for more details for PCA.

# the PCs 
hpv_pca_res$x[1:5, 1:5]

seq_depth <- colSums(hpv_mat)

cor(hpv_pca_res$x[, 1], seq_depth) %>%
  abs()
# A correlation of 0.92 for PC1!!

cor(hpv_pca_res$x[, 2], seq_depth) %>%
  abs()
#PC2 is not correlated with sequencing depth. As we can see PC2 separates the cancer types




hpv_meta$seq_depth <- seq_depth

autoplot(hpv_pca_res, data = hpv_meta, color = "seq_depth") +
  scale_color_viridis_b() +
  ggtitle("TCGA HNSC VS CESC")



# STEP 3: NORMALIZE THE COUNTS 
# PCA after CPM (counts oere million) normalization 

#convert the raw counts to log(cpm+1) using edgeR
hpv_cpm <- edgeR::cpm(hpv_mat, log= TRUE, prior.count = 1)

# select the top 1000 most variable genes 
hpv_gene_idx2<- order(rowVars(hpv_cpm), decreasing = TRUE)[1:1000]

hpv_cpm_sub <- hpv_cpm[hpv_gene_idx2, ]

hpv_pca_res2 <- prcomp(t(hpv_cpm_sub), scale. = TRUE)

autoplot(hpv_pca_res2, data = hpv_meta,  color = "cancer_type") +
  scale_color_manual(values = c("blue", "red")) + 
  ggtitle("TCGA HNSC VS CESC")










