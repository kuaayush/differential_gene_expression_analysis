

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)

# Load The Cancer Genome Atlas (TCGA) data (example : LUAD - Lung Adenocarcinoma)
LUAD_query <- GDCquery(project = "TCGA-LUAD",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification",
                 workflow.type = "STAR - Counts")

GDCdownload(LUAD_query)


# This returns a summarizedExperiment object
TCGA_LUAD_data <- GDCprepare(LUAD_query)
saveRDS(TCGA_LUAD_data, "/Users/aayushojha/Documents/DGE/tcga/blog_data/TCGA_LUAD_SummarizedExperiment.rds")

# Read the saved RDS file
TCGA_LUAD_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/blog_data/TCGA_LUAD_SummarizedExperiment.rds")
 
#There are many metadata for each sample 
colData(TCGA_LUAD_data) %>%
  colnames() %>%
  tail()

# CIMP methylation subtypes
colData(TCGA_LUAD_data) %>%
  as.data.frame() %>%
  janitor::tabyl('paper_CIMP.methylation.signature.')


# save the raw count matrix
TCGA_LUAD_mat <- assay(TCGA_LUAD_data) 
TCGA_LUAD_mat[1:5, 1:5]
dim(TCGA_LUAD_mat)





# We will now download the TCGA_LUSC data in same way  for Lung Squamous Cell Carcinoma
# Load The Cancer Genome Atlas (TCGA) data (example : LUAD - Lung Adenocarcinoma)
LUSC_query <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")

GDCdownload(LUSC_query)

# This returns a summarizedExperiment object
TCGA_LUSC_data <- GDCprepare(LUSC_query)
saveRDS(TCGA_LUSC_data, "/Users/aayushojha/Documents/DGE/tcga/blog_data/TCGA_LUSC_SummarizedExperiment.rds")

# Read the saved RDS file
TCGA_LUSC_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/blog_data/TCGA_LUSC_SummarizedExperiment.rds")

#There are many metadata for each sample 
colData(TCGA_LUSC_data) %>%
  colnames() %>%
  tail()

# CIMP methylation subtypes
colData(TCGA_LUSC_data) %>%
  as.data.frame() %>%
  janitor::tabyl('paper_Expression.Subtype')

# save the raw count matrix
TCGA_LUSC_mat <- assay(TCGA_LUSC_data) 
TCGA_LUSC_mat[1:5, 1:5]
dim(TCGA_LUSC_mat)





# Convert ensembl id to gene symbols 
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# remove the version number in the end of the ENSEMBL id
TCGA_LUAD_genes<- rownames(TCGA_LUAD_mat) %>%
  tibble::enframe() %>%
  mutate(ENSEMBL =stringr::str_replace(value, "\\.[0-9]+", ""))

head(TCGA_LUAD_genes)


# There are duplicated gene symbols for different ENSEMBL id
# Check to see if there are any duplicated ENSEMBL
clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
  janitor::get_dupes(SYMBOL) %>%
  head()



# keep only one of the duplicated Ids 
TCGA_LUAD_gene_map<- clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL,.keep_all = TRUE)


TCGA_LUAD_gene_map<- TCGA_LUAD_gene_map %>%
  left_join(TCGA_LUAD_genes)

head(TCGA_LUAD_gene_map)



# Subset the gene expression matrices and replace the rownames to gene symbol.

TCGA_LUSC_mat <- TCGA_LUSC_mat[TCGA_LUAD_gene_map$value, ]
row.names(TCGA_LUSC_mat) <- TCGA_LUAD_gene_map$SYMBOL
TCGA_LUSC_mat[1:5, 1:5]

TCGA_LUAD_mat<- TCGA_LUAD_mat[TCGA_LUAD_gene_map$value, ]
row.names(TCGA_LUAD_mat)<- TCGA_LUAD_gene_map$SYMBOL
TCGA_LUAD_mat[1:5, 1:5]

dim(TCGA_LUSC_mat)
dim(TCGA_LUAD_mat)



# Combine TCGA LUSC and LUAD
# PCA plot use top variable genes.

# double check the genes are the same
all.equal(rownames(TCGA_LUAD_mat), rownames(TCGA_LUSC_mat))

# 
TCGA_lung_mat <- cbind(TCGA_LUSC_mat, TCGA_LUAD_mat)
TCGA_lung_meta <- data.frame(cancer_type = c(rep( "LUSC", ncol(TCGA_LUSC_mat)), 
                                             rep("LUAD", ncol(TCGA_LUAD_mat))))

dim(TCGA_lung_mat)



# PCA plot using ggplot2

library(ggplot2)
library(ggfortify)

# select the top 1000 most variable genes 
TCGA_gene_idx<- order(rowVars(TCGA_lung_mat), decreasing = TRUE)[1:1000]

TCGA_lung_mat_sub <- TCGA_lung_mat[TCGA_gene_idx, ]

TCGA_pca_res <- prcomp(t(TCGA_lung_mat_sub), scale. = TRUE)

autoplot(TCGA_pca_res, data = TCGA_lung_meta , color ="cancer_type") +
  scale_color_manual(values = c("blue", "red")) +
  ggtitle("TCGA NSCLC")


# LUAD and LUSC samples are not separated in the first PC.
# Let’s see how sequencing depth is correlated with first PC
# read https://divingintogeneticsandgenomics.com/post/pca-in-action/ for more details for PCA.

# the PCs 
TCGA_pca_res$x[1:5, 1:5]

seq_depth <- colSums(TCGA_lung_mat)

cor(TCGA_pca_res$x[, 1], seq_depth) %>%
  abs()
# A correlation of 0.98 for PC1!!

cor(TCGA_pca_res$x[, 2], seq_depth) %>%
  abs()
#PC2 is not correlated with sequencing depth. As we can see PC2 separates the cancer types



TCGA_lung_meta$seq_depth <- seq_depth

autoplot(TCGA_pca_res, data = TCGA_lung_meta, color = "seq_depth") +
  scale_color_viridis_b() +
  ggtitle("TCGA NSCLC")

# PCA after CPM (counts oere million) normalization 


#convert the raw counts to log(cmp+1) using edgeR
TCGA_lung_cpm <- edgeR::cpm(TCGA_lung_mat, log= TRUE, prior.count = 1)

# select the top 1000 most variable genes 
TCGA_gene_idx2<- order(rowVars(TCGA_lung_cpm), decreasing = TRUE)[1:1000]

TCGA_lung_cpm_sub <- TCGA_lung_cpm[TCGA_gene_idx2, ]

TCGA_pca_res2 <- prcomp(t(TCGA_lung_cpm_sub), scale. = TRUE)

autoplot(TCGA_pca_res2, data = TCGA_lung_meta,  color = "cancer_type") +
  scale_color_manual(values = c("blue", "red")) + 
  ggtitle("TCGA NSCLC")

# Now we see the first PC separates cancer type: LUSC vs LUAD, which makes sense!
# It is interesting to see some of the LUSC samples mix with the LUAD.
# It will be interesting to further investigate those samples.



set.seed(123)
library(ComplexHeatmap)

TCGA_ha <- HeatmapAnnotation(df = TCGA_lung_meta,
                             col = list(
                               cancer_type = c("LUSC" = "red", "LUAD" = "blue")
                             ))

Heatmap(t(scale(t(TCGA_lung_cpm_sub))),
        name = "TCGA lung RNA",
        show_column_names = FALSE,
        show_row_names = FALSE,
        show_row_dend = FALSE,
        # column_split = TCGA_lung_meta$cancer_type,
        top_annotation = TCGA_ha,
        row_names_gp = gpar(fontsize = 3),
        border = TRUE,
        row_km = 3
        )



  