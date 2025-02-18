

library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)
library(tibble)
library(dplyr)


# Step 1 : Prepare data for DESeq2
    # 1 a.  Read the saved RDS file 
            HNSC_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/blog_data2/HNSC_SummarizedExperiment.rds")
            CESC_data <- readRDS("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/blog_data2/CESC_SummarizedExperiment.rds")

    # 1 b. Create expression matrix for each cancer type that contains raw counts
            # Here, raw counts refers to th enumber of reads mapped to each gene in a given sample.  
              HNSC_mat <- assay(HNSC_data)
              HNSC_mat[1:5, 1:5]
              dim(HNSC_mat)
              
              CESC_mat <- assay(CESC_data)
              CESC_mat[1:5, 1:5]
              dim(CESC_mat)
    
            # Remove the version number in the end of the ENSEMBL id
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
        
            # Keep only one of the duplicated Ids 
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
        
    # 1 c. Prepare the count matrix and metadata
          # Combine the expression matrices into single count matrix. 
              all.equal(rownames(HNSC_mat), rownames(CESC_mat)) # Double checking if  the genes are the same in two matrices
              combined_counts <- cbind(CESC_mat, HNSC_mat)
              dim(combined_counts)
      
          # Create metadata for DESeq2
            sample_info <- data.frame(
              sample_id = colnames(combined_counts),
              cancer_type = c(rep("CESC", ncol(CESC_mat)), rep("HNSC", ncol(HNSC_mat)))
            )
            head(sample_info)
        
          # Set rownames
            rownames(sample_info) <- sample_info$sample_id
            head(sample_info)

  

#  Step 2 : Filter low-expressed genes
    # Keep genes with at least 10 counts in at least 3 samples
      keep_genes <- rowSums(combined_counts >=10) >=3
      filtered_counts <- combined_counts[keep_genes, ]
      dim(filtered_counts)

# Step 3: Create DESeq2 dataset and normalize
    # 3 a. Create DESEq2 dataset
      dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                    colData = sample_info,
                                    design = ~ cancer_type)
      
    # 3 b. Normalize and run DESEq2  
      dds <- DESeq(dds)
    
    # 3 c. Extract results 
          res <- results(dds, contrast = c("cancer_type", "HNSC", "CESC"), alpha = 0.05)
          res <- as.data.frame(res)
          # Save results
          write.csv(res, "/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/DEG_HNSC_vs_CESC.csv")

  # Step 4 :VISUALIZATION 
    # 4 a. MA plot
          res <- results(dds, contrast = c("cancer_type", "HNSC", "CESC"), alpha = 0.05)
          # changing res from datafram eto original deseq format before MA plot 
          plotMA(res, ylim = c(-8, 8))
          
    
    # 4 b. Volcano Plot
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = "log2FoldChange",
                    y = "pvalue",
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    title = "Using Bulk RNA-seq data from TCGA for DGEA",
                    subtitle = "Volcano Plot showing differentialy expressed genes in Head and Neck Squamous Cell carcinoma (HNSC),\nwhen compared to Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (CESC)",
                    caption = "log2FC >1, p < 0.05",
                    legendPosition = "right")
      ggsave("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/figures/volcano_plot_hnsc_vs_cesc.png", width = 16, height = 12, dpi = 300)
    
  # Step 5: Data interpretation
    # 5a. Sort Differential expressed genes by significance     
        res <- as.data.frame(res)
        res_sorted <- res[order(res$padj), ]
        head(res_sorted, 20) 
        # This gives us the most significantly differentially expressed genes between HNSC and CESC.
        
    # 5b. Identify up-regulated and down-regulated genes
        
        # Define significant DEGs (padj < 0.05 and absolute log2FC > 1)
        significant_genes <- res_sorted[res_sorted$padj < 0.05 & abs(res_sorted$log2FoldChange) > 1, ]
        
        # Upregulated genes
        # Upregulated genes in HNSC → Higher expression in Head and Neck Cancer vs. Cervical Cancer.
        upregulated_genes <- significant_genes[significant_genes$log2FoldChange > 1, ]
        head(upregulated_genes, 10)

        # Down-regulated genes
        # Downregulated genes in HNSC → Lower expression in HNSC, meaning they are more active in CESC. 
        downregulated_genes <- significant_genes[significant_genes$log2FoldChange < 1, ]
        head(downregulated_genes, 10)
        
        
        
        
  
  # Step 6 : Functional Enrichment Analysis : Investigating Biological functions
       
    # BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA"))
    # Load libraries
      library(clusterProfiler)
      library(org.Hs.eg.db)
      library(enrichplot)
      library(ReactomePA)
      library(ggplot2)
        
    # 6a. Prepare the list  of ENTREZ IDs of the genes that are significantly regulated
        # Convert gene symbols to ENTREZ ID
          gene_symbols <- rownames(significant_genes)
          gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        # Merge IDs with DEGs
          significant_genes$ENTREZID <- gene_entrez$ENTREZID[match(rownames(significant_genes), gene_entrez$SYMBOL)]
        # Remove NA values
          significant_genes <- significant_genes[!is.na(significant_genes$ENTREZID), ]
        # Extract ENtrez IDs for enrichment analysis 
          gene_list <- gene_entrez$ENTREZID
    
    # 6b. Perform  and Visualize GO enrichment analysis
          go_enrichment <- enrichGO(gene = gene_list,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENTREZID",
                                    ont = "all", # all 3 ontologies; write "BP" for Biological Process only
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
          
          dotplot(go_enrichment, split="ONTOLOGY", showCategory = 5) +
            ggtitle("GO enrichment analysis") + 
            theme(axis.text.y = element_text(size = rel(0.75))) +
            facet_grid(ONTOLOGY~., scale="free")
          ggsave("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/figures/enrichment_go_hnsc_vs_cesc.png", dpi = 300)
          
          
    # 6d. Perform and Visualize KEGG pathway analysis
          kegg_enrichment <- enrichKEGG(gene = gene_list,
                                        organism = "hsa",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05)
          
          dotplot(kegg_enrichment, showCategory = 10) + 
            ggtitle("KEGG pathway enrichment analysis")
          ggsave("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/figures/enrichment_kegg_hnsc_vs_cesc.png", dpi = 300)
          
          
    # Perform and Visualize Reactome Pathway Analysis
          reactome_enrichment <- enrichPathway(gene = gene_list, 
                                               organism = "human", 
                                               pAdjustMethod = "BH",
                                               pvalueCutoff = 0.05)
          
          dotplot(reactome_enrichment, showCategory = 10) +
            ggtitle("Reactome pathway enrichment analysis")
          ggsave("/Users/aayushojha/Documents/DGE/tcga/headneck_vs_cervical/figures/enrichment_reactome_hnsc_vs_cesc.png", dpi = 300)
          
          
          
          

                    
        
          
          
         
           
          
      
        