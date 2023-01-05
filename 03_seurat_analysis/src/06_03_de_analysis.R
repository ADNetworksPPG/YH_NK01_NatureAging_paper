#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/",
                  "results/snRNA_mm10/analysis_hapoe_chr/",
                  "05_seurat_analysis_filter_percent_mt_noGFAPCre/data")
indir <- paste0(basedir,
                "/04_clustering_noS2/")
outdir <- paste0(basedir,
                 "/06_analysis_of_cluster_of_interest/")
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#1. For cluster 15, determine the uniquely enriched genes and pathways of... 
#cluster 15 versus other two oligodendrocyte clusters (1 and 2).
de_cluster15_vs_1and2 <- FindMarkers(object = dat,
                                     ident.1 = 15,
                                     ident.2 = c(1,2),
                                     assay = "SCT",
                                     slot = "data",
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
de_cluster15_vs_1and2 <- cbind(gene=rownames(de_cluster15_vs_1and2),
                               de_cluster15_vs_1and2)
write.csv(de_cluster15_vs_1and2,
          file = paste0("enriched_genes/",
                        "de_genes_cluster15_vs_clusters1and2_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#2. For cluster 15, determine the uniquely enriched genes and pathways of...
#cluster 15 versus cluster 1.
de_cluster15_vs_1 <- FindMarkers(object = dat,
                                 ident.1 = 15,
                                 ident.2 = 1,
                                 assay = "SCT",
                                 slot = "data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0.1)
de_cluster15_vs_1 <- cbind(gene=rownames(de_cluster15_vs_1),
                           de_cluster15_vs_1)
write.csv(de_cluster15_vs_1,
          file = paste0("enriched_genes/",
                        "de_genes_cluster15_vs_cluster1_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#3. For cluster 7, determine the uniquely enriched genes and pathways of...
#cluster 7 versus other excitatory neuron clusters...
#(clusters 4, 5, 6, 9, 17, 19, 20, 21, 22, 23, 26, 28, 30, 32)
de_cluster7_vs_excitatoryneurons <- FindMarkers(object = dat,
                                                ident.1 = 7,
                                                ident.2 = c(4, 5, 6, 9, 17, 19, 
                                                            20, 21, 22, 23, 26, 
                                                            28, 30, 32),
                                                assay = "SCT",
                                                slot = "data",
                                                test.use = "wilcox",
                                                logfc.threshold = 0.1)
de_cluster7_vs_excitatoryneurons <-
  cbind(gene = rownames(de_cluster7_vs_excitatoryneurons),
        de_cluster7_vs_excitatoryneurons)
write.csv(de_cluster7_vs_excitatoryneurons,
          file = paste0("enriched_genes/",
                        "de_genes_cluster7_vs_allexcitatoryneurons_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#4. For cluster 18, determine the uniquely enriched genes and pathways of...
#cluster 18 versus other excitatory neuron clusters...
#(clusters 4, 5, 6, 9, 17, 19, 20, 21, 22, 23, 26, 28, 30, 32)
de_cluster18_vs_excitatoryneurons <- FindMarkers(object = dat,
                                                 ident.1 = 18,
                                                 ident.2 = c(4, 5, 6, 9, 17, 19, 
                                                             20, 21, 22, 23, 26, 
                                                             28, 30, 32),
                                                 assay = "SCT",
                                                 slot = "data",
                                                 test.use = "wilcox",
                                                 logfc.threshold = 0.1)
de_cluster18_vs_excitatoryneurons <-
  cbind(gene = rownames(de_cluster18_vs_excitatoryneurons),
        de_cluster18_vs_excitatoryneurons)
write.csv(de_cluster18_vs_excitatoryneurons,
          file = paste0("enriched_genes/",
                        "de_genes_cluster18_vs_allexcitatoryneurons_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#5. find the list of background genes
all_data <- GetAssayData(dat, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("enriched_genes/",
                        "nonzero_bakcground_gene_list_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#6. DE gene list for PC15/Res0.7 cluster 12 astrocyte (NOT sub-cluster):
#PS19-fE4/Syn1-Cre vs PS19-fE4 in cluster 12
de_cluster12_syn1_vs_e4 <- FindMarkers(object = dat,
                                       ident.1 = "PS19-fE4 Syn1-Cre", 
                                       ident.2 = "PS19-fE4",
                                       verbose = TRUE, 
                                       group.by="genotype",
                                       subset.ident = 12,
                                       assay = "SCT", 
                                       slot = "data", 
                                       test.use = "wilcox",
                                       logfc.threshold = 0.1)
de_cluster12_syn1_vs_e4 <- cbind(gene=rownames(de_cluster12_syn1_vs_e4), 
                                 de_cluster12_syn1_vs_e4)
write.csv(de_cluster12_syn1_vs_e4,
          file = paste0("enriched_genes/",
                        "de_genes_cluster12_syn1-cre_vs_e4_",
                        "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)


print("********** Script completed! **********")


################## END ################## 

