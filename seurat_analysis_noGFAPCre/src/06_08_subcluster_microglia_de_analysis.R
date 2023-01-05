#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/",
                  "results/snRNA_mm10/analysis_hapoe_chr/",
                  "05_seurat_analysis_filter_percent_mt_noGFAPCre/data/",
                  "06_analysis_of_cluster_of_interest")
indir <- paste0(basedir,
                "/sub_cluster_microglia_clusters_14_25_29/")
outdir <- paste0(basedir,
                 "/enriched_genes/")
setwd(outdir)

#load the Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_micro <- readRDS(paste0(indir,
                      "microglia_data_post_cluster_noGFAPCre_noS2_pcadims_",
                      pca_dim_micro,
                      "_res_",
                      cluster_res_micro,
                      ".rds"))

# 1. DE gene for:
## a. Cluster 4 vs all other microglia
## b. Cluster 6 vs all other microglia
## c. Cluster 8 vs all other microglia
clusters_for_de <- c(4,6,8)
other_clusters <- c(1,2,3,5,7,9,10,11,12,13,14,15)
for(de_clus in clusters_for_de){
  de_list <- FindMarkers(object = dat_micro, 
                         ident.1 = de_clus,
                         ident.2 = other_clusters,
                         assay = "SCT", 
                         slot = "data", 
                         test.use = "wilcox",
                         logfc.threshold = 0.1)
  de_list <- cbind(gene=rownames(de_list), 
                   de_list)
  write.csv(de_list,
            file = paste0("de_genes_microglia_subcluster",
                          de_clus,
                          "_vs_all_other_microglia_",
                          "noGFAPCre_noS2_pcadims_",
                          pca_dim_micro,
                          "_res_",
                          cluster_res_micro,
                          ".csv"),
            row.names = FALSE)
}



#2. find the list of background genes
all_data <- GetAssayData(dat_micro, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "microglia_data_noGFAPCre_noS2_pcadims_",
                        pca_dim_micro,
                        "_res_",
                        cluster_res_micro,
                        ".csv"),
          row.names = FALSE)

print("********** Script completed! **********")

################## END ################## 

