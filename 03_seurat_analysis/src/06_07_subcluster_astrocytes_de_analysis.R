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
                "/sub_cluster_astrocytes_cluster_12/")
outdir <- paste0(basedir,
                 "/enriched_genes/")
setwd(outdir)

#load the Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astro <- readRDS(paste0(indir,
                      "astrocyte_data_post_cluster_noGFAPCre_noS2_pcadims_",
                      pca_dim_astro,
                      "_res_",
                      cluster_res_astro,
                      ".rds"))

# 1. DE gene for:
## a. Cluster 1 vs all other astrocytes
## b. Cluster 3 vs all other astrocytes
## c. Cluster 5 vs all other astrocytes
clusters_for_de <- c(1,3,5)
other_clusters <- c(2,4,6,7,8,9,10,11,12,13,14,15)
for(de_clus in clusters_for_de){
  de_list <- FindMarkers(object = dat_astro, 
                         ident.1 = de_clus,
                         ident.2 = other_clusters,
                         assay = "SCT", 
                         slot = "data", 
                         test.use = "wilcox",
                         logfc.threshold = 0.1)
  de_list <- cbind(gene=rownames(de_list), 
                   de_list)
  write.csv(de_list,
            file = paste0("de_genes_astrocytes_subcluster",
                          de_clus,
                          "_vs_other_astrocytes_",
                          "noGFAPCre_noS2_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)
}



#2. find the list of background genes
all_data <- GetAssayData(dat_astro, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "astrocyte_data_noGFAPCre_noS2_pcadims_",
                        pca_dim_astro,
                        "_res_",
                        cluster_res_astro,
                        ".csv"),
          row.names = FALSE)

print("********** Script completed! **********")

################## END ################## 

