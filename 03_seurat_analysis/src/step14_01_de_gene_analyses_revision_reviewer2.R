#!/usr/bin/env Rscript

#load required packages
library(Seurat)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                  "results/snRNA_mm10/analysis_hapoe_chr/",
                  "05_seurat_analysis_filter_percent_mt_noGFAPCre/data")
indir <- paste0(basedir,
                "/04_clustering_noS2/")
outdir <- paste0(basedir,
                 "/Revision_enriched_genes/")
if(!dir.exists(outdir)){
  dir.create(outdir)
}
setwd(outdir)


############
#load the data
############
#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))
#load the astrocytes Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astrocytes <- readRDS(paste0(basedir,
                                 "/06_analysis_of_cluster_of_interest/",
                                 "sub_cluster_astrocytes_cluster_12/",
                                 "astrocyte_data_post_cluster_noGFAPCre_",
                                 "noS2_pcadims_",
                                 pca_dim_astro,
                                 "_res_",
                                 cluster_res_astro,
                                 ".rds"))
#load the microglia Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_microglia <- readRDS(paste0(basedir,
                                "/06_analysis_of_cluster_of_interest/",
                                "sub_cluster_microglia_clusters_14_25_29/",
                                "microglia_data_post_cluster_noGFAPCre_",
                                "noS2_pcadims_",
                                pca_dim_micro,
                                "_res_",
                                cluster_res_micro,
                                ".rds"))


############
#perform DE analysis
#for the entire Seurat object
############
clusters_for_de <- c(7,15,18)
for(cluster in clusters_for_de){
  #de genes for PS19-fE4-Syn1-Cre vs PS19-fE4 
  de_cluster_syn1_vs_e4 <- FindMarkers(object = dat,
                                       ident.1 = "PS19-fE4 Syn1-Cre", 
                                       ident.2 = "PS19-fE4",
                                       verbose = TRUE, 
                                       group.by="genotype",
                                       subset.ident = cluster,
                                       assay = "SCT", 
                                       slot = "data", 
                                       test.use = "wilcox",
                                       logfc.threshold = 0.1)
  #gene names are by default row names
  #add gene names as a column in the results data frame
  de_cluster_syn1_vs_e4 <- cbind(gene=rownames(de_cluster_syn1_vs_e4), 
                                 de_cluster_syn1_vs_e4)
  #export the results in a csv file
  write.csv(de_cluster_syn1_vs_e4,
            file = paste0("de_genes_cluster",
                          cluster,
                          "_syn1-cre_vs_e4_",
                          "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  
  #de genes for PS19-fE3 vs PS19-fE4 
  de_cluster_e3_vs_e4 <- FindMarkers(object = dat,
                                     ident.1 = "PS19-fE3", 
                                     ident.2 = "PS19-fE4",
                                     verbose = TRUE, 
                                     group.by="genotype",
                                     subset.ident = cluster,
                                     assay = "SCT", 
                                     slot = "data", 
                                     test.use = "wilcox",
                                     logfc.threshold = 0.1)
  #gene names are by default row names
  #add gene names as a column in the results data frame
  de_cluster_e3_vs_e4 <- cbind(gene=rownames(de_cluster_e3_vs_e4), 
                               de_cluster_e3_vs_e4)
  #export the results in a csv file
  write.csv(de_cluster_e3_vs_e4,
            file = paste0("de_genes_cluster",
                          cluster,
                          "_e3_vs_e4_",
                          "post_clustering_sct_noGFAPCre_noS2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
}



############
#perform DE analysis
#for the astrocyte Seurat object
############
astro_subclusters_for_de <- c(1,5)
for(subcluster in astro_subclusters_for_de){
  #de genes for PS19-fE4-Syn1-Cre vs PS19-fE4
  #only for subcluster 1 as subcluster 5 does not have any cells from PS19-fE4-Syn1-Cre genotype
  if(subcluster != 5){
    de_subcluster_syn1_vs_e4 <- FindMarkers(object = dat_astrocytes,
                                            ident.1 = "PS19-fE4 Syn1-Cre", 
                                            ident.2 = "PS19-fE4",
                                            verbose = TRUE, 
                                            group.by="genotype",
                                            subset.ident = subcluster,
                                            assay = "SCT", 
                                            slot = "data", 
                                            test.use = "wilcox",
                                            logfc.threshold = 0.1)
    #gene names are by default row names
    #add gene names as a column in the results data frame
    de_subcluster_syn1_vs_e4 <- cbind(gene=rownames(de_subcluster_syn1_vs_e4), 
                                      de_subcluster_syn1_vs_e4)
    #export the results in a csv file
    write.csv(de_subcluster_syn1_vs_e4,
              file = paste0("de_genes_astrocytes_subcluster",
                            subcluster,
                            "_syn1-cre_vs_e4_",
                            "noGFAPCre_noS2_pcadims_",
                            pca_dim_astro,
                            "_res_",
                            cluster_res_astro,
                            ".csv"),
              row.names = FALSE)
  }
  
  #de genes for PS19-fE3 vs PS19-fE4 
  de_subcluster_e3_vs_e4 <- FindMarkers(object = dat_astrocytes,
                                        ident.1 = "PS19-fE3", 
                                        ident.2 = "PS19-fE4",
                                        verbose = TRUE, 
                                        group.by="genotype",
                                        subset.ident = subcluster,
                                        assay = "SCT", 
                                        slot = "data", 
                                        test.use = "wilcox",
                                        logfc.threshold = 0.1)
  #gene names are by default row names
  #add gene names as a column in the results data frame
  de_subcluster_e3_vs_e4 <- cbind(gene=rownames(de_subcluster_e3_vs_e4), 
                                  de_subcluster_e3_vs_e4)
  #export the results in a csv file
  write.csv(de_subcluster_e3_vs_e4,
            file = paste0("de_genes_astrocytes_subcluster",
                          subcluster,
                          "_e3_vs_e4_",
                          "noGFAPCre_noS2_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)
}

############
#perform DE analysis
#for the microglia Seurat object
############
micro_subclusters_for_de <- c(4,6,8)
for(subcluster in micro_subclusters_for_de){
  #de genes for PS19-fE4-Syn1-Cre vs PS19-fE4 
  #only for subcluster 4 and 8 as subcluster 6 has <3 cells from PS19-fE4-Syn1-Cre genotype
  if(subcluster != 6){
    de_subcluster_syn1_vs_e4 <- FindMarkers(object = dat_microglia,
                                            ident.1 = "PS19-fE4 Syn1-Cre", 
                                            ident.2 = "PS19-fE4",
                                            verbose = TRUE, 
                                            group.by="genotype",
                                            subset.ident = subcluster,
                                            assay = "SCT", 
                                            slot = "data", 
                                            test.use = "wilcox",
                                            logfc.threshold = 0.1)
    #gene names are by default row names
    #add gene names as a column in the results data frame
    de_subcluster_syn1_vs_e4 <- cbind(gene=rownames(de_subcluster_syn1_vs_e4), 
                                      de_subcluster_syn1_vs_e4)
    #export the results in a csv file
    write.csv(de_subcluster_syn1_vs_e4,
              file = paste0("de_genes_microglia_subcluster",
                            subcluster,
                            "_syn1-cre_vs_e4_",
                            "noGFAPCre_noS2_pcadims_",
                            pca_dim_micro,
                            "_res_",
                            cluster_res_micro,
                            ".csv"),
              row.names = FALSE)
  }
  
  #de genes for PS19-fE3 vs PS19-fE4 
  de_subcluster_e3_vs_e4 <- FindMarkers(object = dat_microglia,
                                        ident.1 = "PS19-fE3", 
                                        ident.2 = "PS19-fE4",
                                        verbose = TRUE, 
                                        group.by="genotype",
                                        subset.ident = subcluster,
                                        assay = "SCT", 
                                        slot = "data", 
                                        test.use = "wilcox",
                                        logfc.threshold = 0.1)
  #gene names are by default row names
  #add gene names as a column in the results data frame
  de_subcluster_e3_vs_e4 <- cbind(gene=rownames(de_subcluster_e3_vs_e4), 
                                  de_subcluster_e3_vs_e4)
  #export the results in a csv file
  write.csv(de_subcluster_e3_vs_e4,
            file = paste0("de_genes_microglia_subcluster",
                          subcluster,
                          "_e3_vs_e4_",
                          "noGFAPCre_noS2_pcadims_",
                          pca_dim_micro,
                          "_res_",
                          cluster_res_micro,
                          ".csv"),
            row.names = FALSE)
}



print("********** Script completed! **********")


################## END ################## 
