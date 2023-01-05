#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/",
                  "results/snRNA_mm10/analysis_hapoe_chr/",
                  "05_seurat_analysis_filter_percent_mt_noGFAPCre")
indir <- paste0(basedir,
                "/data/04_clustering_noS2/")
outdir <- basedir
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

#sub cluster the astrocytes cluster 12
dat_astrocytes <- subset(x = dat, idents = 12)
dat_astrocytes <- SCTransform(dat_astrocytes, 
                              method="glmGamPoi")
dat_astrocytes <- RunPCA(dat_astrocytes, 
                         assay = "SCT", 
                         verbose = FALSE, 
                         approx=FALSE)
pdf(paste0(outdir,
           "/plot/06_analysis_of_cluster_of_interest/",
           "sub_cluster_astrocytes_cluster_12/",
           "astrocyte_data_noGFAPCre_noS2_pcaplot.pdf"))
print(ElbowPlot(dat_astrocytes))
print(ElbowPlot(dat_astrocytes, ndims = 50))
dev.off()

pca_dims <- seq(10,30,5)
cluster_res <- seq(0.4,0.9,0.1)
marker_genes_astrocytes <- read.csv(paste0("/gladstone/bioinformatics/",
                                           "adnetworksppg/Project_1/NK01/",
                                           "rawdata/snRNA_mm10/",
                                           "Hippocampus_cellcluster_",
                                           "MarkerGenes.csv"
)
)
marker_genes_astrocytes <- 
  marker_genes_astrocytes$marker_gene_symbol[
    marker_genes_astrocytes$hippocampal_region == "Astrocytes"]
marker_genes_astrocytes <- c(marker_genes_astrocytes,
                             "Apoe","hapoE-transgene",
                             "Mapt","Human-MAPT")

for(i in 1:length(pca_dims)){
  dat_astrocytes_prepca <- NULL
  dat_astrocytes_prepca <- FindNeighbors(dat_astrocytes, 
                                          dims = 1:pca_dims[i])
  dat_astrocytes_prepca <- RunUMAP(dat_astrocytes_prepca, 
                                    dims = 1:pca_dims[i])
  for(j in 1:length(cluster_res)){
    dat_astrocytes_postpca <- NULL
    dat_astrocytes_postpca <- FindClusters(dat_astrocytes_prepca, 
                                           resolution =  cluster_res[j])
    
    #rename the cluster ids to start from 1 instead of 0
    dat_astrocytes_postpca$orig_seurat_clusters <- 
      dat_astrocytes_postpca$seurat_clusters
    dat_astrocytes_postpca$seurat_clusters <- factor(as.numeric(
      as.character(dat_astrocytes_postpca$orig_seurat_clusters))+1)
    Idents(dat_astrocytes_postpca) <- dat_astrocytes_postpca$seurat_clusters
    
    #visualizations
    #i. #Generate UMAP without cluster labels
    DimPlot(dat_astrocytes_postpca, 
            pt.size = 1,
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            reduction = "umap") %>%
      ggsave(file = file.path(outdir, 
                              paste0("plot/06_analysis_of_cluster_of_interest/",
                                     "sub_cluster_astrocytes_cluster_12/",
                                     "astrocyte_data_noGFAPCre_noS2", 
                                     "_pcadims_",
                                     pca_dims[i],
                                     "_res_",
                                     cluster_res[j],
                                     "_umap.pdf")),
             plot = .,
             width = 12,
             height = 7)
    
    #ii. Generate UMAP with cluster labels
    DimPlot(dat_astrocytes_postpca, 
            pt.size = 1,
            raster = FALSE, 
            order = TRUE, 
            label = TRUE, 
            reduction = "umap") %>%
      ggsave(file = file.path(outdir, 
                              paste0("plot/06_analysis_of_cluster_of_interest/",
                                     "sub_cluster_astrocytes_cluster_12/",
                                     "astrocyte_data_noGFAPCre_noS2", 
                                     "_pcadims_",
                                     pca_dims[i],
                                     "_res_",
                                     cluster_res[j],
                                     "_labeled_umap.pdf")),
             plot = .,
             width = 12,
             height = 7)
    
    #iii. Generate UMAP split by genotypes
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_astrocytes_cluster_12/",
                      "astrocyte_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_split_genotype_umap.pdf"),
        width = 15,
        height = 7)
    print(DimPlot(dat_astrocytes_postpca, 
                  group.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_astrocytes_postpca, 
                  split.by = "genotype", 
                  group.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_astrocytes_postpca, 
                  split.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_astrocytes_postpca, 
                  split.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE,
                  label = TRUE))
    dev.off()
    
    #iv. Generate featureplot of all astrocyte marker genes
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_astrocytes_cluster_12/",
                      "astrocyte_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_featureplot.pdf"),
        width = 14,
        height = 7)
    for(f in 1:length(marker_genes_astrocytes)){
      print(FeaturePlot(dat_astrocytes_postpca, 
                        features = marker_genes_astrocytes[f],
                        pt.size = 1.5,
                        raster=FALSE,
                        order = TRUE))
    }
    dev.off()
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_astrocytes_cluster_12/",
                      "astrocyte_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_labeled_featureplot.pdf"),
        width = 14,
        height = 7)
    for(f in 1:length(marker_genes_astrocytes)){
      print(FeaturePlot(dat_astrocytes_postpca, 
                        features = marker_genes_astrocytes[f],
                        pt.size = 1.5,
                        raster=FALSE,
                        order = TRUE,
                        label = TRUE))
    }
    dev.off()
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_astrocytes_cluster_12/",
                      "astrocyte_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_split_genotype_featureplot.pdf"),
        width = 25,
        height = 7)
    for(f in 1:length(marker_genes_astrocytes)){
      print(FeaturePlot(dat_astrocytes_postpca, 
                        features= marker_genes_astrocytes[f], 
                        pt.size = 1.5,
                        split.by = "genotype",
                        raster = FALSE, 
                        order = TRUE, 
                        label = FALSE))
      print(FeaturePlot(dat_astrocytes_postpca, 
                        features= marker_genes_astrocytes[f],
                        pt.size = 1.5,
                        split.by = "genotype",
                        raster = FALSE, 
                        order = TRUE, 
                        label = TRUE))
    }
    dev.off()
    
    #v. Generate dotplot of all astrocyte marker genes
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_astrocytes_cluster_12/",
                      "astrocyte_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_dotplot.pdf"),
        width = 14,
        height = 7)
    print(DotPlot(dat_astrocytes_postpca, 
                  features = marker_genes_astrocytes
    ))
    print(DotPlot(dat_astrocytes_postpca, 
                  features = marker_genes_astrocytes, 
                  group.by = "genotype"
    ))
    dev.off()
    
    print(paste0("***********",pca_dims[i]," and ",
                 cluster_res[j],"***********"))
    print(table(dat_astrocytes_postpca$genotype, 
                dat_astrocytes_postpca$seurat_clusters))
    saveRDS(dat_astrocytes_postpca, 
            file = paste0(outdir,
                          "/data/06_analysis_of_cluster_of_interest/",
                          "sub_cluster_astrocytes_cluster_12/",
                          "astrocyte_data_post_cluster_noGFAPCre_noS2_pcadims_",
                          pca_dims[i],
                          "_res_",
                          cluster_res[j],
                          ".rds"))
  }
}

print("*********** Script completed! ***********")

############### END ###############

