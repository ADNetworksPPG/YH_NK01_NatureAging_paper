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

#sub cluster the microglia clusters 14, 25 and 29
dat_microglia <- subset(x = dat, 
                        idents = c(14,25,29))
dat_microglia <- SCTransform(dat_microglia, 
                             method="glmGamPoi")
dat_microglia <- RunPCA(dat_microglia, 
                        assay = "SCT", 
                        verbose = FALSE, 
                        approx=FALSE)
pdf(paste0(outdir,
           "/plot/06_analysis_of_cluster_of_interest/",
           "sub_cluster_microglia_clusters_14_25_29/",
           "microglia_data_noGFAPCre_noS2_pcaplot.pdf"))
print(ElbowPlot(dat_microglia))
print(ElbowPlot(dat_microglia, ndims = 50))
dev.off()

pca_dims <- seq(10,30,5)
cluster_res <- seq(0.4,0.9,0.1)
marker_genes_microglia <- read.csv(paste0("/gladstone/bioinformatics/",
                                          "adnetworksppg/Project_1/NK01/",
                                          "rawdata/snRNA_mm10/",
                                          "Hippocampus_cellcluster_",
                                          "MarkerGenes.csv"))
marker_genes_microglia <- 
  marker_genes_microglia$marker_gene_symbol[
    marker_genes_microglia$hippocampal_region == "Microglia"]
marker_genes_microglia <- c(marker_genes_microglia,
                            "Apoe","hapoE-transgene",
                            "Mapt","Human-MAPT")

marker_genes_homeostatic_microglia <- c("P2ry12","Csf1r","Hexb","Cst3",
                                        "Cx3cr1","Siglech","Tgfbr1","Selplg",
                                        "Mef2a","Serinc3")

marker_genes_dam <- c("Cd9","Fth1","Plp1")

marker_genes_erbb_pathway <- c("Nrg1","Nrg2","Nrg3","Camk2a",
                               "Akt3","Ptk2","Mapk8","Erbb4")

marker_genes_alzheimers_disease <- c("Grin1","Grin2a","Grin2b","Itpr1",
                                     "Itpr2","Ryr2","Ryr3","Plcb1",
                                     "Plcb4","Cacna1c","Ppp3ca")
marker_genes_microglia <- unique(c(marker_genes_microglia,
                                   marker_genes_homeostatic_microglia,
                                   marker_genes_dam,
                                   marker_genes_erbb_pathway,
                                   marker_genes_alzheimers_disease))


for(i in 1:length(pca_dims)){
  dat_microglia_prepca <- NULL
  dat_microglia_prepca <- FindNeighbors(dat_microglia, 
                                        dims = 1:pca_dims[i])
  dat_microglia_prepca <- RunUMAP(dat_microglia_prepca, 
                                  dims = 1:pca_dims[i])
  for(j in 1:length(cluster_res)){
    dat_microglia_postpca <- NULL
    dat_microglia_postpca <- FindClusters(dat_microglia_prepca, 
                                          resolution =  cluster_res[j])
    
    #rename the cluster ids to start from 1 instead of 0
    dat_microglia_postpca$orig_seurat_clusters <- 
      dat_microglia_postpca$seurat_clusters
    dat_microglia_postpca$seurat_clusters <- factor(as.numeric(
      as.character(dat_microglia_postpca$orig_seurat_clusters))+1)
    Idents(dat_microglia_postpca) <- dat_microglia_postpca$seurat_clusters
    
    #visualizations
    #i. #Generate UMAP without cluster labels
    DimPlot(dat_microglia_postpca, 
            pt.size = 1,
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            reduction = "umap") %>%
      ggsave(file = file.path(outdir, 
                              paste0("plot/06_analysis_of_cluster_of_interest/",
                                     "sub_cluster_microglia_clusters_14_25_29/",
                                     "microglia_data_noGFAPCre_noS2", 
                                     "_pcadims_",
                                     pca_dims[i],
                                     "_res_",
                                     cluster_res[j],
                                     "_umap.pdf")),
             plot = .,
             width = 12,
             height = 7)
    
    #ii. Generate UMAP with cluster labels
    DimPlot(dat_microglia_postpca, 
            pt.size = 1,
            raster = FALSE, 
            order = TRUE, 
            label = TRUE, 
            reduction = "umap") %>%
      ggsave(file = file.path(outdir, 
                              paste0("plot/06_analysis_of_cluster_of_interest/",
                                     "sub_cluster_microglia_clusters_14_25_29/",
                                     "microglia_data_noGFAPCre_noS2", 
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
                      "sub_cluster_microglia_clusters_14_25_29/",
                      "microglia_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_split_genotype_umap.pdf"),
        width = 15,
        height = 7)
    print(DimPlot(dat_microglia_postpca, 
                  group.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_microglia_postpca, 
                  split.by = "genotype", 
                  group.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_microglia_postpca, 
                  split.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE))
    print(DimPlot(dat_microglia_postpca, 
                  split.by = "genotype",
                  pt.size = 1,
                  raster=FALSE,
                  order = TRUE,
                  label = TRUE))
    dev.off()
    
    #iv. Generate featureplot of all microglia marker genes
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_microglia_clusters_14_25_29/",
                      "microglia_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_featureplot.pdf"),
        width = 14,
        height = 7)
    for(f in 1:length(marker_genes_microglia)){
      print(FeaturePlot(dat_microglia_postpca, 
                        features = marker_genes_microglia[f],
                        pt.size = 1.5,
                        raster=FALSE,
                        order = TRUE))
    }
    dev.off()
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_microglia_clusters_14_25_29/",
                      "microglia_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_labeled_featureplot.pdf"),
        width = 14,
        height = 7)
    for(f in 1:length(marker_genes_microglia)){
      print(FeaturePlot(dat_microglia_postpca, 
                        features = marker_genes_microglia[f],
                        pt.size = 1.5,
                        raster=FALSE,
                        order = TRUE,
                        label = TRUE))
    }
    dev.off()
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_microglia_clusters_14_25_29/",
                      "microglia_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_split_genotype_featureplot.pdf"),
        width = 25,
        height = 7)
    for(f in 1:length(marker_genes_microglia)){
      print(FeaturePlot(dat_microglia_postpca, 
                        features= marker_genes_microglia[f], 
                        pt.size = 1.5,
                        split.by = "genotype",
                        raster = FALSE, 
                        order = TRUE, 
                        label = FALSE))
      print(FeaturePlot(dat_microglia_postpca, 
                        features= marker_genes_microglia[f],
                        pt.size = 1.5,
                        split.by = "genotype",
                        raster = FALSE, 
                        order = TRUE, 
                        label = TRUE))
    }
    dev.off()
    
    #v. Generate dotplot of all microglia marker genes
    pdf(file = paste0(outdir,
                      "/plot/06_analysis_of_cluster_of_interest/",
                      "sub_cluster_microglia_clusters_14_25_29/",
                      "microglia_data_noGFAPCre_noS2", 
                      "_pcadims_",
                      pca_dims[i],
                      "_res_",
                      cluster_res[j],
                      "_dotplot.pdf"),
        width = 14,
        height = 7)
    print(DotPlot(dat_microglia_postpca, 
                  features = marker_genes_microglia))
    print(DotPlot(dat_microglia_postpca, 
                  features = marker_genes_microglia, 
                  group.by = "genotype"))
    dev.off()
    
    print(paste0("***********",pca_dims[i]," and ",
                 cluster_res[j],"***********"))
    print(table(dat_microglia_postpca$genotype, 
                dat_microglia_postpca$seurat_clusters))
    saveRDS(dat_microglia_postpca, 
            file = paste0(outdir,
                          "/data/06_analysis_of_cluster_of_interest/",
                          "sub_cluster_microglia_clusters_14_25_29/",
                          "microglia_data_post_cluster_noGFAPCre_noS2_pcadims_",
                          pca_dims[i],
                          "_res_",
                          cluster_res[j],
                          ".rds"))
  }
}

print("*********** Script completed! ***********")

############### END ###############

