#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/")
outdir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                 "12_paper_figures_noGFAP/")
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(basedir,
                      "04_clustering_noS2/",
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

############
#figure 5a
############
(DimPlot(dat, 
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap")  +  
   theme(legend.key.size = unit(0.36, "cm"),
         legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size=1),
                              ncol=1) )) %>%
  ggsave(file = paste0(outdir, 
                       "Fig.5a_postcluster_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)


############
#figure 5e
############
#load the astrocytes Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astrocytes <- readRDS(paste0(basedir,
                                 "06_analysis_of_cluster_of_interest/",
                                 "sub_cluster_astrocytes_cluster_12/",
                                 "astrocyte_data_post_cluster_noGFAPCre_",
                                 "noS2_pcadims_",
                                 pca_dim_astro,
                                 "_res_",
                                 cluster_res_astro,
                                 ".rds"))
DimPlot(dat_astrocytes, 
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = paste0(outdir, 
                       "Fig.5e_astrocyte_subcluster_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)


############
#figure 5i
############
#load the microglia Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_microglia <- readRDS(paste0(basedir,
                                "06_analysis_of_cluster_of_interest/",
                                "sub_cluster_microglia_clusters_14_25_29/",
                                "microglia_data_post_cluster_noGFAPCre_",
                                "noS2_pcadims_",
                                pca_dim_micro,
                                "_res_",
                                cluster_res_micro,
                                ".rds"))
DimPlot(dat_microglia, 
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = paste0(outdir, 
                       "Fig.5i_microglia_subcluster_labeled_umap.pdf"),
         plot = .,
         width = 7,
         height = 7)


############
#extended figure 6a
############
features_extended_fig_6a <- c("Syn1","Dgkh","Trhde","Slc17a7",
                              "Man1a","Gad2","Mbp","Vcan",
                              "Acsbg1","Phkg1","Ly86","Cx3cr1",
                              "Hmgb1","Hspa8","Hsp90ab1","Hsp90aa1")
pdf(file = paste0(outdir, 
                  "Extended_Fig.6a_",
                  "dotplot.pdf"),
    width = 12,
    height = 8)
print(DotPlot(dat, 
              features = features_extended_fig_6a, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()


############
#extended figure 6b
############
#extended figure 6b - option #1
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.6b_",
                       "hapoe_maxcutoff_1.5.pdf"),
         plot = .,
         width = 12,
         height = 7)
#extended figure 6b - option #2
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.6b_",
                       "hapoe_maxcutoff_1.pdf"),
         plot = .,
         width = 12,
         height = 7)


############
#extended figure 6c
############
dat$clusters_of_interest <- ifelse(dat$seurat_clusters %in% c(7,15,18),
                                   dat$seurat_clusters,
                                   0)
(DimPlot(dat, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey","#619CFF", "#F8766D","#00BA38"),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster7", 
                                  "Cluster15", "Cluster18"),
                       values = c("grey","#619CFF", "#F8766D","#00BA38"))) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.6c_clusters_of_interest.pdf"),
         plot = .,
         width = 12,
         height = 7)


############
#extended figure 7a
############
dat_astrocytes$clusters_of_interest <- 
  ifelse(dat_astrocytes$seurat_clusters %in% c(1,3,5),
         dat_astrocytes$seurat_clusters,
         0)
(DimPlot(dat_astrocytes, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey","#619CFF", "#F8766D","#00BA38"),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster1", "Cluster3", "Cluster5"),
                       values = c("grey","#619CFF", "#F8766D","#00BA38"))) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.7a_astrocytes.pdf"),
         plot = .,
         width = 12,
         height = 7)


############
#extended figures 7b
############
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE) + 
   theme(legend.position = "right")) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.7b_astrocytes_hAPOE.pdf"),
         plot = .,
         width = 12,
         height = 7)


############
#extended figure 7i
############
dat_microglia$clusters_of_interest <- 
  ifelse(dat_microglia$seurat_clusters %in% c(4,6,8),
         dat_microglia$seurat_clusters,
         0)
(DimPlot(dat_microglia, 
         split.by = "genotype",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey","#619CFF", "#F8766D","#00BA38"),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster4", "Cluster6", "Cluster8"),
                       values = c("grey","#619CFF", "#F8766D","#00BA38"))) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.7i_microglia.pdf"),
         plot = .,
         width = 12,
         height = 7)

############
#extended figures 7j
############
(FeaturePlot(dat_microglia, 
             features = "hapoE-transgene",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0(outdir, 
                       "Extended_Fig.7j_microglia_hAPOE.pdf"),
         plot = .,
         width = 12,
         height = 7)

