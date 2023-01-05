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
                "/data/01_qc/")
outdir <- basedir
setwd(indir)

#get the file paths for all the 12 sample (no GFAP-Cre samples) datasets
file_paths <- list.files(pattern = ".rds")

#based on the QC results, it was decided to remove 
#sample S2 from the downstream analysis
#so remove these samples from the file_paths list
file_paths <- file_paths[!(file_paths == "S2.rds")]
file_names <-  gsub(pattern = "\\.rds$", 
                    replacement = "", 
                    x = file_paths)
#read in all the samples in a list
data_list <- lapply(file_paths, readRDS)
names(data_list) <- file_names 

#merge all the datasets
data <- merge(data_list[[1]], 
              y=unlist(data_list[2:length(data_list)]), 
              add.cell.ids=file_names,
              project = "YH-NK01", 
              merge.data=FALSE)
print(data)
saveRDS(data, 
        file = paste0(outdir,
                      "/data/03_merge_and_normalize_noS2/",
                      "merged_data_noGFAPCre_noS2.rds"))

print("***** Dataset merge completed! *****")


#read in the merged dataset
data <- readRDS(paste0(outdir,
                       "/data/03_merge_and_normalize_noS2/",
                       "merged_data_noGFAPCre_noS2.rds"))

#normalization using SCTransform
data_sct <- SCTransform(data, 
                        method="glmGamPoi")
saveRDS(data_sct, 
        file = paste0(outdir,
                      "/data/03_merge_and_normalize_noS2/",
                      "merged_data_noGFAPCre_noS2_post_sct.rds"))

#Perform and PCA and UMAP on the SCTransformed data
data_sct_processed <- RunPCA(object = data_sct, 
                             assay = "SCT")
pdf(file.path(outdir, 
              paste0("/plot/03_merge_and_normalize_noS2/",
                     "sctransform_merged_data_noGFAPCre_noS2", 
                     "_pca_plot.pdf")),
    width = 12,
    height = 7)
ElbowPlot(data_sct_processed, ndims = 50, reduction = "pca")
ElbowPlot(data_sct_processed, ndims = 20, reduction = "pca")
ElbowPlot(data_sct_processed, ndims = 16, reduction = "pca")
dev.off()

data_sct_processed <- data_sct_processed %>% 
  RunUMAP(., dims = 1:15, verbose = TRUE)
saveRDS(data_sct_processed, 
        file = paste0(outdir,
                      "/data/03_merge_and_normalize_noS2/",
                      "merged_data_noGFAPCre_noS2_post_sct_processed.rds"))

#Generate visualizations using the various metadata
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sample_number", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_sample_number_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sample_number",
        split.by="sample_number", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_sample_number_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 48)
DimPlot(data_sct_processed, 
        raster = FALSE,
        order = TRUE, 
        label = FALSE, 
        group.by="genotype", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_genotype_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="genotype",
        split.by="genotype", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_genotype_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 28)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="mouse_number", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_mouse_number_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="mouse_number",
        split.by="mouse_number", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_mouse_number_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 48)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sex", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_sex_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="sex",
        split.by="sex", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_sex_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="Cre", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_cre_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE,
        group.by="Cre",
        split.by="Cre", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_cre_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="ApoE", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_apoe_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="ApoE",
        split.by="ApoE", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2",
                                 "_apoe_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="dob", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_dob_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="dob",
        split.by="dob", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_dob_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 48)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="age_perfused", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_age_perfused_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="age_perfused",
        split.by="age_perfused", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_age_perfused_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 29)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_perfused",
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_date_perfused_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_perfused",
        split.by="date_perfused", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_date_perfused_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 36)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_nuclear_isolation", 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_date_of_nuclear_isolation_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
DimPlot(data_sct_processed, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by="date_of_nuclear_isolation",
        split.by="date_of_nuclear_isolation", 
        reduction = "umap",
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_date_of_nuclear_isolation_split_umap.pdf")),
         plot = .,
         width = 12,
         height = 15)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_nCount_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_nFeature_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="percent.mt", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_percent.mt_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_merge_and_normalize_noS2/",
                                 "sctransform_merged_data_noGFAPCre_noS2", 
                                 "_nCount_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, paste0("plot/03_merge_and_normalize_noS2/",
                                         "sctransform_merged_data_",
                                         "noGFAPCre_noS2", 
                                         "_nFeature_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

print("***** Normalization completed! *****")


########################## END ########################## 
