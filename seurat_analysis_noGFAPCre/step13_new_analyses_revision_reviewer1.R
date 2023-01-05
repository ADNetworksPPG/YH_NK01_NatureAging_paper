#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(tools)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/")
outdir <- paste0(basedir,
                 "Revision/new_analyses_reviewer1")
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
dat <- readRDS(paste0(basedir,
                      "04_clustering_noS2/",
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

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


############
#Dotplot for all excitatory neuron disease associated markers
############
#excitatory_neuron_clusters <- c(4:7, 9, 10, 17:23, 26, 28, 30, 32)
features_excitatory_neurons <- c("Syn1", "Dgkh", "Trhde", "Slc17a7", "Man1a", 
                                 "Calm1", "Calm2", "Camk2n1", "Fth1", "Nrgn", "Ubb", 
                                 "Unc5d", "Hspa8", "Hsp90aa1", "Hsp90ab1", "Cdh18", 
                                 "Dpp10", "Hs3st4", "Pde4d", "Slc1a2", "Mgat4c", 
                                 "Galntl6", "Gria1", "Grm8", "Pde1a", "Lingo2", 
                                 "Hs6st3", "Epha6", "Kcnip4", "Cacna2d3", "Lrrtm4", 
                                 "Kcnb2", "Kcnq5", "hapoE-transgene")
pdf(file = "dotplot_excitatory_neuron_disease_associated_markers.pdf",
    width = 12,
    height = 14)
print(DotPlot(dat, 
              #idents = excitatory_neuron_clusters,
              features = features_excitatory_neurons, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()


############
#Dotplot for all oligodendrocyte disease associated markers
############
#oligodendrocyte_clusters <- c(1, 2, 15)
features_oligodendrocytes <- c("Mbp", "Plp1", "Mag", "Calm1", "Calm2", "Nrgn", 
                               "Ubb", "Camk1d", "Fth1", "Hspa8", "Hsp90aa1", 
                               "Hsp90ab1", "St18", "Rnf220", "Pde4b", "Prr5l",
                               "Pcdh9", "Kirrel3", "Frmd5", "Serpina3n", "H2-D1", 
                               "H2-K1", "B2m", "C4b", "Klk6", "Cd63", "Il33", 
                               "Sgk1", "Cd9", "Fxyd1", "Plvap", "Fxyd7", "Rnase4", 
                               "hapoE-transgene", "Gstm1", "Snca")
pdf(file = "dotplot_oligodendrocyte_disease_associated_markers.pdf",
    width = 12,
    height = 14)
print(DotPlot(dat, 
              #idents = c(oligodendrocyte_clusters),
              features = features_oligodendrocytes, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()


############
#Dotplot of astrocyte subclusters for disease associated markers
############
features_astrocyte_subclusters <- c("Luzp2", "Slc7a10", "Mfge8", "Gfap", "Id3", "Aqp4", 
                                    "Myoc", "Id1", "Fabp7", "Ctsb", "Vim", "Osmr", 
                                    "Serpina3n", "Gsn", "Ggta1", "Trpm3", "Csmd1", 
                                    "C4b", "Cd9", "Sparcl1", "Plce1", "Sgcd", "Fos", 
                                    "S100a6", "hapoE-transgene", "Clu", "Slc1a2", "Mertk", 
                                    "Gpc5", "Msi2", "Slc1a3", "Npas3", "Prex2", "Kcnip4", 
                                    "Dpp10", "Meg3", "Ptprd", "Ubb", "Tmsb4x", "Calm1", 
                                    "Fth1", "Calm2", "Cst3", "Nrgn", "Camk1d", "Hspa8", 
                                    "Hsp90aa1", "Hsp90ab1")
pdf(file = "dotplot_astrocyte_subclusters_disease_associated_markers.pdf",
    width = 14,
    height = 8)
print(DotPlot(dat_astrocytes,
              features = features_astrocyte_subclusters, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()


############
#Dotplot of microglia subclusters for disease associated markers
############
features_microglia_subclusters <- c("Hexb", "Cst3", "Cx3cr1", "Ctsd", "Csf1r", "Ctss", 
                                    "Sparc", "Tmsb4x", "P2ry12", "C1qa", "C1qb", 
                                    "Tmem119", "Tyrobp", "Ctsb", "hapoE-transgene", 
                                    "B2m", "Fth1", "Lyz2", "Trem2", "Axl", "Cst7", 
                                    "Ctsl", "Lpl", "Cd9", "Csf1", "Ccl6", "Itgax", 
                                    "Clec7a", "Lilr4b", "Timp2", "Pde4b", "Nkain2", 
                                    "St18", "Prr5l", "Pcdh9", "Dpp10", "Kcnip4", "Meg3", 
                                    "Csmd1", "Nrg1", "Nrg3", "Calm1", "Ubb", 
                                    "Calm2", "Nrgn", "Camk1d", "Hspa8", "Hsp90aa1", 
                                    "Hsp90ab1", "Ptprd")
pdf(file = "dotplot_microglia_subclusters_disease_associated_markers.pdf",
    width = 14,
    height = 8)
print(DotPlot(dat_microglia,
              features = features_microglia_subclusters, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()


############
#volcano plots
############
de_gene_lists <- c("de_genes_cluster7_vs_allexcitatoryneurons_post_clustering_sct_noGFAPCre_noS2_pcadims_15_res_0.7.csv",
                   "de_genes_cluster15_vs_clusters1and2_post_clustering_sct_noGFAPCre_noS2_pcadims_15_res_0.7.csv",
                   "de_genes_cluster18_vs_allexcitatoryneurons_post_clustering_sct_noGFAPCre_noS2_pcadims_15_res_0.7.csv",
                   "de_genes_astrocytes_subcluster1_vs_other_astrocytes_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
                   "de_genes_astrocytes_subcluster5_vs_other_astrocytes_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
                   "de_genes_microglia_subcluster4_vs_all_other_microglia_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
                   "de_genes_microglia_subcluster6_vs_all_other_microglia_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
                   "de_genes_microglia_subcluster8_vs_all_other_microglia_noGFAPCre_noS2_pcadims_15_res_0.9.csv")
for(de_file in de_gene_lists){
  #read the de genes file
  de_genes <- read.csv(paste0("../../06_analysis_of_cluster_of_interest/",
                              "enriched_genes/",
                              de_file))
  #set the plot file name
  plot_filename <- paste0("volcano_plot_",
                          tools::file_path_sans_ext(de_file),
                          ".pdf")
  
  #set the height based on the input
  height_selected <- 12
  
  #set the genes to be labeled based on the input
  selected_genes <- NA
  if(grepl("astrocytes", de_file, fixed = TRUE)){
    selected_genes <- features_astrocyte_subclusters
  } else if(grepl("microglia", de_file, fixed = TRUE)){
    selected_genes <- features_microglia_subclusters
  } else if(grepl("cluster15_vs_clusters1and2", de_file, fixed = TRUE)){
    selected_genes <- features_oligodendrocytes
  } else{
    selected_genes <- features_excitatory_neurons
  }
  
  #set the limits for x-axis
  xlim_selected <- c(-3,3)
  if(max(abs(de_genes$avg_log2FC)) <= 2){
    xlim_selected = c(-2,2)
  }
  
  #generate volcano plot
  pdf(plot_filename,
      height = height_selected,
      width = 20)
  print(EnhancedVolcano(de_genes,
                        lab = de_genes$gene,
                        title = tools::file_path_sans_ext(de_file),
                        #titleLabSize = 8,
                        selectLab = selected_genes,
                        x = 'avg_log2FC',
                        y = 'p_val',
                        drawConnectors = TRUE,
                        FCcutoff = 0.4,
                        xlim = xlim_selected,
                        legendLabSize = 10,
                        ylab = bquote(~-Log[10] ~ italic(p-value)),
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value", 
                                         expression(p-value ~ and
                                                    ~ log[2] ~ FC))
  ))
    dev.off()
  
}


############
#Genotype-split (PS19-fE3, PS19-fE4, and PS19-fE4/Syn1-Cre) 
#Syn1 expression feature plot, with 2 or 3 different expression scales
############
#option #1
(FeaturePlot(dat, 
             features = "Syn1",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE) + 
   theme(legend.position = "right")) %>%
  ggsave(file = paste0("featureplot_",
                       "Syn1_nocutoff.pdf"),
         plot = .,
         width = 12,
         height = 7)

#option #2
(FeaturePlot(dat, 
             features = "Syn1",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0("featureplot_",
                       "Syn1_maxcutoff_2.pdf"),
         plot = .,
         width = 12,
         height = 7)
#option #3
(FeaturePlot(dat, 
             features = "Syn1",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0("featureplot_",
                       "Syn1_maxcutoff_1.5.pdf"),
         plot = .,
         width = 12,
         height = 7)
#option #4
(FeaturePlot(dat, 
             features = "Syn1",
             split.by = "genotype",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = paste0("featureplot_",
                       "Syn1_maxcutoff_1.pdf"),
         plot = .,
         width = 12,
         height = 7)

