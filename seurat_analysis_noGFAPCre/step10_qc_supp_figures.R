#Quality control metrics per cluster
#supplementary figure

#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/")
indir <- paste0(basedir, 
                "04_clustering_noS2/")
outdir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                 "12_paper_figures_noGFAP/")

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0("~/Downloads/",
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]

#make a table of the required data
avg_qc_features_per_cluster <- all_metadata %>% group_by(seurat_clusters) %>% 
  summarise(total_cells = n(),
            avg_ngene = mean(nFeature_RNA),
            sd_ngene = sd(nFeature_RNA),
            se_ngene = sd(nFeature_RNA) / sqrt(length(nFeature_RNA)),
            avg_nUMI = mean(nCount_RNA),
            sd_nUMI = sd(nCount_RNA),
            se_nUMI = sd(nCount_RNA) / sqrt(length(nCount_RNA)),
            avg_mt =  mean(percent.mt),
            sd_mt =  sd(percent.mt),
            se_mt = sd(percent.mt) / sqrt(length(percent.mt))) %>% 
  as.data.frame()

############
#ext.fig 4b
############
# Create plot with legend
ggp1_legend <- ggplot(avg_qc_features_per_cluster, 
                      aes(x=seurat_clusters, y=total_cells, 
                          fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  guides(fill=guide_legend(ncol=4)) + 
  scale_fill_discrete(name = "Cluster Identity", 
                      labels = c("1 - Oligodendrocyte",
                                 "2 - Oligodendrocyte",
                                 "3 - In Neuron",
                                 "4 - Ex Neuron",
                                 "5 - Ex Neuron",
                                 "6 - Ex Neuron",
                                 "7 - Ex Neuron",
                                 "8 - In Neuron",
                                 "9 - Ex Neuron",
                                 "10 - Ex Neuron",
                                 "11 - In Neuron",
                                 "12 - Astrocyte",
                                 "13 - In Neuron",
                                 "14 - Microglia",
                                 "15 - Oligodendrocyte",
                                 "16 - OPC",
                                 "17 - Ex Neuron",
                                 "18 - Ex Neuron",
                                 "19 - Ex Neuron",
                                 "20 - Ex Neuron",
                                 "21 - Ex Neuron",
                                 "22 - Ex Neuron",
                                 "23 - Ex Neuron",
                                 "24 - In Neuron",
                                 "25 - Microglia",
                                 "26 - Ex Neuron",
                                 "27 - In Neuron",
                                 "28 - Ex Neuron",
                                 "29 - Microglia",
                                 "30 - Ex Neuron",
                                 "31 - In Neuron",
                                 "32 - Ex Neuron",
                                 "33 - Unknown",
                                 "34 - Unknown")) +
  theme(text = element_text(size = 26))
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(ggp1_legend)

############
#ext.fig 4c
############
p1 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=total_cells, fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Cells per Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Cells")

############
#ext.fig 4d
############
# p2 <- ggplot(avg_qc_features_per_cluster,
#              aes(x=seurat_clusters, y=avg_ngene, fill=seurat_clusters)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_ngene-se_ngene, ymax=avg_ngene+se_ngene),
#                 width=.2,
#                 position=position_dodge(.9)) + 
#   geom_point(data = all_metadata, aes(y=nFeature_RNA),alpha =0.2, size = 0.3) +
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") +
#   ggtitle("Average Genes per Cells by Cluster") +
#   xlab("Cell Cluster") +
#   ylab("Average Number of Genes")
p2 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=nFeature_RNA, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3, position = position_nudge(x = -0.2)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Genes per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Genes per Cell") 

############
#ext.fig 4e
############
# p3 <- ggplot(avg_qc_features_per_cluster,
#              aes(x=seurat_clusters, y=avg_nUMI, fill=seurat_clusters)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_nUMI-se_nUMI, ymax=avg_nUMI+se_nUMI),
#                 width=.2,
#                 position=position_dodge(.9))  + 
#   geom_point(data = all_metadata, aes(y=nCount_RNA),alpha =0.2, size = 0.3) +
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") +
#   ggtitle("Average nUMI per Cells by Cluster") +
#   xlab("Cell Cluster") +
#   ylab("Average nUMI")
p3 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=nCount_RNA, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3, position = position_nudge(x = -0.2)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("nUMI per Cell by Cluster") +
  xlab("Cell Cluster") + 
  ylab("nUMI per Cell") 

############
#ext.fig 4f
############
# p4 <- ggplot(avg_qc_features_per_cluster,
#              aes(x=seurat_clusters, y=avg_mt, fill=seurat_clusters)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_mt-se_mt, ymax=avg_mt+se_mt),
#                 width=.2,
#                 position=position_dodge(.9))  + 
#   geom_point(data = all_metadata, aes(y=percent.mt),alpha =0.2, size = 0.3) +
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") +
#   ggtitle("Average % Mitochondrial Genes per Cells by Cluster") +
#   xlab("Cell Cluster") +
#   ylab("Average % Mito Genes")
p4 <- ggplot(all_metadata, 
             aes(x=seurat_clusters, y=percent.mt, fill=seurat_clusters)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3, position = position_nudge(x = -0.2)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("% Mitochondrial Genes per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("% Mito Genes per Cell") 



# Draw plots with shared legend
#with standard error
pdf(paste0(outdir, "supp_fig4b-f_boxplots_with_nudged_points.pdf"),
    width = 25, height = 26)
grid.arrange(shared_legend,
             arrangeGrob(p1, p2, p3,p4, ncol = 2),
             heights=c(2, 10)) 
dev.off()

print("********** Script completed! **********")

################## END ################## 