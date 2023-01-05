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
outdir <- paste0(basedir,
                 "Revision/new_analyses_reviewer2/")

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_noGFAPCre_noS2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]
all_metadata$mouse_number <- as.factor(all_metadata$mouse_number)

#make a table of the required data
avg_qc_features_per_mouse <- all_metadata %>% group_by(mouse_number) %>% 
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
avg_qc_features_per_mouse$mouse_number <- unique(factor(avg_qc_features_per_mouse$mouse_number, 
                                                        levels = sort(unique(avg_qc_features_per_mouse$mouse_number))))


##########
#Average Number of Genes per mouse
##########
# p1 <- ggplot(avg_qc_features_per_mouse,
#              aes(x=mouse_number, y=avg_ngene, fill=mouse_number)) +
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
#   ggtitle("Average Genes per Cells by individual animal") +
#   xlab("Mouse Number") +
#   ylab("Average Number of Genes")
p1 <- ggplot(all_metadata, 
             aes(x=mouse_number, y=nFeature_RNA, fill=mouse_number)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3, position = position_nudge(x = -0.1)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Genes per Cell by individual animal") +
  xlab("Mouse number") + 
  ylab("Number of Genes per Cell") 

##########
#UMI counts per mouse
##########
# p2 <- ggplot(avg_qc_features_per_mouse,
#              aes(x=mouse_number, y=avg_nUMI, fill=mouse_number)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   geom_errorbar(aes(ymin=avg_nUMI-se_nUMI, ymax=avg_nUMI+se_nUMI),
#                 width=.2,
#                 position=position_dodge(.9)) +
#   geom_point(data = all_metadata, aes(y=nCount_RNA),alpha =0.2, size = 0.3) +
#   scale_x_discrete(expand = c(0.04,0.04))+
#   theme_classic() + theme(
#     panel.border = element_rect(colour = "black", fill=NA, size=2),
#     text = element_text(size = 20),
#     legend.position = "none") +
#   ggtitle("Average nUMI per Cells by individual animal") +
#   xlab("Mouse Number") +
#   ylab("Average nUMI")
p2 <- ggplot(all_metadata, 
             aes(x=mouse_number, y=nCount_RNA, fill=mouse_number)) + 
  geom_boxplot(outlier.color = NA) + stat_boxplot(geom='errorbar') +
  geom_point(alpha =0.2, size = 0.3, position = position_nudge(x = -0.1)) +
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("nUMI per Cell by individual animal") +
  xlab("Mouse number") + 
  ylab("nUMI per Cell")

# Draw plots with shared legend
#with standard error
pdf(paste0(outdir, "supp_fig4g-h_boxplots_with_nudged_points.pdf"),
    width = 25)
grid.arrange(arrangeGrob(p1, p2, ncol = 2)) 
dev.off()

print("********** Script completed! **********")

################## END ################## 