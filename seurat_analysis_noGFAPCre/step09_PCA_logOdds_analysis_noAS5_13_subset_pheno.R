rm(list=ls())
logOddsdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/05_log_odds_calculation/"
logOdds_astro_dir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/06_analysis_of_cluster_of_interest/log_odds_calculation_sub_cluster_astrocytes_cluster_12/"
logOdds_micro_dir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/06_analysis_of_cluster_of_interest/log_odds_calculation_sub_cluster_microglia_clusters_14_25_29/"

histodir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/"

outdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/PCA/"
setwd(outdir)
library(factoextra)
library(FactoMineR)
library(dplyr)
library(magrittr)
library(RColorBrewer)





##logOdds histopathology all clusters
logOdds_histo<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_cluster_subset_pheno.csv"))
logOdds_histo <- logOdds_histo[,c("X", "Cluster7", "Cluster15","Cluster18")]
##logOdds histopathology all astrocyte subclusters
logOdds_histo_astro<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_astrocytes_subset_pheno.csv"))
colnames(logOdds_histo_astro) <- paste("AS", colnames(logOdds_histo_astro), sep = "_")
logOdds_histo_astro <- logOdds_histo_astro[, c(1:5,7:13,16,14)]

##logOdds histopathology all microglia subclusters
logOdds_histo_micro<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_microglia_subset_pheno.csv"))
colnames(logOdds_histo_micro) <- paste("MG", colnames(logOdds_histo_micro), sep = "_")
logOdds_histo_micro <- logOdds_histo_micro[, c(1:6,16,7:15)]


#merge and extract all subclusters clusters - logOdds only
logOdds_histo_clusters_of_interest_all_sub <- merge(merge(logOdds_histo, logOdds_histo_astro, by.x="X", by.y="AS_X") , logOdds_histo_micro,  by.x="X",by.y="MG_X") 
logOdds_histo_clusters_of_interest_all_sub %<>% slice_head(n = 4)
colnames(logOdds_histo_clusters_of_interest_all_sub)[1] <- "logOddsHistopathology"






#all subclusters

t_histo<-t(logOdds_histo_clusters_of_interest_all_sub)

histo_log<-data.frame(cluster=colnames(logOdds_histo_clusters_of_interest_all_sub)[-1], as.numeric(t_histo[-1,1]),as.numeric(t_histo[-1,2]), as.numeric(t_histo[-1,3]), as.numeric(t_histo[-1,4]))
histo_log$cluster<-factor(histo_log$cluster, levels = histo_log$cluster)
colnames(histo_log)<-c("cluster",logOdds_histo_clusters_of_interest_all_sub$logOddsHistopathology)
histo_log_PCA<-histo_log[,-1]
histo_log_PCA<-data.frame(histo_log_PCA)
colnames(histo_log_PCA) <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3" ,     "logOddsRatio_for_unit_prop_AT8_Coverage_Area" , "logOddsRatio_for_unit_prop_MBP_Coverage_Area" ,  "logOddsRatio_for_unit_prop_OPC_Coverage_Area" )

res.pca.norm.histo<- prcomp(histo_log_PCA,  scale=T, center = T, rank.=2)

pca_histo_log<-data.frame(cluster=histo_log$cluster,res.pca.norm.histo$x)



library(RColorBrewer)
#histopatohlogy
#cl17-25 all subcluster
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
nb.cols <- 31

#rownames(res.pca.norm.histo) <- pca_histo_log$cluster
#pdf(paste0(outdir,"PCA_logodds_histopathology_cl17_25_all_subclusters_refined.pdf"))

p <- fviz_pca_ind(res.pca.norm.histo,
             geom.ind = c("point"),
                pointshape = 21,repel = T,
                pointsize = 5.5, title = "", fill.ind = pca_histo_log$cluster, mean.point = FALSE, addEllipses =F, legend.title = "") + 
                guides(x.sec = "axis", y.sec = "axis")+
                scale_fill_manual(values = c("firebrick1", "pink2","forestgreen",getPalette(nb.cols)))+
                #geom_text(aes(label=ifelse(data.frame(res.pca.norm.histo$x)$PC1>2 & data.frame(res.pca.norm.histo$x)$PC2>0,as.character(pca_histo_log$cluster),'')),position = position_dodge(0.9))
                 annotate("text", label = "MG cluster 8", x = 1.5 ,y = 3.4, size = 4) +
                 annotate("text", label = "MG cluster 6", x = -4.8 ,y = -1.75, size = 4) +
                 annotate("text", label = "Cluster 15", x = -4.1 ,y = 0.6, size = 4) +
                 annotate("text", label = "AS cluster 3", x = -1.3 ,y = 1.5, size = 4) +
                 annotate("text", label = "Cluster 7", x = -2.8 ,y = -0.2, size = 4) +
                 annotate("text", label = "Cluster 18", x = -2.3 ,y = 0.9, size = 4)
  ggsave(file = file.path(outdir, paste0("PCA_logodds_histopathology_cl7_15_18_all_subclusters_noAS_5_13_subset_pheno.pdf")),
         plot = p,
         width = 10,
         height = 10)








