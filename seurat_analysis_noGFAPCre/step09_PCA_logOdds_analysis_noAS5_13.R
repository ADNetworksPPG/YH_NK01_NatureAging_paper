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

#logOdds per genotype - all the clusters and sublcusters
logOdds_all_clusters=read.csv(paste0(logOddsdir,"log_odds_ratio_per_genotype_per_cluster.csv"))
logOdds_astro=read.csv(paste0(logOdds_astro_dir,"log_odds_ratio_per_genotype_per_cluster_subcluster_astrocytes.csv"))
logOdds_astro$X <- paste("AS", logOdds_astro$X, sep = "_")

#excluding AS 5 and AS 13
logOdds_astro <- logOdds_astro[-which(logOdds_astro$X == "AS_Cluster5" | logOdds_astro$X == "AS_Cluster13" ),]
logOdds_micro=read.csv(paste0(logOdds_micro_dir,"log_odds_ratio_per_genotype_per_cluster_subcluster_microglia.csv"))
logOdds_micro$X <- paste("MG", logOdds_micro$X, sep = "_")

#merge and extract clusters of interests - logOdds only
logOdds_clusters_of_interest <- bind_rows(logOdds_all_clusters[which(logOdds_all_clusters$X=="Cluster7"  | logOdds_all_clusters$X== "Cluster15" | logOdds_all_clusters$X=="Cluster18" ),],logOdds_astro, logOdds_micro)
logOdds_clusters_of_interest <- logOdds_clusters_of_interest[,c(1:4)]
logOdds_clusters_of_interest$X<-factor(logOdds_clusters_of_interest$X, levels = logOdds_clusters_of_interest$X)
logOdds_all_clusters_PCA<-logOdds_clusters_of_interest[, -1]
res.pca.norm.all<- prcomp(logOdds_all_clusters_PCA,  scale=T, center = T, rank.=3)
pca_logOdds_all_clusters<-data.frame(cluster=logOdds_clusters_of_interest$X,res.pca.norm.all$x)
#genotype



pdf(paste0(outdir,"PCA_logOdds_by_genotype_all_subclusters_noAS5_13.pdf"))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
nb.cols <- 31

fviz_pca_ind(res.pca.norm.all,
             geom.ind = "point",
                pointshape = 21,
                pointsize = 7.5, title = "", fill.ind = pca_logOdds_all_clusters$cluster, mean.point = FALSE, addEllipses =F, legend.title = "") + 
                guides(x.sec = "axis", y.sec = "axis")+
                scale_fill_manual(values = getPalette(nb.cols))
  
  fviz_pca_biplot(res.pca.norm.all, geom.ind = "point",
                pointshape = 21,
                pointsize = 7.5, title = "", fill.ind = pca_logOdds_all_clusters$cluster, mean.point = FALSE, addEllipses =F, legend.title = "") + 
guides(x.sec = "axis", y.sec = "axis")+
                scale_fill_manual(values = getPalette(nb.cols))
   



dev.off()







##logOdds histopathology all clusters
logOdds_histo<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_cluster.csv"))
logOdds_histo <- logOdds_histo[,c("X", "Cluster7", "Cluster15","Cluster18")]
##logOdds histopathology all astrocyte subclusters
logOdds_histo_astro<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_astrocytes.csv"))
colnames(logOdds_histo_astro) <- paste("AS", colnames(logOdds_histo_astro), sep = "_")
logOdds_histo_astro <- logOdds_histo_astro[, c(1:5,7:13,16,14)]

##logOdds histopathology all microglia subclusters
logOdds_histo_micro<-read.csv(paste0(histodir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_microglia.csv"))
colnames(logOdds_histo_micro) <- paste("MG", colnames(logOdds_histo_micro), sep = "_")
logOdds_histo_micro <- logOdds_histo_micro[, c(1:6,16,7:15)]


#merge and extract all subclusters clusters - logOdds only
logOdds_histo_clusters_of_interest_all_sub <- merge(merge(logOdds_histo, logOdds_histo_astro, by.x="X", by.y="AS_X") , logOdds_histo_micro,  by.x="X",by.y="MG_X") 
logOdds_histo_clusters_of_interest_all_sub %<>% slice_head(n = 8)
colnames(logOdds_histo_clusters_of_interest_all_sub)[1] <- "logOddsHistopathology"






#all subclusters

t_histo<-t(logOdds_histo_clusters_of_interest_all_sub)

histo_log<-data.frame(cluster=colnames(logOdds_histo_clusters_of_interest_all_sub)[-1], as.numeric(t_histo[-1,1]),as.numeric(t_histo[-1,2]), as.numeric(t_histo[-1,3]), as.numeric(t_histo[-1,4]), as.numeric(t_histo[-1,5]), as.numeric(t_histo[-1,6]),as.numeric(t_histo[-1,7]),as.numeric(t_histo[-1,8]))
histo_log$cluster<-factor(histo_log$cluster, levels = histo_log$cluster)
colnames(histo_log)<-c("cluster",logOdds_histo_clusters_of_interest_all_sub$X)
histo_log_PCA<-histo_log[,-1]
histo_log_PCA<-data.frame(histo_log_PCA)
colnames(histo_log_PCA) <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3" ,     "logOddsRatio_for_unit_prop_AT8_Coverage_Area"  , "logOddsRatio_for_unit_prop_CD68_Coverage_Area",  "logOddsRatio_for_unit_prop_GFAP_Coverage_Area",  "logOddsRatio_for_unit_prop_IBA1_Coverage_Area" , "logOddsRatio_for_unit_prop_MBP_Coverage_Area" ,  "logOddsRatio_for_unit_prop_OPC_Coverage_Area" ,  "logOddsRatio_for_unit_prop_S100B_Coverage_Area")

res.pca.norm.histo<- prcomp(histo_log_PCA,  scale=T, center = T, rank.=5)

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
                 annotate("text", label = "MG cluster 8", x = -1.3 ,y = 3.5, size = 4) +
                 annotate("text", label = "MG cluster 6", x = -5.2 ,y = -3.5, size = 4) +
                 annotate("text", label = "Cluster 15", x = -5.8 ,y = -0.7, size = 4) +
                 annotate("text", label = "AS cluster 3", x = -3 ,y = 1.35, size = 4) +
                 annotate("text", label = "Cluster 7", x = -2.8 ,y = -0.72, size = 4) +
                 annotate("text", label = "Cluster 18", x = -4.3 ,y = 0.2, size = 4)
  ggsave(file = file.path(outdir, paste0("PCA_logodds_histopathology_cl7_15_18_all_subclusters_noAS_5_13.pdf")),
         plot = p,
         width = 10,
         height = 10)








