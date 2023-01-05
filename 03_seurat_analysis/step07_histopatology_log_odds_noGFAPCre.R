
 #to identify clusters/cell-types whose membership changes with histopathological parameters
rm(list=ls())
library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)

#set all directory paths
sourcedir<-"~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/05_log_odds_calculation/"
histodir <- "~/Dropbox (Gladstone)/YH_NK01-Results/05_histopatology/input/"
outdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/"
plotdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/plots/"

setwd("~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology")


histop_GFAPCre <- read.csv(paste0(histodir,"YH_NK01_snRNA-Seq PS19-fApoE_mouse_List_MT.csv"), header = T, check.names=F)
counts <- read.csv(paste0(sourcedir,"counts_per_sample_per_cluster.csv"), header = T)
histop <- histop_GFAPCre %>% filter(Genotype != "PS19-fE4_GFAP-Cre")
pheno <- counts %>%
  select(sample_id, animal_model, total_numbers_of_cells_per_sample) %>%
  unique()
Clusters <- unique(counts$cluster_id)

histopheno <- merge(histop, pheno, by="sample_id") 
head(histopheno)
dim(histopheno)

##function to estimate the change in the odds of cluster membership from the E4 to the other genotypes
estimateCellStateChange <- function(k, counts, histopheno, optimizer) {
  require(lme4)
  require(gdata)
  print(paste("Cluster", k))
  cluster_counts <- counts %>% 
    filter(cluster_id == k)
  cluster_counts %<>% merge(., histopheno, all.y=TRUE)
  cluster_counts$number_of_cells_per_sample_in_cluster[is.na(cluster_counts$number_of_cells_per_sample_in_cluster)] <- 0
  
  cluster_counts %<>% arrange(animal_model) %>% mutate(proportion=number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample)



  colnames(cluster_counts)[16:23] <- c("Hippocampus_Vol_mm3", "prop_AT8_Coverage_Area", "prop_GFAP_Coverage_Area", "prop_S100B_Coverage_Area", "prop_IBA1_Coverage_Area", "prop_CD68_Coverage_Area", "prop_MBP_Coverage_Area","prop_OPC_Coverage_Area")

if (k %in% c(7, 15, 18)){
for (ph in 16:(ncol(cluster_counts)-1)){
pdf(paste0(plotdir,"proportion_of_cells_per_sample_cluster_",k,"_per_pheno_",colnames(cluster_counts[ph]),".pdf"))

p <- ggplot(cluster_counts, aes(((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample), cluster_counts[,ph]))
print(p + geom_point(aes(colour = factor(animal_model)) ) + scale_x_log10() + geom_smooth(aes(color = animal_model), method = "lm", se = FALSE, size=1) + xlab(paste0("proportion of cells per sample in cluster ", k)) + ylab(paste0(colnames(cluster_counts[ph]))))
dev.off()
}}



FitData <- NULL
  TempRes <- NULL
  
  
  for(i in 1:nrow(cluster_counts)) {
    number_in_cluster <- cluster_counts$number_of_cells_per_sample_in_cluster[i]
    number_not_in_cluster <- cluster_counts$total_numbers_of_cells_per_sample[i] - number_in_cluster
    InCluster=c(rep(1, number_in_cluster), rep(0, number_not_in_cluster))
    animal_model=rep(cluster_counts$animal_model[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    ID=rep(cluster_counts$sample_id[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    sex=rep(cluster_counts$Sex[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    `Hippocampus_Vol_mm3`=rep(cluster_counts$`Hippocampus_Vol_mm3`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    `prop_AT8_Coverage_Area`=rep(cluster_counts$`prop_AT8_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_GFAP_Coverage_Area`=rep(cluster_counts$`prop_GFAP_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_S100B_Coverage_Area`=rep(cluster_counts$`prop_S100B_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_IBA1_Coverage_Area`=rep(cluster_counts$`prop_IBA1_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_CD68_Coverage_Area`=rep(cluster_counts$`prop_CD68_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_MBP_Coverage_Area`=rep(cluster_counts$`prop_MBP_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
	`prop_OPC_Coverage_Area`=rep(cluster_counts$`prop_OPC_Coverage_Area`[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    TempData <- data.frame(ID, InCluster, animal_model, sex, `Hippocampus_Vol_mm3`, `prop_AT8_Coverage_Area`, `prop_GFAP_Coverage_Area`,`prop_S100B_Coverage_Area`,`prop_IBA1_Coverage_Area`, `prop_CD68_Coverage_Area`,`prop_MBP_Coverage_Area`, `prop_OPC_Coverage_Area`,check.names=F)
    FitData <- rbind(FitData, TempData)
  }
  FitData %<>% mutate(animal_model = as.factor(animal_model), InCluster=as.factor(InCluster), sex=as.factor(sex))
  FitData$animal_model <- relevel(FitData$animal_model, ref="PS19-fE4")
  
 # glmerFit <- glmer(InCluster ~ (1|ID) + animal_model + `%_CD68_Coverage_Area` + `%_CD68_Coverage_Area`:animal_model, data=FitData, family = "binomial", control = glmerControl(optimizer = "bobyqa"))

 
  for (j in 5:ncol(FitData)) {
   if(optimizer=="bobyqa"){
   glmerFit <- glmer(InCluster ~  (1|animal_model/ID) + FitData[,j] , data=FitData, family = "binomial", control = glmerControl(optimizer = "bobyqa"))
   }
######test the default optimizer
  if(optimizer=="nloptwrap"){
    glmerFit <- glmer(InCluster ~  (1|animal_model/ID) + FitData[,j] , data=FitData, family = "binomial", control = glmerControl(optimizer = "nloptwrap"))
   }
  sglmerFit1 <- summary(glmerFit)
  TempRes1 <- (sglmerFit1$coefficients[-1,])
  TempRes <- rbind(TempRes, TempRes1)
  }

  rownames(TempRes) <- c("Hippocampus_Vol_mm3", "prop_AT8_Coverage_Area", "prop_GFAP_Coverage_Area", "prop_S100B_Coverage_Area", "prop_IBA1_Coverage_Area", "prop_CD68_Coverage_Area","prop_MBP_Coverage_Area","prop_OPC_Coverage_Area")

print(TempRes)



}

#run the log odds function for all clusters for each histopathological variable
optimizer="bobyqa"
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, histopheno, optimizer)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:16, 25:32)] <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3",
                                         "logOddsRatio_for_unit_prop_AT8_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_GFAP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_S100B_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_IBA1_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_CD68_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_MBP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_OPC_Coverage_Area",
                                         "StdErr_for_unit_Hippocampus_Vol_mm3",
                                         "StdErr_for_unit_prop_AT8_Coverage_Area",
                                         "StdErr_for_unit_prop_GFAP_Coverage_Area",
                                         "StdErr_for_unit_prop_S100B_Coverage_Area",
                                         "StdErr_for_unit_prop_IBA1_Coverage_Area",
                                         "StdErr_for_unit_prop_CD68_Coverage_Area",
                                         "StdErr_for_unit_prop_MBP_Coverage_Area",
                                         "StdErr_for_unit_prop_OPC_Coverage_Area",
  										 "pvalue_for_unit_Hippocampus_Vol_mm3",
                                         "pvalue_for_unit_prop_AT8_Coverage_Area",
                                         "pvalue_for_unit_prop_GFAP_Coverage_Area",
                                         "pvalue_for_unit_prop_S100B_Coverage_Area",
                                         "pvalue_for_unit_prop_IBA1_Coverage_Area",
                                         "pvalue_for_unit_prop_CD68_Coverage_Area",  
                                         "pvalue_for_unit_prop_MBP_Coverage_Area", 
                                         "pvalue_for_unit_prop_OPC_Coverage_Area")
                                         
                                         
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_for_unit_Hippocampus_Vol_mm3`, 
                           ClusterRes$`pvalue_for_unit_prop_AT8_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_GFAP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_S100B_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_IBA1_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_CD68_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_MBP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_OPC_Coverage_Area`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust_for_unit_Hippocampus_Vol_mm3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust_for_unit_prop_AT8_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust_for_unit_prop_GFAP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):(nrow(ClusterRes)*3)]
ClusterRes[,"p.adjust_for_unit_prop_S100B_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*3 + 1):(nrow(ClusterRes)*4)]
ClusterRes[,"p.adjust_for_unit_prop_IBA1_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*4 + 1):(nrow(ClusterRes)*5)]
ClusterRes[,"p.adjust_for_unit_prop_CD68_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*5 + 1):(nrow(ClusterRes)*6)]
ClusterRes[,"p.adjust_for_unit_prop_MBP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*6 + 1):(nrow(ClusterRes)*7)]
ClusterRes[,"p.adjust_for_unit_prop_OPC_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*7 + 1):(length(p.adjust_all))]


##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X17","X18","X19","X20","X21","X22","X23","X24"))]
print(ClusterRes)
write.csv(t(ClusterRes), file = paste0(outdir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_cluster.csv"))











###default optimizer
#run the log odds function for all clusters for each histopathological variable
optimizer="nloptwrap"
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, histopheno, optimizer)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:12, 19:24)] <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3",
                                         "logOddsRatio_for_unit_prop_AT8_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_GFAP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_S100B_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_IBA1_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_CD68_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_MBP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_OPC_Coverage_Area"
                                         "StdErr_for_unit_Hippocampus_Vol_mm3",
                                         "StdErr_for_unit_prop_AT8_Coverage_Area",
                                         "StdErr_for_unit_prop_GFAP_Coverage_Area",
                                         "StdErr_for_unit_prop_S100B_Coverage_Area",
                                         "StdErr_for_unit_prop_IBA1_Coverage_Area",
                                         "StdErr_for_unit_prop_CD68_Coverage_Area",
                                         "StdErr_for_unit_prop_MBP_Coverage_Area",
                                         "StdErr_for_unit_prop_OPC_Coverage_Area",
  										 "pvalue_for_unit_Hippocampus_Vol_mm3",
                                         "pvalue_for_unit_prop_AT8_Coverage_Area",
                                         "pvalue_for_unit_prop_GFAP_Coverage_Area",
                                         "pvalue_for_unit_prop_S100B_Coverage_Area",
                                         "pvalue_for_unit_prop_IBA1_Coverage_Area",
                                         "pvalue_for_unit_prop_CD68_Coverage_Area",
                                         "pvalue_for_unit_prop_MBP_Coverage_Area",
                                         "pvalue_for_unit_prop_OPC_Coverage_Area")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_for_unit_Hippocampus_Vol_mm3`, 
                           ClusterRes$`pvalue_for_unit_prop_AT8_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_GFAP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_S100B_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_IBA1_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_CD68_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_MBP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_OPC_Coverage_Area`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust_for_unit_Hippocampus_Vol_mm3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust_for_unit_prop_AT8_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust_for_unit_prop_GFAP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):(nrow(ClusterRes)*3)]
ClusterRes[,"p.adjust_for_unit_prop_S100B_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*3 + 1):(nrow(ClusterRes)*4)]
ClusterRes[,"p.adjust_for_unit_prop_IBA1_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*4 + 1):(nrow(ClusterRes)*5)]
ClusterRes[,"p.adjust_for_unit_prop_CD68_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*5 + 1):(nrow(ClusterRes)*6)]
ClusterRes[,"p.adjust_for_unit_prop_MBP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*6 + 1):(nrow(ClusterRes)*7)]
ClusterRes[,"p.adjust_for_unit_prop_OPC_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*7 + 1):(length(p.adjust_all))]


##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X17","X18","X19","X20","X21","X22","X23","X24"))]
print(ClusterRes)
write.csv(t(ClusterRes), file = pate0(outdir,"log_odds_ratio_opt_nloptwrap_per_unit_per_histopathology_per_cluster.csv"))


