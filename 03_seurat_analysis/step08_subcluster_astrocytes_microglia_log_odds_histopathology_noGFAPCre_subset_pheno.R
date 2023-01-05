### Between -cluster comparison
#For this analyses, we will use the **lme4** package in R to fit generalized 
#linear mixed effects models. We are going the model the change in the chance 
#(or more formally the odds) of cells from a given mouse belonging to a given 
#cluster from the 4 animal models (genotypes i.e. E4, E3, GFAP-Cre and Syn1-Cre). 
#The random effects part of these models captures the inherent correlation 
#between the cells coming from the same mouse

 #to identify clusters/cell-types whose membership changes with histopathological parameters

library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)

#set all directory paths
sourcedir<-"~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/06_analysis_of_cluster_of_interest/log_odds_calculation_sub_cluster_astrocytes_cluster_12/"
histodir <- "~/Dropbox (Gladstone)/YH_NK01-Results/05_histopatology/input/"
outdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/"
plotdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/plots/"

setwd("~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology")

#astrocytes
histop_GFAPCre <- read.csv(paste0(histodir,"YH_NK01_snRNA-Seq PS19-fApoE_mouse_List_MT.csv"), header = T, check.names=F)
counts <- read.csv(paste0(sourcedir,"counts_per_sample_per_cluster_subcluster_astrocytes.csv"), header = T)
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
  ##plot proportion of cells per genotype
  #pdf(paste0("proportion_of_cells_per_genotype_cluster_",k,".pdf"))
    #print(ggplot(cluster_counts, aes(x=animal_model, y=((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample))) +
        #geom_boxplot() +
        #geom_jitter() +
        #scale_y_log10() +
        #ylab("proportion of cells per genotype"))
    #dev.off()
  
####scatterplot
#5. You can also plot some of the significant associations in the a scatter plot -  
#there would be 16 (or 14?) points in this plot corresponding to each animal in the study, 
#the x-axis would correspond to the proportion of cells (in log scale using scale_x_log10) 
#from a given animal that is present in the cluster under consideration and the y-axis would
# represent the value of the histo path variable for this animal. You could also overlay 4 best 
#fit lines corresponding to each animal model using geom_smooth on the same plot.


#
  colnames(cluster_counts)[16:23] <- c("Hippocampus_Vol_mm3", "prop_AT8_Coverage_Area", "prop_GFAP_Coverage_Area", "prop_S100B_Coverage_Area", "prop_IBA1_Coverage_Area", "prop_CD68_Coverage_Area", "prop_MBP_Coverage_Area","prop_OPC_Coverage_Area")



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
colnames(ClusterRes)[c(1,2,7,8,9,10,15,16, 25,26,31,32)] <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3",
                                         "logOddsRatio_for_unit_prop_AT8_Coverage_Area",
                                           "logOddsRatio_for_unit_prop_MBP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_OPC_Coverage_Area",
                                         "StdErr_for_unit_Hippocampus_Vol_mm3",
                                         "StdErr_for_unit_prop_AT8_Coverage_Area",
                                           "StdErr_for_unit_prop_MBP_Coverage_Area",
                                         "StdErr_for_unit_prop_OPC_Coverage_Area",
  										 "pvalue_for_unit_Hippocampus_Vol_mm3",
                                         "pvalue_for_unit_prop_AT8_Coverage_Area",
                                        "pvalue_for_unit_prop_MBP_Coverage_Area",
                                         "pvalue_for_unit_prop_OPC_Coverage_Area")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_for_unit_Hippocampus_Vol_mm3`, 
                           ClusterRes$`pvalue_for_unit_prop_AT8_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_MBP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_OPC_Coverage_Area`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust_for_unit_Hippocampus_Vol_mm3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust_for_unit_prop_AT8_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust_for_unit_prop_MBP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):(nrow(ClusterRes)*3)]
ClusterRes[,"p.adjust_for_unit_prop_OPC_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*3 + 1):(length(p.adjust_all))]



##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X3","X4","X5","X6","X11","X12","X13","X14","X17","X18","X19","X20","X21","X22","X23","X24","X27","X28","X29","X30"))]
print(ClusterRes)
write.csv(t(ClusterRes), file = paste0(outdir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_astrocytes_subset_pheno.csv"))


#microglia
sourcedir<-"~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/06_analysis_of_cluster_of_interest/log_odds_calculation_sub_cluster_microglia_clusters_14_25_29/"

histop_GFAPCre  <- read.csv(paste0(histodir,"YH_NK01_snRNA-Seq PS19-fApoE_mouse_List_MT.csv"), header = T, check.names=F)
counts <- read.csv(paste0(sourcedir,"counts_per_sample_per_cluster_subcluster_microglia.csv"), header = T)
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
  ##plot proportion of cells per genotype
  #pdf(paste0("proportion_of_cells_per_genotype_cluster_",k,".pdf"))
    #print(ggplot(cluster_counts, aes(x=animal_model, y=((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample))) +
        #geom_boxplot() +
        #geom_jitter() +
        #scale_y_log10() +
        #ylab("proportion of cells per genotype"))
    #dev.off()
  
####scatterplot
#5. You can also plot some of the significant associations in the a scatter plot -  
#there would be 16 (or 14?) points in this plot corresponding to each animal in the study, 
#the x-axis would correspond to the proportion of cells (in log scale using scale_x_log10) 
#from a given animal that is present in the cluster under consideration and the y-axis would
# represent the value of the histo path variable for this animal. You could also overlay 4 best 
#fit lines corresponding to each animal model using geom_smooth on the same plot.


#
  colnames(cluster_counts)[16:23] <- c("Hippocampus_Vol_mm3", "prop_AT8_Coverage_Area", "prop_GFAP_Coverage_Area", "prop_S100B_Coverage_Area", "prop_IBA1_Coverage_Area", "prop_CD68_Coverage_Area", "prop_MBP_Coverage_Area","prop_OPC_Coverage_Area")

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
colnames(ClusterRes)[c(1,2,7,8,9,10,15,16, 25,26,31,32)] <- c("logOddsRatio_for_unit_Hippocampus_Vol_mm3",
                                         "logOddsRatio_for_unit_prop_AT8_Coverage_Area",
                                            "logOddsRatio_for_unit_prop_MBP_Coverage_Area",
                                         "logOddsRatio_for_unit_prop_OPC_Coverage_Area",
                                         "StdErr_for_unit_Hippocampus_Vol_mm3",
                                         "StdErr_for_unit_prop_AT8_Coverage_Area",
                                           "StdErr_for_unit_prop_MBP_Coverage_Area",
                                         "StdErr_for_unit_prop_OPC_Coverage_Area",
  										 "pvalue_for_unit_Hippocampus_Vol_mm3",
                                         "pvalue_for_unit_prop_AT8_Coverage_Area",
                                           "pvalue_for_unit_prop_MBP_Coverage_Area",
                                         "pvalue_for_unit_prop_OPC_Coverage_Area")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue_for_unit_Hippocampus_Vol_mm3`, 
                           ClusterRes$`pvalue_for_unit_prop_AT8_Coverage_Area`,
                             ClusterRes$`pvalue_for_unit_prop_MBP_Coverage_Area`,
                           ClusterRes$`pvalue_for_unit_prop_OPC_Coverage_Area`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust_for_unit_Hippocampus_Vol_mm3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust_for_unit_prop_AT8_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]
ClusterRes[,"p.adjust_for_unit_prop_MBP_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*2 + 1):(nrow(ClusterRes)*3)]
ClusterRes[,"p.adjust_for_unit_prop_OPC_Coverage_Area"] = p.adjust_all[(nrow(ClusterRes)*3 + 1):(length(p.adjust_all))]


##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X3","X4","X5","X6","X11","X12","X13","X14","X17","X18","X19","X20","X21","X22","X23","X24","X27","X28","X29","X30"))]
print(ClusterRes)
write.csv(t(ClusterRes), file = paste0(outdir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_microglia_subset_pheno.csv"))

