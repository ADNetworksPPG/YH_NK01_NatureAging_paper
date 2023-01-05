### Between -cluster comparison
#For this analyses, we will use the **lme4** package in R to fit generalized 
#linear mixed effects models. We are going the model the change in the chance 
#(or more formally the odds) of cells from a given mouse belonging to a given 
#cluster from the 2 genotypes to the PS19-fE4 genotype. The random effects part
#of these models captures the inherent correlation between the cells coming 
#from the same mouse.

#run locally on Ayushi's laptop

library(lme4)
library(dplyr)
require(magrittr)
library(Seurat)
library(ggplot2)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/",
                  "06_analysis_of_cluster_of_interest/")
indir <- paste0(basedir, 
                "sub_cluster_astrocytes_cluster_12/")
outdir <- paste0(basedir, 
                 "log_odds_calculation_sub_cluster_astrocytes_cluster_12/")
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.9
dat <- readRDS(paste0(indir,
                      "astrocyte_data_post_cluster_noGFAPCre_noS2_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]

#make a table of the required data
smp_cluster_counts <- unique(all_metadata %>%
                               group_by(sample_number) %>%
                               mutate(total_numbers_of_cells_per_sample = 
                                        n()) %>%
                               group_by(seurat_clusters, .add=TRUE) %>%
                               mutate(number_of_cells_per_sample_in_cluster = 
                                        n()) %>%
                               select(sample_number,
                                      genotype,
                                      seurat_clusters,
                                      total_numbers_of_cells_per_sample,
                                      number_of_cells_per_sample_in_cluster))
colnames(smp_cluster_counts)[1:3] <- c("sample_id","animal_model","cluster_id")
smp_cluster_counts$sample_id <- factor(smp_cluster_counts$sample_id,
                                       levels = paste0("S",seq(3,16)))
smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$sample_id,
                                               smp_cluster_counts$cluster_id),]
write.csv(smp_cluster_counts, 
          file = "counts_per_sample_per_cluster_subcluster_astrocytes.csv",
          row.names = FALSE)

counts <- read.csv("counts_per_sample_per_cluster_subcluster_astrocytes.csv", 
                   header = T)

pheno <- counts %>%
  select(sample_id, animal_model, total_numbers_of_cells_per_sample) %>%
  unique()
rownames(pheno) <- seq(1,nrow(pheno))
Clusters <- unique(counts$cluster_id)

##function to estimate the change in the odds of cluster membership from... 
##the E4 to the other genotypes
estimateCellStateChange <- function(k, counts, pheno) {
  require(lme4)
  require(gdata)
  print(paste("Cluster", k))
  cluster_counts <- counts %>% 
    filter(cluster_id == k)
  cluster_counts %<>% merge(., pheno, all.y=TRUE)
  cluster_counts[is.na(cluster_counts$number_of_cells_per_sample_in_cluster),
                 "number_of_cells_per_sample_in_cluster"] <- 0
  
  cluster_counts %<>% 
    arrange(animal_model) %>% 
    mutate(proportion=
             number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample)
  ##plot proportion of cells per genotype
  pdf(paste0("proportion_of_cells_per_genotype_cluster_subcluster_astrocytes_",
             k,
             ".pdf"))
  print(ggplot(cluster_counts, 
               aes(x=animal_model, 
                   y=((number_of_cells_per_sample_in_cluster+0.01)/total_numbers_of_cells_per_sample))) +
          geom_boxplot(outlier.color = NA) +
          geom_jitter(position = position_jitter(0.2)) +
          scale_y_log10() +
          ylab("proportion of cells per genotype") +
          ggtitle(paste0("Astrocyte sub cluster ",k)))
  dev.off()
  
  FitData <- NULL
  
  for(i in 1:nrow(cluster_counts)) {
    number_in_cluster <- cluster_counts$number_of_cells_per_sample_in_cluster[i]
    number_not_in_cluster <- cluster_counts$total_numbers_of_cells_per_sample[i] - number_in_cluster
    InCluster=c(rep(1, number_in_cluster), rep(0, number_not_in_cluster))
    animal_model=rep(cluster_counts$animal_model[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    ID=rep(cluster_counts$sample_id[i], cluster_counts$total_numbers_of_cells_per_sample[i])
    TempData <- data.frame(ID, InCluster, animal_model)
    FitData <- rbind(FitData, TempData)
  }
  FitData %<>% mutate(animal_model = as.factor(animal_model), InCluster=as.factor(InCluster))
  FitData$animal_model <- relevel(FitData$animal_model, ref="PS19-fE4")
  
  glmerFit <- glmer(InCluster ~ (1|ID) + animal_model, data=FitData, family = "binomial", control = glmerControl(optimizer = "bobyqa"))
  sglmerFit1 <- summary(glmerFit)
  TempRes1 <- (sglmerFit1$coefficients[-1,])
  
  return(TempRes1)
}

#run the log odds function for all clusters
ClusterRes <- sapply(Clusters, estimateCellStateChange, counts, pheno)
#reformat the results table
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1:4, 7:8)] <- c("logOddsRatio_PS19-fE3_vs_PS19-fE4",
                                       "logOddsRatio_PS19-fE4Syn1-Cre_vs_PS19-fE4",
                                       "standardError_PS19-fE3_vs_PS19-fE4",
                                       "standardError_PS19-fE4Syn1-Cre_vs_PS19-fE4",
                                       "pvalue-PS19-fE3",
                                       "pvalue-PS19-fE4Syn1-Cre")
# make a vector of all p-values and p.adjust for all p-values together
p.adjust_all <- p.adjust(c(ClusterRes$`pvalue-PS19-fE3`, 
                           ClusterRes$`pvalue-PS19-fE4Syn1-Cre`), method = "BH")
##perform multiple-testing correction
ClusterRes[,"p.adjust-PS19-fE3"] = p.adjust_all[1:nrow(ClusterRes)]
ClusterRes[,"p.adjust-PS19-fE4Syn1-Cre"] = p.adjust_all[(nrow(ClusterRes) + 1):(nrow(ClusterRes)*2)]

##output the results
ClusterRes <- ClusterRes[,!(colnames(ClusterRes) %in% c("X5","X6"))]
print(ClusterRes)
write.csv(ClusterRes, 
          file = "log_odds_ratio_per_genotype_per_cluster_subcluster_astrocytes.csv")


#make boxplots for each cluster
#get data in long form
ClusterRes$cluster_id <- rownames(ClusterRes)
x1 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_PS19-fE3_vs_PS19-fE4", 
                    "standardError_PS19-fE3_vs_PS19-fE4",
                    "pvalue-PS19-fE3",
                    "p.adjust-PS19-fE3")]
x1$genotype <- "PS19-fE3"
colnames(x1) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")
x2 <- ClusterRes[,c("cluster_id", 
                    "logOddsRatio_PS19-fE4Syn1-Cre_vs_PS19-fE4", 
                    "standardError_PS19-fE4Syn1-Cre_vs_PS19-fE4",
                    "pvalue-PS19-fE4Syn1-Cre",
                    "p.adjust-PS19-fE4Syn1-Cre")]
x2$genotype <- "PS19-fE4Syn1-Cre"
colnames(x2) <- c("cluster_id","logOddsRatio","standardError","pvalue","p.adjust","genotype")
ClusterRes_plot <- rbind(x1,x2)
rownames(ClusterRes_plot) <- seq(1,nrow(ClusterRes_plot))
ClusterRes_plot$cluster_id <- factor(ClusterRes_plot$cluster_id,
                                     levels = paste0("Cluster",seq(1,max(Clusters))))

#make the box plot for each cluster
#add label for p.adjusted < 0.05
label.df <- ClusterRes_plot[ClusterRes_plot$p.adjust < 0.05,]
if(!(nrow(label.df)==0)){label.df$logOddsRatio <- 4.5}
label.df2 <- ClusterRes_plot[ClusterRes_plot$pvalue < 0.05,]
if(!(nrow(label.df2)==0)){label.df2$logOddsRatio <- 4.1}
pdf("log_odds_barplot_subcluster_astrocytes_with_padj.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()

#clusters 5 and 13 have much higher standard error compared to other clusters
#to zoom into the log odds of other clusters..
#plot the results without clusters 5 and  13
pdf("log_odds_barplot_subcluster_astrocytes_with_padj_nocluster5and13.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot[!(ClusterRes_plot$cluster_id %in% c("Cluster5",
                                                           "Cluster13")),], 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  coord_cartesian(ylim=c(-6,5)) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()

#clusters 5 and 13 have much higher log odds ratio compared to other clusters
#zoom into the y axis to look at the log odds of other clusters..
pdf("log_odds_barplot_subcluster_astrocytes_with_padj_zoomin.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  coord_cartesian(ylim=c(-6,5)) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()


pdf("log_odds_barplot_subcluster_astrocytes_with_pvalue_and_padj.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue < 0.001,], label = "***",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.001 & label.df2$pvalue < 0.01,], label = "**",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.01 & label.df2$pvalue < 0.05,], label = "*",
            col = "gray", #hjust = -0.5,
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()

#clusters 5 and 13 has much higher log odds ratio compared to other clusters
#to zoom into the log odds of other clusters..
#plot the results without clusters 5 and 13
pdf("log_odds_barplot_subcluster_astrocytes_with_pvalue_and_padj_nocluster5and13.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot[!(ClusterRes_plot$cluster_id %in% c("Cluster5",
                                                           "Cluster13")),], 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  coord_cartesian(ylim=c(-6,5)) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue < 0.001,], label = "***",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.001 & label.df2$pvalue < 0.01,], label = "**",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.01 & label.df2$pvalue < 0.05,], label = "*",
            col = "gray", #hjust = -0.5,
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()


#clusters 5 and 13 has much higher log odds ratio compared to other clusters
#zoom into the y axis to look at the log odds of other clusters..
pdf("log_odds_barplot_subcluster_astrocytes_with_pvalue_and_padj_zoomin.pdf",
    height = 20,
    width = 35
)
ggplot(ClusterRes_plot, 
       aes(x=genotype, y=logOddsRatio, fill=genotype)) +
  geom_hline(yintercept=0) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=logOddsRatio-standardError, ymax=logOddsRatio+standardError), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(labels= c("PS19-fE3","PS19-fE4\nSyn1-Cre")) +
  facet_wrap(~ cluster_id) +
  theme(legend.position = "none",
        text = element_text(size = 28)
  ) +
  coord_cartesian(ylim=c(-6,5)) +
  geom_text(data = label.df, label = "*",
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue < 0.001,], label = "***",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.001 & label.df2$pvalue < 0.01,], label = "**",
            col = "gray", #hjust = -0.5,
            size = 15) +
  geom_text(data = label.df2[label.df2$pvalue > 0.01 & label.df2$pvalue < 0.05,], label = "*",
            col = "gray", #hjust = -0.5,
            size = 15) +
  xlab("") +
  ylab("log odds ratio: animal model vs PS19-fE4")
dev.off()

#add session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#################### END ####################