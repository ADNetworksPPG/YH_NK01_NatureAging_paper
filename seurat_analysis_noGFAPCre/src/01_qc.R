#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/")
indir <- paste0(basedir,
                "results/snRNA_mm10/analysis_hapoe_chr/01_cellranger_count/")
outdir <- paste0(basedir,
                 "results/snRNA_mm10/analysis_hapoe_chr/",
                 "05_seurat_analysis_filter_percent_mt_noGFAPCre")
path_postfix_filtered_matrix <- "/outs/filtered_feature_bc_matrix"
setwd(outdir)

#get the list of all samples
samples_list <- list.dirs(indir,full.names=FALSE,recursive=FALSE)
#remove GFAP-Cre samples from the list
samples_list <- samples_list[!(grepl("GFAP-Cre", samples_list))]

#read in the sample metadata
meta.data <- read.csv(paste0(basedir,
                             "rawdata/snRNA_mm10/YH_NK01_snRNAseq_PS19_fApoE_mouse_list.csv"))

#create a data frame to record cell count before and after filtering per sample
cells_per_sample <- data.frame(sample=character(), 
                               filter_stage=character(), 
                               number_of_cells=numeric(),
                               nFeature_RNA_cutoff=numeric(),
                               percent.mt_cutoff=numeric(),
                               stringsAsFactors = FALSE)

#read in data for each sample and perform QC
for(i in 1:length(samples_list)){
  obj.name <- strsplit(samples_list[i],"_")[[1]][1]
  datadir <- paste0(indir,samples_list[i],path_postfix_filtered_matrix)
  obj.data <- Read10X(data.dir = datadir)
  this.obj <- CreateSeuratObject(counts = obj.data)
  cells_pre_filter <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,
                                                   "pre_filter",
                                                   cells_pre_filter,
                                                   0,
                                                   0)
  this.obj[["percent.mt"]] <- PercentageFeatureSet(this.obj, pattern = "^mt-")

  #plot the distribution of QC metrics
  pdf(paste0(outdir,"/plot/01_qc/",obj.name,"_pre_qc_plot.pdf"))
  # Visualize QC metrics as a violin plot
  print(VlnPlot(this.obj, 
                features = c("nFeature_RNA", 
                             "nCount_RNA",
                             "percent.mt"),
                ncol = 3))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "nFeature_RNA"))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "percent.mt"))
  dev.off()
  
  #filter top 1% (99th quantile)
  percent.mt.cutoff <- 0.25
  nfeature.cutoff <- quantile(this.obj$nFeature_RNA, 0.99)
  this.obj <- subset(this.obj,
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < nfeature.cutoff & 
                       percent.mt < percent.mt.cutoff )
  
  #add metadata
  this.obj$sample_number <- strsplit(samples_list[i],"_")[[1]][1]
  pre_mouse_number <- substr(samples_list[i], 
                             nchar(samples_list[i])-3+1, 
                             nchar(samples_list[i]))
  this.obj$mouse_number <- as.numeric(gsub("[^0-9.]", "",  pre_mouse_number))
  this_obj_metadata <- meta.data[(meta.data[,2] == 
                                    unique(this.obj$mouse_number)),]
  this.obj$genotype <- this_obj_metadata$Genotype
  this.obj$sex <- this_obj_metadata$Sex
  this.obj$Cre <- this_obj_metadata$Cre
  this.obj$ApoE <- this_obj_metadata$ApoE
  this.obj$dob <- this_obj_metadata$DOB
  this.obj$age_perfused <- this_obj_metadata$Age.Perfused
  this.obj$date_perfused <- this_obj_metadata$Date.Perfused
  this.obj$date_of_nuclear_isolation <- this_obj_metadata$Date.of.Nuc..Isolation
  #get number of cells post filtering
  cells_post_filter <- ncol(this.obj)
  
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,
                                                   "post_filter",
                                                   cells_post_filter,
                                                   nfeature.cutoff,
                                                   percent.mt.cutoff)
  #save the filtered Seurat object
  saveRDS(this.obj, file = paste0(outdir,"/data/01_qc/",obj.name,".rds"))
  
}

#save the qc metadata
write.csv(cells_per_sample, 
          file = paste0(outdir,"/data/01_qc/qc_per_sample_metadata.csv"), 
          row.names = FALSE)

#create a plot for cell counts before and after filtering for each sample
pdf(paste0(outdir,
           "/plot/01_qc/",
           "qc_cells_per_sample_plot.pdf"),
    width = 10)
print(ggplot(cells_per_sample, 
             aes(factor(sample, levels = paste0("S",seq(1,16))), 
                 as.numeric(number_of_cells), 
                 fill = rev(filter_stage))) + 
        geom_bar(stat="identity", 
                 position = "dodge") + 
        scale_fill_brewer(palette = "Set1",
                          labels = c("Pre-QC", "Post-QC"))+
        guides(fill=guide_legend(title="QC stage")) +
        xlab("Sample id") +
        ylab("Number of cells") +
        ggtitle("Cells per sample - Pre and Post QC") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
  )
dev.off()


########################## END ##########################                      
