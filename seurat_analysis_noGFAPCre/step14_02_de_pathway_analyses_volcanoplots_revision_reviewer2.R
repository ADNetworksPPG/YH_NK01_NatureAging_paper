#run locally on Ayushi's laptop

#paper revision
#Goal 1: KEGG Enrichment Analysis of the enriched genes
#this script uses clusterProfiler::enrichKEGG() function
#Goal 2: volcano plots of the enriched genes using EnhancedVolcano package

#load required packages
library(EnhancedVolcano)
library(tools)
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyr)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/")
indir_background_genes <- paste0(basedir, 
                                 "06_analysis_of_cluster_of_interest/enriched_genes/")
indir_de_genes <- paste0(basedir,
                         "Revision/new_analyses_reviewer2/enriched_genes/")
outdir_kegg <- paste0(basedir,
                      "Revision/new_analyses_reviewer2/KEGG_enriched_pathways/")
outdir_volcano_plot <- paste0(basedir,
                              "Revision/new_analyses_reviewer2/volcano_plot/")
if(!dir.exists(outdir_kegg)){
  dir.create(outdir_kegg)
}

if(!dir.exists(outdir_volcano_plot)){
  dir.create(outdir_volcano_plot)
}


############
#load the data
############
#get the list of all the de gene files
de_files <- list.files(indir_de_genes, pattern = "de_genes")
#load all the background gene lists for the 3 data sets
background_genes_allclusters <- read.table(paste0(indir_background_genes,
                                                  "nonzero_bakcground_gene_list_",
                                                  "post_clustering_sct_noGFAPCre_noS2_pcadims_15_res_0.7.csv"),
                                           sep = ",",
                                           header = TRUE)
background_genes_astrocytes <- read.table( paste0(indir_background_genes,
                                                  "nonzero_bakcground_gene_list_",
                                                  "astrocyte_data_noGFAPCre_noS2_pcadims_15_res_0.9.csv"),
                                           sep = ",",
                                           header = TRUE)
background_genes_microglia <- read.table(paste0(indir_background_genes,
                                                "nonzero_bakcground_gene_list_",
                                                "microglia_data_noGFAPCre_noS2_pcadims_15_res_0.9.csv"),
                                         sep = ",",
                                         header = TRUE)


############
#KEGG analysis and volcano plots
############
for(i in 1:length(de_files)){
  #these files have DE genes that are significant
  signif_res <- read.table(paste0(indir_de_genes, de_files[i]), 
                           sep = ",",header = TRUE)
  colnames(signif_res)[1] <- "Geneid"
  
  #--------------------------
  #select background set of genes
  #--------------------------
  if(grepl("astrocytes", de_files[i], fixed = TRUE)){
    background_genes <- background_genes_astrocytes
  } else if(grepl("microglia", de_files[i], fixed = TRUE)){
    background_genes <- background_genes_microglia
  }else{
    background_genes <- background_genes_allclusters
  }
  entrezIDs <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = background_genes$background_genes%>%
                                       as.character(),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  entrezIDs_signif <- AnnotationDbi::select(org.Mm.eg.db,
                                            keys = signif_res$Geneid %>%
                                              as.character(),
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")

  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- enrichKEGG(
    gene     = entrezIDs_signif$ENTREZID %>% subset(., !is.na(.)),
    universe = entrezIDs$ENTREZID %>% subset(., !is.na(.)),
    organism    = "mmu",
    minGSSize = 10,
    pvalueCutoff = 0.8,
    keyType = "ncbi-geneid"
  )

  #translating gene IDs to human readable symbols
  ekeggbp <- setReadable(ekeggbp, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

  #Visualize
  ## save images
  pdf(paste0(outdir_kegg,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-dot.pdf"),
      height = 8)
  print(dotplot(ekeggbp, showCategory = 20, orderBy="GeneRatio") )
  dev.off()

  #save the list of enriched pathways
  write.csv(ekeggbp,file = paste0(outdir_kegg,
                                  tools::file_path_sans_ext(de_files[i]),
                                  "_ekegg.csv"))
  
  #-----------------------------------
  #Volcano plot
  #------------------------------
  #set the height based on the input
  height_selected <- 12
  
  #set the limits for x-axis
  xlim_selected <- c(-3,3)
  if(max(abs(signif_res$avg_log2FC)) < 2){
    xlim_selected <- c(-2,2)
  }
  
  pdf(paste0(outdir_volcano_plot,
             "volcano_plot_",
             tools::file_path_sans_ext(de_files[i]),
             ".pdf"),
      height = height_selected,
      width = 20)
  print(EnhancedVolcano(signif_res,
                        lab = signif_res$Geneid,
                        title = tools::file_path_sans_ext(de_files[i]),
                        #titleLabSize = 8,
                        x = 'avg_log2FC',
                        y = 'p_val',
                        FCcutoff = 0.4,
                        pCutoff = 0.05,
                        xlim = xlim_selected,
                        drawConnectors = TRUE,
                        legendLabSize = 10,
                        ylab = bquote(~-Log[10] ~ italic(p-value)),
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value", 
                                         expression(p-value ~ and
                                                    ~ log[2] ~ FC))
  ))
  dev.off()
}

#save the session info
writeLines(capture.output(sessionInfo()), 
           paste0(basedir,"Revision/new_analyses_reviewer2/sessionInfo.txt"))

print("********** Script completed! **********")

################## END ################## 
