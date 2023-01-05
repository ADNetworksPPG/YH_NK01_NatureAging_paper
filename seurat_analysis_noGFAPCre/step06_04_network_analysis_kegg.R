#KEGG Enrichment Analysis of the enriched genes
#this script uses clusterProfiler::enrichKEGG() function

#run locally on Ayushi's laptop

#load the required packages
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyr)

#set the working directories
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "11_seurat_analysis_filter_percent_mt_noGFAPCre/",
                  "06_analysis_of_cluster_of_interest/")
indir <- paste0(basedir, 
                "enriched_genes/")
outdir <- paste0(basedir, 
                 "KEGG_enriched_pathways/")
dir.create(outdir, showWarnings = FALSE)
setwd(indir)


de_files <- list.files(pattern = "de_genes")
background_genes_allclusters <- read.table(
  paste0("nonzero_bakcground_gene_list_",
         "post_clustering_sct_noGFAPCre_noS2_pcadims_15_res_0.7.csv"),
  sep = ",",
  header = TRUE)
background_genes_astrocytes <- read.table(
  "nonzero_bakcground_gene_list_astrocyte_data_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
  sep = ",",
  header = TRUE)
background_genes_microglia <- read.table(
  "nonzero_bakcground_gene_list_microglia_data_noGFAPCre_noS2_pcadims_15_res_0.9.csv",
  sep = ",",
  header = TRUE)

for(i in 1:length(de_files)){
  #these files have DE genes that are significant
  signif_res <- read.table(de_files[i], 
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
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-dot.pdf"),
      height = 8)
  print(dotplot(ekeggbp, showCategory = 20, orderBy="GeneRatio") )
  dev.off()
  
  x2 <- enrichplot::pairwise_termsim(ekeggbp)
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-emap.pdf"),
      height = 8)
  print(emapplot(x2, showCategory = 40) )
  dev.off()
  
  #save the list of enriched pathways
  write.csv(ekeggbp,file = paste0(outdir,
                                  tools::file_path_sans_ext(de_files[i]),
                                  "_ekegg.csv"))

}

#save the session info
writeLines(capture.output(sessionInfo()), 
           paste0(outdir,"sessionInfo.txt"))

print("********** Script completed! **********")

################## END ################## 