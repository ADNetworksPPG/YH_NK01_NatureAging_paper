#heatmap for cluster 17 and cluster 25 with and without rowmeans normalization


sourcedir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/07_histopathology/"
outdir <- "~/Dropbox (Gladstone)/YH_NK01-Results/11_seurat_analysis_filter_percent_mt_noGFAPCre/Revision/07_histopathology/heatmaps_previous_method/"
#outdir <- "~/Dropbox (Gladstone)/MT/05_histopatology/020122/Figure4/"
setwd(sourcedir)
library(dplyr)
library(ggplot2)
library(dendsort)
library(pheatmap)
library(RColorBrewer)

histo_all<-read.csv(paste0(sourcedir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_cluster_subset_pheno.csv"))

histo<- histo_all  %>% slice_head(n = 4)

#without rowmeans across clusters normalization
histo_no_norm <- data.frame(histo_parameters=histo[1], histo[2:ncol(histo)])
colnames(histo_no_norm)[1]<-"histo_parameters"


#with rowmeans across clusters normalization
histo_rowMeans_norm <- data.frame(histo_parameters=histo[1], histo[2:ncol(histo)] - rowMeans(histo[2:ncol(histo)]))
colnames(histo_rowMeans_norm)[1]<-"histo_parameters"

#subset with cluster of interest
histo_no_norm_cluster_of_interest <-  histo_no_norm %>% select(-histo_parameters)
histo_rowMeans_norm_cluster_of_interest <-  histo_rowMeans_norm %>% select(-histo_parameters)


#clustering and heatmap
#Without rowmeans normalization

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 3
data <- histo_no_norm_cluster_of_interest
set.seed(123)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss

#plot(1:k.max, wss,
#     type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")





##perform k-means clustering with k=3
kNClust <- 3
kmeanClust <- kmeans(histo_no_norm_cluster_of_interest,kNClust)

row.names(histo_no_norm_cluster_of_interest) <- histo_no_norm$histo_parameters
no_norm_histo_kmeans <- histo_no_norm_cluster_of_interest #[order(kmeanClust$cluster),]


kGaps <- vector(mode = "numeric")
kGaps[1] <- kmeanClust$size[1] + 1
for(i in 2:(kNClust-1)) {
  kGaps[i] <- kGaps[i-1] + kmeanClust$size[i]
}

simpleredbluecols = colorRampPalette(c("blue","white","red"))(400)

df <- data.frame(histo_parameters=histo_no_norm$histo_parameters)
row.names(df) <- histo_no_norm$histo_parameters





paletteLength <- 400
##center the color scale so that while represents 0 in the row-centered matrix
myBreaks <- c(seq(min(histo_no_norm_cluster_of_interest, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(histo_no_norm_cluster_of_interest, na.rm = TRUE)/paletteLength, max(histo_no_norm_cluster_of_interest, na.rm=TRUE), length.out=floor(paletteLength/2)))


#for dendrogram
#mat_cluster_cols <- hclust(dist(t(no_norm_histo_kmeans)))
#sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#cmat_cluster_cols <- sort_hclust(mat_cluster_cols)
#mat_cluster_rows <- sort_hclust(hclust(dist(no_norm_histo_kmeans)))

cl<- data.frame(colnames(no_norm_histo_kmeans))
row.names(cl) <- colnames(no_norm_histo_kmeans)
no_norm_histo_kmeans



pdf(paste0(outdir,"Fig.4e_all_clusters_pathology_34_clusters.pdf"), height=4, width=24)
print(pheatmap(no_norm_histo_kmeans, cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, treeheight_row=T,  gaps_row = F, annot_cols=df,annot_row=cl, color = simpleredbluecols, breaks = myBreaks, display_numbers = TRUE, number_format="%.3f",
         number_color = "black", 
         fontsize_number = 10))
dev.off()







###subclusters 

#astrocytes
histo_all<-read.csv(paste0(sourcedir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_astrocytes_subset_pheno.csv"))

histo_all<- histo_all[,c(1:13,15,16,14)]


histo<- histo_all  %>% slice_head(n = 4)
#hippocampal volume, AT8 area, Iba1 area, CD68 area, GFAP area, and S100b area
histo_ordered<- rbind(histo[c(1,2,3,4),])

#without rowmeans across clusters normalization
histo_no_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)])
colnames(histo_no_norm)[1]<-"histo_parameters"


#with rowmeans across clusters normalization
histo_rowMeans_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)] - rowMeans(histo_ordered[2:ncol(histo_ordered)], na.rm=T))
colnames(histo_rowMeans_norm)[1]<-"histo_parameters"


histo_no_norm_cluster_of_interest <-  histo_no_norm %>% select(-histo_parameters)
histo_rowMeans_norm_cluster_of_interest <-  histo_rowMeans_norm %>% select(-histo_parameters)


#clustering and heatmap
#Without rowmeans normalization

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 3
data <- histo_no_norm_cluster_of_interest
set.seed(123)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss

#plot(1:k.max, wss,
#     type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")





##perform k-means clustering with k=3
kNClust <- 3
kmeanClust <- kmeans(histo_no_norm_cluster_of_interest,kNClust)

row.names(histo_no_norm_cluster_of_interest) <- histo_no_norm$histo_parameters
no_norm_histo_kmeans <- histo_no_norm_cluster_of_interest #[order(kmeanClust$cluster),]


kGaps <- vector(mode = "numeric")
kGaps[1] <- kmeanClust$size[1] + 1
for(i in 2:(kNClust-1)) {
  kGaps[i] <- kGaps[i-1] + kmeanClust$size[i]
}

simpleredbluecols = colorRampPalette(c("blue","white","red"))(400)

df <- data.frame(histo_parameters=histo_no_norm$histo_parameters)
row.names(df) <- histo_no_norm$histo_parameters


histo_no_norm_cluster_of_interest[1,13]<- " -1"
 histo_no_norm_cluster_of_interest$Cluster13<- as.numeric(histo_no_norm_cluster_of_interest$Cluster13)
 summary(histo_no_norm_cluster_of_interest)



paletteLength <- 400
##center the color scale so that while represents 0 in the row-centered matrix
myBreaks <- c(seq(min(histo_no_norm_cluster_of_interest, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(histo_no_norm_cluster_of_interest, na.rm = TRUE)/paletteLength, max(histo_no_norm_cluster_of_interest, na.rm=TRUE), length.out=floor(paletteLength/2)))


#for dendrogram
#mat_cluster_cols <- hclust(dist(t(no_norm_histo_kmeans)))
#sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#cmat_cluster_cols <- sort_hclust(mat_cluster_cols)
#mat_cluster_rows <- sort_hclust(hclust(dist(no_norm_histo_kmeans)))

cl<- data.frame(colnames(no_norm_histo_kmeans))
row.names(cl) <- colnames(no_norm_histo_kmeans)
no_norm_histo_kmeans

no_norm_histo_kmeans[1,13]<- " -1"
 no_norm_histo_kmeans$Cluster13<- as.numeric(no_norm_histo_kmeans$Cluster13)
 summary(no_norm_histo_kmeans)
 
 
pdf(paste0(outdir,"Fig.7g_all_astrocyte_clusters_pathology_15_clusters.pdf"), height=4, width=24)
print(pheatmap(no_norm_histo_kmeans, cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, treeheight_row=T,  gaps_row = F, annot_cols=df,annot_row=cl, color = simpleredbluecols, breaks = myBreaks, display_numbers = TRUE, number_format="%.3f",
         number_color = "black", 
         fontsize_number = 10)) 
dev.off()


###subclusters 
#microglia

histo_all<-read.csv(paste0(sourcedir,"log_odds_ratio_opt_bobyqa_per_unit_per_histopathology_per_subcluster_microglia_subset_pheno.csv"))
histo_all<- histo_all[,c(1:6,16,7:15)]

histo<- histo_all  %>% slice_head(n = 4)
#hippocampal volume, AT8 area, Iba1 area, CD68 area, GFAP area, and S100b area
histo_ordered<- rbind(histo[c(1,2,3,4),])

#without rowmeans across clusters normalization
histo_no_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)])
colnames(histo_no_norm)[1]<-"histo_parameters"


#with rowmeans across clusters normalization
histo_rowMeans_norm <- data.frame(histo_parameters=histo_ordered[1], histo_ordered[2:ncol(histo_ordered)] - rowMeans(histo_ordered[2:ncol(histo_ordered)], na.rm=T))
colnames(histo_rowMeans_norm)[1]<-"histo_parameters"


histo_no_norm_cluster_of_interest <-  histo_no_norm %>% select(-histo_parameters)
histo_rowMeans_norm_cluster_of_interest <-  histo_rowMeans_norm %>% select(-histo_parameters)


#clustering and heatmap
#Without rowmeans normalization

#Elbow Method for finding the optimal number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 3
data <- histo_no_norm_cluster_of_interest
set.seed(123)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss

#plot(1:k.max, wss,
#     type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")





##perform k-means clustering with k=3
kNClust <- 3
kmeanClust <- kmeans(histo_no_norm_cluster_of_interest,kNClust)

row.names(histo_no_norm_cluster_of_interest) <- histo_no_norm$histo_parameters
no_norm_histo_kmeans <- histo_no_norm_cluster_of_interest #[order(kmeanClust$cluster),]


kGaps <- vector(mode = "numeric")
kGaps[1] <- kmeanClust$size[1] + 1
for(i in 2:(kNClust-1)) {
  kGaps[i] <- kGaps[i-1] + kmeanClust$size[i]
}

simpleredbluecols = colorRampPalette(c("blue","white","red"))(400)

df <- data.frame(histo_parameters=histo_no_norm$histo_parameters)
row.names(df) <- histo_no_norm$histo_parameters





paletteLength <- 400
##center the color scale so that while represents 0 in the row-centered matrix
myBreaks <- c(seq(min(histo_no_norm_cluster_of_interest, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(histo_no_norm_cluster_of_interest, na.rm = TRUE)/paletteLength, max(histo_no_norm_cluster_of_interest, na.rm=TRUE), length.out=floor(paletteLength/2)))


#for dendrogram
#mat_cluster_cols <- hclust(dist(t(no_norm_histo_kmeans)))
#sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
#cmat_cluster_cols <- sort_hclust(mat_cluster_cols)
#mat_cluster_rows <- sort_hclust(hclust(dist(no_norm_histo_kmeans)))

cl<- data.frame(colnames(no_norm_histo_kmeans))
row.names(cl) <- colnames(no_norm_histo_kmeans)
no_norm_histo_kmeans



pdf(paste0(outdir,"Fig.8h_all_microglia_clusters_pathology_15_clusters.pdf"), height=4, width=24)
print(pheatmap(no_norm_histo_kmeans, cluster_rows=F, show_rownames=T, show_colnames=T, cluster_cols=F, treeheight_row=T,  gaps_row = F, annot_cols=df,annot_row=cl, color = simpleredbluecols, breaks = myBreaks, display_numbers = TRUE, number_format="%.3f",
         number_color = "black", 
         fontsize_number = 10))
dev.off()
