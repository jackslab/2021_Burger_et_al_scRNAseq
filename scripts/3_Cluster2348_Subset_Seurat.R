# set up output directory
od <- paste0("Final_Analysis/",Sys.Date(),"/") # this is where the output directory setup will be made
dir.create(od) # run this line if you havent already made the output directory folder
day <- paste0(Sys.Date(),"_") # this will automatically detect the date you are running the code

# Load dependencies
library(patchwork)
library(viridis)

## Cluster4-Cluster8 Subclustering##
dir.create(paste0(od,"Subclustering/"))
od <- paste0(od,"Subclustering/")

lung <- readRDS(file="200918_end_seurat.rds") # end result Seurat object of 1_Preprocessing_DimRed.R

# Subset Data
lung=SetIdent(lung,value='RNA_snn_res.0.7')
lung2 <- subset(lung,idents=c("4","8","2","3")) # subset clusters 2,3,4,8
lung2 <- FindVariableFeatures(lung2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(lung2), 20)

# subclustering
lung2 <- FindNeighbors(lung2, dims = 1:15)
lung2 <- FindClusters(lung2, resolution = 0.5)

# keep old and new clusters
lung2[['cluster']] = lung2@meta.data[["RNA_snn_res.0.7"]]

# dont forget to save
