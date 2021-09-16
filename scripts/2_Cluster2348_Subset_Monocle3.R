# dependencies
library(monocle3)
library(Seurat)
library(dplyr)

# load data
lung <- readRDS(file="200918_2_end_seurat.rds") # from 1_Preprocessing_DimRed.R

# set output directories
od <- "/Users/amandacruz/Dropbox (MIT)/SIIN vs SIY 10X/Final_Analysis/Monocle3/"
day <- Sys.Date()
dir.create(paste0(od,day,"/"))
od <- paste0(od,day,"/")


# format data to monocle3 object
expression.matrix <- lung@assays$RNA@counts # exp matrix
cell.metadata <- data.frame(lung@meta.data) # cell names
gene.annotation <- data.frame(rownames(expression.matrix),row.names=rownames(expression.matrix)) # gene names
gene.annotation[['gene_short_name']] = gene.annotation[['rownames.expression.matrix.']] # short gene name req

cds <- new_cell_data_set(expression.matrix,
                         cell_metadata = cell.metadata,
                         gene_metadata = gene.annotation) # construct cds object

# preprocessing & dim red
cds <- preprocess_cds(cds, num_dim = 30,preprocess_method = "LSI")
cds <- reduce_dimension(cds, reduction_method="UMAP",cores=8,preprocess_method = "LSI",umap.n_neighbors = 30)

# plotting 
pdf(file=paste0(od,day,"umap_seurat_clusters.pdf")) # old Seurat clusters
plot_cells(cds, color_cells_by="RNA_snn_res.0.7",cell_size = 1)
dev.off()

# clustering
cds = cluster_cells(cds, resolution=.001)
pdf(file=paste0(od,day,"umap_monocle_clusters.pdf")) # old Seurat clusters
plot_cells(cds,cell_size = 1)
dev.off()

pdf(file=paste0(od,day,"umap_antigen.pdf")) # antigen
plot_cells(cds,
           color_cells_by='orig.ident',
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 1)
dev.off()

pdf(file=paste0(od,day,"umap_seurat_clusters.pdf")) # old Seurat clusters
plot_cells(cds, color_cells_by="RNA_snn_res.0.4",cell_size = 1, show_trajectory_graph = FALSE)
dev.off()

pdf(file=paste0(od,day,"umap_monocle_clusters.pdf")) # old Seurat clusters
plot_cells(cds,cell_size = 1)
dev.off()

pdf(file=paste0(od,day,"umap_antigen.pdf")) # antigen
plot_cells(cds,
           color_cells_by='orig.ident',
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 1)
dev.off()

# trajectory
cds <- learn_graph(cds,use_partition = FALSE)

saveRDS(cds,paste0(od,day,"cds_monocle3.rds"))

# export coords
low_dim_coords <- reducedDims(cds)[["UMAP"]] # get UMAP coords from Monocle3
write.csv(paste0(od,"201030_subclustered_monocle3_umapCoords.csv"))