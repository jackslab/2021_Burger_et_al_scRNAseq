# Cell state annotation of mouse CD8 T cells using the ProjecTILs pipeline.

# Dependencies:
# r/3.6.0
# Seurat/3.2.2
# ProjecTILs 0.5.1
# Atlases downloaded with ProjecTIL pipeline.
# Query mouse object provided.

library(ProjecTILs)
library(Seurat)
library(ggplot2)

# path.to.atlas <- "ref_TIL_Atlas_mouse_v1.rds"
# path.to.atlas <- "ref_LCMV_Atlas_mouse_v1.rds"
# path.to.query <- "200918_end_seurat.rds
# sample.name <- "mouse_CD8_Tcells"
# outdir <- "results_directory"

dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)
source("utils.R") 

###############################################################################
# Run projection algorithm

ref <- load.reference.map(ref = path.to.atlas)

query.obj <- readRDS(file = path.to.query)

query.projected <- make.projection(query.obj, ref=ref,
                                   filter.cells = T,
                                   query.assay = "RNA",
                                   skip.normalize = TRUE)

query.projected <- cellstate.predict(ref=ref, query=query.projected)

###############################################################################
# Map the subset of cells in projected query back to the original query set of cells
# and assign "NaN" to cells without label.

cells.all <- colnames(query.obj)
cells.kept <- colnames(query.projected)
cells.kept <- gsub("^Q_","",cells.kept)
I <- sapply(cells.kept, function(z) which(cells.all == z))

cell.state.labels <- query.projected@meta.data[["functional.cluster"]]
cell.state.idents <- cells.all
cell.state.idents[I] <- cell.state.labels
cell.state.idents[-I] <- "NA"
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents,
                         col.name = "cell.state.idents")

pred.conf <- query.projected@meta.data[["functional.cluster.conf"]]
conf <- cells.all
conf[I] <- pred.conf
conf[-I] <- "NA"
query.obj <- AddMetaData(object = query.obj, metadata = conf, col.name = "conf")

res <- cbind(cell.state.labels, conf)
write.table(res, file = paste0(outdir, "/", sample.name, "_cell_state_labels_wConf.xls", sep = "\t")

###############################################################################
# Highlight individual cell state labels on original query UMAP (Figure S2A)

create_UMAP_DimPlot(query.obj, with.and.withoutlabels = TRUE,
                    cell.groups = "cell.state.idents",
                    outfile = paste0(outdir, "/top_lvl/", sample.name, "_UMAP_cell_state_idents.pdf"))

dir.create(file.path(paste0(outdir, "/per_cell_state_highlights")), showWarnings = FALSE, recursive = TRUE)

cell.state.idents[-I] <- "Other"
cell.state.idents <- as.factor(cell.state.idents)
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents, col.name = "cell.state.idents")
cell.states <- levels(cell.state.idents)
cell.states <- cell.states[!cell.states %in% "Other"]
group.colors <- c("red", "grey")

for (cellstate in cell.states) {
  levels(query.obj@meta.data[["cell.state.idents"]])[which(levels(query.obj@meta.data[["cell.state.idents"]]) != cellstate)] <- "Other"
  names(group.colors) <- c(cellstate, "Other")
  UMAP_DimPlot(query.obj, with.and.withoutlabels = FALSE,
                      cell.groups = "cell.state.idents", group.colors = group.colors,
                      plot.order = cellstate,
                      outfile = paste0(outdir, "/per_cell_state_highlights/", sample.name, "_UMAP_cell_state_idents_", cellstate, ".pdf"))
  query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents,
                           col.name = "cell.state.idents")
}

###############################################################################
# Highlight Tpex and Tex on original Seurat UMAP (Figure 2D)

# TIL
I <- sort(union(which(cell.state.idents == "CD8_Tpex"), which(cell.state.idents == "CD8_Tex")))
cell.state.idents[-I] <- "Other"
table(cell.state.idents)
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents, col.name = "cell.labels")
table(query.obj$cell.labels)
ccols <- c("orange","purple","gray") # order of table of idents: CD8_Tex, CD8_Tpex, Other
pdf(file = paste0(outdir, "/", sample.name, "_TIL_Tpex_Tex_UMAP_highlight_plots_unordered_UMAP.pdf")
DimPlot(query.obj, reduction = "umap", group.by = "cell.labels", cols = ccols)
dev.off()

# LCMV
# I <- sort(union(which(cell.state.idents == "Tpex"), which(cell.state.idents == "Tex")))
# cell.state.idents[-I] <- "Other"
# table(cell.state.idents)
# query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents,
#                        col.name = "cell.labels")
# table(query.obj$cell.labels)
# ccols <- c("orange","purple","gray") # order of table of idents: CD8_Tex, CD8_Tpex, Other
# pdf(file = paste0(outdir, "/", sample.name, "_LCMV_Tpex_Tex_UMAP_highlight_plots_unordered_UMAP.pdf")
# DimPlot(query.obj, reduction = "umap", group.by = "cell.labels", cols = ccols)
# dev.off()

###############################################################################
# Predicted states overlay on UMAP - showing only labels w/ pred.conf > 0.5.

I <- which(pred.conf < 0.5)
pred.conf[I] <- NaN
print(sprintf("Number of cells with a confidence score < 0.05: %d", length(which(is.nan(pred.conf)))))

cell.state.labels <- query.projected@meta.data[["functional.cluster"]]
cell.state.labels[I] <- "NA"

# Map the subset of cells in projected query back to the original query set of cells
# and assign "NaN" to cells without a confidence score.
cells.all <- colnames(query.obj)
cells.kept <- colnames(query.projected)
cells.kept <- gsub("^Q_","",cells.kept)
I.kept <- sapply(cells.kept, function(z) which(cells.all == z)) # named integer
query.pred.conf <- cells.all
query.pred.conf[I.kept] <- pred.conf
query.pred.conf[-I.kept] <- "NaN"
query.pred.conf <- as.numeric(query.pred.conf)

query.obj <- AddMetaData(query.obj, metadata=query.pred.conf, col.name = "functional.cluster.conf")

outfile = paste0(outdir, "/", sample.name, "_", "functional_cluster_conf_UMAP_highlight.pdf")
pdf(outfile)
print( FeaturePlot(query.obj, reduction = "umap",
            features = "functional.cluster.conf") &
  scale_colour_gradientn(colours = c("#0571B0", "#92C5DE", "#D3D3D3", "#F4A582","#CA0020")) )
dev.off()

cell.state.idents <- cells.all
cell.state.idents[I.kept] <- cell.state.labels

cell.state.idents[-I.kept] <- "filtered"
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents, col.name = "cell.state.idents")
create_UMAP_DimPlot(query.obj, with.and.withoutlabels = TRUE,
                    cell.groups = "cell.state.idents",
                    outfile = paste0(outdir, "/top_lvl/", sample.name, "_UMAP_cell_state_idents_NAs_notlumped.pdf"))

cell.state.idents[-I.kept] <- "NA"
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents,
                         col.name = "cell.state.idents")
create_UMAP_DimPlot(query.obj, with.and.withoutlabels = TRUE,
                    cell.groups = "cell.state.idents",
                    outfile = paste0(outdir, "/top_lvl/", sample.name, "_UMAP_cell_state_idents_NAs_lumped.pdf"))

dir.create(file.path(paste0(outdir, "/per_cell_state_highlights")), showWarnings = FALSE, recursive = TRUE)
cell.state.idents[-I] <- "Other"
cell.state.idents <- as.factor(cell.state.idents)
query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents, col.name = "cell.state.idents")
cell.states <- levels(cell.state.idents)
cell.states <- cell.states[!cell.states %in% "Other"]
group.colors <- c("red", "grey")
for (cellstate in cell.states) {
  levels(query.obj@meta.data[["cell.state.idents"]])[which(levels(query.obj@meta.data[["cell.state.idents"]]) != cellstate)] <- "Other"
  names(group.colors) <- c(cellstate, "Other")
  create_UMAP_DimPlot(query.obj, with.and.withoutlabels = FALSE,
                      cell.groups = "cell.state.idents", group.colors = group.colors,
                      plot.order = cellstate,
                      outfile = paste0(outdir, "/per_cell_state_highlights/",
                                       sample.name, "_UMAP_cell_state_idents_", cellstate, ".pdf"))
  query.obj <- AddMetaData(object = query.obj, metadata = cell.state.idents,
                           col.name = "cell.state.idents")
  
}
