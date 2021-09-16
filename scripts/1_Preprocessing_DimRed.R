### --- Setup --- ###

# Load in package dependencies:
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# set up output directory
od <- paste0("Final_Analysis/",Sys.Date(),"/") # this is where the output directory setup will be made
dir.create(od) # run this line if you havent already made the output directory folder
day <- paste0(Sys.Date(),"_") # this will automatically detect the date you are running the code

# input filepaths (to download data, see GEO accession in README)
siy.path <- "Data/raw_siy_filteredbcmatrix"
siin.path <- "Data/raw_siin_filteredbcmatrix"

# plotting aesthetics
antigens <- c("SIIN", "SIY")
an_colors = c(siin = "#F8766D",siy = "#00BFC4")

# Load Data 
siy.data <- Read10X(data.dir = siy.path)
siin.data <- Read10X(data.dir = siin.path)

# Create Seurat Objects
siin <- CreateSeuratObject(counts = siin.data, 
                           project = "siin",
                           min.cells = 3,
                           min.features = 5)
siy <- CreateSeuratObject(counts = siy.data, 
                          project = "siy", 
                          min.cells = 3, 
                          min.features = 5)

### --- Individual Library QC --- ###
# This will calculate various QC metrics for each library individually. This information is utilized to compare the differences between each library. 
# For informational purposes only

# define output subdirectory
dir.create(paste0(od,"Pre-QC/"))
dir.create(paste0(od,"Pre-QC/Individual_libraries/"))
od.qc.indiv <- paste0(od,"Pre-QC/Individual_libraries/") 

# quantify mitochondrial genes
siin[["percent.mt"]] <- PercentageFeatureSet(siin, pattern = "^mt-") 
siy[["percent.mt"]] <- PercentageFeatureSet(siy, pattern = "^mt-") 

# violin plots of QC metrics (preQC)
qc.indiv.vp <- list()
qc.indiv.vp[[1]] <- VlnPlot(siin, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            ncol = 3, pt.size = 0.25,
                            cols = an_colors[1])
qc.indiv.vp[[2]] <- VlnPlot(siy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            ncol = 3, 
                            pt.size = 0.25, 
                            cols = an_colors[2])

# scatter plots of features (preQC)
qc.indiv.sp <- list()
qc.indiv.sp[[1]] <- FeatureScatter(siin, 
                                   feature1 = "nCount_RNA", 
                                   feature2 = "nFeature_RNA", 
                                   cols = an_colors[1]) + ggtitle("Reads vs Genes Detected")
qc.indiv.sp[[2]] <- FeatureScatter(siy, 
                                   feature1 = "nCount_RNA", 
                                   feature2 = "nFeature_RNA", 
                                   cols = an_colors[2]) + ggtitle("Reads vs Genes Detected")

# save plots
pdf(file = paste0(od.qc.indiv, day, "individual_qc_violinPlots.pdf"), height = 14, width = 9)
wrap_plots(qc.indiv.vp, ncol = 1, nrow = 2)
dev.off()

pdf(file = paste0(od.qc.indiv, day, "individual_qc_scatterPlots_counts-vs-genes.pdf"), height = 14, width = 9)
wrap_plots(qc.indiv.sp, ncol = 1, nrow = 2)
dev.off()

# filter out cells that do not pass QC thresholds
siin <- subset(siin, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 20000)
siy <- subset(siy, subset = nFeature_RNA > 100 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 20000)


### --- Preliminary Individual Library Analysis --- ###
# for informational purposes only

# normalize data 
siin <- NormalizeData(siin, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000, 
                      vars.to.regress = "percent.mt")
siy <- NormalizeData(siy, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     vars.to.regress = "percent.mt")

# variable feature detection 
siin <- FindVariableFeatures(siin, 
                             selection.method = "vst", 
                             nfeatures = 2000)
siy <- FindVariableFeatures(siy, 
                            selection.method = "vst", 
                            nfeatures = 2000)

top10.siin <- head(VariableFeatures(siin), 10)
top10.siy <- head(VariableFeatures(siy), 10)


# plot variable features with and without labels
vfp <- VariableFeaturePlot(siin)
vfp.l <- LabelPoints(plot = vfp, 
                     points = top10.siin, 
                     repel = TRUE)
pdf(file = paste0(od.qc.indiv,day,"siin_variablegenes.pdf"), width = 9, height = 7)
vfp + vfp.l
dev.off()

vfp2 <- VariableFeaturePlot(siy)
vfp2.l <- LabelPoints(plot = vfp2, points = top10.siy, repel = TRUE)
pdf(file = paste0(od.qc.indiv, day, "siy_variablegenes.pdf"), width = 9, height = 7)
vfp2 + vfp2.l
dev.off()

### --- Merge Libraries --- ###
# combine into one seurat object
siin <- ScaleData(siin, vars.to.regress = "percent.mt")
siy <- ScaleData(siy, vars.to.regress = "percent.mt")

cd8 <- merge(siin, y = siy, 
             add.cell.ids = c("siin", "siy"), 
             project = "kplung_cd8", 
             merge.data = TRUE)

rm(siin.data, siy.data) # clear  from memory - you won't need this anymore

### --- Merge Library QC & Filtering --- ###

# output subdirectory
dir.create(paste0(od,"Post-QC/"))
dir.create(paste0(od,"Post-QC/merged/"))
od.qc.m <- paste0(od,"Post-QC/merged/") 

# save file
pdf(file = paste0(od.qc.m,day,"merge_postQC_violin.pdf"), width = 9, height = 7)
VlnPlot(cd8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25)
dev.off()

# Add in VDJ & Hashing Data
# At this point, as part of the QC filtering, we would like to filter out empirically determined doublets. 
# These doublets are detected using TotalSeqC antibody hashing. 
# The premise is that cells with a hashing ID that come from more than one mouse are likely doublets. 
# This processing uses a combination of Python and R, broken up in multiple notebooks.

# Export cell IDs to CSV:
dir.create(paste0(od,"Metadata/"))
od.md <- paste0(od,"Metadata/")
write.csv(cd8@meta.data, file = paste0(od.md,day,"md_prefilter.csv"))

## -- Merge Library QC & Filtering -- ##

# At this point, we will use python to fill in information about cells to add to the metadata, 
# add in VDJ information, and add in Hashing. 
# The Hashing analysis/processing is covered in a separate notebook. 
# This can be done in R, it was just easier to do it this way at the time we did this analysis. 
# Please see file "AddHashingVDJ.ipynb"

# Load metadata: 
md <- read.csv(file = paste0(od.md, day, "md_posthash.csv"))
cd8[['demux']] = md[['demux']]
md[['top_clonotypes']] = factor(md[['top_clonotypes']])
cd8 = RenameCells(cd8, new.names = md[['barcode']]) 

# Remove Doublets
cd8 = SetIdent(cd8, value='demux')
not_doublets <- WhichCells(object = cd8, idents = 'Doublet', invert = TRUE)
cd8 <- cd8[,not_doublets]

write.csv(cd8@meta.data,paste0(od.md,day,"md_doublets_out.csv"))
# At this point, we will use another python script (see Clonotype-FinalMetadataProcessing.ipynb notebook)
# to import ProjectTIL metadata and perform clonotype calculations after filtering out the doublets.
# This will be the final metadata that doesn't include calculations made by Seurat and other downstream analyses. 

# Load in final metadata:
md <- read.csv(file = paste0(od.md, day, "md_final.csv"))
md = data.frame(md)

# add as metadata; this was glitchy so its not done in a standard way
for (i in colnames(md)) {
  cd8[[i]] = md[[i]]
}

# Save
saveRDS(cd8,paste0(od,day,"preprocessed_norm_seurat.rds"))

### -- DimRed -- ###

# Feature Selection & Scaling
all.genes <- rownames(cd8)
cd8 <- FindVariableFeatures(cd8,
                            selection.method = "vst",
                            nfeatures = 2000)
cd8 <- ScaleData(cd8, features = all.genes)

# Load in final metadata:
md <- read.csv(file = paste0(od.md, day, "md_final.csv"))
md = data.frame(md)

# add as metadata; this was glitchy so its not done in a standard way
for (i in colnames(md)) {
  cd8[[i]] = md[[i]]
}

# PCA
cd8 <- RunPCA(cd8,  features = VariableFeatures(object = cd8))
dir.create(paste0(od,"dimensionality_reduction/"))
od.dr <- paste0(od,"dimensionality_reduction/")

# plot first couple PCs and visualize dimensions
pdf(file=paste0(od.dr,day,"merge_PCA1.pdf"), width=9, height=7)
DimPlot(cd8, reduction = "pca")
dev.off()

pdf(file=paste0(od.dr,day,"PCA1_heatmap.pdf"), width=9, height=7)
DimHeatmap(cd8, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(file=paste0(od.dr,day,"PCAplus_heatmap.pdf"), width=9, height=7)
DimHeatmap(cd8, dims = 1:5, cells = 500, balanced = TRUE)
dev.off()

# plot PCs in the upper range of the predicted # of PCs by JackStraw/Elbow plots for downstream dimensionality reductions 
pdf(file=paste0(od.dr,day,"PCA25-35_heatmap.pdf"), width=9, height=7)
DimHeatmap(cd8, dims = 25:35, cells = 500, balanced = TRUE)
dev.off()

### -- Cluster & Viz -- ###

cd8 <- FindNeighbors(cd8, dims = 1:30)
cd8 <- FindClusters(cd8, resolution = 0.7)
cd8 <- RunUMAP(cd8, dims = 1:30, seed.use = 42)
DimPlot(cd8, reduction = "umap")

# dont forget to save!

