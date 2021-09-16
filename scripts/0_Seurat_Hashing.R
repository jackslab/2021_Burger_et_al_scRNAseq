library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)


# read in the data
hash.siin.data <- Read10X(data.dir = "cellRanger/Hash_new/siin/filtered_feature_bc_matrix") # See GEO
hash.siy.data <- Read10X(data.dir = "cellRanger/Hash_new/siy/filtered_feature_bc_matrix") # See GEO

# create Seurat objects with appropriate assay
hash.siin <- CreateSeuratObject(counts = hash.siin.data, project = "siin", assay = "HTO")
hash.siy <- CreateSeuratObject(counts = hash.siy.data, project = "siy", assay = "HTO")

# normalize data internally
hash.siin <- NormalizeData(hash.siin, assay="HTO", normalization.method="CLR")
hash.siy <- NormalizeData(hash.siy, assay="HTO", normalization.method="CLR")

# demultiplex
hash.siin <- HTODemux(hash.siin, assay="HTO", positive.quantile=0.99)
hash.siy <- HTODemux(hash.siy, assay="HTO", positive.quantile=0.99)

# return global classification results & summaries
write.csv(x = table(hash.siin$HTO_classification.global), paste0(od,day,"siin_hash_summary.csv"))
write.csv(x = table(hash.siy$HTO_classification.global), paste0(od,day,"siy_hash_summary.csv"))

write.csv(x = hash.siin$HTO_classification.global, paste0(od,day,"siin_hash_doubletcalls.csv"))
write.csv(x = hash.siin$hash.ID, paste0(od, day, "siin_hash_hashids.csv"))

write.csv(x = hash.siy$HTO_classification.global, paste0(od,day,"siy_hash_doubletcalls.csv"))
write.csv(x = hash.siy$hash.ID, paste0(od, day, "siy_hash_hashids.csv"))

# ridgeplots
Idents(hash.siin) <- "HTO_maxID"
pdf(file = paste0(od, day, "siin_hashing_ridge.pdf"), height = 40)
RidgePlot(hash.siin, assay="HTO", features = rownames(hash.siin), ncol = 1)
dev.off()

Idents(hash.siy) <- "HTO_maxID"
pdf(file = paste0(od, day, "siy_hashing_ridge.pdf"), height = 40)
RidgePlot(hash.siy, assay="HTO", features = rownames(hash.siy), ncol = 1)
dev.off()

# HTO signal pairing
pdf(file = paste0(od, day, "siin_hashing_scatter.pdf"))
FeatureScatter(hash.siin, feature1 = "hto_MB9452", feature2 = "hto_MB9458")
dev.off()

pdf(file = paste0(od, day, "siy_hashing_scatter.pdf"))
FeatureScatter(hash.siy, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
dev.off()
                    