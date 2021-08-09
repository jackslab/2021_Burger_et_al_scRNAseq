# Utility functions for producing multiple plots

###############################################################################
UMAP_DimPlot <- function (Seurat.obj, with.and.withoutlabels = FALSE,
                                 group.colors = NULL, cell.groups = NULL, 
                                 plot.order = NULL, extra.params = NULL, outfile) {
  # extra.params: any other parameters that DimPlot (Seurat) takes
  
  # Create two pdfs - one without labels and one with labels.
  if (with.and.withoutlabels == TRUE)  {
    dir.create(dirname(outfile), showWarnings = FALSE)
    pdf(outfile)
    print( do.call(DimPlot, c(list(Seurat.obj, label = FALSE, reduction = "umap",
                                    cols = group.colors, group.by = cell.groups,
                                    order = plot.order), extra.params) ) )
    print( do.call(DimPlot, c(list(Seurat.obj, label = TRUE, reduction = "umap",
                                   cols = group.colors, group.by = cell.groups,
                                   order = plot.order), extra.params)) )
    dev.off()
    
  # Create a single pdf - without labels.
  } else { 
    dir.create(dirname(outfile), showWarnings = FALSE)
    pdf(outfile)
    print( do.call(DimPlot, c(list(Seurat.obj, label = FALSE, reduction = "umap",
                                   cols = group.colors, group.by = cell.groups,
                                   order = plot.order), extra.params) ) )
    dev.off()
  }
  
}

###############################################################################
tSNE_DimPlot <- function (Seurat.obj, with.and.withoutlabels = FALSE,
                                 group.colors = NULL, cell.groups = NULL,
                                 plot.order = NULL, extra.params = NULL, outfile) {
  # extra.params: any other parameters that DimPlot (Seurat) takes
  
  # Create two pdfs - one without labels and one with labels.
  if (with.and.withoutlabels == TRUE) {
    dir.create(dirname(outfile), showWarnings = FALSE)
    pdf(outfile)
    print( do.call(DimPlot, c(list(Seurat.obj, label = FALSE, reduction = "tsne", 
                           cols = group.colors, group.by = cell.groups, 
                           order = plot.order), extra.params)) )
    print( do.call(DimPlot, c(list(Seurat.obj, label = TRUE, reduction = "tsne", 
                                   cols = group.colors, group.by = cell.groups, 
                                   order = plot.order), extra.params)) )
    dev.off()
    
  # Create a single pdf - without labels.
  } else { 
    dir.create(dirname(outfile), showWarnings = FALSE)
    pdf(outfile)
    print( do.call(DimPlot, c(list(Seurat.obj, label = FALSE, reduction = "tsne", 
                                   cols = group.colors, group.by = cell.groups, 
                                   order = plot.order), extra.params)) )
    dev.off()
  }
  
}

