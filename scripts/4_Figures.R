# Dependencies
library(Seurat)
library(ggplot2)
library(grid)
library(ArchR)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(philentropy)
library(plyr)
library(tidyr)
library(dplyr)
library(cowplot)
library(gridExtra)
library(stringr)
library(monocle3)
library(dendextend)

# load data
lung <- readRDS(file="200918_2_end_seurat.rds") # from 1_Preprocessing_DimRed.R
lung2 <- readRDS(file="201009_subclustered_rds.rds") # from 2_Cluster2348_Subset_Seurat.R

# change UMAPs for subplot
lung2 = SetIdent(lung2, value = "RNA_snn_res.0.7")
m3_umap = read.csv(file="data/201030_subclustered_monocle3_umapCoords.csv")
names(m3_umap)[names(m3_umap) == "V1"] <- "UMAP_1"
names(m3_umap)[names(m3_umap) == "V2"] <- "UMAP_2"
as.matrix(m3_umap)
lung2@reductions$umap@cell.embeddings[,1] = m3_umap$UMAP_1
lung2@reductions$umap@cell.embeddings[,2] = m3_umap$UMAP_2

# add sigs
TIL <- read.csv(file = "mouse_samples_per_cell_attributes.csv") # formatted as CSV
rownames(TIL) = TIL$Cell.ID
lung@meta.data$TIL_miller = TIL$query.projected.functional.cluster
lung@meta.data$TIL_TILPRED = TIL$query.projected.TILPRED


# set output directories
od <- "Final_Analysis/Figures/"
day <- Sys.Date()
dir.create(paste0(od,day,"/"))
od <- paste0(od,day,"/")

# aesthetics
sc.cols = c('0' = "#56251a",
            '1' = "#DC267F",
            '2' = "#e62325",
            '3' = "#0072B2",
            '4' = "#648FFF",
            '5' =  "#44AA99",
            '6' = "#FFB000",
            '7' = "#fc835c",
            '8' =  "#a91560",
            '9' = "#830ebe",
            '10' = "#293352")
featcmap <- c("gray75","gray87","gray89","palegoldenrod","lightgoldenrod","goldenrod","darkorange3","brown")
an_colors = c("#F8766D","#00BFC4")


#### FIGURE 2A ####
lung <- SetIdent(lung, value = 'orig.ident')
p <- DimPlot(lung, cols = an_colors, pt.size = 0.05) + theme_void()
ggsave2(filename = paste0(od,"fig2a_umap",".pdf"),
        plot = (p + theme(legend.position = 'none')),
        device = "pdf",
        units = "cm",
        width = 5.5,
        height = 5.5,
        scale = 1)

#### FIGURE 2B ####
lung <- SetIdent(lung, value = 'RNA_snn_res.0.7')

p <- DimPlot(lung, cols = sc.cols, pt.size = 0.05) +
  theme_void() +
  theme(legend.text = element_text(size = 8, face = 'bold'),
        legend.key.width = unit(0.25, units = 'cm'),
        legend.key.height = unit(0.4, units = 'cm'))

ggsave2(filename = paste0(od,"fig2b_umap",".pdf"),
        plot = (p + theme(legend.position = 'none')),
        device = "pdf",
        units = "cm",
        width = 5.5,
        height = 5.5,
        scale = 1)

# legend
pl.l <- cowplot::get_legend(p)
pdf(file = paste0(od,"fig2b_umap_legend",".pdf"), width = 0.5, height = 2)
grid.newpage()
grid.draw(pl.l)
dev.off()

#### FIGURE S2B ####
# Miller
signatures <- read.csv(file = "data/200922_miller.csv") # formatted as CSV
signatures.id <- colnames(signatures)
for (i in signatures.id) {
  signatures[[i]] <- str_to_title(signatures[[i]])} # change to title case for mouse gene compatibility
exp.gene <- rownames(lung)

# Miller
for (i in signatures.id) {
  #sig.id <- signatures.id[i] # get signature name 
  signatures.gene <- signatures[i] # get gene lit subset
  signatures.gene <- signatures.gene[!apply(signatures.gene == "", 1, all), ] # get rid of empty variables
  signatures.gene <- intersect(signatures.gene, exp.gene)
  #signatures.gene <- subset(signatures.gene, signatures.gene = exp.gene)
  #signatures.gene <- list(signatures.gene)
  lung = AddModuleScore(
    object = lung,
    features = list(signatures.gene),
    name = paste0(i,"_score"))

   p <- FeaturePlot(object=lung, reduction="umap", features=paste0(i,"_score1"), order = T, pt.size = 0.25) + ggtitle(paste0("Signature:",i))   &
   theme_void() &
     theme(legend.position = 'right') &
     labs(title = NULL) &
     theme(plot.tag = NULL, plot.subtitle = NULL, strip.text.x = element_text(size=8, face = "bold"),
           legend.title = element_blank(),
           legend.key.height = unit(0.25,"cm"),
           legend.text = element_text(size=6, face = "bold"),
           legend.key.width = unit(0.5,"line"),
           legend.direction = "horizontal") &
   scale_colour_gradientn(colours =
                            c("gray75","palegoldenrod","goldenrod","darkorange3","firebrick"),
                          labels = c(0, 0.5, 1),
                          breaks = c(0, 0.5 ,1),
                          limits = c(0,1))
   ggsave2(filename = paste0(od,"figS2B_umap_",i,".pdf"),
           plot = (p + theme(legend.position = 'none')),
           device = "pdf",
           units = "cm",
           width = 3.5,
           height = 3.5,
           scale = 1)
  # legend
   pl.l <- cowplot::get_legend(p)
   pdf(file = paste0(od,"figS2B_umap_legend_",i,".pdf"), width = 1, height = 1.5)
   grid.newpage()
   grid.draw(pl.l)
   dev.off()
}

# stats
scores = list(Tumor_Progenitor = lung@meta.data$tumor_progenitor_score1,
              Tumor_Terminal = lung@meta.data$tumor_terminal_score1,
              LCMV_Progenitor = lung@meta.data$lcmv_progenitor_score1,
              LCMV_Terminal = lung@meta.data$lcmv_terminal_score1)
scores_names = names(scores)
scores = as_tibble(scores)

scores = scores %>% mutate(zscore_Tumor_Progenitor = (Tumor_Progenitor- mean(Tumor_Progenitor))/sd(Tumor_Progenitor))
scores = scores %>% mutate(zscore_Tumor_Terminal = (Tumor_Terminal- mean(Tumor_Terminal))/sd(Tumor_Terminal))
scores = scores %>% mutate(zscore_LCMV_Progenitor = (LCMV_Progenitor- mean(LCMV_Progenitor))/sd(LCMV_Progenitor))
scores = scores %>% mutate(zscore_LCMV_Terminal = (LCMV_Terminal- mean(LCMV_Terminal))/sd(LCMV_Terminal))

scores = scores %>% mutate(assign_Tumor_Progenitor = replace(zscore_Tumor_Progenitor, zscore_Tumor_Progenitor > 0.5, "Tumor Progenitor"))
scores = scores %>% mutate(assign_Tumor_Progenitor = replace(assign_Tumor_Progenitor, zscore_Tumor_Progenitor < 0.5, "False"))

scores = scores %>% mutate(assign_Tumor_Terminal = replace(zscore_Tumor_Terminal, zscore_Tumor_Terminal > 0.5, "Tumor Terminal"))
scores = scores %>% mutate(assign_Tumor_Terminal = replace(assign_Tumor_Terminal, zscore_Tumor_Terminal < 0.5, "False"))

scores = scores %>% mutate(assign_LCMV_Progenitor = replace(zscore_LCMV_Progenitor, zscore_LCMV_Progenitor > 0.5, "LCMV Progenitor"))
scores = scores %>% mutate(assign_LCMV_Progenitor = replace(assign_LCMV_Progenitor, zscore_LCMV_Progenitor < 0.5, "False"))

scores = scores %>% mutate(assign_LCMV_Terminal = replace(zscore_LCMV_Terminal, zscore_LCMV_Terminal > 0.5, "LCMV Terminal"))
scores = scores %>% mutate(assign_LCMV_Terminal = replace(assign_LCMV_Terminal, zscore_LCMV_Terminal < 0.5, "False"))

lung@meta.data$tumor_progenitor = scores$assign_Tumor_Progenitor
lung@meta.data$tumor_terminal = scores$assign_Tumor_Terminal
lung@meta.data$lcmv_progenitor = scores$assign_LCMV_Progenitor
lung@meta.data$lcmv_terminal = scores$assign_LCMV_Terminal

scores_names = c("tumor_progenitor",
                 "tumor_terminal",
                 "lcmv_progenitor",
                 "lcmv_terminal")

lung = SetIdent(lung,value='orig.ident')
md.qa = table(Idents(lung))
p.val = array(dim=c(2,4))
row.names(p.val) = row.names(md.qa)
colnames(p.val) = scores_names
for (i in seq(1,4)) {
  for (p in seq(1,2)) {
    name = scores_names[i]
    md.q = (table(Idents(lung), lung@meta.data[,name]))
    md.qp = prop.table(table(Idents(lung), lung@meta.data[,name]))
    x = md.q[p,2] # number of cells for a given antigen assigned to a given cluster
    m = md.qa[p] # number of cells for a given antigen
    k = sum(md.q[,2]) # number of cells assigned to a given cluster
    n = sum(md.qa) - m # number of cells assigned to another antigen
    p.val[p,i] = phyper(q=x-1,m=m,n=n,k=k,lower.tail=F)
  }
}

write.csv(p.val,file="210614_signature-antigen_hypergeotest.csv")


#### FIGURE 2E ####
fig2e = c('Tcf7','Il7r','Gzmb','Havcr2')
for (gene in fig2e) {
  p <- FeaturePlot(lung, features = gene, order = T, pt.size = 0.005) &
    theme_void() &
    theme(legend.position = 'right') &
    scale_colour_gradientn(colours = 
                             c("gray75","palegoldenrod","goldenrod","darkorange3","firebrick")) &
    labs(title = NULL) &
    theme(legend.text = element_text(size = 6, face = 'bold'),
          legend.key.width = unit(0.25, units = 'cm'),
          legend.key.height = unit(0.25, units = 'cm'))
  
  ggsave2(filename = paste0(od,"fig2e_umap_",gene,".pdf"),
          plot = (p + theme(legend.position = 'none')),
          device = "pdf",
          units = "cm",
          width = 3.5,
          height = 3.5,
          scale = 1)
  # legend
  pl.l <- cowplot::get_legend(p)
  pdf(file = paste0(od,"fig2e_umap_legend_",gene,".pdf"), width = 0.5, height = 1.5)
  grid.newpage()
  grid.draw(pl.l)
  dev.off()
}

#### FIGURE S2C ####
lung <- SetIdent(lung,value='RNA_snn_res.0.7')
cluster.markers.old <- FindAllMarkers(lung, min.pct=0.25)

colMeta_df.f <- data.frame(Cluster = lung@meta.data$RNA_snn_res.0.7)
rownames(colMeta_df.f) = colnames(lung)
top10 <- cluster.markers.old %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
exp = GetAssayData(lung,slot='scale.data')

genes = c('Rps19', 'Bcl2', 'Cdkn1a', 'Vim', 'Id2', 'Ctla2a', 'Ccl5', 'Cd52', 'S100a6', 'Cxcr6', 'Gzmk', 'Rgs16', 'Tnfrsf9', 'Lag3', 'Nr4a2', 'Xcl1', 'Tnfrsf4', 'Ramp1', 'Tmem176b', 'Asns', 'Mthfd2', 'Ddit3', 'Ahnak', 'Ptprc', 'Itgb1', 'AY036118', 'Trbv13-1', 'S1pr1', 'Klf2', 'Tcf7', 'IL7r', 'Ifit1', 'Isg15', 'Top2a', 'Tuba1b')
exp <- as.matrix(exp[rownames(exp) %in% top10$gene, ])
ha.cl <- HeatmapAnnotation(Cluster = colMeta_df.f$Cluster,
                           col = list(Cluster = sc.cols),
                           annotation_name_side = "right",
                           show_legend = F,
                           annotation_name_gp= gpar(fontsize = 6),
                           #width = unit(0.2,"mm"),
                           #gap = unit(0,"mm"),
                           annotation_label = c("Cluster"))
hal.cl = Legend(legend_gp = gpar(fill = sc.cols),
                title = "Cluster",
                at = colMeta_df.f$Cluster,
                labels = names(sc.cols),
                grid_width = unit(0.25, "cm"),
                grid_height = unit(0.2,"cm"),
                title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 6))
row.ord <- rownames(exp) %in% top10$gene
row.ord <- unique(top10$gene)
gene.loc <- match(genes, rownames(exp))
ha.gene <- rowAnnotation(Gene = anno_mark(at = gene.loc, labels = genes, which = "row", padding = unit(1, "mm"),
                                          labels_gp = gpar(fontsize = 6, fontface = 'bold'),extend = unit(0, "mm"),
                                          link_width = unit(4, "mm")),
                         annotation_name_side = "top",
                         show_legend = F,
                         annotation_name_gp= gpar(fontsize = 6))

p2 <- ComplexHeatmap::Heatmap(
  matrix = (scale((as.matrix(exp)))),
  column_title = "Cell",
  column_title_side = "bottom",
  name = "Expression",
  border = F,
  row_title = "Gene",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  show_column_dend = F,
  show_row_dend = F,
  cluster_rows = F,
  show_row_names = F,
  show_column_names = F,
  cluster_columns = F,
  row_order = row.ord,
  height = unit(140,"mm"),
  width = unit(75, "mm"),
  use_raster = T,
  column_split = colMeta_df.f$Cluster,
  column_gap = unit(1, "mm"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, face = "bold"),
                              labels_gp = gpar(fontsize = 8), face='bold',
                              legend_direction = "vertical",
                              legend_height = unit(2.5,"cm")),
  top_annotation = ha.cl,
  right_annotation = ha.gene
)
# legend
heat.legend <- Legend(col_fun =colorRamp2(breaks = c(-4,-2,0,2,4), colors = c("#0000FFFF", "#9267FCFF", "#EEEEEEFF", "#FF7C5CFF", "#FF0000FF")),
                      title = "Expression",
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      at = c(-4,-2,0,2,4),
                      labels = c(-4,-2,0,2,4),
                      labels_gp = gpar(fontsize = 6),
                      direction = "vertical",
                      legend_width = unit(2,"cm"),
                      grid_height = unit(0.25, "cm"),)

hal <- packLegend(hal.cl, heat.legend)
pdf(file=paste0(od,"FigureS2C_heat.pdf"),width = 5, height = 7)
draw((p2),show_heatmap_legend = F)
dev.off()
pdf(file=paste0(od,"FigureS2C_heatlegend.pdf"),width = 1, height = 3)
draw(hal)
dev.off()

#### FIG S2E ####
figS2E = c('AY036118', 'Cd44', 'Cd69')
for (gene in figS2E) {
  p <- FeaturePlot(lung, features = gene, order = T, pt.size = 0.005) &
    theme_void() &
    theme(legend.position = 'bottom') &
    scale_colour_gradientn(colours = 
                             c("gray75","palegoldenrod","goldenrod","darkorange3","firebrick")) &
    labs(title = NULL) &
    theme(legend.text = element_text(size = 6, face = 'bold'),
          legend.key.width = unit(0.25, units = 'cm'),
          legend.key.height = unit(0.25, units = 'cm'))
  
  ggsave2(filename = paste0(od,"figS2E_umap_",gene,".pdf"),
          plot = (p + theme(legend.position = 'none')),
          device = "pdf",
          units = "cm",
          width = 3.5,
          height = 3.5,
          scale = 1)
  # legend
  pl.l <- cowplot::get_legend(p)
  pdf(file = paste0(od,"figS2E_umap_legend_",gene,".pdf"), width = 0.5, height = 1.5)
  grid.newpage()
  grid.draw(pl.l)
  dev.off()
}

#### FIGURE 2C (STATISTICAL TEST ONLY) ####

## FIGURE 2C QUANT
lung = SetIdent(lung,value='orig.ident')
md.q = (table(Idents(lung), lung$RNA_snn_res.0.7))
md.qp = prop.table(table(Idents(lung), lung$RNA_snn_res.0.7))
md.qa = table(Idents(lung))
p.val = array(dim=c(2,11))
row.names(p.val) = row.names(md.qa)
colnames(p.val) = seq(0,10)
for (i in seq(1,11)) {
  for (p in seq(1,2)) {
    x = md.q[p,i] # number of cells for a given antigen assigned to a given cluster
    m = md.qa[p] # number of cells for a given antigen
    k = sum(md.q[,i]) # number of cells assigned to a given cluster
    n = sum(md.qa) - m # number of cells assigned to another antigen
    p.val[p,i] = phyper(q=x-1,m=m,n=n,k=k,lower.tail=F)
  }
}

write.csv(p.val,file="210104_cluster-antigen_hypergeotest.csv")
write.csv(md.q,file="210104_cluster-antigen_hypergeotest-md_q.csv")
write.csv(md.qa,file="210104_cluster-antigen_hypergeotest-md_qa.csv")
write.csv(md.qp,file="210104_cluster-antigen_proportions_mdqp.csv")

#### FIGURE S2C ####
exprs_mat <- t(as.matrix(lung@assays$RNA@data)) # get expression data for just that 
exprs_mat <- reshape2::melt(exprs_mat) # make matrix into long format rather than wide
colnames(exprs_mat) <- c("Cell", "Gene", "Expression") # define colnames of matrix
exprs_mat$Gene <- as.character(exprs_mat$Gene) # data type for gene name conversion
cell_group <- lung@meta.data$RNA_snn_res.0.7 # define cell groups by seurat clusters
names(cell_group) = names(lung@meta.data$orig.ident) # carry over cluster names and cell names
exprs_mat$Group <- cell_group[exprs_mat$Cell] # populate matrix with cluster ID for each cell
exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) ==
                                          FALSE) # remove NA values (shouldnt be any)
ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
  dplyr::summarize(
    mean = mean(log(Expression + pseudocount)),
    percentage = sum(Expression > lower_threshold) / length(Expression)
  ) # take summary for each cluster and calculate mean expression as defined by mean(log(expression +1))
# defines percentage as the sum of expression values that pass lower_threshold (by default zero)
ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min,
                      ExpVal$mean) # scaling
ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max,
                      ExpVal$mean) # scaling

inf = reshape::cast(ExpVal, Gene~Group, value = 'mean')

# setup
siin = SetIdent(siin,value = "raw_clonotype_id")
siy = SetIdent(siy, value = "raw_clonotype_id")
# get table of clonotype sizes
md.qa.siin = table(Idents(siin))
md.qa.siy = table(Idents(siy))
# remove counts for cells with no clonotype assignment
md.qa.siin = md.qa.siin[-1]
md.qa.siin = md.qa.siin[-1]

cM.siin = as_tibble(cM.siin)
cM.ft.siy = as_tibble(cM.ft.siy)

# flatten table for plotting
cM.siin.ecdf = melt(cM.siin)
cM.ft.siy.ecdf = melt(cM.ft.siy %>% select(c('C4C8', 'nC4C8')))

# add antigen to variable name
cM.ft.siy.ecdf$Antigen = paste0(cM.ft.siy.ecdf$variable, " SIY")
cM.siin.ecdf$Antigen = paste0(cM.siin.ecdf$variable, " SIIN")

# combine into one dataframe for plotting
cM.ecdf = rbind(cM.siin.ecdf, cM.ft.siy.ecdf)
cM.ecdf2 = cM.ecdf %>% filter(variable != 'nC4C8')

# plot
ggplot(cM.ecdf2, aes(value, color = Antigen)) +
  stat_ecdf(geom = "line") + stat_ecdf(geom = "point") +
  scale_color_manual(values =  c("#F8766D", "#00BFC4"),
                     labels = c("SIIN", "SIY")) +
  xlab("Proportion of Cells in an Individual Clonotype") +
  ylab("Proportion of Clonotypes In C4 and C8") +
  labs(color = "") +
  theme_classic() +
  xlim(0, 1)

#### FIGURE S2D ####
prog <- c('Slamf6','Xcl1','Sell','Ccr7','Id3','Lef1')
modules = c(prog)
for (i in modules) {
  for (gene in i) {
    p <- FeaturePlot(object=lung, reduction="umap", features=gene, order = T, pt.size = 0.05) &
      theme_void() &
      labs(title = NULL) &
      theme(plot.tag = NULL, plot.subtitle = NULL, strip.text.x = element_text(size=8, face = "bold"),
            legend.title = element_blank(),
            legend.key.height = unit(0.25,"cm"),
            legend.text = element_text(size=6, face = "bold"),
            legend.key.width = unit(0.5,"line"),
            legend.direction = "horizontal",
            plot.title = element_text(size=8,face="bold",hjust=0.5)) &
      ggtitle(gene) &
      scale_colour_gradientn(colours = 
                               c("gray75","palegoldenrod","goldenrod","darkorange3","firebrick"))
    
    ggsave2(filename = paste0(od,"figS2D_umap_",gene,".pdf"),
            plot = (p + theme(legend.position = 'none')),
            device = "pdf",
            units = "cm",
            width = 3.5,
            height = 3.5,
            scale = 1)
    # legend
    pl.l <- cowplot::get_legend(p)
    pdf(file = paste0(od,"fig2E_umap_legend_",i,".pdf"), width = 1, height = 1.5)
    grid.newpage()
    grid.draw(pl.l)
    dev.off()
  }
}


#### FIGURE 2F ####
bad_clonos <- read.csv(file='Data/200922_multiplemouse_clonotypes.csv') # From Clonotype-FinalMetadataProcessing.ipynb

## confusion matrix -- all clonotypes ###
cluster.clonotype <- table(lung$raw_clonotype_id) # table of number of cells in each clonotype
cM <- confusionMatrix(paste0(lung$raw_clonotype_id), paste0(lung$seurat_clusters)) # confusion matrix
cM <- cM[!rownames(cM) %in% bad_clonos$X0,] # filter out clonotypes with multiple mice
total = Matrix::rowSums(cM) # calculates total number of cells per clonotype
total.c4c8 = cM[,6] + cM[,7] # calculates total number of cells belonging to C3 and C2 (old 4,8)

## plotting conditions ##
an_colors = list(
  Antigen = c(siin = "#F8766D", siy = "#00BFC4"),
  Clonotype_Cluster_All = c(
    '1' = "#222222",
    '2' = "#F6222E",
    '3' = "#FE00FA",
    '4' = "#16FF32",
    '5' = "#3283FE",
    '6' = "#FEAF16",
    '7' = "#B00068",
    '8' = "#1CFFCE",
    '9' = "#90AD1C",
    '10' = "#2ED9FF",
    '11' = "#DEA0FD",
    '12' = "#AA0DFE"
  ),
  Mouse = c(
    "MB9452" = "#AA0DFE",
    "MB9457" = "#F6222E",
    "MB9458" = "#FE00FA",
    "MB9459" = "#16FF32",
    "MB9509" = "#3283FE",
    "MB9561" = "#FEAF16",
    "MB9619" = "#B00068",
    "MB9621" = "#DEA0FD",
    "MB9622" = "darkorange4",
    "MB9625" = "tan"
  ),
  Size = colorRamp2(c(5,10,88),colors=c("grey","goldenrod","firebrick")))
genes <- c("Gzmb",
           "Cx3cr1",
           "Il17a",
           "Lag3",
           "Tcf7",
           "Ccr6")
clonotypes = unique(lung$raw_clonotype_id)
clonotypes[0] = NULL
lung = SetIdent(lung,value="raw_clonotype_id")
exp = lung@assays$RNA@data[genes,]

sizecols = colorRamp2(c(1,10,88),
                      colors = c("palegoldenrod","lightgoldenrod","darkorange"))(seq(1,max(colMeta_df.f$Size)))

names(sizecols) = seq(1,length(sizecols))


# reorder confusion matrix to be in order of cluster number
phcl.order <- list(seq(0,10,by=1))
phcl.order <- as.character(phcl.order[[1]])

### calculations ###
filter.f = total > 5 # cell number minimum 
cM.f = cM[filter.f,] # filter out clonotypes that dont meet cell number cutoff
filter.f2 = cM.f[,6] > 1 | cM.f[,7] > 1 # include only clusters that have a cell in Cluster 4 or cluster 8
cM.f = cM.f[filter.f2,] # filter out non cluster 4 or 8 containing clonotypes
cM.ftot = Matrix::rowSums(cM.f)
cM.f <- cM.f / Matrix::rowSums(cM.f) # normalize
cM.f = cM.f[,phcl.order] # put in order

ec.d <- dist(cM.f)
phcl.hc <- hclust(ec.d)
phcl.hc2 <- cutree(phcl.hc, h = 0.5) # assigns clonotypes to clusters

write.csv(phcl.hc2, file=paste0(od,"pheatmap_clonotypes.csv"))

clmd <- read.csv(file="data/C4C8_clonotypes_metadata.csv")

colMeta_df.f <- data.frame(Antigen = sub(".*_", "", rownames(cM.f)),
                           #'Clonotype Cluster' = as.character(phcl.hc2),
                           Mouse = clmd$Mouse,
                           Size = cM.ftot,
                           stringsAsFactors = F,
                           row.names = rownames(cM.f)) # provide confusion matrix metadata with new clusters

### plotting ###
# heatmap annotations
ha.size <- row_anno_barplot(colMeta_df.f$Size,
                            gp = gpar(fill = sizecols))

ha.ant <- rowAnnotation(Antigen = colMeta_df.f$Antigen,
                        Mouse = colMeta_df.f$Mouse,
                        # Size = ha.size,
                        Size = colMeta_df.f$Size,
                        col = list(Antigen = (an_colors$Antigen),
                                   Mouse = (an_colors$Mouse),
                                   Size = (sizecols)),
                        annotation_name_side = "bottom",
                        show_legend = F,
                        gp = gpar(fontsize = 6, col = "white",lwd = 0.25),
                        annotation_name_gp= gpar(fontsize = 6, fontface = "bold"),
                        #annotation_width = unit(0.2,"mm"),
                        gap = unit(0,"mm"),
                        simple_anno_size_adjust = TRUE, width = unit(0.5, "cm"))

hal.ant = Legend(legend_gp = gpar(fill = an_colors$Antigen),
                 title = "Antigen",
                 at = colMeta_df.f$Antigen,
                 labels = c("SIIN","SIY"),
                 grid_width = unit(0.25, "cm"),
                 grid_height = unit(0.2,"cm"),
                 title_gp = gpar(fontsize = 8, fontface = "bold"),
                 labels_gp = gpar(fontsize = 6,fontface = "bold"))

hal.mouse = Legend(title = "Mouse",
                   at = names(an_colors$Mouse),
                   legend_gp = gpar(fill = an_colors$Mouse),
                   grid_width = unit(2.5, "mm"),
                   grid_height = unit(0.2,"cm"),
                   title_gp = gpar(fontsize = 8, fontface = "bold"),
                   labels_gp = gpar(fontsize = 6,fontface = "bold"))

hal.size = Legend(title = "Clonotype Size",
                  col_fun = an_colors$Size,
                  grid_height = unit(0.25, "cm"),
                  title_gp = gpar(fontsize = 8, fontface = "bold"),
                  labels_gp = gpar(fontsize = 6, fontface = "bold"),
                  at = c(1, 50, 100),
                  labels = c(1,50, "100+"),
                  legend_width = unit(1.25, "cm"),
                  direction = "horizontal")
### heatmap

#### column dendrogram
jaccard_dist <- function(x, y) {
  #Discretize x,y by clonotypeID (index used)
  floor.val=0.0
  for (i in 1:length(x)) {if (x[i] > floor.val) {x <- replace(x,i,i)} else {x <- replace(x,i,0)}}
  for (i in 1:length(y)) {if (y[i] > floor.val) {y <- replace(y,i,i)} else {y <- replace(y,i,0)}}
  x<-x[x != 0 | 1]; y<-y[y != 0 | 1] #Drop 0s (no clonotype)
  1 - length(intersect(x, y))/length(union(x, y))
}

## calculate column clustering
# set up matrix
nr = nrow(t(cM.f))
mat2 = matrix(NA, nrow = nr, ncol = nr)
rownames(mat2) = colnames(mat2) = rownames(t(cM.f))
# calculate jaccard distance
for(i in 2:nr) {
  for(j in 1:(nr-1)) {
    mat2[i, j] = jaccard_dist(x = t(cM.f)[i, ], y = t(cM.f)[j, ])
  }
}
# plotting/formatting
dend.dist = as.dist(mat2) # convert to distance function
dend = hclust(dend.dist, method = "ward.D") # hclust
dend = cutree(dend, k = 4) # cutoffs

#### row dendrogram
ec.d <- dist(cM.f) # calculate euclidean distance
phcl.hc <- hclust(ec.d) # hclust
phcl.hc2 <- cutree(phcl.hc, k = 4) # cutoff

cf = colorRamp2(c(0,0.0001,0.01,0.1,1),c("grey95",color[1],color[2],color[3],color[7]))

# plot heatmap
p2 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix((cM.f)),
  col = cf,
  column_title = "Cluster",
  column_title_side = "bottom",
  name = "Proportion of Cells",
  rect_gp = gpar(col = "grey", lwd = 0.25),
  border = F,
  row_title = "Clonotype",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  heatmap_width = unit(2,"in"),
  heatmap_height = unit(3.55,"in"),
  show_column_dend = T,
  show_row_dend = T,
  show_row_names = F,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, face = "bold"),
                              labels_gp = gpar(fontsize = 8), face='bold',
                              legend_direction = "vertical",
                              legend_height = unit(1.25,"cm"),
                              legend_width = unit(1.25,"cm")),
  clustering_distance_columns = jaccard_dist,
  clustering_method_columns = "ward.D",
  clustering_method_rows = "ward.D",
  row_split = 4,
  row_gap = unit(2,"mm"),
  row_dend_width = unit(0.5,"cm"),
  #cluster_columns = cluster_within_group((as.matrix(cM.f)), dend),
  #cluster_rows = cluster_within_group(t(as.matrix(cM.f)), phcl.hc2), 
  column_dend_height = unit(0.5,"cm")
)


# heatmap legend
heat.legend <- Legend(cf,
                      title = "Proportion of Cells",
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      at = c(0,0.5,1),
                      labels = c("0","0.5","1"),
                      labels_gp = gpar(fontsize = 6,fontface = "bold"),
                      direction = "horizontal",
                      legend_width = unit(2,"cm"),
                      grid_height = unit(0.25, "cm"),)

hal <- list(hal.ant,hal.mouse,hal.size, heat.legend)

pdf(paste0(od,"fig2f_heatmap.pdf"),width = 4.3, height = 4.55) # plot final pheatmap with clonotype clustering info
draw((p2 + ha.ant),show_heatmap_legend = F, annotation_legend_list = hal, annotation_legend_side = "right",show_annotation_legend = T)
dev.off()



#### FIGURE S4A ####
bad_clonos <- read.csv(file='200922_multiplemouse_clonotypes.csv')  # From Clonotype_MetadataProcessingForHeatmap.ipynb
clmd = read.csv(file="2020-12-03_all_clonotypes_metadata.csv")  # From Clonotype_MetadataProcessingForHeatmap.ipynb

# setup
cluster.clonotype <- table(lung$raw_clonotype_id) # table of number of cells in each clonotype
### calculations ###
filter.f = total > 5 # cell number minimum 
cM.f = cM[filter.f,] # filter out clonotypes that dont meet cell number cutoff
cM.ftot = Matrix::rowSums(cM.f)
write.csv(cM.f,"clonotypesAcrossClustersRaw.csv")

colnames(cM.f) = paste0('C',colnames(cM))
cM.f = cM.f / cM.ftot


# perform hierarchial clustering
ec.da <- dist(cM.f)
phcl.hca <- hclust(ec.da)
phcl.hca2 <- cutree(phcl.hca, h = 0.6) # assigns clonotypes to clusters
#write.csv(phcl.hca2, file=paste0(od,"heatmap_clonotypes.csv")) #export to python

## follow Clonotype_MetadataProcessingForHeatmap.ipynb!!!
## then come back!! 

clmd = read.csv(file="data/2020-12-31_all_clonotypes_metadata.csv") # Clonotype_MetadataProcessingForHeatmap.ipynb
rownames(clmd) = clmd$clonotype
# update metadata
colMeta_df.f <- data.frame(Antigen = sub(".*_", "", rownames(cM.f)),
                           'Clonotype Cluster' = as.character(phcl.hca2),
                           Mouse = clmd$Mouse,
                           Size = cM.ftot,
                           #Sparsity = gini,
                           stringsAsFactors = F,
                           row.names = rownames(cM.f)) # provide confusion matrix metadata with new clusters


### plotting ###
# plotting conditions #
an_colors = list(
  Antigen = c(siin = "#F8766D", siy = "#00BFC4"),
  Clonotype_Cluster_All = c(
    '1' = "#222222",
    '2' = "#F6222E",
    '3' = "#FE00FA",
    '4' = "#16FF32",
    '5' = "#3283FE",
    '6' = "#FEAF16",
    '7' = "#B00068",
    '8' = "#1CFFCE",
    '9' = "#90AD1C",
    '10' = "#2ED9FF",
    '11' = "#DEA0FD",
    '12' = "#AA0DFE"
  ),
  Mouse = c(
    "MB9452" = "#AA0DFE",
    "MB9457" = "#F6222E",
    "MB9458" = "#FE00FA",
    "MB9459" = "#16FF32",
    "MB9509" = "#3283FE",
    "MB9561" = "#FEAF16",
    "MB9619" = "#B00068",
    "MB9621" = "#DEA0FD",
    "MB9622" = "darkorange4",
    "MB9625" = "tan"
  ),
  Size = colorRamp2(c(5,10,88),colors=c("grey","goldenrod","firebrick")))
genes <- c("Gzmb",
           "Cx3cr1",
           "Il17a",
           "Lag3",
           "Tcf7",
           "Ccr6",
           "Havcr2")
clonotypes = unique(clmd$clonotype)
lung = SetIdent(lung,value="raw_clonotype_id")

# calculate # of cells exp gene; this takes some time so you can run it once and reload with commented script below
for (cl in clonotypes){
  cells = WhichCells(lung,idents=cl)
  l = exp['Havcr2',cells ] > .05
  tot[cl,'Havcr2'] = length(l[l == TRUE])
}
tot = read.csv('/Users/amandacruz/210426_geneCount.csv')
write.csv(tot,"geneCount.csv")

#tot = read.csv('/Users/amandacruz/Dropbox (MIT)/SIIN vs SIY 10X/Final_Analysis/210601_geneCount.csv')
#rownames(tot) = tot$X
#tot$X = NULL
#tot$clonotypes = NULL
#tot = tot[clmd$clonotype,]
#exp = lung@assays$RNA@data[genes,]


rownames(tot) = tot$X
tot$X = NULL
tot$clonotypes = NULL
tot = tot[clmd$clonotype,]
exp = lung@assays$RNA@data[genes,]

cM.ftot = rowSums(tot)
tot2 = tot/cM.ftot

sizecols = colorRamp2(c(1,10,88),
                      colors = c("palegoldenrod",
                                 "lightgoldenrod",
                                 "darkorange")
)(seq(1, max(colMeta_df.f$Size)))

names(sizecols) = seq(1, length(sizecols))

# row annotations
ha.ant <- rowAnnotation(Antigen = colMeta_df.f$Antigen,
                        Mouse = colMeta_df.f$Mouse,
                        Size = colMeta_df.f$Size,
                        #Sparsity = colMeta_df.f$Sparsity,
                        Clonotype_Cluster_All = colMeta_df.f$Clonotype.Cluster,
                        col = list(Antigen = (an_colors$Antigen),
                                   Mouse = (an_colors$Mouse),
                                   Size = (sizecols),
                                   #Sparsity =  colorRamp2(c(min(colMeta_df.f$Sparsity),
                                   #                         median(colMeta_df.f$Sparsity),
                                   #                        max(colMeta_df.f$Sparsity)
                                   #                       ),
                                   #                    c("#E7E7E7",
                                   #                     "#9B9B9B",
                                   #                    "#3B3B3B")
                                   #                 ),
                                   Clonotype_Cluster_All = an_colors$Clonotype_Cluster_All),
                        annotation_name_side = "bottom",
                        show_legend = F,
                        gp = gpar(fontsize = 6, col = "white",lwd = 0.5),
                        annotation_name_gp= gpar(fontsize = 6, 
                                                 fontface = "bold"),
                        gap = unit(0,"mm"),
                        simple_anno_size_adjust = TRUE,
                        width = unit(1, "cm"),
                        annotation_label = c("Antigen",
                                             "Mouse",
                                             "Clonotype Size",
                                             #"Sparsity",
                                             "Clonotype Cluster"))

# legends                        
hal.ant = Legend(legend_gp = gpar(fill = an_colors$Antigen),
                 title = "Antigen",
                 at = colMeta_df.f$Antigen,
                 labels = c("SIIN","SIY"),
                 grid_width = unit(0.25, "cm"),
                 grid_height = unit(0.2,"cm"),
                 title_gp = gpar(fontsize = 8,
                                 fontface = "bold"),
                 labels_gp = gpar(fontsize = 6, 
                                  fontface = "bold"))

hal.mouse = Legend(title = "Mouse",
                   at = names(an_colors$Mouse),
                   legend_gp = gpar(fill = an_colors$Mouse),
                   grid_width = unit(2.5, "mm"),
                   grid_height = unit(0.2,"cm"),
                   title_gp = gpar(fontsize = 8, 
                                   fontface = "bold"),
                   labels_gp = gpar(fontsize = 6, 
                                    fontface = "bold"),
                   nrow = 5)

hal.size = Legend(title = "Clonotype Size",
                  col_fun = an_colors$Size,
                  grid_height = unit(0.25, "cm"),
                  title_gp = gpar(fontsize = 8, 
                                  fontface = "bold"),
                  labels_gp = gpar(fontsize = 6, 
                                   fontface = "bold"),
                  at = c(1, 50, 100),
                  labels = c(1,50, "100+"),
                  legend_width = unit(1.25, "cm"),
                  direction = "horizontal")

# hal.sparsity = Legend(title = "Sparsity",
#                   col_fun = colorRamp2(c(min(colMeta_df.f$Sparsity),
#                                          median(colMeta_df.f$Sparsity),
#                                          max(colMeta_df.f$Sparsity)),
#                                        c("#E7E7E7",
#                                          "#9B9B9B",
#                                          "#3B3B3B")
#                                        ),
#                   grid_height = unit(0.25, "cm"),
#                   title_gp = gpar(fontsize = 8, 
#                                   fontface = "bold"),
#                   labels_gp = gpar(fontsize = 6, 
#                                    fontface = "bold"),
#                   at = c(0, 1),
#                   labels = c(0, 1),
#                   legend_width = unit(1.25, "cm"),
#                   direction = "horizontal")

hal.clonoclust = Legend(legend_gp = gpar(fill = an_colors$Clonotype_Cluster_All),
                        title = "Clonotype Cluster",
                        at = colMeta_df.f$Clonotype.Cluster,
                        labels = names(an_colors$Clonotype_Cluster_All),
                        grid_width = unit(0.25, "cm"),
                        grid_height = unit(0.2,"cm"),
                        nrow = 4,
                        title_gp = gpar(fontsize = 8, 
                                        fontface = "bold"),
                        labels_gp = gpar(fontsize = 6, 
                                         fontface = "bold"))


# gene enrichment color schemes
cfG <- colorRamp2(c(min(tot2$Gzmb),
                    max(tot2$Gzmb)),
                  c("grey99",
                    "chartreuse4"))
cfCx <- colorRamp2(c(min(tot2$Cx3cr1),
                     max(tot2$Cx3cr1)),
                   c("grey99",
                     "burlywood4"))
cfI <- colorRamp2(c(min(tot2$Il17a),
                    max(tot2$Il17a)),
                  c("grey99",
                    "cyan4"))
cfH <- colorRamp2(c(min(tot2$Havcr2),
                    max(tot2$Havcr2)),
                  c("grey99",
                    "lightpink3"))
cfT <- colorRamp2(c(min(tot2$Tcf7),
                    max(tot2$Tcf7)),
                  c("grey99",
                    "goldenrod2"))
cfC <- colorRamp2(c(min(tot2$Ccr6),
                    max(tot2$Ccr6)),
                  c("grey99",
                    "purple4"))

# row annotations for gene enrichment
sph <- rowAnnotation(Gzmb = tot2$Gzmb,
                     Havcr2 = tot2$Havcr2,
                     Cx3cr1 = tot2$Cx3cr1,
                     Tcf7 = tot2$Tcf7,
                     Ccr6 = tot2$Ccr6,
                     Il17a = tot2$Il17a,
                     col = list(Gzmb = cfG,
                                Havcr2 = cfH,
                                Cx3cr1 = cfCx,
                                Tcf7 = cfT,
                                Ccr6 = cfC,
                                Il17a = cfI),
                     annotation_label = c("Gzmb",
                                          "Havcr2",
                                          "Cx3cr1",
                                          "Tcf7",
                                          "Ccr6",
                                          "Il17a"),
                     gp = gpar(fontsize = 6, col = "grey",lwd = 0.1),
                     annotation_name_gp= gpar(fontsize = 6, fontface = "bold"),
                     simple_anno_size_adjust = TRUE, width = unit(1, "cm"),
                     gap = unit(0,"mm"))

# gene enrichment legends
cfl.G = Legend(title = "%Gzmb+",
               col_fun = cfG,
               grid_height = unit(0.25, "cm"),
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               at = c(0,0.5,1),
               labels = c(0,0.5,1),
               legend_width = unit(1.25, "cm"),
               direction = "horizontal")
cfl.Cx = Legend(title = "%Cx3cr1+",
                col_fun = cfCx,
                grid_height = unit(0.25, "cm"),
                title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 6, fontface = "bold"),
                at = c(0,0.5,1),
                labels = c(0,0.5,1),
                legend_width = unit(1.25, "cm"),
                direction = "horizontal")
cfl.I = Legend(title = "%Il17a+",
               col_fun = cfI,
               grid_height = unit(0.25, "cm"),
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               at = c(0,0.5,1),
               labels = c(0,0.5,1),
               legend_width = unit(1.25, "cm"),
               direction = "horizontal")
cfl.H = Legend(title = "%Havcr2+",
               col_fun = cfH,
               grid_height = unit(0.25, "cm"),
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               at = c(0,0.5,1),
               labels = c(0,0.5,1),
               legend_width = unit(1.25, "cm"),
               direction = "horizontal")
cfl.T = Legend(title = "%Tcf7+",
               col_fun = cfT,
               grid_height = unit(0.25, "cm"),
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               at = c(0,0.5,1),
               labels = c(0,0.5,1),
               legend_width = unit(1.25, "cm"),
               direction = "horizontal")
cfl.C = Legend(title = "%Ccr6+",
               col_fun = cfC,
               grid_height = unit(0.25, "cm"),
               title_gp = gpar(fontsize = 8, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               at = c(0,0.5,1),
               labels = c(0,0.5,1),
               legend_width = unit(1.25, "cm"),
               direction = "horizontal")

# heatmap ordering
row_dend = as.dendrogram(cluster_within_group(t(as.matrix(cM.f)), phcl.hca2))

# heatmap color scheme
color = brewer.pal(7,"Reds")
cf = colorRamp2(c(0,
                  0.0001,
                  0.01,
                  0.1,
                  1),
                c("grey95",
                  color[1],
                  color[2],
                  color[3],
                  color[7])
)
colnames(cM.f) = colnames(cM)

# heatmap
p2 <- ComplexHeatmap::Heatmap(
  matrix = (as.matrix((cM.f))),
  col = cf,
  column_title = "Cluster",
  column_title_side = "bottom",
  name = "Proportion of Cells",
  rect_gp = gpar(col = "grey", lwd = 0.2),
  border = F,
  row_title = "Clonotype",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  heatmap_width = unit(2,"in"),
  heatmap_height = unit(8,"in"),
  show_column_dend = F,
  show_row_dend = T,
  cluster_rows = row_dend,
  show_row_names = F,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, face = "bold"),
                              labels_gp = gpar(fontsize = 8), face='bold',
                              legend_direction = "vertical",
                              legend_height = unit(2.5,"cm"))
)

# legend
heat.legend <- Legend(col_fun = cf,
                      title = "Proportion of Cells",
                      #title_position="leftcenter-rot",
                      title_gp = gpar(fontsize = 8, fontface = "bold"),
                      at = c(0,0.5,1),
                      labels = c("0","0.5","1"),
                      labels_gp = gpar(fontsize = 6, fontface = "bold"),
                      direction = "horizontal",
                      legend_width = unit(2,"cm"),
                      grid_height = unit(0.25, "cm"),)

hal.l <- packLegend(hal.ant, hal.mouse, hal.size, 
                    #hal.sparsity,
                    hal.clonoclust, heat.legend, cfl.G, cfl.T, cfl.H, cfl.Cx, cfl.C, cfl.I)#max_height = unit(10,"cm"))
hal <- list(hal.l, heat.legend)

pdf(paste0(od,"figS4a_heatmap.pdf"),width = 5, height = 8.5) # plot final pheatmap with clonotype clustering info
draw((p2 + ha.ant+sph),show_heatmap_legend = F, annotation_legend_list = hal.l, annotation_legend_side = "right",show_annotation_legend = T)
dev.off()


### FIGURE S4B ###
size = read.csv('md_AddClonotypeSize.csv') # from AddClonotypeSize_md.ipynb
size$X = NULL
lung@meta.data$clonotype_size = size$clonotype_size
lung = SetIdent(lung,value = "orig.ident")
siin = subset(lung, idents = "siin")
siy = subset(lung, idents = "siy")

# plot and save all data together
p <- FeaturePlot(lung, 
                 features = 'clonotype_size', 
                 order = T,
                 pt.size = 0.005) &
  theme_void() &
  theme(legend.position = 'right') &
  scale_colour_gradientn(colours = 
                           c("gray75",
                             "palegoldenrod",
                             "goldenrod",
                             "darkorange3",
                             "firebrick")) &
  labs(title = NULL) &
  theme(legend.text = 
          element_text(size = 6, face = 'bold'),
        legend.key.width =
          unit(0.25, units = 'cm'),
        legend.key.height =
          unit(0.25, units = 'cm'))

ggsave(filename = paste0(od,"figS4B_all_clonotypesize.pdf"),
       plot = (p + theme(legend.position = 'none')),
       device = "pdf",
       units = "cm",
       width = 3.5,
       height = 3.5,
       scale = 1)



# plot and save siin clonotype size
p <- FeaturePlot(siin, features = 'clonotype_size', 
                 order = T, 
                 pt.size = 0.005) &
  theme_void() &
  theme(legend.position = 'right') &
  scale_colour_gradientn(colours = 
                           c("gray75",
                             "palegoldenrod",
                             "goldenrod",
                             "darkorange3",
                             "firebrick")) &
  labs(title = NULL) &
  theme(legend.text = element_text(size = 6, face = 'bold'),
        legend.key.width = unit(0.25, units = 'cm'),
        legend.key.height = unit(0.25, units = 'cm'))

ggsave(filename = paste0(od,"figS4D_siin_clonotypesize.pdf"),
       plot = (p + theme(legend.position = 'none')),
       device = "pdf",
       units = "cm",
       width = 3.5,
       height = 3.5,
       scale = 1)

# legend
pl.l <- cowplot::get_legend(p)
pdf(file = paste0(od,"figS4B_umap_siy_clonotypeSize_legend.pdf"), width = 0.5, height = 1.5)
grid.newpage()
grid.draw(pl.l)
dev.off()

# plot
p <- FeaturePlot(siy, features = 'clonotype_size',
                 order = T,
                 pt.size = 0.005) &
  theme_void() &
  theme(legend.position = 'right') &
  scale_colour_gradientn(colours = 
                           c("gray75",
                             "palegoldenrod",
                             "goldenrod",
                             "darkorange3",
                             "firebrick")) &
  labs(title = NULL) &
  theme(legend.text = 
          element_text(size = 6, face = 'bold'),
        legend.key.width = 
          unit(0.25, units = 'cm'),
        legend.key.height = 
          unit(0.25, units = 'cm'))

ggsave(filename = paste0(od,"figS4D_siy_clonotypesize.pdf"),
       plot = (p + theme(legend.position = 'none')),
       device = "pdf",
       units = "cm",
       width = 3.5,
       height = 3.5,
       scale = 1)
# legend
pl.l <- cowplot::get_legend(p)
pdf(file = paste0(od,"figS4B_umap_siy_clonotypeSize_legend.pdf"), width = 0.5, height = 1.5)
grid.newpage()
grid.draw(pl.l)
dev.off()


genes <- c('Lag3','Pdcd1','Tigit','Tox','Havcr2','Icos','Tnfrsf9','Tnfrsf4','Cd160','Ifng','Pfn1','Cx3cr1','Gzmb','Gzma','Gzmk','Eomes','Cxcr6','Ccl5','Itgae','Il7r','Tcf7','Xcl1','Il17a','Rorc','Tmem176a','Ramp1')

genes2 <- c('Ccl3','Ccl4','Tnfrsf9','Rgs16','Tox','Ifitm2','Tmem176a','Lgals3','Rsrp1','Gzmk','Cd5','Ctla2a','Cd6','Itgae','Ifitm1','Ahnak','Malat1','Klf6','Bcl2','Itgae','Fgl2','Itgb7','Ccl5','Cxcr6','Mt1','Atf5','Bcl2','Cdkn1a','Slc7a5','Ramp1','Cd81','Tmem176a','Tmem176b','Ccr6','Rorc','Il17a','Prf1','Lag3','Nfatc1','Havcr2','Cd27','Gzmb','Cx3cr1','Klf2','Ly6c2','Nfe2l2','Klf12','Dnmt3b','Xcl1','Egr1','Il7r','Nr4a1','Tcf7','Cd83','Cd27','Tox','Cd160','Rgs16','Stmn1','Birc5','Rrm2','Hmgb2','Tubb5')
fig5b <- list(genes2) # make list of lists to iterate through

# gene set names
fig5b.names <- c("Markers")

# heatmap calculation parameters
pseudocount = 1
lower_threshold = 0
scale_min = -3
scale_max = 3
major_axis <- 1
minor_axis <- 2
max.size = 10

# aethetic color scheme for heatmap
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(min(res2, na.rm = T), 0,
                       max(res2, na.rm = T)),
                     c("blue", "grey99", "red"))
cds <- readRDS('Monocle3/201006/201006_total_cds.rds')
cds@colData$clonotype_cluster_all = lung@meta.data$clonotype_cluster_all
cds2 = cds[module,names(lung3@active.ident)]
for (i in 1:length(fig5b)) {
  
  module = fig5b[[i]] # extract list of genes for each geneset
  exprs_mat <- t(as.matrix(lung3@assays$RNA@data[genes2,])) # get expression data for just that 
  exprs_mat <- reshape2::melt(exprs_mat) # make matrix into long format rather than wide
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression") # define colnames of matrix
  exprs_mat$Gene <- as.character(exprs_mat$Gene) # data type for gene name conversion
  cell_group <- lung3@meta.data$clonotype_cluster_all # define cell groups by seurat clusters
  names(cell_group) = names(lung3@meta.data$orig.ident) # carry over cluster names and cell names
  exprs_mat$Group <- cell_group[exprs_mat$Cell] # populate matrix with cluster ID for each cell
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) ==
                                            FALSE) # remove NA values (shouldnt be any)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(
      mean = mean(log(Expression + pseudocount)),
      percentage = sum(Expression > lower_threshold) / length(Expression)
    ) # take summary for each cluster and calculate mean expression as defined by mean(log(expression +1))
  # defines percentage as the sum of expression values that pass lower_threshold (by default zero)
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min,
                        ExpVal$mean) # scaling
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max,
                        ExpVal$mean) # scaling
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, "gene_short_name"] # fData allows you to access cds rowData
  res <-
    reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 +
                                                                                major_axis]) #reshape
  # heatmap annotation (cluster ID) definition
  ha.ant <- HeatmapAnnotation(
    Cluster = as.character(seq(1,12,1)),
    col = list(Cluster = an_colors$Clonotype_Cluster_All),
    annotation_name_side = "left",
    annotation_label = "",
    annotation_height = unit(0.1,"mm"),
    show_legend = F,
    gp = gpar(
      fontsize = 6,
      col = "white",
      lwd = 0.5
    )
  )
  ha.ant <- rowAnnotation(
    Cluster = as.character(seq(1,12,1)),
    col = list(Cluster = an_colors$Clonotype_Cluster_All),
    #annotation_name_side = "left",
    annotation_label = "",
    simple_anno_size_adjust = TRUE, width = unit(1, "mm"),
    show_legend = F,
    gp = gpar(
      fontsize = 6,
      col = "white",
      lwd = 0.5
    )
  )
  ## my code from here on out
  # format matrix
  res2 = res 
  rownames(res2) <- res$Group
  res2$Group <- NULL
  # scale data by *cluster* (not overall)
  res2 = t(scale(as.matrix(res2)))
  
  # plot
  ht = Heatmap(
    t(res2),
    show_column_dend = F,
    show_row_dend = F,
    show_row_names = T,
    cluster_rows = T,
    cluster_columns = T,
    #row_order = as.vector(module),
    #column_order = as.character(seq(1,12,1)),
    row_names_max_width = unit(0.75,"cm"),
    col = col_fun,
    #rect_gp = gpar(col = "black", lwd = 0.5), # for borders on cells
    heatmap_width = unit(3.5,"in"), # size
    heatmap_height = unit(1.25, "in"), # size
    column_names_rot = 90, # rotate column titles
    column_names_gp = gpar(fontsize = 5), # format col titles
    row_names_gp = gpar(fontsize = 6, fontface = "bold"), # format row titles
    show_heatmap_legend = F, # exclude legend
    left_annotation = ha.ant
    #top_annotation = ha.ant # add cluster annotations and put it on the top
  )
  
  # assign variable names to plots through iterations
  nam <- paste("ht", i, sep = "") 
  assign(nam, ht)
  
  # save plot
  pdf(file = paste0(od,"figS4_heatmap_",fig5b.names[i],".pdf"), height = 2.5, width = 4)
  draw(ht)
  dev.off()
  
}

htg.lc = c(sc.cols[3], sc.cols[4], sc.cols[5], sc.cols[9])
# make cluster legend
htg.l = Legend(legend_gp = gpar(fill = an_colors$Clonotype_Cluster_All),
               title = "Cluster",
               at = seq(1,12,1),
               grid_width = unit(2.5, "mm"),
               grid_height = unit(2.5,"mm"),
               title_gp = gpar(fontsize = 6, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               direction = "vertical",
               ncol = 1,
               title_position = "topcenter",
               legend_height = unit(2,"cm"))
# plot legend & save
pdf(paste0(od,"fig5b_legend.pdf"), height = 2, width = 0.5)
draw(htg.l)
dev.off()

# scale bar
htg.s = Legend(col_fun = col_fun,
               title = "Log10(Mean Expression)",
               title_gp = gpar(fontsize = 6, fontface = "bold"),
               legend_height = unit(2, "cm"), labels_gp = gpar(fontsize = 6, fontface = "bold"),
               direction = "horizontal",
               grid_width = unit(2,"mm"))

# save scalebar
pdf(paste0(od,"fig5b_scale.pdf"), height = 1, width = 2)
draw(htg.s)
dev.off()


#### FIGURE S4C ####

# we will test C4C8 vs All other clusters and then compare SIIN and SIY.
cM.ft.siin = cM.ft.siin %>% dplyr::mutate(nC2 = rowSums(cM.ft.siin[,(names(cM.ft.siin) != 'C2' & names(cM.ft.siin) != 'C4' & names(cM.ft.siin) != 'C4C8' & names(cM.ft.siin) != 'nC4C8')]))
cM.ft.siy = cM.ft.siy %>% dplyr::mutate(nC2 = rowSums(cM.ft.siy[,(names(cM.ft.siy) != 'C8' & names(cM.ft.siy) != 'C4' & names(cM.ft.siy) != 'C4C8' & names(cM.ft.siy) != 'nC4C8')]))


# flatten table for plotting
cM.ft.siin.ecdf = melt(cM.ft.siin %>% select(c('C2'))) 
cM.ft.siy.ecdf = melt(cM.ft.siy %>% select(c('C2')))

# add antigen to variable name
cM.ft.siy.ecdf$Antigen = paste0(cM.ft.siy.ecdf$variable, " SIY")
cM.ft.siin.ecdf$Antigen = paste0(cM.ft.siin.ecdf$variable, " SIIN")

# combine into one dataframe for plotting
cM.ecdf = rbind(cM.ft.siin.ecdf, cM.ft.siy.ecdf)
#cM.ecdf2 = cM.ecdf %>% filter(variable != 'nC4C8')

# plot
p = ggplot(cM.ecdf, aes(value,color=Antigen)) +
  stat_ecdf(geom = "line") +
  stat_ecdf(geom = "point") +
  scale_color_manual(values=c("#F8766D","#00BFC2","#AA5049","#006466")) +
  xlab("Proportion of Cells in an Individual Clonotype") +
  ylab("Proportion of Clonotypes Specific for Antigen") +
  labs(color = "") + 
  xlim(0,1)
ggsave2(
  filename = paste0(od, "figS4X_ECDF_C2vsAll_C4C8.pdf"),
  plot = (
    p + theme_classic() + ggtitle(label = paste0('')) +
      theme(
        legend.position = 'none',
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 6, face = "bold")
      ) +
      xlab(NULL) + 
      ylab(NULL)
  ),
  device = "pdf",
  units = "in",
  width = 1.5,
  height = 1.5,
  scale = 1
)

# All Clonotypes 

# flatten table for plotting
cM.siin.ecdf = melt(cM.siin %>% select(c('C2'))) 
cM.siy.ecdf = melt(cM.siy %>% select(c('C2')))

# add antigen to variable name
cM.siy.ecdf$Antigen = paste0(cM.siy.ecdf$variable, " SIY")
cM.siin.ecdf$Antigen = paste0(cM.siin.ecdf$variable, " SIIN")

# combine into one dataframe for plotting
cM.ecdf = rbind(cM.siin.ecdf, cM.siy.ecdf)
#cM.ecdf2 = cM.ecdf %>% filter(variable != 'nC4C8')

# plot
p = ggplot(cM.ecdf, aes(value,color=Antigen)) +
  stat_ecdf(geom = "line") +
  stat_ecdf(geom = "point") +
  scale_color_manual(values=c("#F8766D","#00BFC2","#AA5049","#006466")) +
  xlab("Proportion of Cells in an Individual Clonotype") +
  ylab("Proportion of Clonotypes Specific for Antigen") +
  labs(color = "") + 
  xlim(0,1)
ggsave2(
  filename = paste0(od, "figS4X_ECDF_C2vsAll_All.pdf"),
  plot = (
    p + theme_classic() + ggtitle(label = paste0('')) +
      theme(
        legend.position = 'none',
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 6, face = "bold")
      ) +
      xlab(NULL) + 
      ylab(NULL)
  ),
  device = "pdf",
  units = "in",
  width = 1.5,
  height = 1.5,
  scale = 1
)

#### S4B #### 
# assigns clonotype cluster to each cell
lung@meta.data$clonotype_cluster_all = revalue(lung@meta.data$raw_clonotype_id,
                                               phcl.hca2)
# make sure its not an integer
lung@meta.data$clonotype_cluster_all = as.character(lung@meta.data$clonotype_cluster_all)


# UMAP plot all of the clonotype cluster cells, aggregated by clonotype cluster
lung <- SetIdent(lung,
                 value = 'clonotype_cluster_all')
# for for loop
ph.clusters.max <- max(phcl.hca2)

# for all clonotype clusters plot
for (i in 1:ph.clusters.max) {
  # setup
  cl <- WhichCells(lung,
                   idents = i)
  # plot
  pl = DimPlot(
    object = lung,
    reduction = "umap",
    cells.highlight = cl,
    group.by = 'clonotype_cluster_all',
    sizes.highlight = 0.05) +
    guides(fill = guide_legend(title = NULL)) +
    ggtitle(label = paste0("Cluster ", as.character(i))) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, face = "bold")) +
    scale_color_manual(
      labels = c("other", paste0("Cluster ", as.character(i))), 
      values = c("grey", an_colors$Clonotype_Cluster_All[i])
    )
  # save
  ggsave2(
    filename = paste0(od, "figS4C_umap_clonoCluster", i, ".pdf"),
    plot = (
      pl + theme_void() + ggtitle(label = paste0("Clonotype Cluster ", as.character(i))) +
        theme(
          legend.position = 'none',
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
        )
    ),
    device = "pdf",
    units = "in",
    width = 1.5,
    height = 1.5,
    scale = 1
  )
  
}

```

```{r}
lung <- SetIdent(lung,value="clonotype_cluster_all")
lung3 <- subset(lung,idents=c("1","2","3","4","5","6","7","8","9","10","11","12"))
test = FindAllMarkers(lung3)
write.csv(test,paste0(od,"clonotype_cluster_markers.csv"))
```

```{r}
#### STATISTICAL TEST FOR SPARSITY IN C4C8 SIINvSIY ####
filter.f = total > 5 # cell number minimum 
#cM.f = cM[filter.f,] # filter out clonotypes that dont meet cell number cutoff
filter.f2 = cM.f[,6] > 0 | cM.f[,7] > 0 # include only clusters that have a cell in Cluster 4 or cluster 8
cM.f2 = cM.f[filter.f2,]
colMeta_df.f2 = colMeta_df.f[rownames(cM.f2),]
gini = apply(cM.f2, 1, function(x) Gini(x))
colMeta_df.f2$Sparsity = gini
# subset gini index calculations by antigen
cmdf2.siin = colMeta_df.f2 %>% filter(Antigen == "siin")
cmdf2.siy = colMeta_df.f2 %>% filter(Antigen == "siy")

# ks test
ks.test(as.numeric(cmdf2.siin$Sparsity),as.numeric(cmdf2.siy$Sparsity))
# plot ECDF
k = ggplot(colMeta_df.f2, aes(Sparsity,color=Antigen)) + stat_ecdf(geom = "point") + scale_color_manual(values=c("#F8766D","#00BFC4")) +
  xlab("Gini Index") +
  ylab("Proportion of Clonotypes Specific for Antigen") +
  labs(color = "")

ggsave2(
  filename = paste0(od, "figS4X_ECDF_Sparsity_C4C8.pdf"),
  plot = (
    k + theme_classic() + ggtitle(label = paste0('')) +
      theme(
        legend.position = 'none',
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5)
      ) +
      xlab(NULL) +
      ylab(NULL)
  ),
  device = "pdf",
  units = "in",
  width = 1.5,
  height = 1.5,
  scale = 1
)
```


#### FIGURE 5A ####

### Main Plot ###
# setup
lung <- SetIdent(lung, value = 'RNA_snn_res.0.7')
lung2 <- SetIdent(lung2, value = 'RNA_snn_res.0.7') # subplot

# plot
p <- DimPlot(lung, cols = sc.cols, pt.size = 0.05) + theme_void()
pl.l <- cowplot::get_legend(p)

# save main plot
ggsave2(
  filename = paste0(od, "fig5a_umap", ".pdf"),
  plot = (p + theme(legend.position = 'none')),
  device = "pdf",
  units = "cm",
  width = 5.5,
  height = 5.5,
  scale = 1
)

# legend
pdf(
  file = paste0(od, "fig5a_umap_legend2", ".pdf"),
  width = 3,
  height = 0.5
)
grid.newpage()
grid.draw(pl.l)
dev.off()


# subplot
p <- DimPlot(lung2,
             cols = sc.cols,
             pt.size = 0.05) +
  theme_void() +
  theme(
    legend.text =
      element_text(size = 8, face = 'bold'),
    legend.key.width =
      unit(0.4, units = 'cm'),
    legend.key.height =
      unit(0.25, units = 'cm'),
    legend.direction =
      "horizontal"
  )

# save subplot
ggsave2(
  filename = paste0(od, "fig5a_subplot_umap.pdf"),
  plot = (p + theme(legend.position = 'none')),
  device = "pdf",
  units = "cm",
  width = 3,
  height = 3,
  scale = 1)

### gene plots ###

# genes to plot
fig5a <- c("Tcf7",
           "Havcr2")

# plot geneplots and save
for (gene in fig5a) {
  # plot cells and theming
  p <-
    plot_cells(
      cds,
      reduction_method = "UMAP",
      genes = gene,
      show_trajectory_graph = F,
      label_groups_by_cluster = FALSE,
      label_branch_points = F,
      label_leaves = F,
      label_roots = F,
      label_cell_groups = F,
      #trajectory_graph_segment_size = .25,
      cell_size = 0.5
    ) +
    theme_void() +
    theme(
      plot.tag = NULL,
      plot.subtitle = NULL,
      strip.text.x = element_text(size = 8, face = "bold"),
      legend.title = element_blank(),
      legend.key.height = unit(0.25, "cm"),
      legend.text = element_text(size = 6,
                                 face = "bold"),
      legend.key.width = unit(0.5, "line"),
      legend.direction = "horizontal"
    )  +
    scale_colour_gradientn(
      colours = c(
        "gray75",
        "palegoldenrod",
        "goldenrod",
        "darkorange3",
        "firebrick"
      ),
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)
    )
  # save plot (no legend)
  ggsave2(filename = paste0(od,"fig5a_umap_",gene,".pdf"),
          plot = (p + theme(legend.position = 'none')),
          device = "pdf",
          units = "cm",
          width = 3.5,
          height = 3.5,
          scale = 1)
  # extract legend from plot and save
  pl.l <- cowplot::get_legend(p)
  pdf(file = paste0(od,"fig5a_umap_legend_",gene,".pdf"), width = 0.5, height = 1.5)
  grid.newpage()
  grid.draw(pl.l)
  dev.off()
}

#### FIGURE S5A ####

signatures <- read.csv(
  file = "Signatures/CD8_published_tolerance_signatures_for_module_analysis.csv") # formatted as CSV
signatures.id <- colnames(signatures)

# change to title case for mouse gene compatibility
for (i in signatures.id) {
  signatures[[i]] <- str_to_title(signatures[[i]])}
exp.gene <- rownames(lung)

# get signatures and add their scores to metadata
plots <- list()
for (i in signatures.id) {
  
  # get gene lit subset
  signatures.gene <- signatures[i] 
  
  # get rid of empty variables
  signatures.gene <- signatures.gene[!apply(signatures.gene == "", 1, all), ] 
  signatures.gene <- intersect(signatures.gene, exp.gene)
  
  # get scores
  lung.ms = AddModuleScore(
    object = lung.ms,
    features = list(signatures.gene),
    name = paste0(i,"_score"))
  
  # plot
  p <- FeaturePlot(object=lung.ms, 
                   reduction="umap", 
                   features=paste0(i,"_score1"), 
                   order = T, 
                   pt.size = 0.25) + 
    ggtitle(paste0("Signature:",i))   &
    theme_void() &
    theme(legend.position = 'right') &
    scale_colour_gradientn(
      colours =
        c(
          "gray75",
          "gray87",
          "gray89",
          "palegoldenrod",
          "lightgoldenrod",
          "goldenrod",
          "darkorange3",
          "brown"
        ),
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)) &
    labs(title = NULL) &
    theme(
      plot.tag = NULL,
      plot.subtitle = NULL,
      strip.text.x = element_text(size = 8, face = "bold"),
      legend.title = element_blank(),
      legend.key.height = unit(0.25, "cm"),
      legend.text = element_text(size = 6, face = "bold"),
      legend.key.width = unit(0.5, "line"),
      legend.direction = "horizontal"
    )
  pl.l <- cowplot::get_legend(p)
  
  # save
  ggsave2(filename = paste0(od,"figS5A_umap_",i,".pdf"),
          plot = (p + theme(legend.position = 'none')),
          device = "pdf",
          units = "cm",
          width = 3.5,
          height = 3.5,
          scale = 1)
  
  # legend
  pdf(file = paste0(od,"figS5A_umap_legend_",i,".pdf"), width = 1, height = 1.5)
  grid.newpage()
  grid.draw(pl.l)
  dev.off()
}

### Tc17 signature ###
# plot
p <- FeaturePlot(lung.ms,
                 features = 'Mm_Imiquimod_Tc17_v_Tc1_sig1',
                 order = T,
                 pt.size = 0.25)  &
  theme_void() &
  theme(legend.position = 'right') &
  scale_colour_gradientn(
    colours =
      c(
        "gray75",
        "gray87",
        "gray89",
        "palegoldenrod",
        "lightgoldenrod",
        "goldenrod",
        "darkorange3",
        "brown"
      ),
    labels = c(0, 0.5, 1),
    breaks = c(0, 0.5, 1)
  ) &
  labs(title = NULL) &
  theme(
    plot.tag = NULL,
    plot.subtitle = NULL,
    strip.text.x = element_text(size = 8, face = "bold"),
    legend.title = element_blank(),
    legend.key.height = unit(0.25, "cm"),
    legend.text = element_text(size = 6, face = "bold"),
    legend.key.width = unit(0.5, "line"),
    legend.direction = "horizontal"
  )

# legend
pl.l <- cowplot::get_legend(p)

# save
ggsave2(filename = paste0(od,"figS5_umap_Tc17_immiquimod",".pdf"),
        plot = (p + theme(legend.position = 'none')),
        device = "pdf",
        units = "cm",
        width = 3.5,
        height = 3.5,
        scale = 1)

# legend
pdf(file = paste0(od,"fig2c_umap_legend_imiquimod",".pdf"), width = 1, height = 1.5)
grid.newpage()
grid.draw(pl.l)
dev.off()


#### FIGURE 5B #### 

## dependencies/methodology
# written using Monocle3 object and based off of the plot_genes_by_group function in Monocle3.
# transcription factor lists
pvse.tfs <- c('Eomes','Prdm1','Id2','Tox','Tcf7','Lef1','Id3','Runx1')
irs <- c('Pdcd1','Lag3','Havcr2','Entpd1','Tigit','Ctla4','Cd244','Cd160')
traf <- c('Cxcr3','Cx3cr1','Ccl5','Cxcr6','S1pr1','Klf2','Ccr7','Xcl1')
tol <- c('Dapl1','Nrn1','Cd200','Cd83','Rgs16','Nr4a2','Egr2','Nfatc1')
eff <- c('Gzmb','Gzma','Gzmk','Prf1','Ifng','Tnf','Fasl','Il2')
srs <- c('Tnfrsf4', 'Tnfrsf9', 'Icos','Vsir','Il7r','Cd27','Cd28', 'Sell')
tc17 <- c('Ccr6','Ramp1','Tmem176a','Igfbp7','Il23r','Il17a','Rorc','Maf')

fig5b <- list(pvse.tfs,irs,traf,tol,eff,srs,tc17) # make list of lists to iterate through

# gene set names
fig5b.names <- c("Progenitor vs Exhausted Transcription Factors",
                 "Inhibitory Receptors",
                 "Trafficking",
                 "Tolerance",
                 "Effector Molecules",
                 "Cell Surface Receptors",
                 "Tc17")

# heatmap calculation parameters
pseudocount = 1
lower_threshold = 0
scale_min = -3
scale_max = 3
major_axis <- 1
minor_axis <- 2
max.size = 10

# aethetic color scheme for heatmap
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(min(res2, na.rm = T), 1,
                       max(res2, na.rm = T)),
                     c("purple", "pink", "chocolate"))

for (i in 1:length(fig5b)) {
  
  module = fig5b[[i]] # extract list of genes for each geneset
  exprs_mat <- t(as.matrix(exprs(cds)[module,])) # get expression data for just that 
  exprs_mat <- reshape2::melt(exprs_mat) # make matrix into long format rather than wide
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression") # define colnames of matrix
  exprs_mat$Gene <- as.character(exprs_mat$Gene) # data type for gene name conversion
  cell_group <- colData(cds)[, "RNA_snn_res.0.7"] # define cell groups by seurat clusters
  names(cell_group) = colnames(cds) # carry over cluster names and cell names
  exprs_mat$Group <- cell_group[exprs_mat$Cell] # populate matrix with cluster ID for each cell
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) ==
                                            FALSE) # remove NA values (shouldnt be any)
  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(
      mean = mean(log(Expression + pseudocount)),
      percentage = sum(Expression > lower_threshold) / length(Expression)
    ) # take summary for each cluster and calculate mean expression as defined by mean(log(expression +1))
  # defines percentage as the sum of expression values that pass lower_threshold (by default zero)
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min,
                        ExpVal$mean) # scaling
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max,
                        ExpVal$mean) # scaling
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, "gene_short_name"] # fData allows you to access cds rowData
  res <-
    reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 +
                                                                                major_axis]) #reshape
  # heatmap annotation (cluster ID) definition
  ha.ant <- HeatmapAnnotation(
    Cluster = c("2", "3", "4", "8"),
    col = list(Cluster = sc.cols),
    annotation_name_side = "left",
    annotation_label = "",
    annotation_height = unit(0.1,"mm"),
    show_legend = F,
    gp = gpar(
      fontsize = 6,
      col = "white",
      lwd = 0.5
    )
  )
  ## my code from here on out
  # format matrix
  res2 = res 
  rownames(res2) <- res2$Group
  res2$Group <- NULL
  # scale data by *cluster* (not overall)
  res2 = t(scale(as.matrix(res2)))
  
  # plot
  ht = Heatmap(
    res2,
    show_column_dend = F,
    show_row_dend = F,
    show_row_names = T,
    cluster_rows = F,
    cluster_columns = F,
    row_order = module,
    row_names_max_width = unit(0.75,"cm"),
    col = col_fun,
    rect_gp = gpar(col = "black", lwd = 0.5), # for borders on cells
    heatmap_width = unit(1.75,"cm"), # size
    heatmap_height = unit(1, "in"), # size
    column_names_rot = 0, # rotate column titles
    column_names_gp = gpar(fontsize = 6, fontface = "bold"), # format col titles
    row_names_gp = gpar(fontsize = 6), # format row titles
    show_heatmap_legend = F, # exclude legend
    top_annotation = ha.ant # add cluster annotations and put it on the top
  )
  
  # assign variable names to plots through iterations
  nam <- paste("ht", i, sep = "") 
  assign(nam, ht)
  
  # save plot
  pdf(file = paste0(od,"fig5b_heatmap_",fig5b.names[i],".pdf"), height = 1.25, width = 1)
  draw(ht)
  dev.off()
  
}

htg.lc = c(sc.cols[3], sc.cols[4], sc.cols[5], sc.cols[9])
# make cluster legend
htg.l = Legend(legend_gp = gpar(fill = htg.lc),
               title = "Cluster",
               at = c(2,3,4,8),
               grid_width = unit(2.5, "mm"),
               grid_height = unit(2.5,"mm"),
               title_gp = gpar(fontsize = 6, fontface = "bold"),
               labels_gp = gpar(fontsize = 6, fontface = "bold"),
               direction = "vertical",
               ncol = 1,
               title_position = "topcenter",
               legend_height = unit(2,"cm"))
# plot legend & save
pdf(paste0(od,"fig5b_legend.pdf"), height = 2, width = 0.5)
draw(htg.l)
dev.off()

# scale bar
htg.s = Legend(col_fun = col_fun, title = NULL, at = c(-1, 2), labels = c("min","max"),
               legend_height = unit(2, "cm"), labels_gp = gpar(fontsize = 6, fontface = "bold"),
               grid_width = unit(2.5,"mm"))

# save scalebar
pdf(paste0(od,"fig5b_scale.pdf"), height = 2, width = 0.5)
draw(htg.s)
dev.off()


#### FIGURE S5B ####
prog <- c('S1pr1','Klf2','Id3','Ccr7')
act <- c('Tox','Icos','Tnfrsf4','Cd160')
tol <- c('Dapl1','Cd83','Egr2','Cd200')
tc17 <- c('Ramp1','Tmem176a','Igfbp7','Maf')
modules = c(prog,act,tol,tc17)
for (i in modules) {
  for (gene in i) {
    p <- FeaturePlot(object=lung.ms, reduction="umap", features=gene, order = T, pt.size = 0.25) &
      theme_void() &
      labs(title = NULL) &
      theme(plot.tag = NULL, plot.subtitle = NULL, strip.text.x = element_text(size=8, face = "bold"),
            legend.title = element_blank(),
            legend.key.height = unit(0.25,"cm"),
            legend.text = element_text(size=6, face = "bold"),
            legend.key.width = unit(0.5,"line"),
            legend.direction = "horizontal",
            plot.title = element_text(size=8,face="bold",hjust=0.5)) &
      ggtitle(gene) &
      scale_colour_gradientn(colours = 
                               c("gray75","gray87","gray89","palegoldenrod","lightgoldenrod","goldenrod","darkorange3","brown"))
    
    ggsave2(filename = paste0(od,"figS5A_umap_",gene,".pdf"),
            plot = (p + theme(legend.position = 'none')),
            device = "pdf",
            units = "cm",
            width = 3.5,
            height = 3.5,
            scale = 1)
    # legend
    pl.l <- cowplot::get_legend(p)
    pdf(file = paste0(od,"figS5A_umap_legend_",i,".pdf"), width = 1, height = 1.5)
    grid.newpage()
    grid.draw(pl.l)
    dev.off()
  }
}


#### FIGURE 5D ####

cds <- readRDS(file="201009_cds_monocle3.rds")

fig5d <- c("Ccr6","Rorc","Il17a")
for (gene in fig5d) {
  p <-  plot_cells(cds,
                   reduction_method = "UMAP",
                   genes = gene,
                   show_trajectory_graph = TRUE,
                   label_groups_by_cluster = FALSE,
                   label_branch_points = F,
                   label_leaves = F,
                   label_roots = F,
                   label_cell_groups = F,
                   trajectory_graph_segment_size = .25,
                   cell_size = 0.5) + theme_void() +
    theme(plot.tag = NULL, plot.subtitle = NULL, strip.text.x = element_text(size=8, face = "bold"),
          legend.title = element_blank(),
          legend.key.height = unit(0.25,"cm"),
          legend.text = element_text(size=6, face = "bold"),
          legend.key.width = unit(0.5,"line"),
          legend.direction = "horizontal",) +
    scale_colour_gradientn(colours = c("gray75","palegoldenrod","goldenrod","darkorange3","firebrick"),
                           labels = c(0, 0.5, 1),
                           breaks = c(0,0.5,1))
  
  ggsave2(filename = paste0(od,"fig5d_umap_",gene,".pdf"),
          plot = (p + theme(legend.position = 'none')),
          device = "pdf",
          units = "cm",
          width = 3.5,
          height = 3.5,
          scale = 1)
  pl.l <- cowplot::get_legend(p)
  pdf(file = paste0(od,"fig5d_umap_legend_",gene,".pdf"), width = 1.5, height = 0.5)
  grid.newpage()
  grid.draw(pl.l)
  dev.off()
}


#### fIGURE 5E #### 
## this requires having caluclated the clonotype clusters & having those variables in your environment
lung@meta.data$clonotype_cluster_all = revalue(lung@meta.data$raw_clonotype_id, phcl.hca2) # assigns clonotype cluster to each cell
lung@meta.data$clonotype_cluster_all = as.character(lung@meta.data$clonotype_cluster_all) # make sure its not an integer


an_colors = list(Antigen = c(siin = "#F8766D", siy = "#00BFC4"),
                 Clonotype_Cluster = c('1' = "#1F78B4",
                                       '2' = "#B2DF8A",
                                       '3' = "#33A02C",
                                       '4' = "#FB9A99",
                                       '5' = "#E31A1C",
                                       '6' = "#FDBF6F",
                                       '7' = "#FF7F00",
                                       '8' = "#CAB2D6",
                                       '9' = "#6A3D9A",
                                       '10' = "#FFFF99"),
                 Clonotype_Cluster_All = c('1' = "#222222",
                                           '2' = "#F6222E",
                                           '3' = "#FE00FA",
                                           '4' = "#16FF32",
                                           '5' = "#3283FE",
                                           '6' = "#FEAF16",
                                           '7' = "#B00068",
                                           '8' = "#1CFFCE",
                                           '9' = "#90AD1C",
                                           '10' = "#2ED9FF",
                                           '11' = "#DEA0FD",
                                           '12' = "#AA0DFE"))


lung <- SetIdent(lung,value='clonotype_cluster_all')

# plotting 
cl <- WhichCells(lung,idents=c("7"))

p = DimPlot(object=lung,reduction="umap",cells.highlight=cl,group.by='clonotype_cluster_all',pt.size = 0.005,sizes.highlight = 0.1) +
  guides(fill=guide_legend(title=NULL)) +
  theme_void() +
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        legend.key = element_rect(fill = "transparent", colour = NA)) +
  scale_color_manual(labels = c("other", paste0("Cluster ", as.character("7"))),
                     values = c("grey", an_colors$Clonotype_Cluster_All['7']))

ggsave2(filename = paste0(od,"fig5e_umap_clonoCluster7",".pdf"),
        plot = (p+ ggtitle(NULL)),
        device = "pdf",
        units = "in",
        width = 0.5,
        height = 0.5,
        scale = 1.75,
        bg = "transparent")

