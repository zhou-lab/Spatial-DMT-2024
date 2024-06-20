#rm(list=ls())
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(Nebulosa)
library(FigR)
library(BuenColors)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)
library(org.Mm.eg.db)
setwd("/mnt/nas1/Users/Yanxiang/Processed_data/2024/04032024_SpMETRTE14RXT8/SpMETRTE14RXT8_merge_1/output/align/Solo.out/Gene/WNN")
source('./scripts/SpatialPlot_new.R')
## RNA matrix
RNA_dir <- '../RNA/raw'
gene_matrix <- Read10X(data.dir = RNA_dir)
##
## Change spatial barcodes in DNA matrix
bc_df <- data.frame(RNA_bc=colnames(gene_matrix), row.names = colnames(gene_matrix))
bc_df$DNA_bc <- str_replace_all(bc_df$RNA_bc, "C", "T")
rownames(bc_df) <- bc_df$DNA_bc
DNAm_dir <- '../DNAm/DNAm_data.rds'
sampleNames <- 'SpMETRTE14'
DNAm_matrix <- readRDS(DNAm_dir)
rownames(DNAm_matrix) <- sub("\\..*$", "", rownames(DNAm_matrix))
rownames(DNAm_matrix) <- paste0(substr(rownames(DNAm_matrix),9,16), substr(rownames(DNAm_matrix),1,8))
bc_df <- bc_df[rownames(DNAm_matrix),]
all(rownames(bc_df)==rownames(DNAm_matrix))
rownames(DNAm_matrix) <- bc_df$RNA_bc
DNAm_matrix <- t(DNAm_matrix)
##
## Spatial object
data.dir <- '../RNA/'
assay = "RNA"
filter.matrix = TRUE
slice = "slice1"
object <- CreateSeuratObject(counts = gene_matrix, assay = assay)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image
spatial.obj <- object
spatial.obj <- subset(spatial.obj, cells=colnames(DNAm_matrix))
spatial.obj
##
## Spatial RNA analysis
p <- VlnPlot(spatial.obj, features = "nFeature_RNA", pt.size = 0.1, log = TRUE) + NoLegend()
p
median(spatial.obj$nFeature_RNA)
p <- SpatialFeaturePlot(spatial.obj, features = "nFeature_RNA",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p
DefaultAssay(spatial.obj) <- 'RNA'
spatial.obj <- NormalizeData(spatial.obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')
all.genes <- rownames(spatial.obj)
spatial.obj <- ScaleData(spatial.obj, features = all.genes, assay = 'RNA')
spatial.obj <- SCTransform(spatial.obj, assay = "RNA", verbose = FALSE)
spatial.obj <- RunPCA(spatial.obj, assay = "SCT", verbose = FALSE)
spatial.obj <- FindNeighbors(spatial.obj, reduction = "pca", dims = 1:30)
spatial.obj <- FindClusters(spatial.obj, verbose = FALSE, resolution = 0.6)
spatial.obj <- RunUMAP(spatial.obj, reduction = "pca", dims = 1:30)
DefaultAssay(spatial.obj) <- 'SCT'
n_clusters <- length(unique(spatial.obj$SCT_snn_res.0.6))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
spatial.obj$SCT_snn_res.0.6 <- paste0('R', spatial.obj$SCT_snn_res.0.6)
names(cols) <- paste0('R', n_clusters-seq_len(n_clusters))
cols['R12'] <- "#FEE500"
cols['R10'] <- "#D51F26"
cols['R7'] <- "#208A42"

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'SCT_snn_res.0.6', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

DimPlot(spatial.obj, reduction = 'umap', repel = TRUE, group.by = 'SCT_snn_res.0.6', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()

de_markers_RNA <- FindMarkers(spatial.obj, ident.1 = 2, only.pos = TRUE, group.by = 'wsnn_res.1') #ident.1 = 6
de_markers_RNA %>% dplyr::filter(p_val_adj < 0.05) -> de_markers_RNA_sel

p3 <- SpatialFeaturePlot(object = spatial.obj, features = 'Rab3c', alpha = c(0.1, 1), ncol = 3, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) # features = rownames(de_markers_RNA_sel)[1:5]
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3

## Spatial RNA analysis WNN integration
spatial.obj[['DNAm']] <- CreateAssayObject(counts = DNAm_matrix)
DefaultAssay(spatial.obj) <- "DNAm"
VariableFeatures(spatial.obj) <- rownames(spatial.obj[["DNAm"]])
spatial.obj@assays$DNAm@scale.data <- as.matrix(spatial.obj@assays$DNAm@counts)
spatial.obj <- RunPCA(spatial.obj, assay = "DNAm", reduction.name = 'mpca', verbose = FALSE)
ElbowPlot(spatial.obj)
spatial.obj <- FindNeighbors(spatial.obj, reduction = "mpca", dims = 1:10)
spatial.obj <- FindClusters(spatial.obj, verbose = FALSE, resolution = 1)
spatial.obj <- RunUMAP(spatial.obj, reduction = "mpca",  reduction.name = "meth.umap", dims = 1:10)

# plot spatial clusters DNA
n_clusters <- length(unique(spatial.obj$DNAm_snn_res.1))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
spatial.obj$DNAm_snn_res.1 <- paste0('D', spatial.obj$DNAm_snn_res.1)
names(cols) <- paste0('D', n_clusters-seq_len(n_clusters))
cols['D4'] <- "#D51F26"
cols['D6'] <- "#208A42"

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'DNAm_snn_res.1', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

DimPlot(spatial.obj, reduction = 'meth.umap', repel = TRUE, group.by = 'DNAm_snn_res.1', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()

## WNN integration
spatial.obj <- FindMultiModalNeighbors(
  spatial.obj, reduction.list = list("pca", "mpca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "RNA.weight"
)
spatial.obj <- RunUMAP(spatial.obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
spatial.obj <- FindClusters(spatial.obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 1)

# plot spatial clusters WNN
n_clusters <- length(unique(spatial.obj$wsnn_res.1))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
cols
spatial.obj$wsnn_res.1 <- paste0('W', spatial.obj$wsnn_res.1)
names(cols) <- paste0('W', n_clusters-seq_len(n_clusters))
cols['W13'] <- "#C06CAB"
cols['W6'] <- "#D51F26"
cols['W9'] <- "#D8A767"
cols['W5'] <- "#F47D2B"

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'wsnn_res.1', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1
DimPlot(spatial.obj, reduction = 'wnn.umap', repel = TRUE, group.by = 'wsnn_res.1', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()

## RNA marker analysis
DefaultAssay(spatial.obj) <- 'SCT'
Idents(spatial.obj) <- 'wsnn_res.1'
de_markers_RNA <- FindMarkers(spatial.obj, ident.1 = 6, only.pos = TRUE) #6
feature_RNA <- 'Myh7'
p3 <- SpatialFeaturePlot(object = spatial.obj, features = feature_RNA, ncol = 1, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q95", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3
##

## DNAm marker analysis
DefaultAssay(spatial.obj) <- 'DNAm'
Idents(spatial.obj) <- 'wsnn_res.1'
de_markers_DNA <- FindMarkers(spatial.obj, logfc.threshold = 0.1, ident.1 = 6) #6
feature_DNA <- 'chr8.85379165.85382665' #'chr8.85379165.85382665'
p3 <- SpatialFeaturePlot(object = spatial.obj, features = feature_DNA, ncol = 1, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q95", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3
##

## Smooth data
cellkNN <- spatial.obj@neighbors$weighted.nn@nn.idx
dim(cellkNN)
rownames(cellkNN) <- spatial.obj@neighbors$weighted.nn@cell.names
DNAmat.s <- spatial.obj@assays$DNAm@data
DNAmat.s <- DNAmat.s[unique(all_markers_DNA$gene),] # see calculation below
dim(DNAmat.s)
DNAmat.s <- smoothScoresNN(NNmat = cellkNN, mat = DNAmat.s, nCores = 10)
RNAmat.s <- spatial.obj@assays$RNA@data
dim(RNAmat.s) # Genes x Cells
RNAmat.s <- RNAmat.s[Matrix::rowSums(RNAmat.s)!=0,]
dim(RNAmat.s) # Genes x Cells
RNAmat.s <- smoothScoresNN(NNmat = cellkNN,mat = RNAmat.s,nCores = 10)

# create a new assay to store smooth information
DNAm_s_assay <- CreateAssayObject(counts = DNAmat.s)
spatial.obj[["DNAm_smooth"]] <- DNAm_s_assay
DefaultAssay(spatial.obj) <- 'DNAm_smooth'
feature_DNA <- 'chr2.180331404.180333604'
p <- SpatialFeaturePlot(object = spatial.obj, slot = 'counts', features = feature_DNA, ncol = 1, 
                        pt.size.factor = 4.5, image.alpha = 0,  max.cutoff = "q95", stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

RNA_s_assay <- CreateAssayObject(counts = RNAmat.s)
spatial.obj[["RNA_smooth"]] <- RNA_s_assay
##


## Integrative analysis
# map VMR to nearest genes, one cluster
source('./scripts/utils.R')
VMR_gene_df <- Methylation_genes(rownames(de_markers_DNA))
de_markers_DNA_sel <- de_markers_DNA[rownames(VMR_gene_df),]
all(rownames(de_markers_DNA_sel)==rownames(VMR_gene_df))
de_markers_DNA_sel$geneSymbols <- VMR_gene_df$geneSymbols
de_markers_DNA_sel <- de_markers_DNA_sel[order(de_markers_DNA_sel$avg_log2FC, decreasing = FALSE), ]
overlap_genes <- intersect(rownames(de_markers_RNA_sel), de_markers_DNA_sel$geneSymbols)
de_markers_DNA_plot <- de_markers_DNA_sel[de_markers_DNA_sel$geneSymbols %in% overlap_genes, ]
de_markers_RNA_plot <- de_markers_RNA_sel[rownames(de_markers_RNA_sel) %in% overlap_genes, ]
de_markers_DNA_plot <- de_markers_DNA_plot[order(de_markers_DNA_plot$p_val_adj, decreasing = FALSE), ]
de_markers_RNA_plot <- de_markers_RNA_plot[order(de_markers_RNA_plot$p_val_adj, decreasing = FALSE), ]
DoHeatmap(spatial.obj,  assay='DNAm_smooth', slot = 'counts', disp.min=NULL, features=rownames(de_markers_DNA_plot)[1:10]) + 
  scale_fill_gradientn(colors = jdb_palette("solar_extra"))
  #scale_fill_gradientn(colors = c("blue", "white", "red")) # + NoLegend()
DoHeatmap(spatial.obj,  assay='RNA', features=de_markers_DNA_plot$geneSymbols[1:10]) + scale_fill_gradientn(colors = jdb_palette("solar_extra"))
DoHeatmap(spatial.obj,  assay='RNA', features=rownames(de_markers_RNA_plot)[1:10]) + scale_fill_gradientn(colors = jdb_palette("brewer_yes"))
DoHeatmap(spatial.obj,  assay='DNAm_smooth', disp.min=NULL, 
          features=rownames(de_markers_DNA_plot[de_markers_DNA_plot$geneSymbols %in% rownames(de_markers_RNA_plot)[1:10], ])) + 
          scale_fill_gradientn(colors = jdb_palette("solar_extra"))

# Methylation level comparison
meth_dir <- "../DNAm/methylation_fractions.csv.gz"
meth_frac <- read.csv(meth_dir, header = TRUE, sep = ",", quote = "\"", row.names = 1, comment.char = "", stringsAsFactors = FALSE)
rownames(meth_frac) <- sub("\\..*$", "", rownames(meth_frac))
rownames(meth_frac) <- paste0(substr(rownames(meth_frac),9,16), substr(rownames(meth_frac),1,8))
meth_frac <- meth_frac[rownames(bc_df),]
all(rownames(bc_df)==rownames(meth_frac))
rownames(meth_frac) <- bc_df$RNA_bc
Idents(spatial.obj) <- 'wsnn_res.1'
cells_g1 <- Cells(subset(x = spatial.obj, idents = "6"))
cells_g2 <- Cells(subset(x = spatial.obj, idents = "6", invert = TRUE))
meth_frac_g1 <- meth_frac[cells_g1, features, drop=FALSE]
dim(meth_frac_g1)
meth_frac_g1 <- na.omit(meth_frac_g1)
colnames(meth_frac_g1) <- 'feature'
dim(meth_frac_g1)
meth_frac_g1$group <- 'cluster1'
meth_frac_g2 <- meth_frac[cells_g2, features, drop=FALSE]
dim(meth_frac_g2)
meth_frac_g2 <- na.omit(meth_frac_g2)
colnames(meth_frac_g2) <- 'feature'
dim(meth_frac_g2)
meth_frac_g2$group <- 'cluster2'
meth_frac_plot <- rbind(meth_frac_g1, meth_frac_g2)
p <- ggplot(meth_frac_plot, aes(x=group, y=feature, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white", outlier.size = -1)+
  labs(y = "Methylation level")
p <- p + theme_classic() + #scale_x_discrete can be used to change the order of items
  scale_x_discrete(limits=unique(meth_frac_plot$group), labels=labels) +
  theme(axis.title.x=element_blank(), legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size = 20))
p

p <- ggplot(meth_frac_plot, aes(x = group, y = feature, fill = group))+
            geom_boxplot()+
            geom_jitter(height = 0, width = .1)+
            scale_x_discrete(name = "group") + # note the x-axis is discrete
            scale_y_continuous(name = "Methylation level")+
            scale_fill_discrete(guide = FALSE) # this suppresses the legend because we don't need it
p <- p + theme_classic() + #scale_x_discrete can be used to change the order of items
  scale_x_discrete(limits=unique(meth_frac_plot$group), labels=labels) +
  theme(axis.title.x=element_blank(), legend.position = "none",
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size = 20))
p

DoHeatmap(spatial.obj, features = de_markers_RNA) + NoLegend()

