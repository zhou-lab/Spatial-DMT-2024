#rm(list=ls())

library(ggplot2)
library(dplyr)
library(Seurat)
#library(ArchR)
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

# png(filename = 'nFrags.png', width = 1200, height = 1200, res = 300)
# p
# dev.off()

p <- SpatialFeaturePlot(spatial.obj, features = "nFeature_RNA",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p
# png(filename = paste0(sampleNames, '_nFrags_spatial.png'), width = 2400, height = 2400, res = 300)
# p
# dev.off()
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

#cols <- c('10'="#a1d99b",'9'= "#636363",'0'="#FEE100",'8'= "#89288F",'7'= "#E44CAB",'6'= "#D51F26",'5'= "#F47D2B",'4'= "#272E6A",'3'= "#208A42",'2'= "#89C75F", '1'="#90D5E4")

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'SCT_snn_res.0.6', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

DimPlot(spatial.obj, reduction = 'umap', repel = TRUE, group.by = 'SCT_snn_res.0.6', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()

# png(filename = paste0(sampleNames, '_clusters_spatial.png'), width = 2400, height = 2400, res = 300)
# p1
# dev.off()

de_markers_RNA <- FindMarkers(spatial.obj, ident.1 = 2, only.pos = TRUE, group.by = 'wsnn_res.1') #ident.1 = 6
de_markers_RNA %>% dplyr::filter(p_val_adj < 0.05) -> de_markers_RNA_sel
#de_markers_RNA_sel$pct_diff <- de_markers_RNA_sel$pct.1 - de_markers_RNA_sel$pct.2

p3 <- SpatialFeaturePlot(object = spatial.obj, features = 'Rab3c', alpha = c(0.1, 1), ncol = 3, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) # features = rownames(de_markers_RNA_sel)[1:5]
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3


## Spatial RNA analysis WNN integration
spatial.obj[['DNAm']] <- CreateAssayObject(counts = DNAm_matrix)

DefaultAssay(spatial.obj) <- "DNAm"

VariableFeatures(spatial.obj) <- rownames(spatial.obj[["DNAm"]])
# spatial.obj@assays$Spatial@var.features <- rownames(spatial.obj)
# spatial.obj <- FindVariableFeatures(spatial.obj, selection.method = "vst", nfeatures = 3000)
# VariableFeaturePlot(spatial.obj)

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

#cols <- c('10'="#a1d99b",'9'= "#636363",'0'="#FEE100",'8'= "#89288F",'7'= "#E44CAB",'6'= "#D51F26",'5'= "#F47D2B",'4'= "#272E6A",'3'= "#208A42",'2'= "#89C75F", '1'="#90D5E4")

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'DNAm_snn_res.1', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

DimPlot(spatial.obj, reduction = 'meth.umap', repel = TRUE, group.by = 'DNAm_snn_res.1', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()
##


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

#cols <- c('10'="#a1d99b",'9'= "#636363",'0'="#FEE100",'8'= "#89288F",'7'= "#E44CAB",'6'= "#D51F26",'5'= "#F47D2B",'4'= "#272E6A",'3'= "#208A42",'2'= "#89C75F", '1'="#90D5E4")

p1 <- SpatialDimPlot(spatial.obj, label = FALSE, group.by = 'wsnn_res.1', label.size = 3,  pt.size.factor = 4.5, image.alpha = 0, stroke = 0, cols = cols) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) + DarkTheme() + NoAxes() + NoGrid()
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

#DimPlot(spatial.obj, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5)
DimPlot(spatial.obj, reduction = 'wnn.umap', repel = TRUE, group.by = 'wsnn_res.1', cols = cols, pt.size = 2) + DarkTheme() + NoGrid()

# VlnPlot(spatial.obj, features = "SCT.weight", group.by = 'wsnn_res.1', sort = TRUE, pt.size = 0.1) + NoLegend()
# VlnPlot(spatial.obj, features = "DNAm.weight", group.by = 'wsnn_res.1', sort = TRUE, pt.size = 0.1) + NoLegend()
##

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
#DNAmat.s <- DNAmat.s[rownames(de_markers),]
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
  

# map VMR to nearest genes, all clusters
# DefaultAssay(spatial.obj) <- 'DNAm'
# Idents(spatial.obj) <- 'wsnn_res.1'
# 
# all_markers_DNA <- FindAllMarkers(spatial.obj)
# all_markers_DNA %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC < -5) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10_DNA




# Correlation between meth and RNA
FeatureScatter(spatial.obj, feature1 = "dnamsmooth_chr1.172014428.172016628", feature2 = "rnasmooth_Vangl2")


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


# Heatmap

DoHeatmap(spatial.obj, features = de_markers_RNA) + NoLegend()

##
















































########################################## old code for reference
class(de_markers)
#save(de_markers,file="SpnbRNA53.Rdata")
load("SpnbRNA53.Rdata")
spatial.obj <- FindSpatiallyVariableFeatures(spatial.obj, assay = "SCT", features = VariableFeatures(spatial.obj)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(spatial.obj, selection.method = "moransi"), 1)



markers_cluster <- as.data.frame(de_markers_6)
markers_cluster

features_spatial <- c('Adgrl3', 'Ddx39b', 'Cd24a', 'Cdh2', 'Sept3', 'Zmiz1', 'Celf2', 'Serf2', 'G3bp2', 'H3f3a', 'Syncrip', 'Pcdh9', 
'Eef1b2', 'Kif5b', 'Zfp462', 'Nrxn1', 'Apc', 'Ubn2', 'Tubb3', 'Zfp618', 'Whsc1', 'Bptf', 'Tnrc6c', 'Rps28', 'Gria2', 'Hmgcs1', 'Gm37899',
'Mex3a', 'Rbfox2', 'Mllt3', 'Kmt2c', 'Elavl3', 'Ogt', 'Raph1', 'Sorbs2', 'Meis1', 'Crmp1', 'Rufy3', 'Calm2', 'Tcaf1', 'Msi2', 'Ube2d3',
'Epha5', 'Rps18', 'Eif4g3', 'Nrxn2', 'Slc1a2', 'Cbfa2t3', 'Elavl2', 'Rtn1', 'Plxna2', 'Rps14', 'Olfm1', 'Kif5c', 'Ank2', 'Bcl11b', 'Robo3',
'Rmst', 'Bcl11a', 'Celf4', 'Cyp26b1', 'Chl1', 'Stmn2', 'Hbb-bt', 'Igf1', 'Ptprd', 'Map2',
'Nrxn3', 'Runx1t1', 'Meg3', 'Mapt', 'Gap43', 'Meis2', 'Ank3', 'Zfhx3', 'Tuba1a', 'Dcx', 'Dpysl3', 'Sox11', 'Ttn', 'Snhg11', 'Ina')
features_spatial <- c('Gm42418', 'mt-Rnr1', 'mt-Rnr2', 'H19', 'Col1a2', 'Igf2', 'Sox11', 'Col4a1', 'Col3a1', 'Malat1', 'Hba-a1', 'Lars2', 'Rbp1', 
                      'Atp1a2', 'Arhgap29', 'Rps8', 'Gja1', 'Akap12', 'Lgals1', 'Fstl1', 'Col4a2', 'Meis2', 'Hba-x', 'Rpl13', 'Rpl41', 'Sparc', 'Vcan',
                      'Rps12', 'Rpl36', 'Col1a1', 'Rplp1', 'Ina', 'Ncam1', 'Akap9', 'Snhg11', 'Rrbp1', 'Map1b', 'Rpl37a', 'Cdkn1c', 'Map2', 'Ank3', 'Rps20',
                      'Rps26', 'Vim', 'Auts2', 'Nrxn3', 'Atxn7l3b', 'Tubb2b', 'Dcx', 'Tuba1a', 'Mir124-2hg', 'Igfbp5', 'Rmst', 'Sept3', 'Bcl11a', 'Kif1b', 'Rtn1',
                      'Zfhx3', 'Meg3', 'Basp1', 'Dpysl3', 'Mapt', 'Stmn2',
                      'Gsk3b', 'Ank2', 'Bcl11b', 'Nnat', 'Celf4', 'Elavl3', 'Ebf1', 'Syt11', 'Runx1t1', 'Cyp26b1')
features_spatial <- c('Mir124-2hg', 'Ina', 'Map2', 'Sox11', 'Dcx', 'Snhg11', 'Rtn1', 'Gria2', 'Rmst', 'Chl1', 'Igfbpl1', 'Dpysl3', 'Ank3',
                      'Elavl3', 'Nrxn3', 'Ncam1', 'Celf4', 'Tuba1a', 'Map1b', 'Scn3a', 'Thsd7a', 'Tubb2b', 'Nrxn2', 'Stmn2', 'mt-Rnr2', 'Nrxn1', 'Sept3',
                      'Slc1a2', 'Mapt', 'Bcl11b', 'Myt1l', 'Slc17a6', 'Pou2f2', 'Grin2b', 'Rbfox1', 'mt-Rnr1', 'Robo3', 'Crmp1', 'Meg3', 'Malat1', 'Tubb3', 
                      'Runx1t1', 'Rnf165', 'Rab3c', 'Snhg14', 'H19', 'Rufy3', 'Basp1', 'Elavl4', 'Ank2', 'Cdk5r1', 'Apc', 'Plxna2', 'Cbfa2t3', 'Olfm1', 
                      'Sv2a', 'Elavl2', 'Kif1b', 'Syt1', 'Dusp8', 'Nnat', 'Bcl11a', 'Gnao1', 'Zfhx3', 'Kif5c', 'Csrnp3', 'Epha5', 'Draxin', 'Onecut2',
                      'Tfap2b', 'Nova2', 'Igf2', 'Col1a2', 'Lrp8', 'Meis2', 'Nfasc', 'Gdpd1', 'Pcdh9', 'Rbfox3', '9330159F19Rik', '9330162G02Rik', 'Gad2', 
                      'Kcnq2', 'Slco5a1', 'Akap9', 'Nrcam', 'Pcsk1n', 'Hecw1', 'Gdap1', 'Nsg2', 'Cdh2', 'Nsg1', 'Celf2', 'Sgip1', 'Rims2', 'Ephb1', 'Srrm3', 
                      'Cmip', 'Agap1', 'Celf5', 'Pid1', 'Tmem178b', '5330434G04Rik', 'Nbea', 'Scn2a1', 'Jph4', 'Hbb-y', 'Pak7', 'Soga3', 'St6galnac5', 
                      'Cntn2', 'Cacna1b', 'Cacnb4', 'Dclk1', 'Plxna4', 'Reep1', 'Dpysl2', 'A230103L15Rik', 'Mab21l2', 'Sbk1', 'Meis1', 'Klf7', 'Phyhipl',
                      'Pclo', 'Gabrg2', 'Rab6b', 'Kidins220', 'Cbln2', 'Aplp1', 'Sorbs2', 'Car10', 'Dnm3', 'Sox4', 'Ccnd2', 'Unc79', 'St8sia4', 'Adgrl3',
                      'Map6', 'Pax8', 'Dner', 'B230334C09Rik', 'Nav1', 'Gm43889', 'Fnbp1l', 'Xpr1', 'Gap43', 'Pou3f1', '2900027M19Rik', 'Jakmip2', 'Ncan',
                      'Cep170', 'Fstl1', 'Mir124a-1hg', 'Lsamp', 'Sparc', 'Zfhx4', 'Tshz3', 'Tulp4', 'Fat3', 'Ap3b2', 'Kif1a', 'Zfhx2', 'Hdgfrp3', 'Gm44799',
                      'Tenm2', 'Shox2', 'Cacna2d2', 'Ctnna2', 'Myt1', 'Syt11', 'Atp2b2', 'Gabrb3', 'Zfp462', 'Stmn1', 'B3gat1', 'Cacng2', 'Ccdc177',
                      'C1qtnf4', 'Sptbn2', 'Cxadr', 'Enc1', 'Erc2', 'St8sia2', 'March1', 'Apc2', 'Xkr6', 'Pde10a', 'Dpf1', 'Add2', 'Zcchc18', 'Nefl',
                      'Rasal2', 'Rnf152', 'Cyp26b1', 'Cntnap2', 'Magi2', 'Gng2', 'L1cam', '1500004A13Rik', 'Gsk3b', 'Atp8a1', 'Cnr1', 'Astn1', 'Cadm1',
                      'Hba-a1', 'Robo2', 'Zbtb20', 'Ankrd12', 'Gnaq', 'Hmgcs1', 'Ebf1', 'Dcc', 'Bach2', 'Scg3', 'Bicd1', 'Nova1', 'Sez6l2', 'Ppp3ca', 
                      'Kif21b', 'Ttc3', 'Mtus2', 'Gm2694', 'Arpp21', 'Akap6', 'Nap1l5', 'Mapk8ip2', 'Ppm1l', 'Npdc1', 'Neurod2', 'Camk2n1', 'Dpysl4', 
                      'Cd24a', 'Fry', 'Ubn2', 'Cux2', 'Gm44433', 'Tubb2a', 'Lin7a', 'Atp1a3', 'Lamp5', 'Clcn4', 'Dpysl5', 'Kif21a', 'Cyfip2', 'Cask', 
                      'Gria4', 'Phox2b', 'Mycbp2', 'Tcaf1', 'Camsap1', 'Mcf2l', '6330415B21Rik', 'Gm44435', 'Igfbp5', 'Uncx', 'Kalrn', 'Mapk10', 'Kif5a', 
                      'Dclk2', 'Tspan13', 'Nav3', 'Raph1', 'Lhx1', 'Mllt3', 'Fbn2', 'Herc1', 'Atp1b1', 'Cpe', 'Ptpn5', 'Gpm6a', 'Cxxc4', 'Tmeff1', 'Anks1b',
                      'Syt4', 'Ntrk3', 'Col3a1', 'Ahi1', 'Vezt', 'Hoxb3', 'Hba-a2', 'Whsc1', 'Trim2', 'Fgd4', 'C2cd5', 'Celsr2', 'Asxl3', 'Plppr3', 'Atcay',
                      'Smpd3', 'Gm44235', 'Kifap3', 'Nxph1', 'Tanc2', 'Ntm', 'Ppip5k2', 'Thra', 'Scn3b', 'Zfp57', 'Adgrb3', 'Camta1', 'Stmn3', 'Mast1',
                      'Tmcc1', 'Mapre2', 'Lhfpl4', 'Gm26871', 'Mfsd6', 'D430041D05Rik', 'Tnrc6c', 'Abi2', 'Cacna2d1', 'E230020A03Rik', 'Arhgef12', 'Usp29',
                      'Samd14', 'Lars2', 'Ckb', 'Fam171b', 'Efr3b', 'Gatsl2', 'Ntng1', 'Strbp', 'Atp13a2', 'Sptan1', 'Pfkp', 'Tmod2', 'Prrc2b', 'Gm3764', 
                      'Pbx3', 'Ebf3', 'Gm38340', 'Phactr1', 'Wdfy1', 'Cacna1a', 'Zyg11b', 'Dopey1', 'Cpeb4', 'Ptprd', 'Rpl41', 'Zfp811', 'Rusc2', 'Tmem44',
                      'Abr', 'Lrba', 'Pgm2l1', 'Mllt11', 'Pcdh7', 'Polb', 'Dync1i2', 'Hspa12a', 'Ywhag', 'Nhlh2', 'Spag9', 'Dmxl2', 'Sestd1', 'Sorl1', 
                      'Trp53i11', 'Tubb5', 'Pak3', 'App', 'Gm38042', 'Mex3a', 'Ptbp2', 'Ttc28', 'Miat', '0610010F05Rik', 'Csnk1e', 'G3bp2', 'Klc1', 'Cdkn1c',
                      'Inpp5f', 'Zfp60', 'Fzd3', 'Tnik', 'Rps8', 'Gabbr1', 'Dlgap4', 'Pcp4', 'Adgrl1', 'Auts2', 'Cers6', 'Fam117b', 'Rbm4b', 'Madd',
                      'E130308A19Rik', 'Dock7', 'C530008M17Rik', 'RP23-127M7.1', 'C130071C03Rik', 'Srgap3', 'Cacng4', 'Ppp2r3a', 'D10Wsu102e', 'Sbf2',
                      'Rtn4', 'Pik3r3', 'Nedd4l', 'Ago1', 'Zswim6', 'Napb', 'Dnm3os', 'Ssh2', 'Kif3a', 'Actb', 'Hook3', 'Hn1', 'Kmt2e', 'Plekha1', 'Rnd3',
                      'Tmem57', 'Sh3bp5', 'Trp53bp1', 'Npepps', '4930402H24Rik', 'Smarcd1', 'Gpc3', 'Krit1', 'Rtn3', 'Serinc1', 'Scn8a', 'Tmx4', 'Arhgap21', 
                      'Tshz2', 'Gdi1', 'Gm42992', 'Slc38a1', 'Clcn3', 'Ccdc88a', 'Magi1', 'Pnmal2', 'Tspan5', 'Tia1', 'Cdh4', 'Acvr2a', 'Enah', 'Bcr', 'Mdk',
                      'Macf1', 'Frmd4a', 'Myo5a', 'Atxn7l3b', 'Zfp618', 'Ptpn1', 'Pkia', 'Fsd1l', 'Uba1', 'Kdm7a', 'Ndn', 'Abcc5', 'Atp9a', 'Usp22', 'Idh1', 
                      'Pabpn1', 'Arhgef2', 'Srgap1', 'Dbn1', 'Cpsf6', 'Lcorl', 'Rbfox2', 'Socs2', 'Trim8', 'Dach1', 'Sugp2', '1700025G04Rik', 'Uchl1', 'Maf',
                      'Mir99ahg', 'Mapre1', 'Top2b', 'Rc3h2', 'Acvr2b', 'Hipk2', 'Fbrsl1', 'Zmiz1', 'Dynlt1c', 'Mex3b', 'Kdm6b', 'Zfp292', 'Nefm', 'Afap1', 
                      'Wsb1', 'Tbc1d16', 'Zfp638', 'Mapk8', 'Rpl13', 'Sobp', 'Kmt2c', 'Zfand5', 'Fign', 'Exoc5', 'Bzw2', 'Baz2b', 'Eif4g3', 'Zmynd8', 'Nav2',
                      'Lpgat1', 'Dennd5b', 'Zcchc14', 'Mecp2', 'Znrf1', 'Ttn', 'Zfp157', 'Itsn1', 'Zfp329', 'Ntrk2', 'Rb1cc1', 'Gapvd1', 'Peg3', 'Fyn', 'Sms', 'Jarid2', 'Wdfy3', 'Tbc1d24', 'Adgrl2', 'Ptprs', 'Actg1', 'Clasp1', 'Col2a1', 'Epc1', 'Pygo1', 'R3hdm1', 'Dhx36', 'Trove2', 'Rfx7', 'Calm2', 'Myo9a', 'Rnasel', 'Plekha5', 'Dnajc5', 'Fryl', 'Stox2', 'Fam168a', 'Vangl2', 'Zfp451', 'Zic1', 'H3f3a', 'Akap11', 'Khdrbs1', 'Arid2', 'Dync1h1', 'Phip', 'Rab6a', 'Zfp280d', 'Mau2', 'Rbm5', 'Kcnq1ot1', 'Calm1', 'Brwd1', 'Ssbp3', 'Atp2b1', 'Gm37899', '2410089E03Rik', 'Col4a1', 'Msi2', 'Rps19', 'Crabp1', 'Csnk2a1', 'Mtf2', 'Rbm28', 'H3f3b', 'Kdm5b', 'Ddx5', 'Tmem2', 'Ubb', 'Kmt2a', 'Map4k4', 'Rev3l', 'Rrbp1', 'Srrm2', 'Lmo4', 'Myef2', 'Tra2a', 'Col11a1', 'Pdgfra', 'Ilf2', 'Unc5c', 'Arglu1', 'Top2a', 'Ash1l', 'Pabpc1', 'Rn18s-rs5', 'Hspa5', 'Igf1', 'Rplp1', 'Chd4', 'Igfbp4', 'Lamb1', 'Rps20', 'Rpsa', 'Prrx1', 'Vcan', 'Ahnak', 'Serpinh1', 
                      'Hba-x', 'Hmga2', 'Tpm1', 'Hsp90ab1', 'Rps24', 'Tgfb2', 'Ptn', 'Grb10', 'Mest', 'Foxp1', 'Col27a1', 'Zeb2', 'Rpl37a', 'Hsp90b1', 'Lgals1',
                      'Hmcn1', 'Mki67', 'Hbb-bt', 'Rps27a', 'Gm26917', 'Rps3a1', 'Rpl17', 'Rpl36a', 'Rplp2', 'Tcf4', 'Rps16', 'Rps12')
features_spatial <- c('Dnm3os', 'Ccnd2', 'Crym', 'Lef1', 'Col2a1', 'Tenm4', 'Meg3', 'Cped1', 'Hba-a1', 'Hmga2', 'Sox9', 'Pou3f4', 'Ina', 'Six1', 'Col11a1', 'Mecom', 'Hbb-y', 'Hist1h2ap', 'Sfrp1', 'Map1b', 'Tubb2b', 'Fbn2', 'Col23a1', 'Tuba1a', 'Foxp1', 'Mki67', 'Col27a1', 'Top2a', 'Hba-a2', 'Ank3', 'Dcx', 'Tcf7l2', 'Ccnd1', 'Dpysl3', 'Snhg11', 'Mapt', 'Elavl3', 'Rtn1', 'Tmsb4x', 'Foxp2', 'Tubb3', 'Mmp14', 'Hjurp', 'Fbln2', 'Stmn2', 'Wwp2', 'Sema5a', 'Prrx1', 'Plod2', 'Sema6a', '2810417H13Rik', 'Gpc6', 'Vcan', 'Hba-x', 'Kif1b', 'Celf2', 'Hmcn1', 'Elavl2', 'Ncam1', 'Cntn2', 'Map2', 'Crmp1', 'Hnrnpr', 'Fstl1', 'Igsf3', 'Zbtb20', 'Calm1', 'Nrxn3', 'Gap43', 'Auts2', 'Nefl', 'Ttn', 'Gria2', 'Nucks1', 'Celf4', 'Rmst', 'Runx1t1', 'Sept3', 'Syt11', 'Nefm', 'Rian', 'Nnat', 'Sox11', 'Igfbp5', 'Hbb-bt', 'Chl1', 'Cyp26b1')
features_spatial <- c('Hba-a2', 'Hbb-y', 'Hba-a1', 'Hbb-bt', 'Hba-x', 'Hbb-bs', 'Hbb-bh1', 'Slc4a1', 'Sox11', 'Ttn', 'Nudt4', 'App', 'Auts2', 'Zfp462', 'Celf4', 'Meis2', 'Tubb2b', 'Mir99ahg', 'Prrc2b', 'Enah', 'Kif1b', 'Bcl11b', 'Akap9', 'Matr3', 'Map2', 'Thsd7a', 'Rtn1', 'Dcx', 'Runx1t1', 'Cyp26b1', 'Gria2', 'Pbx1', 'Dpysl3', 'Igf1', 'Nrxn3', 'Zfhx3', 'Ebf1', 'Nnat', 'Zbtb20', 'Col2a1', 'Bcl11a', 'Tcf4', 'Snhg11', 'Elavl3', 'Ncam1', 'Mapt', 'Ank3', 'Ina', 'Map1b')
features_spatial <- c('Kcnq1ot1', 'Igf1', 'Dnm3os', 'Peg3', 'Neb', 'Grem1', 'Igf2', 'Gas1', 'Col1a2', 'Rian', 'Map1b', 'Col3a1', 'Cdkn1c', 'Tnnt1', 'Clcn5', 'Palld', 'Pdgfra', 'Lpar1', 'C430049B03Rik', 'H19', 'Myh3', 'Igfbp5', 'Adamts9', 'Nrk', 'Snhg11', 'Ttn', 'Ina', 'Gm37899', 'Prrx1', 'Dlk1', 'Igfbp4', 'Fam101b', '9430076C15Rik', 'Tuba1a', 'Fzd2', 'Plagl1', 'Lamb1', 'Postn', 'Nfib', 'Ahnak', 'Parm1', 'Mdk', 'Ptn', 'Mef2c', 'Dpysl3', 'Tubb2b', 'Cdh11', 'Pbx1', 'Mylpf', 'Ncam1', 'Fstl1', 'Timp3', 'Rtn1', 'Stmn2', 'Fus', 'Basp1', 'Tpm1', 'Tgfb2', 'Dcx', 'Map2', 'Tnnc1', 'Mapt', 'Cxcl12', 'Nnat', 'Gpc3', 'Actc1', 'Sparc', 'Rmst', 'Tpm2', 'Gm26917', 'H2afy', 'Fbn2', 'Actb', 'Hnrnpr', 'Sorbs2', 'Fn1', 'Sox11', 'Elavl3', 'Timp2', 'Tcf4', 'Rrbp1', 'Rbms1', 'Top2a', 'Celf4', 'Chl1', 'Itgb1', 'Nrxn3', 'Thsd7a', 'Stmn1', 'Slc1a2', 'mt-Nd1', 'Mest', 'Elavl4', 'Nefl', 'Gria2', 'Mir124-2hg', 'Gap43', 'Epha5', 'Nrxn1', 'Tubb3', 'Zbtb20', 'Syt11', 'Cyp26b1', 'Cntn2', 'Kif1b', 'Hbb-y', 'Sept3', 'Hba-a2', 'Bcl11a', 'Col2a1', 'Hba-a1', 'Zfhx3', 'Col11a1')
features_spatial <- c('Gm42418', 'Map1b', 'Sox11', 'Lhx1', 'Cacng4', 'H19', 'Igf2', 'Ina', 'Mapt', 'Tubb2b', 'Igf2r', 'Rtn1', 'Tuba1a', 'Serpinh1', 'Nrxn3', 'Dcx', 'Mdk', 'Celf4', 'mt-Rnr1', 'Pdgfra', 'Mest', 'Grb10', 'Ahnak', 'Calu', 'Cdkn1c', 'Ncam1', 'Itgb1', 'Peg3', 'Stmn2', 'Nfib', 'Rpl8', 'Ttn', '2810417H13Rik', 'Ebf3', 'Runx1t1', 'Tgfb2', 'Vim', 'Kcnq1ot1', 'Rps18', 'Igf1', 'Sept3', 'Col3a1', 'Dnm3os', 'Fbn2', 'Map2', 'Igfbp5', 'Col1a2', 'Fstl1', 'Hbb-bt', 'Rian', 'Hba-x', 'Col11a1', 'Hba-a2', 'Hba-a1', 'Hbb-y')
features_spatial <- c('Hopx', 'Pax3', 'Ttyh1', 'Adgrv1', 'Fabp7', 'Sox2', 'Rfx4', '2900037B21Rik', 'Gm3764', 'Wnt7b', 'Gm37608', 'Fgf15', 'Notch1', 'Hes5', 'Fgfr3', 'A330074H02Rik', 'Ascl1', 'Mir9-3hg', 'Cyp26b1', 'AW047730', 'Gm37812', 'Nckap5', 'Celsr2', 'Ednrb', 'Ptprz1', 'Lrp2', 'C130071C03Rik', 'Tenm2', 'Kcnq3', 'Pla2g7', 'Sox2ot', 'E130114P18Rik', 'Kcnk10', 'Ildr2', 'Pou3f2', 'Pea15a', 'Adcyap1r1', 'Vcam1', 'Fzd10', 'Qk', 'Irx2', 'Gm37457', 'Gpm6b', 'Kif21a', 'Itgb8', 'Epha3', 'Jam2', 'Adgra1', 'Nfib', 'Prex1', 'Adgrg1', 'Wscd1', 'Wwc1', 'Msx3', 'Ndnf', 'Ncald', 'Tcf4', 'Phyhipl', 'Nfia', 'Cecr2', 'Pantr1', 'Ckb', 'Pax6', 'Msi2', 'Ptn', 'Slc1a3', 'Lhfp', 'Pgpep1', 'Grid2', 'Tox3', 'Slc35f1', 'Meg3', 'Msi1', 'Sox9', 'Btg2', 'Sox21', 'Chd7', 'Zic1', 'Dtx4', 'Fzd3', 'Paqr8', 'Syt11', 'Marcks', 'Rassf4', 'Vangl2', 'Vim', 'Phgdh', 'Dbi', 'Pdpn', 'Wnt5a', 'Gm44645', 'Zbtb20', 'Npas3', 'Bcl2', 'Tcf12', 'Draxin', 'Gm37699', 'Agrn', 'Igf2', 'Ccnd2', 'Celf2', 'Kbtbd11', 'Mir99ahg', 'Prox1', 'mt-Rnr2', 'Ddr1', 'Epha4', 'Trp53inp1', 'Zic3', 'Cspg5', 'Cachd1', 'Plpp3', 'H19', 'Id1', 'Rian', 'Tubb2b', 'Zfp462', 'Col1a2', 'Serpine2', 'Zic4', 'Gnai2', 'Stt3b', 'Pou3f3', 'Rbl1', 'Hspa4l', 'Kif1b', 'Ddah1', 'Sema5b', 'Fbn2', 'Gm42418', 'Rnd2', 'Zfp536', 'Rragd', 'Sox11', 'Nkd1', 'Nras', 'Skp1a', 'Pou2f1', 'Cdh2', 'C130075A20Rik', 'Sall3', 'Col3a1', 'Erbb4', 'Setbp1', 'Akirin1', 'Akap6', 'Nfix', 'Sept3', 'Mllt4', 'Meis1', 'Zfp24', 'Rgma', 'Satb1', 'Cst3', 'Cenpf', 'Akap9', 'Pebp1', 'Lrp8', 'Aplp2', 'mt-Rnr1', 'Dnm3os', 'Trim59', 'Ttf1', 'Gli3', 'Spire1', 'Zic2', 'Fam120a', 'Gpc1', 'Golm1', 'Auts2', 'Ptprf', '4632427E13Rik', 'Cdc14a', 'Arl6ip1', 'Gsk3b', 'Igf2bp3', 'Rap2b', 'Miat', 'Hmgcs1', 'Srgap3', 'Zfp646', 'Mib1', 'Pgam1', 'Aspm', 'Hbb-y', 'Hes6', 'Limch1', 'Tmem47', 'Lrrc58', 'Map2', 'Zic5', 'Tshz1', 'Vezf1', 'Atp8a1', 'Hmgb2', 'Pkm', 'Gapdh', 'Fgfr1', 'Cdkn1b', 'Ybx1', 'Spred1', 'Hcfc1', 'Bcl7a', 'Rsf1', 'Nfasc', 'Slc38a1', 'Scd2', 'Brwd1', 'Smad4', 'Fbxw7', 'Rnf165', 'Dach1', 'Pogk', 'Tfdp2', 'Gpc3', 'Slc1a2', 'Lpcat1', 'Gm37298', 'Tmem178b', 'Fat3', 'C530008M17Rik', 'Bcl11a', 'Ywhaz', 'Rgs12', 'Trp53i11', 'Igsf8', 'Atp5e', 'Uhrf2', 'Jmjd1c', 'Dctn2', 'Odf2', 'Stip1', 'Gm37606', 'Meis2', 'Kmt2e', 'Hsph1', 'Arid1a', 'Sec63', 'Snrpf', 'Sox5', 'Ttn', 'Igfbp2', 'Zmiz1', 'Fat4', 'Cdh4', 'Dek', 'Nipbl', 'Arvcf', 'Dclk1', 'Rasal2', 'Rsrc2', 'Tspan18', 'Psmc3', 'Igfbp5', 'Cbl', 'Eif4g3', 'H3f3a', 'H2afy', 'Rbmx', 'Eif4g2', 'Erbb2ip', 'Zfp811', 'Atf2', 'Zfhx3', 'Pom121', 'Nr2f2', 'Vps36', 'Tuba1a', 'Tfap2b', 'Zfp608', 'Rbx1', 'Top2a', 'Tsn', 'Plxnb1', 'BC005561', 'Ebf1', 'Socs2', 'Arhgap5', 'Hnrnpd', 'Nr2f1', 'Tmsb4x', 'Prrc2b', 'Ubxn4', 'Hmgb1', 'Tead2', 'Cdkn1c', 'Pak2', 'Arhgap21', 'Ubr5', 'Zfp871', 'Id3', 'Mbtd1', 'Ptprs', 'Trim24', 'Impad1', 'Ube2d3', 'Smarcb1', 'Reep3', 'Calm1', 'Smarca5', 'Amn1', 'Rmst', 'Postn', 'Atxn7l3b', 'Whsc1', 'Hist1h2ap', 'Hspd1', 'Set', 'Csde1', 'Trrap', 'Marcksl1', 'Hnrnpk', 'Jarid2', 'Cbx5', 'Dpysl2', 'Rc3h2', 'Igfbp4', 'Serbp1', 'C430049B03Rik', 'Peg3', 'Khdrbs1', 'Psip1', 'Tmem2', 'Clmp', 'Top1', 'Caprin1', 'Tead1', 'Cpsf6', 'Sptbn1', 'Ktn1', 'Tra2a', 'Igf1', 'Hba-a1', 'Hba-a2', 'Tnrc6b', 'Spag9', 'Canx', 'Ddx5', 'Sept11', 'Tia1', 'Cspp1', 'Rps11', 'Rps2', 'Ilf2', 'Pdgfra', 'Flrt2', 'Lars2', 'Kmt2a', 'Kctd12', 'Wdr26', 'Usp9x', 'Csnk1a1', 'Phip', 'Prrx1', 'Col11a1', 'Hsp90b1', 'Chd3', 'Plxna4', 'Ptp4a2', 'Dync1h1', 'Macf1', 'Cnot6', 'Mat2a', '2610203C20Rik', 'Hnrnpu', 'Maf', 'Fn1', 'Fstl1', 'Hmcn1', 'Rbms1', 'Foxp1', 'Gm26917', 'Hba-x', 'Tshz2', 'Tenm4', 'Ank3', 'Snhg11', 'Hbb-bt', 'Col4a1', 'Col2a1', 'Ncam1', 'Runx1t1', 'Mapt')

features_spatial <- markerGenes
feature <- features_spatial[10]
feature <- 'Thy-1'

p4 <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0,min.cutoff = "q5", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22) # set spots to square shape
p4

png(filename = paste0('./markers_plot_1/', feature, '_upReg_spatial.pdf'), width = 1200, height = 1200, res = 300)

plotPDF(p4, name = paste0(feature, '_clusters_spatial.pdf'),  addDOC = FALSE, width = 5, height = 5)
p4
dev.off()


##density plot

pd <- plot_density(spatial.obj, rownames(de_markers)[1], reduction = 'umap')
pd_data <- pd$data
pd_data <- pd_data[, c(3), drop=FALSE]

all(Cells(spatial.obj) == row.names(pd_data))
spatial.obj <- AddMetaData(object = spatial.obj, metadata = pd_data)




p <- SpatialFeaturePlot(spatial.obj, features = "feature",  pt.size.factor = 4.5, image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p


##
SpatialFeaturePlot(object = spatial.obj, features = 'Sox2', alpha = c(0.1, 1), ncol = 3, pt.size.factor = 4.5,  max.cutoff = "q90", image.alpha = 1, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))

SpatialFeaturePlot(object = spatial.obj, features = 'Car10', alpha = c(0.1, 1), ncol = 1, pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 1, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))


##################### Ignore the code below
p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4.5, image.alpha = 0, stroke = 0)
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

png(filename = paste0(sampleNames, '_clusters_spatial.png'), width = 3600, height = 3000, res = 300)
p1
dev.off()



# png(filename = paste0(sampleNames, '_clusters_umap.png'), width = 3600, height = 3000, res = 300)
# p2
# dev.off()

markers_cluster <- as.data.frame(markerList_pos$C3)
markers_cluster$name[1:5]

features_spatial <- c('Cd3e', 'Cd4', 'Prom1', 'Pdgfra', 'Ptprc', 'Cd34', 'Thy1')


features_spatial <- as.vector(row.names(de_markers_6))

#features_spatial <- markerGenes
# feature <- features_spatial[1]
#feature <- 'Fgfr3'

# p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
#   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
# p
#
# png(filename = paste0('./markers_plot/', feature, '_upReg_spatial.png'), width = 1200, height = 1200, res = 300)
# p
# dev.off()


# plot list
# features_spatial <- markerGenes
#features_spatial <- as.data.frame(markerList_pos$C7)$name[11:13]

plot_features <- function(feature){
  p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

# source('SpatialPlot_new_set_limit.R')
# plot_features <- function(feature){
#   p <- SpatialPlot_new_set_limit(spatial.obj, features = feature, pt.size.factor = 4, image.alpha = 0, stroke = 0, limits=c(0.32,0.55)) +
#     theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
#   p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
#   p
# }
library(grid)

ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

lapply(seq(length(ggList)),
       function(x)ggsave(filename=paste0('./markers_plot/clusterR7/', features_spatial[x],".png"), plot=ggList[[x]]))
1


##density plot
pd <- plot_density(spatial.obj, features, reduction = 'umap')
pd_data <- pd$data
pd_data <- pd_data[, c(3), drop=FALSE]

all(Cells(spatial.obj) == row.names(pd_data))
spatial.obj <- AddMetaData(object = spatial.obj, metadata = pd_data)

p <- SpatialFeaturePlot(spatial.obj, features = "feature",  pt.size.factor = 4.5, image.alpha = 0, stroke = 0) +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p
##
