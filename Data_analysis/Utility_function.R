### Give a list of VMR, map them to overlapping genes
Methylation_genes <- function(vmr_list){
  coordinates <- t(sapply(vmr_list, function(x) {
    parts <- unlist(strsplit(x, split = "[.]"))
    return(list(chr = parts[1], start = as.integer(parts[2]), end = as.integer(parts[3])))}))
  coordinates_df <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  coordinates_df$start <- as.integer(coordinates_df$start)
  coordinates_df$end <- as.integer(coordinates_df$end)
  coordRanges <- GRanges(
    seqnames = Rle(as.character(coordinates_df$chr)),
    ranges = IRanges(start = coordinates_df$start, end = coordinates_df$end)
  )
  genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
  overlap_GenesIndices <- findOverlaps(coordRanges, genes, select="first")
  coordinates_df$overlap_GenesIndices <- overlap_GenesIndices
  coordinates_df_sel <- na.omit(coordinates_df)
  proximalGenes_id <- genes$gene_id[coordinates_df_sel$overlap_GenesIndices]
  coordinates_df_sel$proximalGenes_id <- proximalGenes_id
  geneSymbols <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = proximalGenes_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  coordinates_df_sel$geneSymbols <- geneSymbols
  coordinates_df_sel[, c('geneSymbols'), drop=FALSE]
}

### Iterative PCA for VMR imputation
imputeColMean <- function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- colMeans(mtx, na.rm=TRUE)[k[,2]]
  mtx
}
prcomp_iterative <- function(x, n=10, n_iter=100, min_gain=0.01, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  #x <- imputeColMean(x) # average is our first guess for VMR
  x[na_loc] = 0  # 0 is our first guess for residual
  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = T, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  list(pr,x)
}

### Plot spatial_map
spatial_plot <- function(partition,spatial_barcode,continuous=FALSE,add_NA=FALSE){
  set.seed(123)
  group <- partition
  max <- max(as.numeric(group))
  spatial_coordinate <- read.table("~/SpM/spatial_barcodes.txt",sep="\t")
  idx <- match(spatial_barcode,spatial_coordinate$V4)
  spatial_coordinate_filter <- spatial_coordinate[idx,]
  spatial_heatmap <- cbind(spatial_coordinate_filter$V2,spatial_coordinate_filter$V3,group)
  rownames(spatial_heatmap) <- spatial_barcode
  empty_barcode_idx <- which(!spatial_coordinate$V4 %in% rownames(spatial_heatmap))
  empty_barcode <- cbind(spatial_coordinate[empty_barcode_idx,]$V2,spatial_coordinate[empty_barcode_idx,]$V3,rep(NA,length(empty_barcode_idx)))
  spatial_heatmap <- as.data.frame(rbind(spatial_heatmap,empty_barcode))
  spatial_cluster <- cbind(rownames(spatial_heatmap),spatial_heatmap$group)
  colnames(spatial_cluster) <- c("barcode","group")
  spatial_cluster <- as.data.frame(spatial_cluster)
  spatial_cluster$group <- as.numeric(spatial_cluster$group)
  spatial_heatmap <- as.data.frame(spatial_heatmap)
  spatial_heatmap$V1 <- as.factor(spatial_heatmap$V1)
  spatial_heatmap$V2 <- as.factor(spatial_heatmap$V2)
  if(continuous == FALSE){
    spatial_heatmap$group <- as.factor(spatial_heatmap$group)
    spatial_heatmap$group <- ifelse(is.na(spatial_heatmap$group),spatial_heatmap$group,paste0("C",spatial_heatmap$group))
    spatial_heatmap$group <- factor(spatial_heatmap$group, levels = paste0("C", c(1:max)))
    spatial <- ggplot(spatial_heatmap, aes(V2, V1))+
      geom_point(aes(color=group),size = 1.5) +
      theme_bw() + xlab("") + ylab("")+scale_x_discrete(limits = rev(levels(spatial_heatmap$V2)))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
      theme(legend.title=element_blank())+
      coord_fixed(ratio=1)+ scale_y_discrete(limits=rev)+scale_color_discrete(na.value="white",na.translate = FALSE)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
  }
  else{
    spatial_heatmap$group <- as.numeric(spatial_heatmap$group)
    if(add_NA == TRUE){
      spatial <- ggplot(spatial_heatmap, aes(V2, V1,color=group))+
        geom_point(size = 1.5) +
        theme_bw() + xlab("") + ylab("")+scale_x_discrete(limits = rev(levels(spatial_heatmap$V2)))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
        theme(legend.title=element_blank())+
        coord_fixed(ratio=1)+ scale_y_discrete(limits=rev)+scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral")),na.value = "gray") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
              axis.ticks.y = element_blank())
    }
    else{
      spatial <- ggplot(spatial_heatmap, aes(V2, V1,fill=group))+
        geom_tile(color = "black", size = 0.1,width = 0.8, height = 0.8) +
        theme_bw() + xlab("") + ylab("")+scale_x_discrete(limits = rev(levels(spatial_heatmap$V2)))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())+
        coord_fixed(ratio=1)+ scale_y_discrete(limits=rev)+scale_color_gradientn(colors = rev(brewer.pal(n = 11, name = "Spectral")),na.value = "white") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
              axis.ticks.y = element_blank())+geom_tile(data = subset(spatial_heatmap, is.na(group)), color = "white", size = 0.3)
    }
  }
  return(spatial)
}
