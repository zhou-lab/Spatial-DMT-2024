#!/usr/bin/env Rscript

# Extract the small pieces needed to audit FigR smoothing from the deposited
# E11 Seurat object.  Run on a high-memory node because the GEO RDS is large.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: extract_published_e11_wnn.R INPUT_RDS_GZ OUTPUT_DIR")
}

suppressPackageStartupMessages(library(Seurat))

input <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", input)
# GEO applied an outer gzip layer to an RDS that was itself saved with R's
# default gzip compression, so two decompression layers are required.
outer <- gzfile(input, open = "rb")
inner <- gzcon(outer)
object <- readRDS(inner)
close(inner)
cells <- colnames(object)

summary_file <- file.path(outdir, "object_summary.txt")
sink(summary_file)
cat("class\t", paste(class(object), collapse = ","), "\n", sep = "")
cat("version\t", as.character(object@version), "\n", sep = "")
cat("cells\t", length(cells), "\n", sep = "")
cat("assays\t", paste(names(object@assays), collapse = ","), "\n", sep = "")
for (name in names(object@assays)) {
  assay <- object@assays[[name]]
  cat("assay\t", name, "\t", paste(dim(assay), collapse = "x"),
      "\tlayers=", paste(Layers(assay), collapse = ","), "\n", sep = "")
}
cat("reductions\t", paste(names(object@reductions), collapse = ","), "\n", sep = "")
for (name in names(object@reductions)) {
  cat("reduction\t", name, "\t",
      paste(dim(Embeddings(object, reduction = name)), collapse = "x"), "\n", sep = "")
}
cat("neighbors\t", paste(names(object@neighbors), collapse = ","), "\n", sep = "")
for (name in names(object@neighbors)) {
  nn <- object@neighbors[[name]]
  cat("neighbor\t", name, "\t", paste(dim(nn@nn.idx), collapse = "x"), "\n", sep = "")
}
cat("graphs\t", paste(names(object@graphs), collapse = ","), "\n", sep = "")
cat("images\t", paste(names(object@images), collapse = ","), "\n", sep = "")
cat("metadata\t", paste(colnames(object[[]]), collapse = ","), "\n", sep = "")
cat("cell_head\t", paste(head(cells), collapse = ","), "\n", sep = "")
print(object)
sink()

metadata <- object[[]]
metadata$cell <- rownames(metadata)
write.table(metadata, file.path(outdir, "cell_metadata.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

for (name in names(object@images)) {
  image <- object@images[[name]]
  coords <- image@coordinates
  coords$cell <- rownames(coords)
  write.table(coords, file.path(outdir, paste0("coordinates_", name, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

if ("weighted.nn" %in% names(object@neighbors)) {
  nn <- object@neighbors[["weighted.nn"]]
  idx <- nn@nn.idx
  dist <- nn@nn.dist
  rownames(idx) <- nn@cell.names
  rownames(dist) <- nn@cell.names
  colnames(idx) <- paste0("neighbor_", seq_len(ncol(idx)))
  colnames(dist) <- paste0("distance_", seq_len(ncol(dist)))
  write.table(data.frame(cell = rownames(idx), idx, check.names = FALSE),
              file.path(outdir, "weighted_nn_indices.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  write.table(data.frame(cell = rownames(dist), dist, check.names = FALSE),
              file.path(outdir, "weighted_nn_distances.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
}

# Features explicitly named in the repository script.  Exporting these from
# the raw and already-smoothed deposited assays provides a compact, exact
# visual check of the authors' saved FigR result.
feature_candidates <- c(
  "chr2.180331404.180333604",
  "chr8.85379165.85382665"
)
if (all(c("DNAm_frac", "DNAm_frac_smooth") %in% names(object@assays))) {
  common <- Reduce(intersect, list(
    feature_candidates,
    rownames(object[["DNAm_frac"]]),
    rownames(object[["DNAm_frac_smooth"]])
  ))
  if (length(common) > 0L) {
    raw <- LayerData(object, assay = "DNAm_frac", layer = "counts", features = common)
    smooth <- LayerData(object, assay = "DNAm_frac_smooth", layer = "counts", features = common)
    feature_values <- data.frame(cell = colnames(raw), check.names = FALSE)
    for (feature in common) {
      feature_values[[paste0(feature, "__raw")]] <- as.numeric(raw[feature, ])
      feature_values[[paste0(feature, "__smooth")]] <- as.numeric(smooth[feature, ])
    }
    write.table(feature_values, file.path(outdir, "script_feature_values.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

message("Wrote compact audit data to ", outdir)
