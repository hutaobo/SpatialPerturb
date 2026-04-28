args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_seurat_to_tables.R <input_rds_or_rds_gz> <output_dir>")
}

input_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

read_seurat_object <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(readRDS(con))
  }
  readRDS(path)
}

get_counts_matrix <- function(object) {
  assay_name <- DefaultAssay(object)
  suppressWarnings({
    counts <- tryCatch(
      GetAssayData(object = object, assay = assay_name, slot = "counts"),
      error = function(...) NULL
    )
  })
  if (!is.null(counts)) {
    return(counts)
  }
  counts <- tryCatch(
    LayerData(object = object, assay = assay_name, layer = "counts"),
    error = function(...) NULL
  )
  if (is.null(counts)) {
    stop("Could not extract a counts matrix from the Seurat object.")
  }
  counts
}

object <- read_seurat_object(input_path)
counts <- get_counts_matrix(object)
meta <- object[[]]

feature_names <- rownames(counts)
if (is.null(feature_names)) {
  stop("Counts matrix is missing feature names.")
}
cell_names <- colnames(counts)
if (is.null(cell_names)) {
  stop("Counts matrix is missing cell barcodes.")
}

meta <- meta[cell_names, , drop = FALSE]
meta$cell_id <- rownames(meta)

var <- data.frame(
  feature_id = feature_names,
  gene = feature_names,
  row.names = feature_names,
  stringsAsFactors = FALSE
)

Matrix::writeMM(obj = counts, file = file.path(output_dir, "matrix.mtx"))
utils::write.csv(meta, file = file.path(output_dir, "obs.csv"), row.names = TRUE)
utils::write.csv(var, file = file.path(output_dir, "var.csv"), row.names = TRUE)
