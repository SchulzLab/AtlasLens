#!/usr/bin/env Rscript
# =============================================================================
# Build a Seurat object for AtlasLens from raw input files
# (matrix + features + barcodes, OR a dense count table, + optional metadata)
#
#   AtlasLens currently accepts a Seurat object (.rds). This script turns the
#   common upstream file formats into that object. Pick ONE input mode below,
#   set the paths, run:  Rscript build_seurat_from_files.R
# =============================================================================
#
# -----------------------------------------------------------------------------
# REQUIRED FILE FORMATS  — your files MUST match one of these shapes
# -----------------------------------------------------------------------------
#
# MODE "mtx"  — the 10x / Matrix-Market triplet
#   * counts matrix  : .mtx or .mtx.gz, sparse Matrix-Market format.
#                      Orientation: GENES in rows, CELLS in columns.
#                      Values: raw integer UMI/read counts (NOT normalized).
#   * features file  : plain text (.tsv/.txt, optionally .gz), ONE GENE PER LINE.
#                      Number of lines == number of matrix ROWS.
#                      If it is a 10x features.tsv with several columns
#                      (EnsemblID <tab> Symbol <tab> Type), set FEATURE_COLUMN
#                      to the column holding the gene SYMBOL (usually 2).
#   * barcodes file  : plain text, ONE CELL BARCODE PER LINE.
#                      Number of lines == number of matrix COLUMNS.
#
# MODE "dense" — a single rectangular table (.csv/.tsv, optionally .gz)
#   * First column  = gene symbols   (these become the row names).
#   * Header row    = cell barcodes  (these become the column names).
#   * Body          = raw integer counts. Orientation: GENES rows, CELLS columns.
#       gene,    CellA, CellB, CellC
#       Actb,    5,     0,     2
#       Gapdh,   3,     7,     1
#
# METADATA file (optional but recommended for AtlasLens)  — .csv/.tsv(.gz)
#   * ONE ROW PER CELL.
#   * One column must hold the cell barcode that EXACTLY matches the matrix
#     column names (set BARCODE_COLUMN). If omitted, rows are matched by order
#     and must be in the SAME order/length as the barcodes file.
#   * Include the annotation columns you want to explore in AtlasLens, e.g.:
#       cell_type, sample/dataset, condition, sex, tissue, and — for the
#       temporal feature — an ORDERED time column (e.g. day0,day3,...).
#
# GENERAL RULES
#   * Cell barcodes must be UNIQUE.
#   * Gene names cannot contain "_" — Seurat will replace them with "-".
#   * Counts must be the RAW matrix; normalization is done here, not before.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
set.seed(42)

# =============================  CONFIG  ======================================
# 1) Choose ONE input mode:
INPUT_MODE <- "mtx"          # "mtx"  or  "dense"

# 2a) Paths for MODE "mtx"
MTX_FILE       <- "path/to/matrix.mtx.gz"
FEATURES_FILE  <- "path/to/features.tsv.gz"
BARCODES_FILE  <- "path/to/barcodes.tsv.gz"
FEATURE_COLUMN <- 1          # column in features file holding the gene SYMBOL

# 2b) Path for MODE "dense"
DENSE_FILE     <- "path/to/counts.csv.gz"
DENSE_SEP      <- ","        # "," for CSV, "\t" for TSV

# 3) Metadata (set to NULL to skip)
META_FILE      <- "path/to/metadata.csv"
META_SEP       <- ","
BARCODE_COLUMN <- "barcode"  # column in metadata matching matrix barcodes;
                             # set to NULL to match by row order instead.

# 4) Output + processing
PROJECT_NAME   <- "AtlasLens_dataset"
OUT_RDS        <- "AtlasLens_input.rds"
RUN_PROCESSING <- TRUE       # NormalizeData -> PCA -> UMAP -> clusters
N_PCS          <- 30
N_HVG          <- 2000
# =============================================================================


# ----------------------------- helpers ---------------------------------------
.read_lines_col <- function(path, col = 1) {
  df <- read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
                   header = FALSE, stringsAsFactors = FALSE)
  df[[col]]
}

load_mtx <- function(mtx, features, barcodes, feature.column = 1) {
  m     <- Matrix::readMM(if (grepl("\\.gz$", mtx)) gzfile(mtx) else mtx)
  genes <- .read_lines_col(features, feature.column)
  bcs   <- .read_lines_col(barcodes, 1)
  if (nrow(m) == length(genes) && ncol(m) == length(bcs)) {
    # genes x cells (expected)
  } else if (nrow(m) == length(bcs) && ncol(m) == length(genes)) {
    m <- Matrix::t(m)                        # transpose cells x genes -> genes x cells
  } else {
    stop(sprintf("Dimension mismatch: mtx is %dx%d but features=%d, barcodes=%d.\n%s",
                 nrow(m), ncol(m), length(genes), length(bcs),
                 "Check that the matrix is genes-in-rows and the files line up."))
  }
  rownames(m) <- make.unique(as.character(genes))
  colnames(m) <- as.character(bcs)
  as(m, "CsparseMatrix")
}

load_dense <- function(path, sep = ",") {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  df  <- read.delim(con, sep = sep, header = TRUE, row.names = 1,
                    check.names = FALSE, stringsAsFactors = FALSE)
  m <- as(as.matrix(df), "CsparseMatrix")    # genes (rows) x cells (cols)
  rownames(m) <- make.unique(rownames(m))
  m
}

# ----------------------------- 1. counts -------------------------------------
counts <- switch(INPUT_MODE,
  mtx   = load_mtx(MTX_FILE, FEATURES_FILE, BARCODES_FILE, FEATURE_COLUMN),
  dense = load_dense(DENSE_FILE, DENSE_SEP),
  stop("INPUT_MODE must be 'mtx' or 'dense'."))

if (any(duplicated(colnames(counts))))
  stop("Cell barcodes are not unique — fix the barcodes/header before continuing.")
message("Counts: ", nrow(counts), " genes x ", ncol(counts), " cells")

# ----------------------------- 2. metadata -----------------------------------
meta <- NULL
if (!is.null(META_FILE)) {
  con  <- if (grepl("\\.gz$", META_FILE)) gzfile(META_FILE) else META_FILE
  meta <- read.delim(con, sep = META_SEP, header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE)
  message("Metadata columns: ", paste(colnames(meta), collapse = ", "))

  if (!is.null(BARCODE_COLUMN)) {
    if (!BARCODE_COLUMN %in% colnames(meta))
      stop("BARCODE_COLUMN '", BARCODE_COLUMN, "' not found in metadata.")
    rownames(meta) <- as.character(meta[[BARCODE_COLUMN]])
    common <- intersect(colnames(counts), rownames(meta))
    if (length(common) == 0)
      stop("No barcodes in metadata match the matrix columns. ",
           "Check BARCODE_COLUMN and that the barcodes are identical strings.")
    if (length(common) < ncol(counts))
      warning(ncol(counts) - length(common),
              " cells have no metadata row and will be dropped.")
    counts <- counts[, common, drop = FALSE]
    meta   <- meta[common, , drop = FALSE]
  } else {
    if (nrow(meta) != ncol(counts))
      stop("Matching by row order, but metadata rows (", nrow(meta),
           ") != number of cells (", ncol(counts), ").")
    rownames(meta) <- colnames(counts)
  }
}

# ----------------------------- 3. build --------------------------------------
obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = PROJECT_NAME)
message("Seurat object: ", ncol(obj), " cells, ", nrow(obj), " genes")

# ----------------------------- 4. processing ---------------------------------
# AtlasLens expects a processed (and, for multi-sample atlases, integrated)
# object. This is a standard single-batch pass; if your data spans batches you
# should integrate (Harmony / Seurat anchors) instead of / after this block.
if (RUN_PROCESSING) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = N_PCS, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:N_PCS, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:N_PCS, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
}

# ----------------------------- 5. save ---------------------------------------
saveRDS(obj, file = OUT_RDS)
message("Saved: ", normalizePath(OUT_RDS))

# ----------------------------- 6. verify -------------------------------------
cat("\n================ VERIFICATION ================\n")
cat("Cells x genes :", ncol(obj), "x", nrow(obj), "\n")
cat("Assays        :", paste(Assays(obj), collapse = ", "), "\n")
cat("Reductions    :", paste(Reductions(obj), collapse = ", "), "\n")
cat("Metadata cols :\n  ", paste(colnames(obj@meta.data), collapse = ", "), "\n")
cat("==============================================\n")
