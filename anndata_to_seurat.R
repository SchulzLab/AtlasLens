#!/usr/bin/env Rscript
# anndata_to_seurat.R ---------------------------------------------------------
# Convert a scanpy AnnData (.h5ad) file into a Seurat .rds object for AtlasLens.
#
# This runs in YOUR OWN R, as a one-time pre-processing step. It is NOT part of
# the AtlasLens app — the app keeps loading a plain Seurat .rds as before, so no
# new dependencies are added to the app's Docker image / conda env.
#
# It uses anndataR's NATIVE R reader (rhdf5). No Python / reticulate is involved,
# which is why it works where the older Python-wrapping converters did not.
#
# It deliberately builds a CLASSIC (v3) Seurat assay rather than letting anndataR
# emit a v5 Assay5. A v3 assay stores plain @counts/@data matrices with no
# per-layer cell/feature maps, so the .rds loads under ANY Seurat version. (A v5
# object written here can fail to load in an app running a different Seurat
# version with: invalid class "Assay5" object: All layers must have a record in
# the cells map. Building v3 avoids that entirely.)
#
# One-time setup (in your R, not the app):
#   - R >= 4.5.0                         # anndataR requires this
#   - install.packages("BiocManager")
#   - BiocManager::install(c("anndataR", "rhdf5"))
#   - Seurat / SeuratObject              # you most likely already have these
#
# Usage:
#   Rscript anndata_to_seurat.R  input.h5ad  output.rds
#
# Expected .h5ad structure (standard scanpy):
#   layers['counts'] or X  -> expression. RAW COUNTS are preferred (AtlasLens's
#                             DEA needs counts). The script uses layers['counts']
#                             if present, otherwise X.
#   obs                    -> per-cell metadata (cell type, condition, timepoint)
#                             becomes Seurat meta.data.
#   var (index)            -> gene symbols / Ensembl IDs -> rownames.
#   obsm['X_umap'], ['X_pca'] -> become 'umap' / 'pca' reductions (optional;
#                             if absent, AtlasLens computes a UMAP on first load).
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript anndata_to_seurat.R input.h5ad output.rds", call. = FALSE)
in_h5ad <- args[[1]]
out_rds <- args[[2]]

if (!file.exists(in_h5ad))
  stop("Input file not found: ", in_h5ad, call. = FALSE)
if (getRversion() < "4.5.0")
  stop("anndataR requires R >= 4.5.0; this R is ", getRversion(),
       ". Upgrade R (or convert on a machine that has R >= 4.5).", call. = FALSE)

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Missing package '", pkg, "'. Install with: BiocManager::install('",
         pkg, "')", call. = FALSE)
}
invisible(lapply(c("anndataR", "rhdf5", "Seurat", "Matrix"), need))
suppressPackageStartupMessages({ library(anndataR); library(Seurat); library(Matrix) })

# Read the NATIVE AnnData object (no Seurat conversion yet) so we control how the
# Seurat object is assembled.
message("Reading ", in_h5ad, " ...")
adata <- read_h5ad(in_h5ad)

# --- expression matrix: prefer raw counts, fall back to X --------------------
layer_names <- tryCatch(names(adata$layers), error = function(e) character(0))
if ("counts" %in% layer_names) {
  message("Using layers[['counts']] as raw counts.")
  mat <- adata$layers[["counts"]]
} else {
  message("No 'counts' layer found - using X. ",
          "(If X is normalised, AtlasLens detects that and won't re-normalise.)")
  mat <- adata$X
}
if (is.null(mat))
  stop("No expression matrix found (both layers['counts'] and X are empty).",
       call. = FALSE)

mat <- Matrix::t(mat)                       # AnnData is cells x genes -> genes x cells
mat <- as(mat, "CsparseMatrix")
rownames(mat) <- make.unique(as.character(adata$var_names))
colnames(mat) <- as.character(adata$obs_names)

# --- cell metadata ----------------------------------------------------------
meta <- as.data.frame(adata$obs)
rownames(meta) <- colnames(mat)

# --- build a CLEAN, classic (v3) Seurat object ------------------------------
# A classic v3 "Assay" stores plain @counts/@data matrices with NO per-layer
# cell/feature maps, so the .rds loads under ANY Seurat version. A v5 "Assay5"
# written under one Seurat version can be rejected by another with:
#   invalid class "Assay5" ... All layers must have a record in the cells map.
options(Seurat.object.assay.version = "v3")
seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta)

# Guarantee v3 even if this Seurat ignores the option above. CreateAssayObject()
# ALWAYS builds a v3 assay regardless of the global option, so rebuild the RNA
# assay from its counts whenever CreateSeuratObject produced a v5 Assay5.
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  message("CreateSeuratObject produced a v5 assay; rebuilding it as classic v3 ...")
  seurat_obj[["RNA"]] <- CreateAssayObject(
    counts = GetAssayData(seurat_obj, assay = "RNA", layer = "counts"))
}
DefaultAssay(seurat_obj) <- "RNA"

# --- carry over embeddings (obsm: X_umap -> umap, X_pca -> pca, ...) ---------
obsm_names <- tryCatch(names(adata$obsm), error = function(e) character(0))
message("obsm keys in file: ",
        if (length(obsm_names)) paste(obsm_names, collapse = ", ") else "(none)")
for (key in obsm_names) {
  ok <- tryCatch({
    emb <- as.matrix(adata$obsm[[key]])
    if (nrow(emb) != ncol(seurat_obj))
      stop("rows (", nrow(emb), ") != cells (", ncol(seurat_obj), ")")
    if (ncol(emb) < 1) stop("no columns")
    rownames(emb) <- colnames(seurat_obj)
    red <- tolower(gsub("[^A-Za-z0-9]", "", sub("^X_", "", key)))   # X_umap -> umap
    if (!nzchar(red)) red <- "reduction"
    if (grepl("^[0-9]", red)) red <- paste0("r", red)               # key must start a-z
    colnames(emb) <- paste0(red, "_", seq_len(ncol(emb)))
    seurat_obj[[red]] <- CreateDimReducObject(embeddings = emb,
                                              key = paste0(red, "_"), assay = "RNA")
    red
  }, error = function(e) { message("  Skipped obsm['", key, "']: ",
                                   conditionMessage(e)); NA_character_ })
  if (!is.na(ok)) message("  Added '", ok, "' reduction from obsm['", key, "'].")
}

# NOTE: deliberately NO UpdateSeuratObject() here - in Seurat v5 it can re-upgrade
# the v3 assay back to a v5 Assay5, which is exactly what we are avoiding.

assay_class <- class(seurat_obj[["RNA"]])[1]
cat(sprintf("\nConverted: %d genes x %d cells\n", nrow(seurat_obj), ncol(seurat_obj)))
cat("RNA assay class: ", assay_class,
    if (identical(assay_class, "Assay")) "  (v3 - loads under any Seurat version)"
    else "  *** EXPECTED 'Assay' (v3) - tell the developer if this says Assay5 ***",
    "\n", sep = "")
cat("Metadata columns (used for celltype / condition / timepoint):\n  ",
    paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n", sep = "")
cat("Reductions: ",
    if (length(seurat_obj@reductions))
      paste(names(seurat_obj@reductions), collapse = ", ")
    else "(none found in obsm - AtlasLens will compute a UMAP on first load)",
    "\n", sep = "")

saveRDS(seurat_obj, out_rds)
cat("\nWrote ", normalizePath(out_rds),
    "\nPoint AtlasLens at it via the DATASET_PATH environment variable.\n",
    sep = "")
