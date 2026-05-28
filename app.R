# =============================================================================
# AtlasLens
# =============================================================================

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(msigdbr)
library(plotly)
suppressWarnings(library(geneCOCOA)) # Suppress annotate replacement warning
library(shinycssloaders)
library(viridis)
library(ggrepel)
library(digest)
library(future)
library(promises)
library(DT)
library(shinyjs)
library(qs)
library(RColorBrewer)
# patchwork comes in as a Seurat dependency; we use the `/` operator to
# stack the three QC violin plots for the combined download.
suppressWarnings(library(patchwork))

# Defensive null-coalescing operator. Shiny exports its own %||% in recent
# versions, but defining it here makes the file portable across Shiny
# versions and avoids "could not find function" errors at startup.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}


# === CONFIGURATION (DYNAMIC) ===

# 1. Determine Dataset Path
# Deployers point AtlasLens at their dataset via the DATASET_PATH environment
# variable, e.g.  -e DATASET_PATH="/mnt/data/file.rds"  in `docker run`, or
# `export DATASET_PATH=...` before launching R. When the variable is unset
# DEFAULT_DATA_PATH stays NULL and the Introduction tab surfaces a clear
# user-facing notice instead of trying to load from a hard-coded location.
env_data_path <- Sys.getenv("DATASET_PATH")

if (nzchar(env_data_path)) {
  DEFAULT_DATA_PATH <- env_data_path
  message(paste("Config: Using dataset from environment variable:", DEFAULT_DATA_PATH))
} else {
  DEFAULT_DATA_PATH <- NULL
  message("Config: DATASET_PATH is not set. Set the DATASET_PATH environment ",
          "variable to the absolute path of your .rds Seurat object before ",
          "launching AtlasLens.")
}

# 2. Determine Cache Directory (Ephemeral)
# ~/tmp inside the container
CACHE_DIR <- path.expand("~/tmp")

if (!dir.exists(CACHE_DIR)) {
  tryCatch({
    dir.create(CACHE_DIR, recursive = TRUE)
    message(paste("Config: Created ephemeral cache directory at", CACHE_DIR))
  }, error = function(e) {
    message("Config Warning: Could not create ~/tmp. Using local 'cache' directory.")
    CACHE_DIR <<- "cache"
    dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  })
} else {
  message(paste("Config: Ephemeral cache directory exists at", CACHE_DIR))
}

SAMPLESIZE <- 10
PSEUDOCOUNT <- 0.1
N_SIMULATIONS <- 100
MAX_UPLOAD_SIZE <- 20000 
COCOA_MAX_CELLS <- 500 # Limit cells per group for COCOA speed (Downsampling threshold)

options(shiny.maxRequestSize = MAX_UPLOAD_SIZE * 1024^2)
options(future.globals.maxSize = 800000 * 1024^2) 
plan(multicore, workers = 4)

# === LOGGING & LOCKING ===

get_lock_file <- function(filename) {
  file.path(CACHE_DIR, paste0("lock_", filename, ".txt"))
}

is_locked <- function(filename) {
  file.exists(get_lock_file(filename))
}

set_lock <- function(filename) {
  file.create(get_lock_file(filename))
}

remove_lock <- function(filename) {
  f <- get_lock_file(filename)
  if (file.exists(f)) unlink(f)
}

# === HELPER FUNCTIONS ===

info_icon <- function(text) {
  tags$span(class = "custom-tooltip",
            icon("info-circle"),
            tags$span(class = "tooltiptext", text)
  )
}



get_valid_metadata_columns <- function(seurat_obj, show_hidden = FALSE) {
  # Returns the metadata columns offered to the user for grouping / filtering /
  # colouring. By default we hide three classes of non-useful columns; pass
  # show_hidden = TRUE to bypass the filter (the Dataset Overview tab exposes
  # a toggle so reviewers can inspect everything).
  #
  # Hidden by default:
  #   * Technical blacklist: per-cell barcode and raw QC counters that the
  #     app never uses as grouping factors (n_genes, n_counts,
  #     observation_joinid, nCount_RNA, nFeature_RNA, cell).
  #   * Single-valued columns: any column that takes one value across the
  #     whole object - no analysis can compare across a single group.
  #   * Non-atomic columns (e.g. sf 'geometry' list-columns) - they would
  #     break the UI controls outright.
  #
  # Non-atomic columns are hidden in both modes because the UI cannot
  # render them.
  meta <- seurat_obj@meta.data
  if (ncol(meta) == 0) return(character())
  cols <- names(meta)
  is_atomic_col <- vapply(cols, function(col) is.atomic(meta[[col]]), logical(1))
  
  if (show_hidden) return(cols[is_atomic_col])
  
  blacklist <- c("n_genes", "n_counts", "observation_joinid",
                 "nCount_RNA", "nFeature_RNA", "cell")
  is_blacklisted <- cols %in% blacklist
  is_single_valued <- vapply(cols, function(col) {
    if (!is.atomic(meta[[col]])) return(FALSE)
    v <- meta[[col]]; v <- v[!is.na(v)]
    length(unique(v)) <= 1L
  }, logical(1))
  cols[is_atomic_col & !is_blacklisted & !is_single_valued]
}

# Companion: returns the metadata columns hidden by get_valid_metadata_columns
# grouped by the reason. Used by the Dataset Overview UI to be transparent
# about what is filtered out and let the user toggle the filter off.
get_hidden_metadata_columns <- function(seurat_obj) {
  meta <- seurat_obj@meta.data
  empty <- list(technical = character(), single_valued = character(),
                non_atomic = character())
  if (ncol(meta) == 0) return(empty)
  cols <- names(meta)
  is_atomic_col <- vapply(cols, function(col) is.atomic(meta[[col]]), logical(1))
  blacklist <- c("n_genes", "n_counts", "observation_joinid",
                 "nCount_RNA", "nFeature_RNA", "cell")
  technical    <- cols[is_atomic_col & cols %in% blacklist]
  single_valued <- cols[is_atomic_col & !(cols %in% blacklist) &
                          vapply(cols, function(col) {
                            if (!is.atomic(meta[[col]])) return(FALSE)
                            v <- meta[[col]]; v <- v[!is.na(v)]
                            length(unique(v)) <= 1L
                          }, logical(1))]
  non_atomic <- cols[!is_atomic_col]
  list(technical = technical, single_valued = single_valued,
       non_atomic = non_atomic)
}

# Helper to generate enough distinct colors
# Heuristically detect the metadata columns that play standard biological
# roles. Each role is matched by an ordered list of regular expressions and
# the first matching column name wins. Returns NULL for any role the dataset
# does not expose, so callers must check before use.
detect_role_columns <- function(meta) {
  cols <- colnames(meta)
  if (length(cols) == 0) return(list())

  # A timepoint column must not only be NAMED like a timepoint, its VALUES
  # must actually look like ordered timepoints. This rejects columns whose
  # values are ontology identifiers (e.g. "PARO:0000461", "CL:0000236") or
  # arbitrary free-text labels, which would otherwise yield a meaningless
  # temporal analysis. Accepts pure numbers, an embedded number
  # (Day_3, D7, 24h, 7dpi), or a recognised ordinal vocabulary
  # (baseline/acute/chronic, pre/during/post, ...).
  is_timepoint_like <- function(values) {
    v <- as.character(values)
    v <- v[!is.na(v) & nzchar(v)]
    v <- v[!tolower(trimws(v)) %in%
             c("no_data", "no data", "none", "n/a", "na", "unknown", "-")]
    v <- unique(v)
    if (length(v) == 0) return(FALSE)
    # Reject ontology-style identifiers outright (PREFIX:digits, e.g. PARO:0000461).
    if (any(grepl("^[A-Za-z]+:[0-9]+$", v))) return(FALSE)
    lt <- tolower(trimws(v))
    # Stage words that legitimately appear on a time axis without a number
    # (so a "Control"/"Baseline" reference level alongside Day1/Day3 still counts).
    ordinal <- c("baseline", "acute", "subacute", "chronic", "control", "ctrl",
                 "naive", "sham", "mock", "untreated", "pre", "during", "post",
                 "early", "intermediate", "mid", "middle", "late", "terminal")
    looks_timed <- function(x) {
      if (!is.na(suppressWarnings(as.numeric(x)))) return(TRUE)              # 7, 0.5
      if (!is.na(suppressWarnings(as.numeric(                                # Day_3, 24h, 7dpi
        sub(".*?(-?[0-9]*\\.?[0-9]+).*", "\\1", x))))) return(TRUE)
      x %in% ordinal                                                        # baseline/acute/...
    }
    # Require a strong majority so categorical columns (cell types, sex, ...)
    # are still rejected, while tolerating the odd stray label in an otherwise
    # temporal column.
    mean(vapply(lt, looks_timed, logical(1))) >= 0.8
  }

  pick <- function(patterns, validate = NULL) {
    for (p in patterns) {
      for (h in grep(p, cols, ignore.case = TRUE, value = TRUE)) {
        if (is.null(validate) || isTRUE(validate(meta[[h]]))) return(h)
      }
    }
    NULL
  }
  list(
    timepoint = pick(c("^timeseriesinfo$", "^time_series", "^timepoint$",
                       "^time$", "^day$", "^dpi$", "^hpi$",
                       "_dpi$", "_hpi$", "_hour", "timepoint"),
                     validate = is_timepoint_like),
    condition = pick(c("^condition$", "^treatment$", "^group$", "^status$",
                       "^genotype$", "^phenotype$", "^disease", "condition")),
    celltype  = pick(c("^celltype_annotated$", "^annotated_celltype$",
                       "^celltype$", "^cell_type$", "^cell\\.type$",
                       "^cluster_label$", "celltype")),
    dataset   = pick(c("^datasetid$", "^dataset_id$", "^dataset$",
                       "^study$", "^batch$", "^sample$", "^donor$",
                       "dataset", "batch"))
  )
}

# Mask cells that have a usable timepoint annotation. Cells with NA, an
# empty string, or any of the common placeholder labels (no_data, none,
# n/a, unknown) are dropped. Used by every Time Series observer so the
# tab works regardless of the placeholder convention a dataset uses.
ts_has_timepoint_mask <- function(values) {
  if (is.null(values)) return(logical(0))
  v <- as.character(values)
  bad <- is.na(v) | v == "" |
    tolower(v) %in% c("no_data", "no data", "none", "n/a", "na", "unknown", "-")
  !bad
}

# Order timepoint labels chronologically. Handles pure numbers, labels with an
# embedded number (Day_3, D7, 24h, week2), and known ordinal vocabularies
# (baseline/acute/chronic, pre/during/post, ...). Falls back to input order so
# string-labelled timepoints never collapse into alphabetical order.
sort_timepoints <- function(tps) {
  tps <- unique(as.character(tps))
  if (length(tps) <= 1) return(tps)
  num <- suppressWarnings(as.numeric(tps))
  if (!anyNA(num)) return(tps[order(num)])
  emb <- suppressWarnings(as.numeric(sub(".*?(-?[0-9]*\\.?[0-9]+).*", "\\1", tps)))
  if (!anyNA(emb)) return(tps[order(emb)])
  ladders <- list(
    c("baseline", "acute", "subacute", "chronic"),
    c("control", "ctrl", "naive", "early", "intermediate", "mid", "middle", "late", "terminal"),
    c("pre", "during", "post")
  )
  lt <- tolower(trimws(tps))
  for (lad in ladders) {
    idx <- match(lt, lad)
    if (!anyNA(idx)) return(tps[order(idx)])
  }
  tps
}

# Infer the species of a Seurat object from gene-symbol casing. Human
# symbols are almost always all-caps (TNF, IL6, GAPDH); mouse symbols are
# title-case (Tnf, Il6, Gapdh). Returns one of "Homo sapiens" or
# "Mus musculus" — the two species expected by msigdbr::msigdbr().
detect_species <- function(seurat_obj) {
  genes <- rownames(seurat_obj)
  if (length(genes) == 0) return("Mus musculus")
  sample_genes <- if (length(genes) > 2000) sample(genes, 2000) else genes
  # Ensembl gene IDs carry an unambiguous species prefix (ENSMUSG... = mouse,
  # ENSG... = human). Check that first — the all-caps casing heuristic below
  # would otherwise misread every Ensembl ID as a human symbol.
  if (mean(grepl("^ENSMUS", sample_genes)) > 0.5) return("Mus musculus")
  if (mean(grepl("^ENSG[0-9]", sample_genes)) > 0.5) return("Homo sapiens")
  alphabetic <- sample_genes[grepl("^[A-Za-z]", sample_genes)]
  if (length(alphabetic) == 0) return("Mus musculus")
  upper_frac <- mean(alphabetic == toupper(alphabetic))
  if (upper_frac > 0.7) "Homo sapiens" else "Mus musculus"
}

# Detect plausible "control" condition values within a vector of condition
# labels. Recognises the common scRNA-seq controls (Control, control, Ctrl,
# Sham, WT, wild-type, baseline, untreated) plus a numeric "0" used in some
# bench datasets as a t = 0 baseline.
detect_control_conditions <- function(conds) {
  if (length(conds) == 0) return(character(0))
  pat <- "^(control|ctrl|sham|wt|wild[-_ ]?type|untreated|baseline|0)(/.*)?$"
  conds[grepl(pat, conds, ignore.case = TRUE)]
}

get_expanded_palette <- function(n) {
  # A manually curated list of 55 highly distinct colors to avoid "muddy" interpolation
  distinct_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", 
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", 
    "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", 
    "#7F7F7F", "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5", "#1f77b4"
  )
  
  if (n <= length(distinct_colors)) {
    return(distinct_colors[1:n])
  } else {
    # Only fallback to interpolation if we have a massive number of groups (>55)
    return(colorRampPalette(distinct_colors)(n))
  }
}

# === WORKER: LOAD & PROCESS ===
load_and_process_seurat_worker <- function(file_path, cache_dir) {
  filename <- basename(file_path)
  tryCatch({
    cached_path <- file.path(cache_dir, paste0("processed_", filename, ".qs"))
    
    if (file.exists(cached_path)) {
      tryCatch({
        seurat_obj <- qs::qread(cached_path) 
        
        # Check for either umap or scVI
        if ("umap" %in% names(seurat_obj@reductions) || "scVI" %in% names(seurat_obj@reductions)) {
          dataset_hash <- digest(cached_path, algo="md5")
          return(list(data = seurat_obj, error = NULL, hash = dataset_hash, msg = "Loaded cached processed data."))
        }
      }, error = function(e) unlink(cached_path))
    }
    
    if (!file.exists(file_path)) return(list(data = NULL, error = paste("File not found:", file_path), hash = NULL))
    
    if (grepl("\\.qs$", file_path)) { seurat_obj <- qs::qread(file_path) } else { seurat_obj <- readRDS(file_path) }
    gc() 
    
    if (ncol(seurat_obj@meta.data) > 0) {
      seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) {
        if (is.factor(x)) levels(x) <- iconv(levels(x), "UTF-8", "UTF-8", sub = "_")
        if (is.character(x)) x <- iconv(x, "UTF-8", "UTF-8", sub = "_")
        return(x)
      })
    }
    
    
    if (!"RNA" %in% names(seurat_obj@assays)) {
      default_assay <- DefaultAssay(seurat_obj)
      seurat_obj[["RNA"]] <- seurat_obj[[default_assay]]
    }
    DefaultAssay(seurat_obj) <- "RNA"
    if (inherits(seurat_obj[["RNA"]], "Assay5")) suppressWarnings({ seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]]) })
    
    was_processed <- FALSE
    # If neither umap nor scVI exists, run standard pipeline
    if (!("umap" %in% names(seurat_obj@reductions)) && !("scVI" %in% names(seurat_obj@reductions))) {
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE); gc()
      seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE); gc()
      seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE); gc()
      seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE); gc()
      suppressWarnings({ seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE, n.neighbors = 15, min.dist = 0.5) })
      gc()
      was_processed <- TRUE
    }
    
    msg <- "Loaded successfully."
    if (was_processed || !file.exists(cached_path)) {
      qs::qsave(seurat_obj, cached_path)
      msg <- "Data processed and cached! Next load will be instant."
    }
    
    # Make sure the QC metadata columns (nFeature_RNA / nCount_RNA / percent.mt)
    # are available so the QC Violins subtab can render even on objects that
    # were saved without them.
    seurat_obj <- ensure_qc_columns(seurat_obj)
    
    dataset_hash <- digest(paste(ncol(seurat_obj), nrow(seurat_obj), colnames(seurat_obj@meta.data), collapse = "_"), algo = "md5")
    return(list(data = seurat_obj, error = NULL, hash = dataset_hash, msg = msg))
  }, error = function(e) {
    return(list(data = NULL, error = conditionMessage(e), hash = NULL))
  })
}

# --- WORKER FUNCTIONS ---

# =============================================================================
# COCOA worker 
# =============================================================================

run_cocoa_worker <- function(counts_matrix, gene, gene_sets, samplesize, n_sims,
                             log_transform = FALSE) {
  # Calculate true library sizes BEFORE subsetting genes to prevent 0-sum NaN errors
  lib_sizes <- colSums(counts_matrix)
  lib_sizes[lib_sizes == 0] <- 1 # Fallback for completely empty cells
  
  # Limit to top 500 expressed genes (plus target gene) for speed
  if (nrow(counts_matrix) > 500) {
    gene_sums <- rowSums(counts_matrix)
    top_genes <- names(sort(gene_sums, decreasing = TRUE))[1:500]
    # Ensure target gene is included
    if (gene %in% rownames(counts_matrix) && !(gene %in% top_genes)) {
      top_genes[500] <- gene
    }
    counts_matrix <- counts_matrix[top_genes, , drop = FALSE]
  }
  
  cpm <- t(t(as.matrix(counts_matrix)) / lib_sizes * 1e6) + 0.1
  if (isTRUE(log_transform)) cpm <- log2(cpm)
  expr_df <- as.data.frame(cpm)
  expr_info <- geneCOCOA::get_expr_info(expr = expr_df, GOI = gene)
  if (is.null(expr_info) || nrow(expr_info$expr_df) == 0) return(NULL)
  available <- rownames(expr_info$expr_df)
  
  # === OPTIMIZATION: Max 500 genes per group ===
  if (length(available) > 500) {
    means <- rowMeans(expr_info$expr_df[available, , drop = FALSE])
    available <- names(sort(means, decreasing = TRUE))[1:500]
    expr_info$expr_df <- expr_info$expr_df[available, , drop = FALSE]
  }
  
  runtime_sets <- lapply(gene_sets, function(s) intersect(s, available))
  runtime_sets <- runtime_sets[sapply(runtime_sets, length) > samplesize]
  if (length(runtime_sets) == 0) return(NULL)
  suppressWarnings({ result <- geneCOCOA::get_stats(geneset_collection = runtime_sets, GOI = gene, GOI_expr = expr_info$GOI_expr, expr_df = expr_info$expr_df, nsims = n_sims, samplesize = samplesize) })
  gc()
  if (!is.null(result) && !is.null(result$p_value_df)) return(result$p_value_df) else return(NULL)
}

# =============================================================================
# DEA worker 
# =============================================================================
run_dea_worker <- function(subset_counts, subset_meta, meta_col, group1, group2) {
  # === Wilcoxon + Presto using Minimal Object ===
  mini_seurat <- Seurat::CreateSeuratObject(counts = subset_counts, meta.data = subset_meta)
  suppressWarnings({ mini_seurat <- Seurat::NormalizeData(mini_seurat, verbose = FALSE) })
  gc()
  
  markers <- NULL
  suppressMessages({
    suppressWarnings({
      Seurat::Idents(mini_seurat) <- mini_seurat@meta.data[[meta_col]]
      markers <- Seurat::FindMarkers(mini_seurat, ident.1 = group1, ident.2 = group2,
                                     test.use = "wilcox", use_presto = TRUE,
                                     logfc.threshold = 0.585, min.pct = 0.1, verbose = FALSE,
                                     assay = "RNA", slot = "data", recorrect_umi = FALSE)
    })
  })
  gc()
  
  if (!is.null(markers) && nrow(markers) > 0) {
    markers$gene <- rownames(markers)
    markers$group1 <- group1
    markers$group2 <- group2
  }
  return(markers)
}

# === GO ENRICHMENT WORKER ===
# Runs clusterProfiler::enrichGO inside a future worker. Returns a list with
# `mode`, the enrichment data frames (all_df / up_df / down_df), and the
# rrvgo semantic-similarity reductions (all_rrvgo / up_rrvgo / down_rrvgo)
# that drive the treemap and scatter plots. Plain data frames + matrices
# serialise cleanly back from the future worker.
#
# org_db / ontology are passed through to rrvgo because the reduction step
# needs them to look up GO term annotations.
run_go_worker <- function(gene_list, species, ontology, p_cutoff, q_cutoff,
                          mode, key_type) {
  org_db <- if (species == "Homo sapiens") "org.Hs.eg.db" else "org.Mm.eg.db"
  
  # Run clusterProfiler::enrichGO on one gene set and return a plain data
  # frame (plain data frames serialise cleanly back from the future worker).
  enrich_df <- function(genes) {
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (length(genes) < 1) return(NULL)
    ego <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org_db,
      keyType       = key_type,
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    )
    if (is.null(ego)) return(NULL)
    df <- as.data.frame(ego)
    if (nrow(df) == 0) return(NULL)
    df
  }
  
  # Build the rrvgo semantic-similarity reduction for an enrichGO data frame.
  # Returns NULL if rrvgo is unavailable, too few terms, or the similarity
  # matrix cannot be built (e.g. all GO IDs are obsolete in the OrgDb release).
  # Output: list(sim_matrix, reduced) where `reduced` is the data frame
  # returned by rrvgo::reduceSimMatrix() - that's what the treemap + scatter
  # both consume.
  rrvgo_reduce <- function(df) {
    if (is.null(df) || nrow(df) < 3) return(NULL)
    if (!requireNamespace("rrvgo", quietly = TRUE)) return(NULL)
    # rrvgo's similarity matrix is only defined for BP / MF / CC ontologies.
    if (!ontology %in% c("BP", "MF", "CC")) return(NULL)
    tryCatch({
      go_ids <- df$ID
      sim_mat <- rrvgo::calculateSimMatrix(go_ids, orgdb = org_db,
                                           ont = ontology, method = "Rel")
      if (is.null(sim_mat) || nrow(sim_mat) < 2) return(NULL)
      scores <- stats::setNames(-log10(df$p.adjust), df$ID)
      scores <- scores[rownames(sim_mat)]
      reduced <- rrvgo::reduceSimMatrix(sim_mat, scores,
                                        threshold = 0.7, orgdb = org_db)
      list(sim_matrix = sim_mat, reduced = reduced)
    }, error = function(e) NULL)
  }
  
  result <- list(mode = mode,
                 all_df = NULL, up_df = NULL, down_df = NULL,
                 all_rrvgo = NULL, up_rrvgo = NULL, down_rrvgo = NULL)
  if (mode == "compare") {
    # Up- and down-regulated genes enriched separately, shown side by side.
    result$up_df     <- enrich_df(gene_list$Upregulated)
    result$down_df   <- enrich_df(gene_list$Downregulated)
    result$up_rrvgo   <- rrvgo_reduce(result$up_df)
    result$down_rrvgo <- rrvgo_reduce(result$down_df)
  } else {
    # All DE genes pooled into a single enrichment.
    result$all_df    <- enrich_df(gene_list[[1]])
    result$all_rrvgo <- rrvgo_reduce(result$all_df)
  }
  return(result)
}

# GO dot plot: GeneRatio on the x-axis, GO terms on the y-axis, one facet per
# gene set. `df_list` is a named list of enrichGO data frames.
create_go_dotplot <- function(df_list, top_n = 15) {
  df_list <- df_list[!vapply(df_list, is.null, logical(1))]
  if (length(df_list) == 0) return(NULL)
  parts <- lapply(names(df_list), function(nm) {
    d <- df_list[[nm]]
    d <- d[order(d$p.adjust), , drop = FALSE]
    d <- utils::head(d, top_n)
    ratio <- vapply(strsplit(as.character(d$GeneRatio), "/"),
                    function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    data.frame(Description = d$Description, GeneRatio = ratio,
               Count = d$Count, p.adjust = d$p.adjust,
               Panel = nm, stringsAsFactors = FALSE)
  })
  plot_df <- do.call(rbind, parts)
  if (is.null(plot_df) || nrow(plot_df) == 0) return(NULL)
  plot_df$Description <- factor(
    plot_df$Description,
    levels = unique(plot_df$Description[order(plot_df$GeneRatio)]))
  p <- ggplot(plot_df, aes(x = GeneRatio, y = Description,
                           size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#e74c3c", high = "#3498db", name = "p.adjust") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    labs(title = "GO Enrichment", x = "GeneRatio", y = NULL) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 9))
  if (length(df_list) > 1) p <- p + facet_wrap(~ Panel, scales = "free_y")
  p
}

# GO treemap: collapses redundant GO terms via rrvgo and tiles each parent
# term as a coloured rectangle whose area scales with the score (-log10 p).
# `rrvgo_out` is the named list returned by rrvgo_reduce(): list(sim_matrix,
# reduced). Returns NULL when no reduction was computed (rrvgo missing,
# too few terms, ontology not BP/MF/CC).
#
# `vp` is an optional grid viewport. rrvgo::treemapPlot() draws via grid
# (treemap::treemap), so to place several treemaps side-by-side (compare
# mode: Up vs Down) the caller pushes a grid layout and passes each panel's
# viewport here; treemapPlot forwards it to treemap::treemap() through `...`.
# When vp is NULL treemap draws to the whole current device (single panel).
create_go_treemap <- function(rrvgo_out, panel_title = NULL, vp = NULL) {
  if (is.null(rrvgo_out) || is.null(rrvgo_out$reduced)) return(NULL)
  if (!requireNamespace("rrvgo",   quietly = TRUE)) return(NULL)
  if (!requireNamespace("treemap", quietly = TRUE)) return(NULL)
  tryCatch({
    rrvgo::treemapPlot(rrvgo_out$reduced,
                       title = panel_title %||% "GO term semantic clusters",
                       vp = vp)
  }, error = function(e) NULL)
}

# GO scatter: 2D MDS projection of the GO term similarity matrix, points
# coloured by parent term and sized by score. Same inputs as the treemap
# helper; returns a ggplot object (or NULL on failure / no rrvgo data).
create_go_scatter <- function(rrvgo_out, panel_title = NULL) {
  if (is.null(rrvgo_out) || is.null(rrvgo_out$sim_matrix) ||
      is.null(rrvgo_out$reduced)) return(NULL)
  if (!requireNamespace("rrvgo", quietly = TRUE)) return(NULL)
  tryCatch({
    p <- rrvgo::scatterPlot(rrvgo_out$sim_matrix, rrvgo_out$reduced,
                            labelSize = 3)
    if (!is.null(panel_title)) p <- p + ggplot2::ggtitle(panel_title)
    p
  }, error = function(e) NULL)
}

# === ENSEMBL -> SYMBOL CONVERSION (biomaRt) ===
# TRUE if the object's gene identifiers are predominantly Ensembl IDs.
is_ensembl_object <- function(seurat_obj) {
  if (is.null(seurat_obj)) return(FALSE)
  genes <- rownames(seurat_obj)
  if (length(genes) == 0) return(FALSE)
  sample_genes <- if (length(genes) > 2000) sample(genes, 2000) else genes
  mean(grepl("^ENS", sample_genes)) > 0.5
}

# Query Ensembl BioMart for Ensembl-gene-ID -> symbol mappings. Runs inside a
# future worker; only the (small) ID vector is sent, never the Seurat object.
run_biomart_query <- function(ens_ids, species) {
  dataset <- if (species == "Homo sapiens") "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"
  # The "dataset not valid" error users hit is almost always a transient
  # Ensembl outage on the default host, not a genuinely missing dataset.
  # Try the main site then the regional mirrors before giving up, and run
  # the whole query (connect + getBM) per mirror so a mid-query failure also
  # rolls over to the next one.
  mirrors  <- c("www", "useast", "asia")
  last_err <- NULL
  for (mr in mirrors) {
    res <- tryCatch({
      mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset, mirror = mr)
      bm <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters    = "ensembl_gene_id",
        values     = ens_ids,
        mart       = mart
      )
      bm[!is.na(bm$external_gene_name) & nzchar(bm$external_gene_name), , drop = FALSE]
    }, error = function(e) { last_err <<- e; NULL })
    if (!is.null(res)) return(res)
  }
  stop("Could not reach Ensembl BioMart for dataset '", dataset, "' after trying ",
       "mirrors (", paste(mirrors, collapse = ", "), "). This is usually a temporary ",
       "Ensembl outage or a network/firewall block - please try again shortly. ",
       "Last error: ",
       if (!is.null(last_err)) conditionMessage(last_err) else "unknown",
       call. = FALSE)
}

# Rebuild a Seurat object's RNA assay under new gene names. Cell-level data
# (metadata, UMAP) is untouched; only the gene identifiers change.
rename_genes_in_object <- function(seurat_obj, new_names) {
  cnt <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  dat <- tryCatch(GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
                  error = function(e) NULL)
  rownames(cnt) <- new_names
  new_assay <- suppressWarnings(CreateAssayObject(counts = cnt))
  seurat_obj[["RNA"]] <- new_assay
  DefaultAssay(seurat_obj) <- "RNA"
  if (!is.null(dat) && nrow(dat) == length(new_names)) {
    rownames(dat) <- new_names
    seurat_obj <- tryCatch(
      SetAssayData(seurat_obj, assay = "RNA", layer = "data", new.data = dat),
      error = function(e) tryCatch(
        SetAssayData(seurat_obj, assay = "RNA", slot = "data", new.data = dat),
        error = function(e2) seurat_obj))
  }
  seurat_obj
}

# === PLOTTING & CACHE HELPERS ===
create_dea_volcano_plot <- function(dea_results, highlight_gene = NULL, p_threshold = 0.05, logfc_threshold = 0.25, filter_list = NULL, show_significant = FALSE) {
  if (is.null(dea_results) || nrow(dea_results) == 0) return(NULL)
  
  if (show_significant) {
    dea_results <- dea_results %>% filter(p_val_adj < p_threshold & abs(avg_log2FC) > logfc_threshold)
  }
  
  plot_df <- dea_results %>% mutate(NegLogP = -log10(p_val_adj), NegLogP = ifelse(is.infinite(NegLogP), 300, NegLogP), Significant = p_val_adj < p_threshold & abs(avg_log2FC) > logfc_threshold, IsHighlight = if (!is.null(highlight_gene) && highlight_gene != "") gene == highlight_gene else FALSE) %>% arrange(IsHighlight)
  top_genes <- head(plot_df$gene, 10)
  label_genes <- if (!is.null(highlight_gene) && highlight_gene != "") unique(c(top_genes, highlight_gene)) else top_genes
  plot_df$Label <- ifelse(plot_df$gene %in% label_genes, plot_df$gene, NA)
  filter_str <- "Global Filters: None"
  if (!is.null(filter_list) && length(filter_list) > 0) { f_parts <- sapply(filter_list, function(f) { val_str <- paste(head(f$vals, 2), collapse=","); if(length(f$vals)>2) val_str <- paste0(val_str, "..."); paste0(f$col, "=(", val_str, ")") }); filter_str <- paste("Global Filters:", paste(f_parts, collapse="; ")) }
  
  # UPDATED: Use linewidth instead of size
  p <- ggplot(plot_df, aes(x = avg_log2FC, y = NegLogP)) + 
    geom_point(aes(color = Significant), alpha = 0.5, size = 1.5, na.rm = TRUE) + 
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "#95a5a6", alpha = 0.6, linewidth = 0.8) + 
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#95a5a6", alpha = 0.6, linewidth = 0.8) + 
    scale_color_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#bdc3c7"), guide = "none") + 
    geom_text_repel(aes(label = Label), max.overlaps = 20, box.padding = 0.5, na.rm = TRUE) + 
    labs(title = paste0("Volcano Plot: ", dea_results$group1[1], " vs ", dea_results$group2[1]), subtitle = paste0(if(!is.null(highlight_gene) && highlight_gene != "") paste("Highlighted Gene:", highlight_gene) else "All Genes", "\n", filter_str), x = "Average Log2 Fold Change", y = "-log10(Adj. P-Value)") + 
    theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold"))
  
  if (!is.null(highlight_gene) && highlight_gene != "") p <- p + geom_point(data = subset(plot_df, IsHighlight), color = "black", fill = "yellow", shape = 21, size = 5, stroke = 1.5, na.rm = TRUE)
  return(p)
}

create_gene_violin_plot <- function(seurat_obj, gene, meta_col, group1, group2, filter_list = NULL) {
  if (is.null(gene) || gene == "") return(NULL)
  cells_mask <- rep(TRUE, ncol(seurat_obj)); if (!is.null(filter_list) && length(filter_list) > 0) { masks <- lapply(filter_list, function(f) seurat_obj@meta.data[[f$col]] %in% f$vals); cells_mask <- Reduce("&", masks) }
  sub_obj <- seurat_obj[, cells_mask]; keep_cells <- sub_obj@meta.data[[meta_col]] %in% c(group1, group2); sub_obj <- sub_obj[, keep_cells]
  expr <- GetAssayData(sub_obj, layer = "data")[gene, ]; df <- data.frame(Expression = expr, Group = sub_obj@meta.data[[meta_col]])
  ggplot(df, aes(x = Group, y = Expression, fill = Group)) + geom_violin(alpha = 0.6, scale = "width", trim = FALSE) + geom_jitter(width = 0.2, size = 0.5, alpha = 0.4) + geom_boxplot(width = 0.1, fill = "white", alpha = 0.8, outlier.shape = NA) + scale_fill_manual(values = c("#3498db", "#e74c3c")) + labs(title = paste("Expression Distribution:", gene), y = "Log-Normalized Expression", x = NULL) + theme_minimal(base_size = 14) + theme(legend.position = "none")
}

# Fast point layer for on-screen UMAP rendering. scattermore rasterises the
# points in C, which makes large atlases (>100k cells) draw in a fraction of a
# second instead of the ~20s that the vector geom_point path takes. We fall
# back to geom_point if scattermore is not installed, so the app still works
# without it (just slower). For high-resolution downloads we deliberately keep
# geom_point (raster = FALSE in the callers below) so exported figures stay
# crisp and vector-clean.
fast_points <- function(pointsize = 3, alpha = 0.9, pixels = c(2000, 2000)) {
  if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(pointsize = pointsize, alpha = alpha,
                                  pixels = pixels, na.rm = TRUE)
  } else {
    geom_point(size = 1.5, alpha = alpha, na.rm = TRUE)
  }
}

# Build a single-gene expression UMAP overlay from a precomputed UMAP dataframe.
# `df` has UMAP_1, UMAP_2 columns; `expr` is a numeric vector aligned to df rows.
# Used by the Metadata UMAP and Coexpression subtabs so they look identical.
# `raster = TRUE` uses the fast scattermore layer (on-screen); set FALSE for
# downloads to get a crisp vector geom_point rendering.
build_expression_umap <- function(df, gene_name, expr, limits = NULL, raster = TRUE) {
  df$expr <- as.numeric(expr)
  df <- df[order(df$expr), ]
  pts <- if (raster) fast_points() else geom_point(size = 1.2, alpha = 0.9, na.rm = TRUE)
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
    pts +
    scale_color_viridis_c(option = "plasma", name = "Expression", limits = limits) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Gene:", gene_name)) +
    theme(plot.title = element_text(face = "bold"))
}

# Parse free-text or uploaded gene list. Genes can be on separate lines, or
# comma/space/tab separated. Returns a deduplicated character vector.
parse_gene_input <- function(text = "", file_path = NULL) {
  raw <- character(0)
  if (!is.null(text) && nzchar(text)) {
    raw <- c(raw, strsplit(text, "[,;\\s]+", perl = TRUE)[[1]])
  }
  if (!is.null(file_path) && file.exists(file_path)) {
    lines <- tryCatch(readLines(file_path, warn = FALSE), error = function(e) character(0))
    raw <- c(raw, unlist(strsplit(lines, "[,;\\s]+", perl = TRUE)))
  }
  raw <- trimws(raw)
  raw <- raw[nzchar(raw)]
  unique(raw)
}

# Ensure standard QC metadata columns exist on a Seurat object.
# Returns the object unchanged if all three are present; otherwise computes
# the missing columns from the counts matrix (or, if counts is unavailable,
# from the normalised data layer). Mitochondrial percentage is computed
# from rows whose gene symbol matches MT-/mt-/Mt- (covers human/mouse and
# legacy annotations). A pre-existing `percent.mito` column is aliased to
# `percent.mt` for compatibility with older preprocessing pipelines.
ensure_qc_columns <- function(seurat_obj) {
  md <- seurat_obj@meta.data
  
  # Alias `percent.mito` -> `percent.mt` if the object came from an older pipeline.
  if (!"percent.mt" %in% colnames(md) && "percent.mito" %in% colnames(md)) {
    md$percent.mt <- md$percent.mito
  }
  
  needs <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  if (all(needs %in% colnames(md))) {
    seurat_obj@meta.data <- md
    return(seurat_obj)
  }
  
  # Prefer raw counts; fall back to the normalised data layer if counts is unavailable.
  counts <- tryCatch(GetAssayData(seurat_obj, layer = "counts"),
                     error = function(e) NULL)
  if (is.null(counts) || ncol(counts) == 0 || nrow(counts) == 0) {
    counts <- tryCatch(GetAssayData(seurat_obj, layer = "data"),
                       error = function(e) NULL)
  }
  if (is.null(counts) || ncol(counts) == 0 || nrow(counts) == 0) {
    # Cannot recompute. Record sentinel zeros for the columns still missing so
    # downstream code does not branch on column existence forever.
    if (!"nCount_RNA"   %in% colnames(md)) md$nCount_RNA   <- 0
    if (!"nFeature_RNA" %in% colnames(md)) md$nFeature_RNA <- 0
    if (!"percent.mt"   %in% colnames(md)) md$percent.mt   <- 0
    seurat_obj@meta.data <- md
    return(seurat_obj)
  }
  
  if (!"nCount_RNA"   %in% colnames(md)) md$nCount_RNA   <- Matrix::colSums(counts)
  if (!"nFeature_RNA" %in% colnames(md)) md$nFeature_RNA <- Matrix::colSums(counts > 0)
  if (!"percent.mt"   %in% colnames(md)) {
    mt_pat <- "^(MT|Mt|mt)[-\\.:]"
    mt_idx <- grep(mt_pat, rownames(counts))
    if (length(mt_idx) > 0) {
      mt_counts     <- Matrix::colSums(counts[mt_idx, , drop = FALSE])
      md$percent.mt <- 100 * mt_counts / pmax(md$nCount_RNA, 1)
    } else {
      md$percent.mt <- 0
    }
  }
  seurat_obj@meta.data <- md
  seurat_obj
}

# Single QC violin plot, optionally grouped by a metadata column. Returns a
# ggplot. log_y applies to nCount / nFeature only (not percent.mt).
create_qc_violin_plot <- function(seurat_obj, qc_col, group_col = NULL,
                                  filter_list = NULL, log_y = FALSE) {
  md <- seurat_obj@meta.data
  if (!qc_col %in% colnames(md)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = paste0("'", qc_col, "' not found."),
                               size = 6, color = "#e67e22") + theme_void())
  }
  cells_mask <- rep(TRUE, nrow(md))
  if (!is.null(filter_list) && length(filter_list) > 0) {
    masks <- lapply(filter_list, function(f) md[[f$col]] %in% f$vals)
    cells_mask <- Reduce("&", masks)
  }
  md <- md[cells_mask, , drop = FALSE]
  if (nrow(md) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "No cells after filtering.",
                               size = 6, color = "#e67e22") + theme_void())
  }
  use_grouping <- !is.null(group_col) && length(group_col) == 1 &&
    nzchar(group_col) && group_col %in% colnames(md)
  df <- if (use_grouping) {
    data.frame(value = md[[qc_col]], group = as.factor(md[[group_col]]))
  } else {
    data.frame(value = md[[qc_col]], group = factor("all cells"))
  }
  ylab <- switch(qc_col,
                 nFeature_RNA = "# detected genes per cell",
                 nCount_RNA   = "# UMIs per cell",
                 percent.mt   = "% mitochondrial UMIs",
                 qc_col)
  # 60-degree rotation + extra bottom margin guarantees long group labels
  # (e.g. "Lymphatic Endothelial Cell") never overlap with each other or
  # bleed into the next plot in the stacked layout.
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_violin(alpha = 0.7, scale = "width", trim = FALSE, na.rm = TRUE) +
    geom_jitter(width = 0.2, size = 0.3, alpha = 0.3, na.rm = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", alpha = 0.85,
                 outlier.shape = NA, na.rm = TRUE) +
    labs(title = qc_col, x = NULL, y = ylab) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x  = element_text(angle = 60, hjust = 1, vjust = 1, size = 11),
          plot.title   = element_text(face = "bold"),
          plot.margin  = margin(t = 12, r = 16, b = 80, l = 16))
  if (log_y && qc_col %in% c("nFeature_RNA", "nCount_RNA")) {
    p <- p + scale_y_log10()
  }
  p
}

# Multi-gene dot plot: mean expression and percent expressing per group.
# Color = z-scored mean expression across groups (Seurat::DotPlot default).
# Built without Seurat::DotPlot so it doesn't depend on a scale.data slot.
create_multi_gene_dotplot <- function(seurat_obj, genes, group_col,
                                      filter_list = NULL,
                                      cluster_genes = TRUE,
                                      scale_expression = TRUE) {
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "None of the requested genes are in this dataset.",
                               size = 6, color = "#e67e22") + theme_void())
  }
  md <- seurat_obj@meta.data
  cells_mask <- rep(TRUE, nrow(md))
  if (!is.null(filter_list) && length(filter_list) > 0) {
    masks <- lapply(filter_list, function(f) md[[f$col]] %in% f$vals)
    cells_mask <- Reduce("&", masks)
  }
  if (!group_col %in% colnames(md)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = paste0("Grouping column '", group_col, "' not found."),
                               size = 6, color = "#e67e22") + theme_void())
  }
  groups <- as.character(md[[group_col]])
  expr_mat <- GetAssayData(seurat_obj, layer = "data")[genes, cells_mask, drop = FALSE]
  groups   <- groups[cells_mask]
  group_levels <- sort(unique(groups))
  if (length(group_levels) < 2) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                               label = "Need at least 2 groups to build a dot plot.",
                               size = 6, color = "#e67e22") + theme_void())
  }
  mean_mat <- sapply(group_levels, function(g) {
    idx <- which(groups == g)
    Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
  })
  pct_mat <- sapply(group_levels, function(g) {
    idx <- which(groups == g)
    Matrix::rowMeans(expr_mat[, idx, drop = FALSE] > 0) * 100
  })
  if (length(genes) == 1) {
    mean_mat <- matrix(mean_mat, nrow = 1, dimnames = list(genes, group_levels))
    pct_mat  <- matrix(pct_mat,  nrow = 1, dimnames = list(genes, group_levels))
  } else {
    rownames(mean_mat) <- genes; rownames(pct_mat) <- genes
  }
  z_mat <- if (isTRUE(scale_expression) && ncol(mean_mat) >= 2) {
    z <- t(scale(t(mean_mat))); z[is.na(z)] <- 0; z
  } else mean_mat
  if (isTRUE(cluster_genes) && nrow(z_mat) >= 3) {
    hc <- tryCatch(hclust(dist(z_mat)), error = function(e) NULL)
    if (!is.null(hc)) genes <- rownames(z_mat)[hc$order]
  }
  long <- expand.grid(Gene = genes, Group = group_levels, stringsAsFactors = FALSE)
  long$AvgExpr <- mapply(function(g, gr) z_mat[g, gr],  long$Gene, long$Group)
  long$PctExpr <- mapply(function(g, gr) pct_mat[g, gr], long$Gene, long$Group)
  long$Gene  <- factor(long$Gene,  levels = rev(genes))
  long$Group <- factor(long$Group, levels = group_levels)
  fill_label <- if (isTRUE(scale_expression)) "Avg expr (z)" else "Avg expr"
  ggplot(long, aes(x = Group, y = Gene, color = AvgExpr, size = PctExpr)) +
    geom_point() +
    scale_color_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
                          midpoint = if (isTRUE(scale_expression)) 0 else
                            mean(long$AvgExpr, na.rm = TRUE),
                          name = fill_label) +
    scale_size_continuous(range = c(0, 8), name = "% cells expr.") +
    labs(title = "Multi-gene expression dot plot",
         x = group_col, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
}

# === CACHE FUNCTIONS ===

create_cache_key <- function(dataset_hash, filter_list, comparison_meta) {
  filter_str <- if(length(filter_list) > 0) { paste(sapply(filter_list, function(f) paste0(f$col, ":", paste(sort(f$vals), collapse = ","))), collapse = "|") } else { "NOFILTER" }
  comp_str <- paste0(comparison_meta$col, ":", paste(sort(comparison_meta$groups), collapse = ","))
  key_str <- paste(dataset_hash, filter_str, comp_str, sep = "_"); digest(key_str, algo = "md5")
}

# =============================================================================
# === R-SCRIPT GENERATORS (reproducibility export) ==========================
# =============================================================================
# Each generator returns a self-contained R script that re-runs the exact
# analysis the user just executed. Scripts depend only on the source Seurat
# object (path supplied via DATASET_PATH at run time) and a handful of CRAN /
# Bioconductor packages. They are deliberately verbose so a reviewer can read
# the analytical steps top-to-bottom without re-launching AtlasLens.
#
# r_literal(): wraps base::deparse() so an arbitrary R value (vector, named
# list, list-of-lists) round-trips into the script as a valid expression.

r_literal <- function(x) {
  paste(deparse(x, width.cutoff = 500L, control = c("keepInteger", "showAttributes")),
        collapse = "\n")
}

# Common prologue: library loading + dataset read + multi-dim filter loop.
# `filter_list` is the same structure stored on vals$active_filters and the
# cache record - a list of list(col=..., vals=...) elements.
.script_prologue <- function(filter_list, comment_header) {
  c(
    "# AtlasLens reproducibility script",
    paste0("# ", comment_header),
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "#",
    "# Set DATASET_PATH to the absolute path of the same .rds Seurat object",
    "# used in the AtlasLens session before running this script:",
    "#   Sys.setenv(DATASET_PATH = '/path/to/your.rds')",
    "",
    "suppressPackageStartupMessages({",
    "  library(Seurat)",
    "  library(SeuratObject)",
    "})",
    "",
    "dataset_path <- Sys.getenv('DATASET_PATH')",
    "if (!nzchar(dataset_path)) stop('Set DATASET_PATH before running.')",
    "seurat_obj <- readRDS(dataset_path)",
    "if (!'RNA' %in% names(seurat_obj@assays)) {",
    "  default_assay <- DefaultAssay(seurat_obj)",
    "  seurat_obj[['RNA']] <- seurat_obj[[default_assay]]",
    "}",
    "DefaultAssay(seurat_obj) <- 'RNA'",
    "if (inherits(seurat_obj[['RNA']], 'Assay5')) {",
    "  seurat_obj[['RNA']] <- SeuratObject::JoinLayers(seurat_obj[['RNA']])",
    "}",
    "",
    "# Multi-dimensional cell filter from the AtlasLens session.",
    paste0("filter_list <- ", r_literal(filter_list %||% list())),
    "keep_cells <- rep(TRUE, ncol(seurat_obj))",
    "for (f in filter_list) {",
    "  keep_cells <- keep_cells & (seurat_obj@meta.data[[f$col]] %in% f$vals)",
    "}",
    "seurat_obj <- seurat_obj[, keep_cells]",
    "message('After filter: ', ncol(seurat_obj), ' cells.')",
    ""
  )
}

# DEA: Wilcoxon (presto-accelerated) FindMarkers between two metadata groups.
generate_dea_script <- function(comparison_meta, filter_list,
                                highlight_gene = NULL,
                                p_thr = 0.05, fc_thr = 0.585) {
  groups <- comparison_meta$groups %||% comparison_meta$g1
  meta_col <- comparison_meta$col %||% ""
  g1 <- if (!is.null(comparison_meta$g1)) comparison_meta$g1 else groups[[1]]
  g2 <- if (!is.null(comparison_meta$g2)) comparison_meta$g2 else groups[[2]]
  body <- c(
    .script_prologue(filter_list,
                     "Differential Expression Analysis (Wilcoxon, presto-accelerated)"),
    paste0("meta_col <- ", r_literal(meta_col)),
    paste0("group1   <- ", r_literal(g1), "  # reference"),
    paste0("group2   <- ", r_literal(g2), "  # comparison"),
    paste0("p_thr    <- ", r_literal(p_thr)),
    paste0("fc_thr   <- ", r_literal(fc_thr)),
    "",
    "Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[meta_col]]",
    "seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)",
    "markers <- Seurat::FindMarkers(seurat_obj,",
    "  ident.1 = group2, ident.2 = group1,",
    "  test.use = 'wilcox', use_presto = TRUE,",
    "  logfc.threshold = fc_thr, min.pct = 0.1, verbose = FALSE,",
    "  assay = 'RNA', slot = 'data', recorrect_umi = FALSE)",
    "markers$gene <- rownames(markers)",
    "markers$Significant <- !is.na(markers$p_val_adj) & markers$p_val_adj < p_thr &",
    "  !is.na(markers$avg_log2FC) & abs(markers$avg_log2FC) >= fc_thr",
    "",
    "write.csv(markers, 'atlaslens_dea_results.csv', row.names = FALSE)",
    paste0("message('Wrote ', nrow(markers), ' rows to atlaslens_dea_results.csv',",
           " ' (', sum(markers$Significant, na.rm = TRUE), ' significant).')")
  )
  if (!is.null(highlight_gene) && nzchar(highlight_gene)) {
    body <- c(body,
              "",
              paste0("# Highlighted gene from the AtlasLens session: ", highlight_gene),
              paste0("print(markers[markers$gene == ", r_literal(highlight_gene), ", , drop = FALSE])"))
  }
  paste(body, collapse = "\n")
}

# COCOA: gene-of-interest co-regulation across the user-defined comparison
# groups. Mirrors the in-app worker (Hallmark gene sets, samplesize/n_sims
# defaults). The downloaded script will need geneCOCOA + msigdbr installed.
generate_cocoa_script <- function(gene, comparison_meta, filter_list,
                                  samplesize = 100L, n_sims = 100L) {
  meta_col <- comparison_meta$col %||% ""
  groups   <- comparison_meta$groups %||% list()
  body <- c(
    .script_prologue(filter_list,
                     "geneCOCOA - context-dependent gene function analysis"),
    "suppressPackageStartupMessages({",
    "  library(geneCOCOA)",
    "  library(msigdbr)",
    "})",
    "",
    paste0("gene_of_interest <- ", r_literal(gene)),
    paste0("meta_col         <- ", r_literal(meta_col)),
    paste0("groups           <- ", r_literal(groups),
           "  # named list of group definitions (each: vector of metadata values)"),
    paste0("samplesize       <- ", r_literal(samplesize)),
    paste0("n_sims           <- ", r_literal(n_sims)),
    "",
    "# Species auto-detected from gene symbol casing / Ensembl prefix.",
    "all_genes <- rownames(seurat_obj)",
    "species   <- if (mean(grepl('^ENSMUS', all_genes)) > 0.5 ||",
    "                  mean(all_genes == toupper(all_genes)) < 0.7)",
    "               'Mus musculus' else 'Homo sapiens'",
    "hallmark <- msigdbr::msigdbr(species = species, category = 'H')",
    "gene_sets <- split(hallmark$gene_symbol, hallmark$gs_name)",
    "",
    "# Run COCOA once per group, on a counts matrix subset to the cells in",
    "# that group. Returns one results frame per group; merge for plotting.",
    "results <- lapply(names(groups), function(g) {",
    "  cells_in_group <- seurat_obj@meta.data[[meta_col]] %in% groups[[g]]",
    "  counts_mat <- Seurat::GetAssayData(seurat_obj[, cells_in_group],",
    "                                      assay = 'RNA', layer = 'counts')",
    "  geneCOCOA::COCOA(counts_matrix = counts_mat, gene = gene_of_interest,",
    "                   gene_sets = gene_sets, samplesize = samplesize,",
    "                   n_sims = n_sims)",
    "})",
    "names(results) <- names(groups)",
    "",
    "saveRDS(results, 'atlaslens_cocoa_results.rds')",
    paste0("message('Saved COCOA results for ', length(results), ' group(s)',",
           " ' to atlaslens_cocoa_results.rds')")
  )
  paste(body, collapse = "\n")
}

# GO enrichment: from a DEA result file (the canonical input) or an uploaded
# gene table, run enrichGO (single-list mode) or compareCluster (compare).
generate_go_script <- function(go_settings, filter_list) {
  body <- c(
    .script_prologue(filter_list,
                     "GO Enrichment (clusterProfiler enrichGO / compareCluster + rrvgo)"),
    "suppressPackageStartupMessages({",
    "  library(clusterProfiler)",
    "  library(rrvgo)",
    "  library(org.Hs.eg.db)",
    "  library(org.Mm.eg.db)",
    "})",
    "",
    paste0("data_source <- ", r_literal(go_settings$data_source %||% "DEA")),
    paste0("mode        <- ", r_literal(go_settings$mode %||% "all"),
           "  # one of: up / down / all / compare"),
    paste0("ontology    <- ", r_literal(go_settings$ontology %||% "BP")),
    paste0("p_cut       <- ", r_literal(go_settings$p_cut    %||% 0.05)),
    paste0("q_cut       <- ", r_literal(go_settings$q_cut    %||% 0.1)),
    paste0("fc_cut      <- ", r_literal(go_settings$logfc_cut %||% 0.585)),
    "",
    "# Load the DEA results table written by atlaslens_dea_results.csv (or",
    "# substitute your own DEA-shaped CSV here).",
    "dea <- read.csv('atlaslens_dea_results.csv', stringsAsFactors = FALSE)",
    "sig <- dea[!is.na(dea$p_val_adj) & dea$p_val_adj < p_cut &",
    "             !is.na(dea$avg_log2FC) & abs(dea$avg_log2FC) >= fc_cut, , drop = FALSE]",
    "",
    "# Species + key type auto-detected from the gene-symbol casing / prefix.",
    "all_genes <- unique(sig$gene)",
    "species  <- if (mean(grepl('^ENSMUS', all_genes)) > 0.5 ||",
    "                 mean(all_genes == toupper(all_genes)) < 0.7)",
    "              'Mus musculus' else 'Homo sapiens'",
    "org_db   <- if (species == 'Homo sapiens') 'org.Hs.eg.db' else 'org.Mm.eg.db'",
    "key_type <- if (mean(grepl('^ENS', all_genes)) > 0.5) 'ENSEMBL' else 'SYMBOL'",
    "",
    "if (mode == 'compare') {",
    "  gene_list <- list(",
    "    Upregulated   = unique(sig$gene[sig$avg_log2FC > 0]),",
    "    Downregulated = unique(sig$gene[sig$avg_log2FC < 0])",
    "  )",
    "  ego <- clusterProfiler::compareCluster(",
    "    geneClusters = gene_list, fun = 'enrichGO',",
    "    OrgDb = org_db, keyType = key_type, ont = ontology,",
    "    pAdjustMethod = 'BH', pvalueCutoff = p_cut, qvalueCutoff = q_cut)",
    "  saveRDS(ego, 'atlaslens_go_compareCluster.rds')",
    "} else {",
    "  genes <- if (mode == 'up')   unique(sig$gene[sig$avg_log2FC > 0])",
    "           else if (mode == 'down') unique(sig$gene[sig$avg_log2FC < 0])",
    "           else unique(sig$gene)",
    "  ego <- clusterProfiler::enrichGO(",
    "    gene = genes, OrgDb = org_db, keyType = key_type, ont = ontology,",
    "    pAdjustMethod = 'BH', pvalueCutoff = p_cut, qvalueCutoff = q_cut)",
    "  saveRDS(ego, 'atlaslens_go_enrichGO.rds')",
    "  # Optional rrvgo semantic-similarity reduction for treemap / scatter.",
    "  if (!is.null(ego) && nrow(as.data.frame(ego)) >= 3) {",
    "    df <- as.data.frame(ego)",
    "    sim <- rrvgo::calculateSimMatrix(df$ID, orgdb = org_db, ont = ontology, method = 'Rel')",
    "    scores <- setNames(-log10(df$p.adjust), df$ID)",
    "    reduced <- rrvgo::reduceSimMatrix(sim, scores, threshold = 0.7, orgdb = org_db)",
    "    saveRDS(list(sim_matrix = sim, reduced = reduced), 'atlaslens_go_rrvgo.rds')",
    "  }",
    "}",
    "message('GO results saved.')"
  )
  paste(body, collapse = "\n")
}

save_dea_cache <- function(dataset_hash, results, filter_list, comparison_meta, highlight_gene) {
  cache_key <- create_cache_key(dataset_hash, filter_list, comparison_meta); cache_file <- file.path(CACHE_DIR, paste0("dea_", cache_key, ".qs"))
  cache_data <- list(dataset_hash = dataset_hash, results = results, filter_list = filter_list, comparison_meta = comparison_meta, timestamp = Sys.time(), gene = highlight_gene, type = "DEA")
  qs::qsave(cache_data, cache_file)
}

save_cocoa_cache <- function(dataset_hash, results, filter_list, comparison_meta, gene) {
  cache_key <- create_cache_key(dataset_hash, filter_list, comparison_meta); cache_file <- file.path(CACHE_DIR, paste0("cocoa_", cache_key, "_", gene, ".qs"))
  cache_data <- list(dataset_hash = dataset_hash, results = results, filter_list = filter_list, comparison_meta = comparison_meta, timestamp = Sys.time(), gene = gene, type = "COCOA")
  qs::qsave(cache_data, cache_file)
}

# Cache a GO Enrichment run for the History tab. `settings` carries the GO
# sidebar inputs so a saved run can be fully restored; comparison_meta is
# filled with GO-shaped values so list_cached_analyses() builds a uniform row.
save_go_cache <- function(dataset_hash, go_results, settings) {
  key_str <- paste(dataset_hash, settings$data_source, settings$mode,
                   settings$ontology, settings$p_cut, settings$q_cut,
                   settings$logfc_cut, sep = "_")
  cache_file <- file.path(CACHE_DIR, paste0("go_", digest(key_str, algo = "md5"), ".qs"))
  mode_label <- switch(settings$mode, up = "Upregulated", down = "Downregulated",
                       all = "All significant genes (combined)",
                       compare = "Compare Up vs Down", settings$mode)
  cache_data <- list(dataset_hash = dataset_hash, results = go_results,
                     go_settings = settings, filter_list = list(),
                     comparison_meta = list(col = settings$ontology, groups = mode_label),
                     timestamp = Sys.time(), gene = "GO Enrichment", type = "GO")
  qs::qsave(cache_data, cache_file)
}

list_cached_analyses <- function() {
  files <- list.files(CACHE_DIR, pattern = "^(dea_|cocoa_|go_).+\\.(rds|qs)$", full.names = TRUE); if (length(files) == 0) return(NULL)
  hist_list <- lapply(files, function(f) { tryCatch({ data <- if(grepl("\\.qs$", f)) qs::qread(f) else readRDS(f); data$cache_file_path <- basename(f); if (is.null(data$gene) || data$gene == "") data$gene <- "All Genes"; if (is.null(data$type)) data$type <- "DEA"; return(data) }, error = function(e) NULL) })
  hist_list <- hist_list[!sapply(hist_list, is.null)]; if (length(hist_list) == 0) return(NULL)
  df <- do.call(rbind, lapply(hist_list, function(h) { data.frame(type = h$type, gene = h$gene, desc = paste(h$comparison_meta$col, paste(h$comparison_meta$groups, collapse = "/"), sep = ": "), cache_file = h$cache_file_path, timestamp = h$timestamp, stringsAsFactors = FALSE) }))
  df <- df[order(df$timestamp, decreasing = TRUE), ]; return(df)
}

# --- Initial Load ---
# This executes once when the R process starts. 
# It pre-loads the default dataset into global memory.
default_load <- list(data = NULL, error = "No data loaded.", hash = NULL)
# DEFAULT_DATA_PATH is NULL when DATASET_PATH is unset; gate file.exists()
# so we don't pass NULL to it and the UI can show a clean error notice.
if (!is.null(DEFAULT_DATA_PATH) && file.exists(DEFAULT_DATA_PATH)) {
  tryCatch({
    seurat_obj <- readRDS(DEFAULT_DATA_PATH)
    if (!"RNA" %in% names(seurat_obj@assays)) { default_assay <- DefaultAssay(seurat_obj); seurat_obj[["RNA"]] <- seurat_obj[[default_assay]] }
    DefaultAssay(seurat_obj) <- "RNA"
    if (inherits(seurat_obj[["RNA"]], "Assay5")) seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
    if (ncol(seurat_obj@meta.data) > 0) { seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) { if (is.factor(x)) levels(x) <- iconv(levels(x), "UTF-8", "UTF-8", sub = "_"); if (is.character(x)) x <- iconv(x, "UTF-8", "UTF-8", sub = "_"); return(x) }) }
    
    # Ensure RNA@data is log-normalised at startup so the visualisation plots
    # (violin, temporal trend, expression UMAP) display real expression values
    # instead of raw counts. The DEA worker handles normalisation internally
    # and is unaffected by this block. SCT-processed objects are skipped
    # because SCT@data is already log-normalised corrected counts.
    #
    # NormalizeData() is not called directly because it raises a
    # meta.features mismatch on certain v5 Assay5 objects. Equivalent
    # arithmetic is applied directly to the sparse @x slot and written back
    # via a SetAssayData fallback chain. Failures here are logged as a
    # warning and load continues; downstream plots fall back to raw counts.
    # The saved .rds on disk is never modified.
    is_sct_loaded <- any(grepl("SCTransform", names(seurat_obj@commands)))
    if (!is_sct_loaded) {
      cnt <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
      dat <- tryCatch(GetAssayData(seurat_obj, assay = "RNA", layer = "data"),
                      error = function(e) NULL)
      data_equals_counts <- !is.null(dat) && length(cnt@x) == length(dat@x) &&
        identical(cnt@x, dat@x)
      # If @counts itself already contains non-integer values the object was
      # delivered with normalised counts (e.g. an scVI- or sctransform-
      # integrated object stored that way). Running LogNormalize on top
      # would double-normalise. Skip the injection in that case.
      counts_already_normalised <- !is.null(cnt) && length(cnt@x) > 0 &&
        any(cnt@x != round(cnt@x))
      if (counts_already_normalised) {
        if (is.null(dat) || nrow(dat) == 0) {
          message("Pre-load: RNA@counts contains non-integer values — treating as already-normalised; copying counts into the empty RNA@data so visualisation plots have an expression layer.")
          new_obj <- tryCatch(SetAssayData(seurat_obj, assay = "RNA", layer = "data", new.data = cnt), error = function(e) NULL)
          if (is.null(new_obj)) new_obj <- tryCatch(SetAssayData(seurat_obj, assay = "RNA", slot = "data", new.data = cnt), error = function(e) NULL)
          if (is.null(new_obj)) new_obj <- tryCatch({
            # v5-aware last-resort: Assay5 stores expression in @layers[["data"]],
            # the legacy v3/v4 Assay class still uses the @data slot.
            if (inherits(seurat_obj@assays$RNA, "Assay5")) {
              seurat_obj@assays$RNA@layers[["data"]] <- cnt
            } else {
              seurat_obj@assays$RNA@data <- cnt
            }
            seurat_obj
          }, error = function(e) NULL)
          if (!is.null(new_obj)) { seurat_obj <- new_obj; message("Pre-load: RNA@data populated from counts.") } else message("Pre-load: WARNING — could not populate RNA@data; visualisation plots may fail.")
        } else {
          message("Pre-load: RNA@counts already-normalised and RNA@data present — leaving RNA@data as-is.")
        }
      } else if (is.null(dat) || nrow(dat) == 0 || data_equals_counts) {
        message("Pre-load: RNA@data missing or equal to raw counts — computing manual LogNormalize and injecting into RNA@data ...")
        
        t_step <- Sys.time()
        flush.console()
        message("Pre-load:   step 1/5 colSums (cells = ", ncol(cnt),
                ", non-zero entries = ", length(cnt@x), ") ...")
        col_sums <- Matrix::colSums(cnt)
        col_sums[col_sums == 0] <- 1
        message("Pre-load:   step 1/5 done in ",
                round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s")
        
        t_step <- Sys.time(); flush.console()
        message("Pre-load:   step 2/5 divide @x by colsum ...")
        norm_mat <- cnt
        norm_mat@x <- norm_mat@x / rep.int(col_sums, diff(norm_mat@p))
        norm_mat@x <- norm_mat@x * 10000
        message("Pre-load:   step 2/5 done in ",
                round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s")
        
        t_step <- Sys.time(); flush.console()
        message("Pre-load:   step 3/5 log1p ...")
        norm_mat@x <- log1p(norm_mat@x)
        message("Pre-load:   step 3/5 done in ",
                round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s")
        
        injected <- FALSE
        t_step <- Sys.time(); flush.console()
        message("Pre-load:   step 4/5 SetAssayData (layer = 'data') ... [v5 path]")
        new_obj <- tryCatch(SetAssayData(seurat_obj, assay = "RNA", layer = "data", new.data = norm_mat),
                            error = function(e) { message("Pre-load:     v5 path errored: ", conditionMessage(e)); NULL })
        if (!is.null(new_obj)) {
          seurat_obj <- new_obj; injected <- TRUE
          message("Pre-load:   step 4/5 done in ",
                  round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s [v5 path]")
        }
        if (!injected) {
          message("Pre-load:   step 4/5 SetAssayData (slot = 'data') ... [v3/v4 path]")
          new_obj <- tryCatch(SetAssayData(seurat_obj, assay = "RNA", slot = "data", new.data = norm_mat),
                              error = function(e) { message("Pre-load:     v3/v4 path errored: ", conditionMessage(e)); NULL })
          if (!is.null(new_obj)) {
            seurat_obj <- new_obj; injected <- TRUE
            message("Pre-load:   step 4/5 done in ",
                    round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s [v3/v4 path]")
          }
        }
        if (!injected) {
          message("Pre-load:   step 4/5 direct slot assignment ... [last-resort path]")
          # v5-aware last-resort: Assay5 stores expression in @layers[["data"]],
          # the legacy v3/v4 Assay class still uses the @data slot.
          new_obj <- tryCatch({
            if (inherits(seurat_obj@assays$RNA, "Assay5")) {
              seurat_obj@assays$RNA@layers[["data"]] <- norm_mat
            } else {
              seurat_obj@assays$RNA@data <- norm_mat
            }
            seurat_obj
          }, error = function(e) { message("Pre-load:     direct path errored: ", conditionMessage(e)); NULL })
          if (!is.null(new_obj)) {
            seurat_obj <- new_obj; injected <- TRUE
            message("Pre-load:   step 4/5 done in ",
                    round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s [direct path]")
          }
        }
        
        t_step <- Sys.time(); flush.console()
        message("Pre-load:   step 5/5 free intermediate matrix ...")
        rm(norm_mat); gc()
        message("Pre-load:   step 5/5 done in ",
                round(as.numeric(difftime(Sys.time(), t_step, units = "secs")), 1), "s")
        
        if (injected) {
          message("Pre-load: log-normalized data injected into RNA@data.")
        } else {
          message("Pre-load: WARNING — could not inject log-normalized data into RNA@data. Visualization plots will display raw counts. DEA is unaffected.")
        }
        flush.console()
      }
    }
    
    if (!("umap" %in% names(seurat_obj@reductions))) {
      message("Pre-load: WARNING — no precomputed UMAP found in object. Running full Seurat UMAP pipeline now; this can take 5-15 minutes for a large dataset. Provide an object with a precomputed 'umap' reduction to skip this step.")
      message("Pre-load:   step 1/5 NormalizeData ...")
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      message("Pre-load:   step 2/5 FindVariableFeatures ...")
      seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
      message("Pre-load:   step 3/5 ScaleData ...")
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
      message("Pre-load:   step 4/5 RunPCA ...")
      seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
      message("Pre-load:   step 5/5 RunUMAP ...")
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
      message("Pre-load: UMAP computation complete.")
    }
    seurat_obj <- ensure_qc_columns(seurat_obj)
    dataset_hash <- digest(paste(ncol(seurat_obj), nrow(seurat_obj), colnames(seurat_obj@meta.data), collapse = "_"), algo = "md5")
    default_load <- list(data = seurat_obj, error = NULL, hash = dataset_hash)
  }, error = function(e) {
    # `<<-` is required so the error propagates to the outer `default_load`
    # binding and is then surfaced in the UI. Local `<-` would write only
    # into this handler's frame and the failure would be silently swallowed.
    msg <- conditionMessage(e)
    message("Pre-load: ERROR during dataset load: ", msg)
    default_load <<- list(data = NULL, error = msg, hash = NULL)
  })
}

# =============================================================================
# UI
# =============================================================================

custom_header <- tags$head(
  tags$script("
    var interval;
    Shiny.addCustomMessageHandler('start_timer', function(message) {
      clearInterval(interval);
      var startTime = Date.now();
      var id_selector = '#' + message.id;
      
      // Reset text immediately
      $(id_selector).text('Time elapsed: 0s');
      
      interval = setInterval(function() {
        var elapsed = Math.floor((Date.now() - startTime) / 1000);
        $(id_selector).text('Time elapsed: ' + elapsed + 's');
      }, 1000);
    });
    
    Shiny.addCustomMessageHandler('stop_timer', function(message) { 
      clearInterval(interval); 
    });
  "),
  tags$style(HTML("
    body { font-family: 'Segoe UI', sans-serif; background-color: #f5f5f5; zoom: 0.65; }
    .irs-line, .irs-grid, .irs-line-mid, .irs-line-left, .irs-line-right { pointer-events: none !important; }
    .irs-handle { pointer-events: auto !important; }
    /* Applied to a sidebar's config wrapper to freeze every control inside it
       (no clicks, no dropdowns) while an analysis runs or while a required
       prerequisite - e.g. Ensembl->symbol conversion - is still pending. */
    .panel-disabled { pointer-events: none; opacity: 0.5; }
    
    /* === PANELS & BOXES === */
    .title-panel { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 25px; margin-bottom: 25px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
    .feature-box { background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 5px solid #667eea; transition: transform 0.2s; }
    .feature-box:hover { transform: translateY(-2px); box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
    
    /* === V127 RESTORED OVERLAYS === */
    .large-status-overlay { position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%); z-index: 99999; background: white; padding: 50px 60px; border-radius: 15px; box-shadow: 0 15px 50px rgba(0,0,0,0.4); min-width: 500px; text-align: center; border: 2px solid #667eea; }
    .status-detail-text { font-family: monospace; color: #d35400; font-weight: bold; background: #ffe4c4; padding: 5px 10px; border-radius: 4px; display: inline-block; margin-top: 10px;}
    .timer-text { font-size: 14px; color: #666; margin-top: 5px; font-weight: bold; }
    
    /* === NAVIGATION HELPERS === */
    .section-header { border-bottom: 2px solid #667eea; padding-bottom: 5px; margin-bottom: 15px; color: #2c3e50; font-weight: bold; margin-top: 20px; font-size: 15px; letter-spacing: 0.5px; }
    .control-label { color: #34495e; font-weight: 600; }
    
    /* === INPUT HIGHLIGHTING (NEW) === */
    .form-control:focus, .selectize-input.focus { border-color: #667eea !important; box-shadow: 0 0 0 0.2rem rgba(102, 126, 234, 0.25) !important; }
    
    /* === CUSTOM TOOLTIPS (NEW - REPLACES NATIVE TITLE) === */
    .custom-tooltip { position: relative; display: inline-block; cursor: help; color: #17a2b8; margin-left: 5px; transition: color 0.2s; }
    .custom-tooltip:hover { color: #2c3e50; }
    .custom-tooltip .tooltiptext { visibility: hidden; width: 220px; background-color: #2c3e50; color: #fff; text-align: center; border-radius: 6px; padding: 8px 10px; position: absolute; z-index: 1000; bottom: 135%; left: 50%; margin-left: -110px; opacity: 0; transition: opacity 0.3s; font-size: 12px; font-weight: normal; box-shadow: 0 4px 10px rgba(0,0,0,0.2); line-height: 1.4; pointer-events: none; }
    .custom-tooltip .tooltiptext::after { content: ''; position: absolute; top: 100%; left: 50%; margin-left: -5px; border-width: 5px; border-style: solid; border-color: #2c3e50 transparent transparent transparent; }
    .custom-tooltip:hover .tooltiptext { visibility: visible; opacity: 1; }

    /* === HISTORY & BADGES === */
    .history-item-row { background: #f8f9fa; padding: 15px; margin: 10px 0; border-left: 5px solid #17a2b8; border-radius: 4px; cursor: pointer; transition: all 0.2s; }
    .history-item-row:hover { background-color: #e9ecef; border-left-width: 8px; transform: translateX(5px); }
    .filter-badge { display: inline-block; background: #667eea; color: white; padding: 4px 10px; border-radius: 12px; margin: 2px; font-size: 12px; }
    .active-filters { margin: 10px 0; padding: 10px; background: #f0f0f0; border-radius: 5px; border-left: 3px solid #667eea; }
    .type-badge { display: inline-block; padding: 2px 6px; border-radius: 4px; font-size: 10px; font-weight: bold; color: white; margin-right: 5px; vertical-align: middle; }
    .type-badge.DEA { background-color: #e67e22; }
    .type-badge.COCOA { background-color: #9b59b6; }
    .type-badge.GO { background-color: #16a085; }
    
    /* === INFO BOXES === */
    .gene-info-box { padding: 10px; background: #d4edda; border: 1px solid #c3e6cb; border-radius: 4px; color: #155724; font-size: 13px; margin-top: 5px; animation: fadeIn 0.5s; }
    .warning-box { padding: 10px; background: #fff3cd; border: 1px solid #ffeeba; border-radius: 4px; color: #856404; font-size: 13px; margin-top: 5px; animation: fadeIn 0.5s; }
    @keyframes fadeIn { from { opacity: 0; transform: translateY(-5px); } to { opacity: 1; transform: translateY(0); } }
    
    /* === FIX CURSOR BLINKING === */
    body { caret-color: transparent; }
    input, textarea, select, .selectize-input { caret-color: auto !important; }

    /* === NOTIFICATIONS: top-centre, wider, errors prominent === */
    /* Move the whole stack from the default bottom-right corner to the top
       centre of the viewport so users see error messages immediately. */
    #shiny-notification-panel {
      position: fixed !important;
      top: 18px !important;
      bottom: auto !important;
      left: 50% !important;
      right: auto !important;
      transform: translateX(-50%);
      width: auto !important;
      max-width: 700px;
      z-index: 99999;
    }
    .shiny-notification {
      width: 640px !important;
      max-width: 90vw;
      padding: 16px 20px !important;
      font-size: 15px !important;
      line-height: 1.45 !important;
      border-radius: 8px !important;
      box-shadow: 0 8px 24px rgba(0, 0, 0, 0.22) !important;
      margin-bottom: 10px !important;
      opacity: 1 !important;
    }
    /* Error notifications: a stronger red border + slightly larger so they
       stand out from messages / warnings, and the close button is bigger
       since errors are persistent (duration = NULL in the server code). */
    .shiny-notification-error {
      background-color: #fdecea !important;
      color: #7b1d13 !important;
      border-left: 6px solid #e74c3c !important;
      font-size: 16px !important;
      font-weight: 500;
    }
    .shiny-notification-warning {
      background-color: #fff8e1 !important;
      color: #7b5d10 !important;
      border-left: 6px solid #f1c40f !important;
    }
    .shiny-notification-message {
      background-color: #e8f4fd !important;
      color: #14517a !important;
      border-left: 6px solid #3498db !important;
    }
    .shiny-notification-close {
      font-size: 22px !important;
      line-height: 1 !important;
      top: 8px !important;
      right: 12px !important;
    }
  "))
)

ui <- navbarPage(
  title = "AtlasLens", id = "main_nav", theme = NULL, header = tagList(custom_header, useShinyjs()),
  
  tabPanel("Introduction", icon = icon("info-circle"),
           div(class = "title-panel", h1("Welcome to AtlasLens")),
           
           # === V127 LOAD SCREEN ===
           div(id = "global_loading_overlay", class = "large-status-overlay", style = "display: none;",
               h3(icon("database"), " Loading Dataset..."), p("Processing large files. Please do not refresh."),
               div(id = "global_status_text", class = "status-detail-text", "Starting load..."),
               div(id = "global_timer", class = "timer-text", "Time elapsed: 0s"), br(),
               div(class = "progress", style = "height: 8px;", div(class = "progress-bar progress-bar-striped active", style = "width: 100%; background-color: #667eea;"))
           ),
           
           fluidRow(
             column(12, 
                    div(class = "feature-box", style = "max-width: 800px; margin: 0 auto; text-align: center;",
                        h3("AtlasLens", style = "text-align: center; color: #2c3e50; font-weight: bold;"), 
                        # Dynamic Check
                        uiOutput("dataset_status_ui"),
                        div(style = "text-align: center; margin: 30px 0;", 
                            actionButton("go_explore", "Start Analysis →", class = "btn-primary btn-lg", icon = icon("play-circle"), style = "padding: 15px 30px; font-size: 18px;"))
                    )
             )
           ),
           
           # Landing content. In the published design this section is
           # populated from the optional JSON configuration file (manuscript
           # §2: personalized "introduction" and "dataset information" shown
           # on the landing page). The JSON loader is not wired up yet, so a
           # neutral placeholder renders here until a config is supplied.
           fluidRow(
             column(12,
                    div(style = "max-width: 800px; margin: 20px auto; color: #555;",
                        uiOutput("landing_content"))
             )
           )
  ),
  # === EXPLORE (VERTICAL UI) ===
  tabPanel("Dataset Exploration", icon = icon("database"),
           div(class = "title-panel", h1("Dataset Exploration")),
           fluidRow(column(12, wellPanel(h3(icon("microscope"), " Dataset Overview"), uiOutput("dataset_info")))),
           fluidRow(column(12, wellPanel(
             h3(icon("map"), " Interactive Visualisations"),
             sidebarLayout(
               sidebarPanel(width = 3,
                            # Shared cell filter. Adding or removing a filter
                            # row updates vals$explore_mask immediately; each
                            # subtab consumes the mask when its own "Show Plot"
                            # button is pressed (or on input change for the
                            # Metadata UMAP subtab).
                            div(class = "section-header", style = "margin-top: 0;", icon("filter"),
                                " Filter cells (shared)"),
                            tags$label("Subset Data:", style = "font-weight: bold;"),
                            info_icon("Pick a column, choose values, then click Add Filter. The filter applies to every subtab on this page."),
                            uiOutput("explore_filter_controls_ui"),
                            uiOutput("explore_active_filters_ui"),
                            hr(),
                            
                            # ---- Metadata UMAP subtab controls ----
                            conditionalPanel(condition = "input.explore_subtab == 'meta' || input.explore_subtab == null",
                                             div(class = "section-header", icon("sliders-h"), " Metadata UMAP"),
                                             tags$label("Color left UMAP by:", style = "font-weight: bold;"),
                                             info_icon("Pick the metadata column used to color the left UMAP."),
                                             selectInput("explore_filter_meta_col", NULL, NULL),
                                             tags$label("Show / hide groups:", style = "font-weight: bold; margin-top: 10px;"),
                                             info_icon("Unchecking values hides those cells from BOTH UMAPs."),
                                             uiOutput("explore_filter_meta_ui"),
                                             hr(),
                                             tags$label("Right UMAP (expression):", style = "font-weight: bold;"),
                                             info_icon("Pick a gene to color the right UMAP."),
                                             selectizeInput("explore_gene", NULL, NULL,
                                                            options = list(placeholder = "Type gene name...",
                                                                           loadThrottle = 500, maxOptions = 100)),
                                             uiOutput("explore_gene_info"),
                                             sliderInput("explore_umap_height", "Plot height (px):",
                                                         min = 500, max = 2400, value = 1000, step = 50),
                                             actionButton("show_umap_plot", "Show Plot",
                                                          icon = icon("eye"), class = "btn-primary btn-block")
                            ),
                            
                            # ---- Coexpression subtab controls ----
                            conditionalPanel(condition = "input.explore_subtab == 'coexp'",
                                             div(class = "section-header", icon("dna"), " Coexpression"),
                                             tags$label("Gene 1:", style = "font-weight: bold;"),
                                             selectizeInput("coexp_gene_a", NULL, NULL,
                                                            options = list(placeholder = "Type gene name...",
                                                                           loadThrottle = 500, maxOptions = 100)),
                                             tags$label("Gene 2:", style = "font-weight: bold; margin-top: 10px;"),
                                             selectizeInput("coexp_gene_b", NULL, NULL,
                                                            options = list(placeholder = "Type gene name...",
                                                                           loadThrottle = 500, maxOptions = 100)),
                                             sliderInput("coexp_plot_height", "Plot height (px):",
                                                         min = 500, max = 1400, value = 1000, step = 50),
                                             actionButton("show_coexp_plot", "Show Plot",
                                                          icon = icon("eye"), class = "btn-primary btn-block")
                            ),
                            
                            # ---- Multi-gene dot plot subtab controls ----
                            conditionalPanel(condition = "input.explore_subtab == 'dotplot'",
                                             div(class = "section-header", icon("braille"), " Multi-gene Dot Plot"),
                                             tags$label("Group cells by:", style = "font-weight: bold;"),
                                             info_icon("X-axis."),
                                             selectInput("dotplot_group_col", NULL, NULL),
                                             tags$label("Genes:", style = "font-weight: bold; margin-top: 10px;"),
                                             info_icon("Type to search."),
                                             selectizeInput("dotplot_genes", NULL, NULL, multiple = TRUE,
                                                            options = list(placeholder = "Type gene name(s)...",
                                                                           loadThrottle = 500, maxOptions = 100)),
                                             tags$label("...or upload a gene list (.txt):", style = "font-weight: bold;"),
                                             info_icon("Plain text. One gene per line, or comma/space-separated. Appended to typed genes."),
                                             fileInput("dotplot_genes_file", NULL,
                                                       accept = c("text/plain", ".txt", ".csv", ".tsv")),
                                             checkboxInput("dotplot_cluster_genes",
                                                           "Cluster genes (hierarchical)", value = TRUE),
                                             checkboxInput("dotplot_scale_expr",
                                                           "Z-score across groups", value = TRUE),
                                             sliderInput("dotplot_plot_height", "Plot height (px):",
                                                         min = 300, max = 1200, value = 700, step = 50),
                                             actionButton("show_dotplot", "Show Plot",
                                                          icon = icon("eye"), class = "btn-primary btn-block")
                            ),
                            
                            # ---- QC violins subtab controls ----
                            conditionalPanel(condition = "input.explore_subtab == 'qc'",
                                             div(class = "section-header", icon("chart-bar"), " QC Violins"),
                                             p(style = "font-size: 12px; color: #555;",
                                               icon("info-circle"),
                                               " Standard scverse QC metrics: nFeature_RNA, nCount_RNA, percent.mt ",
                                               "(computed at startup if missing, MT-/mt- prefix)."),
                                             tags$label("Group cells by:", style = "font-weight: bold;"),
                                             info_icon("Typically a celltype or condition column."),
                                             selectInput("qc_group_col", NULL, NULL),
                                             checkboxInput("qc_log_y",
                                                           "Log-scale y-axis (counts / features only)", value = FALSE),
                                             sliderInput("qc_plot_height", "Plot height per violin (px):",
                                                         min = 300, max = 900, value = 550, step = 25),
                                             actionButton("show_qc_plot", "Show Plot",
                                                          icon = icon("eye"), class = "btn-primary btn-block")
                            )
               ),
               mainPanel(width = 9,
                         tabsetPanel(id = "explore_subtab", type = "tabs",
                                     
                                     # === SUBTAB 1: Metadata UMAP (2-panel view) ===
                                     tabPanel("Metadata UMAP", value = "meta", icon = icon("map"),
                                              br(),
                                              div(style = "display: flex; align-items: center; margin-bottom: 10px;",
                                                  p(icon("search-plus"), " Tip: click and drag on either plot to zoom.",
                                                    style = "color: #555; font-size: 13px; margin: 0; margin-right: 15px;"),
                                                  actionButton("reset_umap_zoom", "Reset Zoom",
                                                               icon = icon("expand"), class = "btn-sm btn-outline-secondary")
                                              ),
                                              fluidRow(
                                                column(6, h4("Metadata UMAP"),
                                                       withSpinner(plotOutput("explore_umap_meta",
                                                                              width = "100%", height = "auto",
                                                                              click = "umap_meta_click",
                                                                              dblclick = "umap_dblclick",
                                                                              brush = brushOpts(id = "umap_brush", resetOnNew = TRUE)),
                                                                   type = 6, color = "#667eea"),
                                                       br(), downloadButton("explore_umap_meta_download",
                                                                            "Download (PNG)", class = "btn-success")),
                                                column(6, h4("Expression UMAP"),
                                                       withSpinner(plotOutput("explore_umap_expr",
                                                                              width = "100%", height = "auto",
                                                                              click = "umap_expr_click",
                                                                              dblclick = "umap_dblclick",
                                                                              brush = brushOpts(id = "umap_brush", resetOnNew = TRUE)),
                                                                   type = 6, color = "#667eea"),
                                                       br(), downloadButton("explore_umap_expr_download",
                                                                            "Download (PNG)", class = "btn-success"))
                                              ),
                                              # Full-width metadata legend rendered as a CSS grid of colour-chip /
                                              # label pairs. The in-plot legend on the Metadata UMAP is suppressed
                                              # (see server) so the two UMAPs stay identically square; every
                                              # category - no matter how many or how long the names - is shown
                                              # here with its full label and exact swatch colour.
                                              fluidRow(column(12, uiOutput("explore_umap_meta_legend"))),
                                              fluidRow(column(12, wellPanel(style = "margin-top: 20px;",
                                                                            h4(icon("hand-pointer"), "Cell Inspector", style = "margin-top: 0;"),
                                                                            p("Click a cell on either UMAP to view its metadata."),
                                                                            htmlOutput("explore_click_info"))))
                                     ),
                                     
                                     # === SUBTAB 2: Coexpression (two side-by-side gene UMAPs) ===
                                     tabPanel("Coexpression", value = "coexp", icon = icon("dna"),
                                              br(),
                                              p(style = "color: #555;", icon("info-circle"),
                                                " Side-by-side expression of two genes on the same UMAP. ",
                                                "Pick genes in the sidebar and click ", tags$b("Show Plot"),
                                                " — plots do not re-render automatically."),
                                              fluidRow(
                                                column(6, h4(textOutput("coexp_g1_title")),
                                                       withSpinner(plotOutput("coexp_umap_g1",
                                                                              width = "100%", height = "auto"),
                                                                   type = 6, color = "#667eea"),
                                                       br(), downloadButton("coexp_g1_download",
                                                                            "Download (PNG)", class = "btn-success")),
                                                column(6, h4(textOutput("coexp_g2_title")),
                                                       withSpinner(plotOutput("coexp_umap_g2",
                                                                              width = "100%", height = "auto"),
                                                                   type = 6, color = "#667eea"),
                                                       br(), downloadButton("coexp_g2_download",
                                                                            "Download (PNG)", class = "btn-success"))
                                              )
                                     ),
                                     
                                     # === SUBTAB 3: Multi-gene Dot Plot ===
                                     tabPanel("Multi-gene Dot Plot", value = "dotplot", icon = icon("braille"),
                                              br(),
                                              p(style = "color: #555;", icon("info-circle"),
                                                " Average expression and the percent of cells expressing each gene, ",
                                                "per group. Color = mean expression z-scored across groups (centred and scaled per gene); dot size = percentage of cells in the group expressing the gene."),
                                              uiOutput("dotplot_status_msg"),
                                              withSpinner(plotOutput("explore_dotplot", height = "700px"),
                                                          type = 6, color = "#3498db"),
                                              br(),
                                              downloadButton("explore_dotplot_download_png", "Download Plot (PNG)", class = "btn-success"),
                                              downloadButton("explore_dotplot_download_pdf", "Download Plot (PDF)", class = "btn-info"),
                                              downloadButton("explore_dotplot_download_csv", "Download Data (CSV)", class = "btn-info")
                                     ),
                                     
                                     # === SUBTAB 4: QC Violins (stacked vertically, not in a row) ===
                                     tabPanel("QC Violins", value = "qc", icon = icon("chart-bar"),
                                              br(),
                                              p(style = "color: #555;", icon("info-circle"),
                                                " Single-cell QC: number of detected genes, total UMI counts, ",
                                                "and percent of mitochondrial UMIs per cell. ",
                                                "Pick a grouping column and click ", tags$b("Show Plot"), "."),
                                              uiOutput("qc_status_msg"),
                                              fluidRow(column(12, h4("nFeature_RNA"),
                                                              withSpinner(plotOutput("qc_violin_nFeature"),
                                                                          type = 6, color = "#27ae60"))),
                                              tags$div(style = "height: 40px;"),
                                              fluidRow(column(12, h4("nCount_RNA"),
                                                              withSpinner(plotOutput("qc_violin_nCount"),
                                                                          type = 6, color = "#27ae60"))),
                                              tags$div(style = "height: 40px;"),
                                              fluidRow(column(12, h4("percent.mt"),
                                                              withSpinner(plotOutput("qc_violin_percentMt"),
                                                                          type = 6, color = "#27ae60"))),
                                              br(),
                                              downloadButton("qc_download_png", "Download all QC violins (PNG)", class = "btn-success"),
                                              downloadButton("qc_download_pdf", "Download all QC violins (PDF)", class = "btn-info")
                                     )
                         )
               )
             )
           )))
  ),
  
  # === COCOA ===
  tabPanel("Gene Function (COCOA)", value = "Analysis", icon = icon("bar-chart"),
           # === V127 LOAD SCREEN ===
           div(id = "cocoa_loading_overlay", class = "large-status-overlay", style = "display: none;",
               h3(icon("spinner", class = "fa-spin"), " Analyzing..."), p("Performing Calculations..."),
               div(id = "cocoa_timer", class = "timer-text", "Time elapsed: 0s"), br(),
               div(class = "progress", style = "height: 8px;", div(class = "progress-bar progress-bar-striped active", style = "width: 100%; background-color: #667eea;"))
           ),
           div(class = "title-panel", h1("Gene Function (COCOA)")),
           sidebarLayout(sidebarPanel(width = 3, div(id = "cocoa_config_wrap", div(class = "section-header", style = "margin-top: 0;", icon("dna"), " 1. Select Gene of Interest"), tags$label("Target Gene:", style = "font-weight: bold;"), info_icon("The gene whose pathway co-regulation you want to analyze."), selectizeInput("analysis_gene", NULL, NULL, options = list(placeholder = "Type to search...", loadThrottle = 500, maxOptions = 100)), uiOutput("analysis_gene_info"), div(class = "section-header", icon("filter"), " 2. Extra Filtering (Optional)"), tags$label("Subset Data:"), info_icon("Apply filters to narrow down the cells"), uiOutput("filter_controls_ui"), uiOutput("active_filters_ui"), div(class = "section-header", icon("columns"), " 3. Group Comparison"), tags$label("Grouping Column:"), info_icon("Select the metadata column defining the groups you want to compare."), selectInput("analysis_meta_col", NULL, choices = NULL), tags$label("Select Groups to Compare:", class = "control-label", style = "font-weight: bold; margin-top: 10px;"), info_icon("Check the specific groups to run analysis on."), uiOutput("analysis_groups_ui"), uiOutput("analysis_cell_count_ui"), hr(), div(class = "section-header", icon("cogs"), " Options"), numericInput("p_val_thresh", "Adjusted P-value threshold:", value = 0.05, min = 0, max = 1, step = 0.001), helpText("Enter a value between 0 and 1."), sliderInput("cocoa_plot_height", "Plot height (px):", min = 400, max = 1400, value = 800, step = 50)), br(), uiOutput("run_btn_ui"), br(), uiOutput("download_ui")), mainPanel(width = 9, h3(icon("chart-bar"), " Analysis Results"), uiOutput("analysis_status_msg"), uiOutput("analysis_results_ui")))
  ),
  
  # === DEA ===
  tabPanel("Differential Expression", value = "DEA", icon = icon("chart-area"),
           # === V127 LOAD SCREEN ===
           div(id = "dea_loading_overlay", class = "large-status-overlay", style = "display: none;",
               h3(icon("spinner", class = "fa-spin"), " Analyzing..."), p("Calculating Genome-Wide DEA..."),
               div(id = "dea_timer", class = "timer-text", "Time elapsed: 0s"), br(),
               div(class = "progress", style = "height: 8px;", div(class = "progress-bar progress-bar-striped active", style = "width: 100%; background-color: #27ae60;"))
           ),
           div(class = "title-panel", h1("Differential Expression Analysis")),
           sidebarLayout(sidebarPanel(width = 3, 
                                      div(id = "dea_config_wrap",
                                      div(class = "section-header", style = "margin-top: 0;", icon("dna"), " 1. Select Gene to Highlight (Optional)"), tags$label("Highlight Gene:", style = "font-weight: bold;"), info_icon("Optional. Select a gene to highlight on the plot."), selectizeInput("dea_gene", NULL, NULL, options = list(placeholder = "Type gene name (optional)...", loadThrottle = 500, maxOptions = 100)), uiOutput("dea_gene_info"), hr(), 
                                      div(class = "section-header", icon("filter"), " 2. Extra Filtering (Optional)"), tags$label("Subset Data:"), info_icon("Narrow down cells for comparison."), uiOutput("dea_filter_controls_ui"), uiOutput("dea_active_filters_ui"), hr(), 
                                      div(class = "section-header", icon("layer-group"), " 3. Select Metadata"), tags$label("Metadata Column:", style = "font-weight: bold;"), info_icon("Choose the column defining conditions."), selectInput("dea_meta_col", NULL, choices = NULL), hr(),
                                      div(class = "section-header", icon("code-branch"), " 4. Choose Two Submetadata"), tags$label("Submetadata 1 (Reference):", style = "font-weight: bold;"), selectInput("dea_group1", NULL, choices = NULL), tags$label("Submetadata 2 (Comparison):", style = "font-weight: bold;", style = "margin-top: 10px;"), selectInput("dea_group2", NULL, choices = NULL), uiOutput("dea_cell_count_ui"), hr(), 
                                      div(class = "section-header", icon("cogs"), " 5. Options"),
                                      checkboxInput("dea_show_significant", "Show significant only", value = FALSE),
                                      tags$label("Adjusted P-value threshold:", style = "font-weight: bold;"),
                                      info_icon("Genes with p_val_adj below this are flagged significant and highlighted in the results table. Type to override."),
                                      numericInput("dea_p_threshold", NULL,
                                                   value = 0.05, min = 0, max = 1, step = 0.001),
                                      helpText("Enter a value between 0 and 1."),
                                      tags$label("Log2FC threshold:", style = "font-weight: bold;"),
                                      info_icon("Effect-size cutoff used together with the p-value cutoff."),
                                      numericInput("dea_logfc_threshold", NULL,
                                                   value = 0.585, min = 0, max = 10, step = 0.05),
                                      helpText("Enter a value of 0 or greater."),
                                      sliderInput("dea_volcano_height", "Volcano plot height (px):",
                                                  min = 400, max = 1100, value = 650, step = 25),
                                      sliderInput("dea_violin_height", "Single-gene violin height (px):",
                                                  min = 300, max = 800, value = 500, step = 25),
                                      hr(),
                                      div(class = "section-header", icon("file-upload"), " 6. Restrict to a Custom Gene List (Optional)"),
                                      tags$label("Upload a CSV or TXT file with one column of gene symbols. The volcano plot, results table, and downloads will be restricted to genes in this list. The Time Series tab's cluster export produces a compatible file."),
                                      info_icon("Accepts a CSV with a 'gene' column or a plain text file with one gene per line. Use the Time Series tab's 'Download cluster genes' button to generate one."),
                                      fileInput("dea_gene_upload", NULL,
                                                accept = c(".csv", ".txt", ".tsv"),
                                                buttonLabel = "Choose file...",
                                                placeholder = "No file uploaded"),
                                      uiOutput("dea_custom_genes_status_ui"),
                                      uiOutput("dea_custom_genes_clear_ui")),
                                      br(),
                                      uiOutput("dea_run_btn_ui")),
                         mainPanel(width = 9, h3(icon("chart-line"), " Analysis Results"), uiOutput("dea_status_msg"), uiOutput("dea_results_ui")))
  ),
  
  # === GO ENRICHMENT ===
  tabPanel("GO Enrichment", value = "GO", icon = icon("project-diagram"),
           # Loading overlay
           div(id = "go_loading_overlay", class = "large-status-overlay", style = "display: none;",
               h3(icon("spinner", class = "fa-spin"), " Analyzing..."), p("Running GO Enrichment..."),
               div(id = "go_timer", class = "timer-text", "Time elapsed: 0s"), br(),
               div(class = "progress", style = "height: 8px;", div(class = "progress-bar progress-bar-striped active", style = "width: 100%; background-color: #8e44ad;"))
           ),
           div(class = "title-panel", h1("GO Enrichment Analysis")),
           sidebarLayout(
             sidebarPanel(width = 3,
                          # --- Data Source ---
                          div(id = "go_config_wrap",
                          div(class = "section-header", style = "margin-top: 0;", icon("database"), " 1. Data Source"),
                          radioButtons("go_data_source", NULL, choices = c("Use DEA Results" = "dea", "Upload gene list" = "upload"), selected = "dea", inline = TRUE),
                          conditionalPanel(condition = "input.go_data_source == 'upload'",
                                           fileInput("go_upload_file", "Upload your significant-gene list:", accept = c(".csv", ".tsv", ".txt")),
                                           helpText("One gene per line, or a table with a 'gene' column. An optional fold-change column (e.g. avg_log2FC) enables the Up-vs-Down comparison.")
                          ),
                          conditionalPanel(condition = "input.go_data_source == 'dea'",
                                           uiOutput("go_dea_status_ui")
                          ),
                          hr(),
                          # --- Mode ---
                          div(class = "section-header", icon("layer-group"), " 2. Analysis Mode"),
                          radioButtons("go_mode", NULL,
                                       choices = c("Compare Up vs Down"      = "compare",
                                                   "All DE genes (combined)" = "all"),
                                       selected = "compare"),
                          helpText("'Compare' enriches the up- and down-regulated genes separately and shows them side by side. 'All DE genes (combined)' pools every DE gene into one enrichment - the only mode that detects pathways perturbed in both directions at once."),
                          hr(),
                          # --- Thresholds ---
                          div(class = "section-header", icon("sliders-h"), " 3. Thresholds"),
                          numericInput("go_p_cutoff", "Significance cutoff (adjusted p):", value = 0.05, min = 0, max = 1, step = 0.01),
                          helpText("Enter a value between 0 and 1."),
                          numericInput("go_logfc_cutoff", "Min |log2FC| (gene selection / up-down split):", value = 0.585, min = 0, max = 10, step = 0.1),
                          helpText("Enter a value of 0 or greater."),
                          numericInput("go_q_cutoff", "Q-value (FDR) cutoff:", value = 0.1, min = 0, max = 1, step = 0.01),
                          helpText("Enter a value between 0 and 1."),
                          hr(),
                          # --- Ontology ---
                          div(class = "section-header", icon("sitemap"), " 4. GO Ontology"),
                          selectInput("go_ontology", "Ontology:", choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"), selected = "BP"),
                          helpText("Species and gene-ID type are detected automatically from your genes."),
                          hr()),
                          # --- Run button ---
                          br(),
                          uiOutput("go_run_btn_ui")
             ),
             mainPanel(width = 9,
                       h3(icon("project-diagram"), " GO Enrichment Results"),
                       uiOutput("go_status_msg"),
                       uiOutput("go_results_ui")
             )
           )
  ),
  
  # === TIME SERIES ===
  tabPanel("Time Series Information", icon = icon("clock"),
           div(class = "title-panel", h1("Time Series Information")),
           sidebarLayout(
             sidebarPanel(width = 3,
                          # Role columns (timepoint / condition / celltype /
                          # dataset) are auto-detected via detect_role_columns()
                          # and stored in vals$ts_roles_detected. The four
                          # selectInputs below are kept in the DOM but hidden
                          # so the existing ts_roles() reactive still resolves
                          # them; they are populated programmatically and not
                          # shown to the user. To re-expose the override UI,
                          # remove the surrounding shinyjs::hidden() wrapper.
                          shinyjs::hidden(
                            tags$div(id = "ts_role_overrides",
                                     selectInput("ts_col_timepoint", "Timepoint column:", choices = NULL),
                                     selectInput("ts_col_condition", "Condition column:", choices = NULL),
                                     selectInput("ts_col_celltype",  "Cell-type column:", choices = NULL),
                                     selectInput("ts_col_dataset",   "Dataset / batch column (optional):", choices = NULL))
                          ),
                          div(id = "ts_config_wrap",
                          div(class = "section-header", style = "margin-top: 0;", icon("database"), " 1. Select Dataset"),
                          uiOutput("ts_dataset_ui"),
                          uiOutput("ts_dataset_info_ui"),
                          hr(),
                          div(class = "section-header", icon("filter"), " 2. Select Cohort"),
                          tags$label("Condition:", style = "font-weight: bold;"),
                          info_icon("Pick a condition. If a control condition  is detected, it is automatically included on every plot for reference."),
                          selectInput("ts_condition", NULL, choices = NULL),
                          uiOutput("ts_control_badge_ui"),
                          tags$label("Celltype:", style = "font-weight: bold; margin-top: 10px;"),
                          selectInput("ts_celltype", NULL, choices = NULL),
                          hr(),
                          div(class = "section-header", icon("dna"), " 3. Select Gene(s)"),
                          tags$label("Gene (single for Trend/Violin, multiple for Heatmap):", style = "font-weight: bold;"),
                          # A sentinel "All genes (heatmap)" entry is the default selection
                          # and is mutually exclusive with individual gene selections; the
                          # server observer enforces the exclusion in both directions.
                          selectizeInput("ts_gene", NULL, NULL, multiple = TRUE,
                                         options = list(placeholder = "Type gene name(s) or keep 'All genes'...",
                                                        loadThrottle = 500, maxOptions = 100)),
                          tags$small(style = "color: #666;",
                                     "Keep 'All genes (heatmap)' for a full transcriptome heatmap, ",
                                     "or type gene names. ",
                                     "Trend and Violin subtabs use the first selected gene."),
                          hr(),
                          div(class = "section-header", icon("clock"), " 4. Select Timepoints"),
                          tags$label("Available Timepoints:", style = "font-weight: bold;"),
                          uiOutput("ts_timepoints_ui"),
                          hr(),
                          div(class = "section-header", icon("ruler-vertical"), " Plot size"),
                          sliderInput("ts_heatmap_height_override", "Heatmap height (px):",
                                      min = 0, max = 3000, value = 400, step = 50),
                          tags$small(style = "color: #666;",
                                     "0 = auto-fit to gene count."),
                          sliderInput("ts_trend_height", "Temporal Trend height (px):",
                                      min = 400, max = 1200, value = 600, step = 25),
                          sliderInput("ts_violin_height", "Violin Plot height (px):",
                                      min = 400, max = 1200, value = 600, step = 25),
                          actionButton("ts_show_plot", "Show Plot", icon = icon("eye"), class = "btn-primary", style = "width: 100%; margin-top: 15px; margin-bottom: 5px;")
                          )
             ),
             mainPanel(width = 9,
                       h3(icon("chart-bar"), " Expression Over Time"),
                       conditionalPanel(
                         condition = "output.ts_timepoint_available != true",
                         div(class = "error-box", style = "margin-top: 15px; font-size: 15px;",
                             icon("circle-info"),
                             tags$b(" Time Series analysis is not available for this dataset."),
                             tags$p(style = "margin: 8px 0 0 0; font-weight: normal;",
                                    "No metadata column with usable timepoint values was detected ",
                                    "(the temporal modules need ordered timepoints such as numeric ",
                                    "times, Day_N / Hour_N labels, or a baseline-acute-chronic ladder). ",
                                    "The controls have been disabled. Every other tab - Dataset ",
                                    "Exploration, Differential Expression, GO Enrichment and Gene ",
                                    "Function (COCOA) - works normally."))
                       ),
                       conditionalPanel(
                         condition = "output.ts_timepoint_available == true",
                       uiOutput("ts_summary_table_ui"),
                       tabsetPanel(
                         tabPanel("Heatmap", icon = icon("th"),
                                  br(),
                                  p(style = "color: #555;", icon("info-circle"), " Select multiple genes (or use all genes). The heatmap shows the mean expression of each gene at each timepoint. Genes are grouped by temporal expression pattern using k-means clustering on the z-scored shape of each gene's trajectory (cluster count = number of timepoints)."),
                                  uiOutput("ts_cluster_summary_ui"),
                                  uiOutput("ts_heatmap_plot_ui"),
                                  br(),
                                  uiOutput("ts_cluster_export_ui"),
                                  uiOutput("ts_heatmap_download_ui")
                         ),
                         tabPanel("Temporal Trend", icon = icon("chart-line"),
                                  br(),
                                  p(style = "color: #555;", icon("info-circle"), " This view shows the mean expression trajectory across ordered timepoints with standard error bars. Select a single gene."),
                                  uiOutput("ts_trend_plot_ui"),
                                  br(),
                                  uiOutput("ts_trend_download_ui")
                         ),
                         tabPanel("Boxplot", icon = icon("square"),
                                  br(),
                                  p(style = "color: #555;", icon("info-circle"), " Per-timepoint expression distribution as boxplots. Box = IQR, line = median, whiskers = 1.5x IQR. Select a single gene."),
                                  uiOutput("ts_boxplot_ui"),
                                  br(),
                                  uiOutput("ts_boxplot_download_ui")
                         ),
                         tabPanel("Violin Plot", icon = icon("chart-bar"),
                                  br(),
                                  p(style = "color: #555;", icon("info-circle"), " Single-gene expression distribution across timepoints. Select a single gene."),
                                  uiOutput("ts_plot_ui"),
                                  br(),
                                  uiOutput("ts_download_ui")
                         )
                       )
                       )
             )
           )
  ),
  
  # === HISTORY ===
  tabPanel("History", icon = icon("clock-rotate-left"),
           div(class = "title-panel", h1("Analysis History"),
               h4("Click an item to load results immediately")),
           wellPanel(actionButton("refresh_history", "Refresh History",
                                  icon = icon("rotate"), class = "btn-primary"),
                     hr(), uiOutput("history_list"),
                     hr(),
                     p(style = "color: #555; font-size: 13px;", icon("info-circle"),
                       " Tip: after restoring a saved analysis the relevant tab",
                       " (DEA / COCOA / GO) hydrates with the cached results.",
                       " From there, click ", strong("Download R script"),
                       " to export a self-contained script that reproduces it",
                       " from the same source dataset."))
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  vals <- reactiveValues(
    data = default_load$data, dataset_hash = default_load$hash, umap_calculated = !is.null(default_load$data),
    global_plot_data = NULL, # Cached DF for plotting (Embeddings + Metadata)
    analysis_res = NULL, analysis_meta = NULL, active_filters = list(), error_msg = NULL, is_analyzing = FALSE,
    cocoa_gene_used = NULL, cocoa_filters_used = NULL, 
    dea_results = NULL, dea_error_msg = NULL, dea_filter_list = list(), dea_is_analyzing = FALSE,
    dea_filters_used = NULL, dea_meta_used = NULL, # Add these for plot isolation
    # Custom gene list uploaded into the DEA tab. When non-NULL the volcano
    # plot, results table, and downloads are restricted to these genes. The
    # Time Series tab's per-cluster CSV export is the canonical producer.
    dea_custom_genes = NULL, dea_custom_genes_filename = NULL,
    explore_mask = NULL, explore_metadata_filter = list(), ts_timepoint_available = FALSE,
    explore_active_subtitle = "", # Track 'locked' subtitle for plot
    dea_history_pending = NULL, # Used to store pending group selections from history load
    cocoa_history_pending = NULL, # Used to store pending group selections from history load
    ts_current_choices = NULL, analysis_current_choices = NULL, explore_current_choices = NULL, # For select/deselect all logic
    ts_active_condition = NULL, ts_active_celltype = NULL, ts_active_tps = NULL, ts_active_gene = NULL, ts_active_dataset = NULL, ts_active_all_genes = FALSE, ts_active_roles = NULL,
    # Dataset Exploration new subtabs
    dotplot_active_genes = NULL, dotplot_plot_cached = NULL, dotplot_data_cached = NULL,
    qc_active = FALSE,
    # Time Series gene picker exclusion bookkeeping
    ts_gene_prev = NULL,
    # Session-level cache of msigdbr Hallmark gene sets keyed by species so
    # the COCOA worker does not re-query msigdbr on every Run-Analysis click.
    hallmark_cache = NULL,
    # GO Enrichment tab state
    go_results = NULL, go_error_msg = NULL, go_is_analyzing = FALSE,
    go_upload_data = NULL,
    # geneCOCOA Ensembl-to-symbol conversion state
    cocoa_is_converting = FALSE,
    # Optional landing-page content from a JSON config (manuscript §2).
    # NULL until a JSON loader is implemented; the Introduction tab shows a
    # placeholder while this is NULL.
    landing_config = NULL
  )
  
  output$dataset_status_ui <- renderUI({
    if (!is.null(vals$data)) {
      div(style="color: #27ae60; font-weight: bold; margin-bottom: 20px;",
          icon("check-circle"), " Dataset Ready")
    } else if (is.null(DEFAULT_DATA_PATH)) {
      # No dataset path was provided by the caller. Tell the user how to
      # supply one rather than dressing it up as a load failure.
      div(style = "margin-bottom: 20px; padding: 14px 16px; background-color: #fff7e6; border-left: 4px solid #f39c12; border-radius: 3px;",
          div(style = "color: #b9770e; font-weight: bold; font-size: 1.05em;",
              icon("exclamation-triangle"), " No dataset configured"),
          p(style = "margin: 8px 0 0; color: #2c3e50;",
            "AtlasLens reads the path to your Seurat object from the ",
            tags$code("DATASET_PATH"), " environment variable."),
          tags$ul(style = "color: #2c3e50; margin-bottom: 0;",
                  tags$li(tags$b("Docker:"), " pass ", tags$code("-e DATASET_PATH=/data/file.rds"), " on ", tags$code("docker run"), "."),
                  tags$li(tags$b("Local R:"), " ", tags$code("export DATASET_PATH=/path/to/file.rds"), " before launching."),
                  tags$li("Restart AtlasLens after setting the variable.")))
    } else {
      div(style = "margin-bottom: 20px; padding: 14px 16px; background-color: #fdecea; border-left: 4px solid #e74c3c; border-radius: 3px;",
          div(style = "color: #c0392b; font-weight: bold; font-size: 1.05em;",
              icon("times-circle"), " Dataset failed to load"),
          p(style = "margin: 8px 0 0; color: #2c3e50;",
            "AtlasLens attempted to load the dataset at ",
            tags$code(DEFAULT_DATA_PATH),
            " but could not read it. Confirm the file exists, is a valid ",
            ".rds Seurat object, and that the process has read permission."))
    }
  })
  
  # Landing-page descriptive content. In the published design this is
  # populated from the optional JSON configuration file (manuscript §2:
  # personalized "introduction" and "dataset information"). vals$landing_config
  # stays NULL until that JSON loader is implemented, so a clearly-marked
  # placeholder renders here in the meantime. When the loader lands, populate
  # vals$landing_config$introduction / $dataset_information and this output
  # will render them automatically.
  output$landing_content <- renderUI({
    cfg <- vals$landing_config
    has_intro   <- !is.null(cfg) && !is.null(cfg$introduction)   && nzchar(cfg$introduction)
    has_dataset <- !is.null(cfg) && !is.null(cfg$dataset_information) && nzchar(cfg$dataset_information)
    if (has_intro || has_dataset) {
      tagList(
        if (has_intro) div(
          h3(icon("circle-info"), " Introduction",
             style = "border-bottom: 2px solid #667eea; padding-bottom: 10px;"),
          p(cfg$introduction)),
        if (has_dataset) div(style = "margin-top: 18px;",
                             h3(icon("database"), " Dataset Information",
                                style = "border-bottom: 2px solid #667eea; padding-bottom: 10px;"),
                             p(cfg$dataset_information))
      )
    } else {
      # Placeholder shown until a JSON config supplies landing content.
      div(style = "border: 1px dashed #b0b8c4; border-radius: 6px; padding: 18px 20px; background-color: #fafbfc;",
          h3(icon("circle-info"), " Introduction & Dataset Information",
             style = "border-bottom: 2px solid #667eea; padding-bottom: 10px; color: #2c3e50;"),
          p(style = "color: #7f8c8d; font-style: italic; margin-bottom: 6px;",
            "Placeholder for landing content."),
          p(style = "color: #555;",
            ""))
    }
  })
  
  observe({ 
    req(vals$data); 
    
    # Pre-calculate global plot data for speed
    # Prioritize 'umap' for visualization coordinates. 'scVI' is the latent representation.
    red_to_use <- if ("umap" %in% names(vals$data@reductions)) "umap" else if ("scVI" %in% names(vals$data@reductions)) "scVI" else "pca"
    emb <- Embeddings(vals$data, red_to_use)
    if (ncol(emb) > 2) emb <- emb[, 1:2]
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    vals$global_plot_data <- cbind(as.data.frame(emb), vals$data@meta.data)
    
    all_genes <- rownames(vals$data); updateSelectizeInput(session, "explore_gene", choices = all_genes, selected = "", server = TRUE); updateSelectizeInput(session, "analysis_gene", choices = all_genes, selected = "", server = TRUE); updateSelectizeInput(session, "dea_gene", choices = all_genes, selected = character(0), server = TRUE); valid_cols <- get_valid_metadata_columns(vals$data, show_hidden = isTRUE(input$show_hidden_metadata));
    dataset_col_explore <- if("Celltype_annotated" %in% valid_cols) "Celltype_annotated" else if("Condition" %in% valid_cols) "Condition" else valid_cols[1];
    dataset_col_others <- if("Condition" %in% valid_cols) "Condition" else valid_cols[1];
    updateSelectInput(session, "explore_meta_col", choices = valid_cols, selected = dataset_col_explore); updateSelectInput(session, "explore_filter_meta_col", choices = valid_cols, selected = dataset_col_explore); updateSelectInput(session, "analysis_meta_col", choices = valid_cols, selected = dataset_col_others); updateSelectInput(session, "dea_meta_col", choices = valid_cols, selected = dataset_col_others); updateSelectInput(session, "dea_replicate_col", choices = c("", valid_cols))
    
    # ---- New Dataset Exploration subtab inputs ----
    updateSelectizeInput(session, "coexp_gene_a", choices = all_genes, selected = "",            server = TRUE)
    updateSelectizeInput(session, "coexp_gene_b", choices = all_genes, selected = "",            server = TRUE)
    updateSelectizeInput(session, "dotplot_genes", choices = all_genes, selected = character(0), server = TRUE)
    celltype_default <- if ("Celltype_annotated" %in% valid_cols) "Celltype_annotated"
    else if (any(grepl("celltype|cell_type", valid_cols, ignore.case = TRUE)))
      valid_cols[grep("celltype|cell_type", valid_cols, ignore.case = TRUE)[1]]
    else valid_cols[1]
    updateSelectInput(session, "dotplot_group_col", choices = valid_cols, selected = celltype_default)
    updateSelectInput(session, "qc_group_col",
                      choices = c("(no grouping)" = "", valid_cols),
                      selected = celltype_default)
    
    # ---- Time Series column mapping ----
    # Auto-detected roles come from detect_role_columns(). Dataset is the
    # only role that is optional; the user can leave it as "(none)" if the
    # object contains only a single dataset / batch. Choices for each role
    # are the FULL metadata column set (not just valid_cols) because some
    # bona-fide role columns (e.g. integer timepoints) are filtered out by
    # get_valid_metadata_columns().
    all_meta_cols <- colnames(vals$data@meta.data)
    roles <- detect_role_columns(vals$data@meta.data)
    # Gate the whole Time Series tab on whether a genuine timepoint column was
    # found. When none exists, the tab shows an explanatory message and its
    # controls are disabled rather than silently analysing an arbitrary column.
    vals$ts_timepoint_available <- !is.null(roles$timepoint)
    role_choices_required <- all_meta_cols
    role_choices_optional <- c("(none)" = "", all_meta_cols)
    updateSelectInput(session, "ts_col_timepoint",
                      choices = role_choices_required,
                      selected = roles$timepoint %||% all_meta_cols[1])
    updateSelectInput(session, "ts_col_condition",
                      choices = role_choices_required,
                      selected = roles$condition %||% all_meta_cols[1])
    updateSelectInput(session, "ts_col_celltype",
                      choices = role_choices_required,
                      selected = roles$celltype  %||% all_meta_cols[1])
    updateSelectInput(session, "ts_col_dataset",
                      choices = role_choices_optional,
                      selected = roles$dataset   %||% "")
  })
  observeEvent(input$go_explore, { updateNavbarPage(session, inputId = "main_nav", selected = "Dataset Exploration") })
  output$dataset_info <- renderUI({
    req(vals$data)
    show_hidden <- isTRUE(input$show_hidden_metadata)
    visible_cols <- get_valid_metadata_columns(vals$data, show_hidden = show_hidden)
    hidden <- get_hidden_metadata_columns(vals$data)
    n_hidden_total <- length(hidden$technical) + length(hidden$single_valued) + length(hidden$non_atomic)
    
    # Compact "hidden columns" summary so reviewers see exactly what is
    # filtered out. Each category is rendered only if non-empty.
    fmt_section <- function(label, cols, tooltip) {
      if (length(cols) == 0) return(NULL)
      div(style = "margin: 4px 0;",
          tags$span(style = "font-weight: 600; color: #2c3e50;", paste0(label, " (", length(cols), "): ")),
          tags$span(style = "color: #7f8c8d; font-style: italic;",
                    title = tooltip, paste(cols, collapse = ", ")))
    }
    
    div(
      p(strong("Cells:"), format(ncol(vals$data), big.mark = ",")),
      p(strong("Genes:"), format(nrow(vals$data), big.mark = ",")),
      p(strong("Metadata columns shown:"),
        paste0(length(visible_cols),
               if (show_hidden && n_hidden_total > 0) paste0(" (all ", length(visible_cols), " atomic columns - filter disabled)")
               else if (n_hidden_total > 0) paste0(" of ", length(visible_cols) + n_hidden_total)
               else "")),
      p(style = "margin-bottom: 4px;", strong("Visible columns: "),
        if (length(visible_cols) == 0) tags$em("(none)") else paste(visible_cols, collapse = ", ")),
      if (!show_hidden && n_hidden_total > 0) {
        div(style = "margin-top: 10px; padding: 10px; background-color: #fdf6e3; border-left: 3px solid #f39c12; border-radius: 3px;",
            div(style = "font-weight: bold; color: #b9770e; margin-bottom: 6px;",
                icon("eye-slash"), sprintf(" %d metadata column(s) are hidden by default", n_hidden_total)),
            fmt_section("Technical / QC", hidden$technical,
                        "Per-cell barcode, raw counts, joinid - never useful as grouping factors."),
            fmt_section("Single-valued", hidden$single_valued,
                        "These columns take exactly one value in this object - no comparison is possible."),
            fmt_section("Non-atomic", hidden$non_atomic,
                        "List-columns / S4 objects that the UI cannot render."))
      },
      div(style = "margin-top: 10px;",
          checkboxInput("show_hidden_metadata",
                        "Show hidden metadata columns in all dropdowns",
                        value = show_hidden))
    )
  })
  
  # When the user flips the "show hidden" toggle (without reloading the
  # dataset) re-emit the metadata choice lists for every selectInput that
  # the on-load observer at line ~1656 wires up. The renderUI-based filter
  # dropdowns refresh themselves automatically because they read the same
  # input reactively, but updateSelectInput-driven controls do not.
  observeEvent(input$show_hidden_metadata, {
    req(vals$data)
    valid_cols <- get_valid_metadata_columns(vals$data,
                                             show_hidden = isTRUE(input$show_hidden_metadata))
    keep_or_first <- function(cur) if (length(cur) && cur %in% valid_cols) cur else valid_cols[1]
    updateSelectInput(session, "explore_meta_col",
                      choices = valid_cols, selected = keep_or_first(input$explore_meta_col))
    updateSelectInput(session, "explore_filter_meta_col",
                      choices = valid_cols, selected = keep_or_first(input$explore_filter_meta_col))
    updateSelectInput(session, "analysis_meta_col",
                      choices = valid_cols, selected = keep_or_first(input$analysis_meta_col))
    updateSelectInput(session, "dea_meta_col",
                      choices = valid_cols, selected = keep_or_first(input$dea_meta_col))
    updateSelectInput(session, "dea_replicate_col",
                      choices = c("", valid_cols),
                      selected = if (!is.null(input$dea_replicate_col) && input$dea_replicate_col %in% valid_cols) input$dea_replicate_col else "")
    # Dotplot / QC grouping selectors live under the Exploration subtabs and
    # also expect a valid metadata column.
    updateSelectInput(session, "dotplot_group_col",
                      choices = valid_cols, selected = keep_or_first(input$dotplot_group_col))
    updateSelectInput(session, "qc_group_col",
                      choices = c("(no grouping)" = "", valid_cols),
                      selected = if (!is.null(input$qc_group_col) && input$qc_group_col %in% valid_cols) input$qc_group_col else "")
  }, ignoreInit = TRUE)
  
  render_gene_info <- function(gene, data) {
    if (gene == "" || is.null(gene) || !(gene %in% rownames(data))) return(NULL)
    expr <- GetAssayData(data, layer = "data")[gene, ]
    n_cells <- sum(expr > 0)
    div(class = "gene-info-box", p(icon("check-circle", style = "color: green;"), strong("Selected")), p(style = "margin: 3px 0;", paste0("• Expressing Cells: ", format(n_cells, big.mark=",")))) 
  }
  
  output$explore_gene_info <- renderUI({ req(input$explore_gene); render_gene_info(input$explore_gene, vals$data) })
  output$analysis_gene_info <- renderUI({ req(input$analysis_gene); render_gene_info(input$analysis_gene, vals$data) })
  output$dea_gene_info <- renderUI({ req(input$dea_gene); render_gene_info(input$dea_gene, vals$data) })
  
  output$explore_filter_meta_ui <- renderUI({ req(vals$data, input$explore_filter_meta_col); col <- input$explore_filter_meta_col; groups <- sort(unique(as.character(vals$data@meta.data[[col]]))); vals$explore_current_choices <- groups; div(div(style = "margin-bottom: 5px; text-align: right;", actionLink("toggle_explore_filter", "Select / Deselect All")), selectizeInput("explore_filter_meta_groups", NULL, choices = groups, selected = groups, multiple = TRUE, width = "100%", options = list(plugins = list("remove_button")))) })
  observeEvent(input$toggle_explore_filter, {
    req(vals$explore_current_choices)
    if (length(input$explore_filter_meta_groups) == length(vals$explore_current_choices)) {
      updateSelectizeInput(session, "explore_filter_meta_groups", selected = character(0))
    } else {
      updateSelectizeInput(session, "explore_filter_meta_groups", selected = vals$explore_current_choices)
    }
  })
  
  output$explore_filter_controls_ui <- renderUI({ req(vals$data); div(fluidRow(column(6, selectInput("explore_filter_col_select", NULL, choices = get_valid_metadata_columns(vals$data, show_hidden = isTRUE(input$show_hidden_metadata)), width = "100%")), column(6, uiOutput("explore_filter_vals_select_ui"))), actionButton("explore_add_filter", "Add Filter", icon = icon("plus"), class = "btn-xs btn-info", style = "width: 100%; margin-top: 5px;")) })
  output$explore_filter_vals_select_ui <- renderUI({ req(input$explore_filter_col_select, vals$data); selectInput("explore_filter_vals_select", NULL, choices = sort(unique(as.character(vals$data@meta.data[[input$explore_filter_col_select]]))), multiple = TRUE, width = "100%") })
  observeEvent(input$explore_add_filter, { col <- input$explore_filter_col_select; val <- input$explore_filter_vals_select; if(is.null(col) || is.null(val)) return(); curr <- vals$explore_metadata_filter; if (!is.list(curr)) curr <- list(); exist <- which(vapply(curr, function(x) x$col == col, logical(1))); if(length(exist)>0) curr[[exist]]$vals <- val else curr[[length(curr)+1]] <- list(col=col, vals=val); vals$explore_metadata_filter <- curr })
  observeEvent(input$explore_remove_filter, { idx <- as.numeric(input$explore_remove_filter); if(length(vals$explore_metadata_filter)>=idx) vals$explore_metadata_filter[[idx]] <- NULL })
  
  output$explore_active_filters_ui <- renderUI({ if(length(vals$explore_metadata_filter)==0) return(NULL); div(class = "active-filters", lapply(seq_along(vals$explore_metadata_filter), function(i) { f <- vals$explore_metadata_filter[[i]]; span(span(class="filter-badge", strong(f$col), ": ", paste(head(f$vals, 2), collapse=","), if(length(f$vals)>2)"..."), actionLink("explore_remove_filter", "✕", onclick = paste0("Shiny.onInputChange('explore_remove_filter', ", i, ");"), style="color: red; margin-left: 5px;")) })) })
  
  # Shared mask. Recomputes automatically whenever the filter list, the
  # Show/Hide Groups checkboxes, or the Color-by column change, so the
  # plots stay in sync without an explicit Update View button.
  observe({
    req(vals$data)
    
    cell_mask  <- rep(TRUE, ncol(vals$data))
    has_filter <- FALSE
    
    # Apply the Show/Hide Groups checkbox set, if present (Metadata UMAP subtab).
    if (!is.null(input$explore_filter_meta_groups) && !is.null(input$explore_filter_meta_col)) {
      all_groups <- sort(unique(as.character(vals$data@meta.data[[input$explore_filter_meta_col]])))
      if (length(input$explore_filter_meta_groups) < length(all_groups)) {
        cell_mask  <- cell_mask & (vals$data@meta.data[[input$explore_filter_meta_col]] %in% input$explore_filter_meta_groups)
        has_filter <- TRUE
      }
    }
    
    # Apply the Add-Filter row list.
    if (length(vals$explore_metadata_filter) > 0) {
      for (f in vals$explore_metadata_filter) {
        cell_mask <- cell_mask & (vals$data@meta.data[[f$col]] %in% f$vals)
      }
      has_filter <- TRUE
    }
    
    if (!has_filter) {
      vals$explore_mask <- NULL
      vals$explore_active_subtitle <- ""
    } else if (sum(cell_mask) > 0) {
      vals$explore_mask <- cell_mask
      vals$explore_active_subtitle <- paste("Filtered:", format(sum(cell_mask), big.mark=","), "cells")
    } else {
      vals$explore_mask <- NULL
      vals$explore_active_subtitle <- ""
      # Not an error: no cell matches the chosen combination, so we simply
      # fall back to the full view. A brief, self-dismissing warning explains
      # what happened without leaving a sticky banner on screen.
      showNotification(
        "No cells match the current filter combination - showing the full, unfiltered view instead.",
        type = "warning", duration = 6)
    }
  })
  
  # === TIME SERIES TAB LOGIC ===
  # Resolve the user-configured role columns once per call. Returns NULL
  # for any role the dataset (or the user) has not provided.
  ts_roles <- reactive({
    list(
      timepoint = if (isTruthy(input$ts_col_timepoint)) input$ts_col_timepoint else NULL,
      condition = if (isTruthy(input$ts_col_condition)) input$ts_col_condition else NULL,
      celltype  = if (isTruthy(input$ts_col_celltype))  input$ts_col_celltype  else NULL,
      dataset   = if (isTruthy(input$ts_col_dataset))   input$ts_col_dataset   else NULL
    )
  })

  # Expose timepoint availability to the UI so the Time Series main panel can
  # swap between the analysis tabs and an explanatory "not available" message
  # via conditionalPanel. suspendWhenHidden = FALSE keeps it evaluated even
  # while the message branch (which itself references it) is the visible one.
  output$ts_timepoint_available <- reactive({ isTRUE(vals$ts_timepoint_available) })
  outputOptions(output, "ts_timepoint_available", suspendWhenHidden = FALSE)

  # Freeze the Time Series sidebar controls when no timepoint column exists.
  observe({
    if (isTRUE(vals$ts_timepoint_available)) {
      shinyjs::removeClass(id = "ts_config_wrap", class = "panel-disabled")
    } else {
      shinyjs::addClass(id = "ts_config_wrap", class = "panel-disabled")
    }
  })
  
  # Dataset selector. If the user mapped a dataset column, render a
  # selectInput listing the unique values that have at least one cell with
  # a usable timepoint label. If no dataset column was mapped, render an
  # empty placeholder so downstream observers can fall through to the
  # all-cells path.
  output$ts_dataset_ui <- renderUI({
    req(vals$data)
    r <- ts_roles()
    if (is.null(r$dataset)) {
      return(tagList(
        tags$small(style = "color: #666;",
                   "No dataset column mapped — treating all cells as one dataset."),
        # Hidden selectInput so input$ts_dataset still resolves to "" downstream.
        tags$div(style = "display:none;",
                 selectInput("ts_dataset", NULL, choices = c("(all cells)" = ""), selected = ""))
      ))
    }
    meta <- vals$data@meta.data
    if (!r$dataset %in% colnames(meta)) return(NULL)
    if (!is.null(r$timepoint) && r$timepoint %in% colnames(meta)) {
      keep <- ts_has_timepoint_mask(meta[[r$timepoint]])
      ts_datasets <- sort(unique(as.character(meta[[r$dataset]][keep])))
    } else {
      ts_datasets <- sort(unique(as.character(meta[[r$dataset]])))
    }
    tagList(
      tags$label("Dataset:", style = "font-weight: bold;"),
      selectInput("ts_dataset", NULL, choices = ts_datasets)
    )
  })
  
  observeEvent(vals$data, {
    req(vals$data)
    all_genes <- rownames(vals$data)
    # The heatmap gene picker offers a sentinel "All genes (heatmap)" entry
    # that is mutually exclusive with any individual gene selection. The
    # sentinel is the default so the heatmap renders the full transcriptome
    # out of the box; picking any concrete gene clears the sentinel.
    ts_choices <- c("All genes (heatmap)" = "__ALL_GENES__", setNames(all_genes, all_genes))
    updateSelectizeInput(session, "ts_gene", choices = ts_choices,
                         selected = "__ALL_GENES__", server = TRUE)
    vals$ts_gene_prev <- "__ALL_GENES__"
  })
  
  # ts_gene mutual exclusion. The previous selection is cached in
  # vals$ts_gene_prev; the diff identifies which side of the conflict was
  # just added so the other can be dropped.
  observeEvent(input$ts_gene, {
    sel <- input$ts_gene
    if (is.null(sel)) sel <- character(0)
    has_sentinel <- "__ALL_GENES__" %in% sel
    others       <- setdiff(sel, "__ALL_GENES__")
    if (has_sentinel && length(others) > 0) {
      prev <- vals$ts_gene_prev
      prev_had_sentinel <- !is.null(prev) && "__ALL_GENES__" %in% prev
      if (prev_had_sentinel) {
        # Sentinel was already there. The latest pick is a concrete gene -> drop sentinel.
        updateSelectizeInput(session, "ts_gene", selected = others)
        vals$ts_gene_prev <- others
      } else {
        # Sentinel was just added; concrete genes were already there -> drop concrete genes.
        updateSelectizeInput(session, "ts_gene", selected = "__ALL_GENES__")
        vals$ts_gene_prev <- "__ALL_GENES__"
      }
      return()
    }
    vals$ts_gene_prev <- sel
  }, ignoreNULL = FALSE)
  
  # ts_dataset_filter returns the boolean mask selecting cells in the
  # currently-active dataset. When no dataset column is mapped, every cell
  # is included; this lets the Time Series tab work on objects that lack a
  # dataset / batch column.
  ts_dataset_filter <- function() {
    req(vals$data)
    meta <- vals$data@meta.data
    r <- ts_roles()
    if (is.null(r$dataset) || !isTruthy(input$ts_dataset))
      return(rep(TRUE, nrow(meta)))
    if (!r$dataset %in% colnames(meta)) return(rep(TRUE, nrow(meta)))
    meta[[r$dataset]] == input$ts_dataset
  }
  
  # Dataset -> Condition cascade
  # The detected control condition is removed from the dropdown because it
  # is automatically included in every Time Series plot. This way the user
  # always picks a "treatment" condition and never has to remember to add
  # the control manually.
  observe({
    req(vals$data)
    r <- ts_roles()
    if (is.null(r$condition)) return()
    meta <- vals$data@meta.data
    if (!r$condition %in% colnames(meta)) return()
    ds_mask <- ts_dataset_filter()
    all_conditions <- sort(unique(as.character(meta[[r$condition]][ds_mask])))
    ctrl <- detect_control_conditions(all_conditions)
    # Remove control from dropdown — it's auto-included on every plot
    display_conditions <- setdiff(all_conditions, ctrl)
    if (length(display_conditions) == 0) display_conditions <- all_conditions
    updateSelectInput(session, "ts_condition", choices = display_conditions)
  })
  
  # Badge showing the auto-detected control condition
  output$ts_control_badge_ui <- renderUI({
    req(vals$data)
    r <- ts_roles()
    if (is.null(r$condition)) return(NULL)
    meta <- vals$data@meta.data
    if (!r$condition %in% colnames(meta)) return(NULL)
    ds_mask <- ts_dataset_filter()
    all_conds <- unique(as.character(meta[[r$condition]][ds_mask]))
    ctrl <- detect_control_conditions(all_conds)
    if (length(ctrl) == 0)
      return(tags$small(style = "color: #999; display: block; margin-bottom: 8px;",
                        icon("info-circle"), " No control condition detected."))
    tags$small(style = "color: #27ae60; display: block; margin-bottom: 8px; font-weight: bold;",
               icon("check-circle"),
               paste0(" Auto-included: ", paste(ctrl, collapse = ", ")))
  })
  
  # Dataset -> Cell-type cascade
  observe({
    req(vals$data)
    r <- ts_roles()
    if (is.null(r$celltype)) return()
    meta <- vals$data@meta.data
    if (!r$celltype %in% colnames(meta)) return()
    ds_mask <- ts_dataset_filter()
    celltypes <- sort(unique(as.character(meta[[r$celltype]][ds_mask])))
    updateSelectInput(session, "ts_celltype", choices = celltypes)
  })
  
  # Dataset info badge. Shows cell count, conditions, the auto-detected
  # control (if any), and the timepoint set. Falls back gracefully on
  # missing role mappings.
  output$ts_dataset_info_ui <- renderUI({
    req(vals$data)
    r <- ts_roles()
    meta <- vals$data@meta.data
    ds_mask <- ts_dataset_filter()
    ds_meta <- meta[ds_mask, , drop = FALSE]
    tps <- if (!is.null(r$timepoint) && r$timepoint %in% colnames(meta)) {
      raw_tps <- as.character(ds_meta[[r$timepoint]])
      raw_tps <- raw_tps[ts_has_timepoint_mask(raw_tps)]
      raw_tps <- unique(raw_tps)
      # numeric sort if possible
      num_tps <- suppressWarnings(as.numeric(raw_tps))
      if (!anyNA(num_tps)) raw_tps[order(num_tps)] else sort(raw_tps)
    } else character(0)
    conds <- if (!is.null(r$condition) && r$condition %in% colnames(meta))
      sort(unique(as.character(ds_meta[[r$condition]]))) else character(0)
    ctrl <- detect_control_conditions(conds)
    ctrl_str <- if (length(ctrl) > 0) ctrl[1] else "None detected"
    div(class = "gene-info-box",
        p(strong("Cells: "), format(nrow(ds_meta), big.mark = ",")),
        if (length(conds) > 0) p(strong("Conditions: "), paste(conds, collapse = ", ")),
        p(strong("Control: "), ctrl_str),
        if (length(tps) > 0) p(strong("Timepoints: "), paste(tps, collapse = ", "))
    )
  })
  
  output$ts_timepoints_ui <- renderUI({
    req(vals$data, input$ts_condition, input$ts_celltype)
    r <- ts_roles()
    if (is.null(r$timepoint))
      return(p(style = "color: red;",
               "No timepoint column mapped. Set one under 'Column mapping'."))
    meta <- vals$data@meta.data
    if (!r$timepoint %in% colnames(meta))
      return(p(style = "color: red;",
               paste("Timepoint column '", r$timepoint, "' not found.")))
    
    # Include timepoints from BOTH the selected condition AND the auto-detected
    # control. Without this, control-only timepoints (e.g. timepoint 0) would
    # never appear in the checkbox list when the user picks a treatment condition.
    mask <- ts_dataset_filter()
    if (!is.null(r$condition) && r$condition %in% colnames(meta)) {
      all_conds_ds <- unique(as.character(meta[[r$condition]][mask]))
      ctrl <- detect_control_conditions(all_conds_ds)
      conditions_for_tps <- unique(c(input$ts_condition, ctrl))
      mask <- mask & meta[[r$condition]] %in% conditions_for_tps
    }
    if (!is.null(r$celltype) && r$celltype %in% colnames(meta))
      mask <- mask & meta[[r$celltype]]  == input$ts_celltype
    
    if (sum(mask) == 0)
      return(p(style = "color: #e67e22;", "No cells found for this combination."))
    
    available_tps <- unique(as.character(meta[[r$timepoint]][mask]))
    available_tps <- available_tps[ts_has_timepoint_mask(available_tps)]
    if (length(available_tps) == 0)
      return(p("No time series information available for these cells."))
    # Chronological sort: numeric, embedded-number, and ordinal labels.
    available_tps <- sort_timepoints(available_tps)
    
    vals$ts_current_choices <- available_tps
    div(
      div(style = "margin-bottom: 5px; text-align: right;", actionLink("toggle_ts_tps", "Select / Deselect All")),
      checkboxGroupInput("ts_selected_tps", NULL, choices = available_tps, selected = available_tps)
    )
  })
  
  observeEvent(input$toggle_ts_tps, {
    req(vals$ts_current_choices)
    if (length(input$ts_selected_tps) == length(vals$ts_current_choices)) {
      updateCheckboxGroupInput(session, "ts_selected_tps", selected = character(0))
    } else {
      updateCheckboxGroupInput(session, "ts_selected_tps", selected = vals$ts_current_choices)
    }
  })
  
  # === TIME SERIES SHOW PLOT OBSERVER ===
  observeEvent(input$ts_show_plot, {
    vals$ts_active_condition <- input$ts_condition
    vals$ts_active_celltype  <- input$ts_celltype
    # Numeric sort so timepoints like 7, 14, 28, 105 order naturally
    # instead of lexically.
    tps <- input$ts_selected_tps
    num_tps <- suppressWarnings(as.numeric(tps))
    vals$ts_active_tps <- if (!anyNA(num_tps)) tps[order(num_tps)] else sort(tps)
    # Translate the "All genes (heatmap)" sentinel into all-genes mode and
    # strip it from the active gene list so downstream code never sees it.
    sel <- input$ts_gene
    vals$ts_active_all_genes <- "__ALL_GENES__" %in% sel
    vals$ts_active_gene <- setdiff(sel, "__ALL_GENES__")
    vals$ts_active_dataset <- input$ts_dataset
    # Snapshot the column mapping so plot reactives stay consistent even if
    # the user changes the mapping after pressing Show Plot.
    vals$ts_active_roles <- ts_roles()
  })
  
  # === SUMMARY TABLE ===
  output$ts_summary_table_ui <- renderUI({
    req(vals$data)
    r <- ts_roles()
    if (is.null(r$timepoint) || is.null(r$condition)) return(NULL)
    meta <- vals$data@meta.data
    if (!(r$timepoint %in% colnames(meta)) || !(r$condition %in% colnames(meta))) return(NULL)
    mask <- ts_dataset_filter() & ts_has_timepoint_mask(meta[[r$timepoint]])
    if (sum(mask) == 0) return(NULL)
    ds_meta <- meta[mask, , drop = FALSE]
    summary_df <- as.data.frame(table(
      Condition = ds_meta[[r$condition]],
      Timepoint = ds_meta[[r$timepoint]]
    ))
    summary_df <- summary_df[summary_df$Freq > 0, ]
    num_tp <- suppressWarnings(as.numeric(as.character(summary_df$Timepoint)))
    summary_df <- summary_df[order(summary_df$Condition,
                                   if (!anyNA(num_tp)) num_tp else summary_df$Timepoint), ]
    colnames(summary_df) <- c(r$condition, r$timepoint, "Cell Count")
    
    title_suffix <- if (isTruthy(input$ts_dataset)) paste("Dataset:", input$ts_dataset)
    else "All cells"
    wellPanel(style = "background-color: #f8f9fa; padding: 10px; margin-bottom: 15px;",
              h4(icon("table"), " Summary — ", title_suffix),
              DT::renderDataTable(
                DT::datatable(summary_df,
                              options = list(pageLength = 20, dom = 'tp', scrollX = TRUE),
                              rownames = FALSE, style = "bootstrap")
              )
    )
  })
  
  # === HEATMAP LOGIC ===
  # One reactive serves both the user-typed gene list and the "All genes
  # (heatmap)" sentinel. The branching only changes (a) which gene set
  # enters the computation, (b) the plot height, (c) whether row labels are
  # rendered, and (d) geom_tile vs geom_raster (raster is a single bitmap
  # so it scales to tens of thousands of rows without ggplot emitting one
  # polygon per cell).
  ts_heatmap_logic <- reactive({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_tps)
    r <- vals$ts_active_roles %||% ts_roles()
    req(!is.null(r$timepoint), !is.null(r$condition), !is.null(r$celltype))
    
    use_all <- isTRUE(vals$ts_active_all_genes)
    if (!use_all) {
      req(vals$ts_active_gene)
      req(length(vals$ts_active_gene) >= 2)
    }
    
    meta <- vals$data@meta.data
    
    # Auto-detect the dataset's control condition so it is always rendered
    # alongside the user-selected condition. The control reference is
    # required for time-course interpretation of treatment-induced
    # trajectories.
    ds_mask <- ts_dataset_filter()
    ds_conditions <- unique(as.character(meta[[r$condition]][ds_mask]))
    ctrl <- detect_control_conditions(ds_conditions)
    conditions_to_use <- vals$ts_active_condition
    if (length(ctrl) > 0 && !ctrl[1] %in% conditions_to_use) {
      conditions_to_use <- c(ctrl[1], conditions_to_use)
    }
    
    mask <- ds_mask &
      meta[[r$condition]] %in% conditions_to_use &
      meta[[r$celltype]]  == vals$ts_active_celltype &
      meta[[r$timepoint]] %in% vals$ts_active_tps
    if (sum(mask) == 0) return(NULL)
    
    sub_obj <- vals$data[, mask]
    
    # All-genes mode takes every gene present in the matrix; otherwise the
    # user-typed list, filtered to names that exist in the dataset.
    valid_genes <- if (use_all) {
      rownames(vals$data)
    } else {
      vals$ts_active_gene[vals$ts_active_gene %in% rownames(vals$data)]
    }
    if (length(valid_genes) < 2) return(NULL)
    
    expr_mat  <- GetAssayData(sub_obj, layer = "data")[valid_genes, , drop = FALSE]
    tp_labels <- as.character(sub_obj@meta.data[[r$timepoint]])
    
    # Mean expression per gene per timepoint, timepoints in chronological order.
    sorted_tps <- sort_timepoints(vals$ts_active_tps)
    mean_mat <- matrix(NA, nrow = length(valid_genes), ncol = length(sorted_tps),
                       dimnames = list(valid_genes, sorted_tps))
    for (tp in sorted_tps) {
      tp_idx <- which(tp_labels == tp)
      if (length(tp_idx) > 0) {
        mean_mat[, tp] <- Matrix::rowMeans(expr_mat[, tp_idx, drop = FALSE])
      }
    }
    
    # Drop rows that are constant or all-NA — there's nothing to learn from
    # them and they'd waste height. This also matters in all-genes mode where
    # most genes won't be expressed in this celltype/condition slice.
    row_var <- apply(mean_mat, 1, var, na.rm = TRUE)
    mean_mat <- mean_mat[!is.na(row_var) & row_var > 0, , drop = FALSE]
    if (nrow(mean_mat) < 2) return(NULL)
    
    # The heatmap displays raw MEAN EXPRESSION per gene per timepoint (no
    # z-scoring), so it has no negative values and matches the manuscript.
    # Genes are still grouped by temporal SHAPE: a z-scored + L2-normalised
    # copy is built purely as the k-means input (correlation-distance
    # clustering), independent of absolute expression level.
    if (ncol(mean_mat) >= 2) {
      clust_mat <- t(scale(t(mean_mat)))
      clust_mat[is.na(clust_mat)] <- 0
      # Number of gene clusters tracks the number of timepoints, so the
      # cluster count reflects the time-course design, not the gene count.
      k <- min(ncol(mean_mat), max(2, floor(nrow(mean_mat) / 2)))
      row_norms <- sqrt(rowSums(clust_mat^2))
      row_norms[row_norms == 0] <- 1
      km_input  <- clust_mat / row_norms
      set.seed(42)
      km <- kmeans(km_input, centers = k, nstart = 10, iter.max = 50)
      km_clusters <- km$cluster
    } else {
      k <- 1
      km_clusters <- rep(1, nrow(mean_mat))
      names(km_clusters) <- rownames(mean_mat)
    }
    fill_label   <- "Mean\nexpression"
    caption_str  <- "Mean expression per gene at each timepoint; genes grouped by temporal pattern (k-means)."
    midpoint_val <- mean(mean_mat, na.rm = TRUE)
    
    # Per-cluster stats: `amplitude` (largest mean-expression swing across
    # timepoints) ranks how dynamic a cluster is; `peak_tp` is the timepoint
    # at which the cluster's mean expression is highest, used to label the
    # cluster so the user sees which timepoint each cluster corresponds to.
    cluster_ids   <- sort(unique(km_clusters))
    cluster_traj  <- lapply(cluster_ids, function(cl)
      colMeans(mean_mat[km_clusters == cl, , drop = FALSE], na.rm = TRUE))
    cluster_stats <- data.frame(
      cluster   = cluster_ids,
      n_genes   = as.integer(table(factor(km_clusters, levels = cluster_ids))),
      amplitude = vapply(cluster_traj, function(traj)
        if (all(is.na(traj))) 0 else diff(range(traj, na.rm = TRUE)), numeric(1)),
      peak_tp   = vapply(cluster_traj, function(traj)
        if (all(is.na(traj))) NA_character_ else sorted_tps[which.max(traj)], character(1)),
      stringsAsFactors = FALSE
    )
    cluster_stats$label <- ifelse(
      is.na(cluster_stats$peak_tp),
      paste0("Cluster ", cluster_stats$cluster),
      paste0("Cluster ", cluster_stats$cluster, " (peak: ", cluster_stats$peak_tp, ")"))
    cluster_stats <- cluster_stats[order(-cluster_stats$amplitude), , drop = FALSE]
    
    # Order rows by k-means cluster so the heatmap shows blocks of similar
    # temporal expression patterns instead of arbitrary alphabetical order.
    gene_order    <- order(km_clusters)
    ordered_genes <- rownames(mean_mat)[gene_order]
    
    plot_mat  <- mean_mat[ordered_genes, , drop = FALSE]
    # Map each gene to its cluster label so the heatmap can be split into
    # labelled per-cluster bands showing which timepoint each cluster peaks at.
    gene_label <- stats::setNames(
      cluster_stats$label[match(km_clusters, cluster_stats$cluster)],
      names(km_clusters))
    long_df   <- expand.grid(Gene = ordered_genes, Timepoint = sorted_tps,
                             stringsAsFactors = FALSE)
    long_df$Value     <- as.vector(plot_mat)
    long_df$Cluster   <- factor(gene_label[long_df$Gene],
                                levels = cluster_stats$label[order(cluster_stats$cluster)])
    long_df$Gene      <- factor(long_df$Gene, levels = rev(ordered_genes))
    long_df$Timepoint <- factor(long_df$Timepoint, levels = sorted_tps)
    
    n_rows <- nrow(plot_mat)
    title_str <- if (k > 1) {
      paste0("Gene Expression Heatmap (", n_rows, " genes, k = ", k, " clusters",
             if (use_all) "; all-genes mode" else "", ")")
    } else {
      paste0("Gene Expression Heatmap (", n_rows, " genes",
             if (use_all) "; all-genes mode" else "", ")")
    }
    
    # I switch from geom_tile to geom_raster once there are enough rows that
    # individual tile borders aren't visible anyway. geom_raster renders the
    # whole grid as one bitmap, which is dramatically faster for big heatmaps.
    use_raster <- n_rows > 200
    
    p <- ggplot(long_df, aes(x = Timepoint, y = Gene, fill = Value))
    p <- if (use_raster) p + geom_raster()
    else            p + geom_tile(color = "white", linewidth = 0.3)
    p <- p +
      scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#e08214",
                           midpoint = midpoint_val, name = fill_label) +
      labs(title = title_str,
           subtitle = paste("Condition:", vals$ts_active_condition,
                            "| Cell type:", vals$ts_active_celltype),
           x = paste("Timepoint",
                     if (!is.null(r$timepoint)) paste0("(", r$timepoint, ")") else ""),
           y = NULL, caption = caption_str) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"),
            panel.grid = element_blank())
    
    # Facet the heatmap into one labelled band per cluster, so the user sees
    # directly which cluster corresponds to which peak timepoint.
    if (k > 1) {
      p <- p + facet_grid(rows = vars(Cluster), scales = "free_y", space = "free_y") +
        theme(strip.text.y     = element_text(angle = 0, hjust = 0, size = 9, face = "bold"),
              strip.background = element_rect(fill = "#ecf0f1", color = NA),
              panel.spacing.y  = grid::unit(3, "pt"))
    }
    
    # Hide row labels once the row count exceeds what is legibly displayable.
    if (n_rows > 150) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    } else {
      p <- p + theme(axis.text.y = element_text(size = max(6, 12 - n_rows / 5)))
    }
    list(plot = p, clusters = km_clusters, k = k, cluster_stats = cluster_stats)
  })
  
  # Helper: how tall should the heatmap be? I keep it scaled to row count but
  # capped so it doesn't hit browser limits in all-genes mode.
  ts_heatmap_height_px <- function(n_rows) {
    if (n_rows <= 0) return(400)
    px_per_row <- if (n_rows <= 50) 18 else if (n_rows <= 500) 8 else 0.4
    h <- round(n_rows * px_per_row + 200)
    max(400, min(h, 8000))
  }
  
  output$ts_heatmap_plot_ui <- renderUI({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype)
    if (is.null(vals$ts_active_tps) || length(vals$ts_active_tps) == 0) {
      return(p(style = "color:#e67e22;", icon("info-circle"),
               " Please select timepoints and click Show Plot."))
    }
    use_all <- isTRUE(vals$ts_active_all_genes)
    if (!use_all && (is.null(vals$ts_active_gene) || length(vals$ts_active_gene) < 2)) {
      return(p(style = "color:#e67e22;", icon("info-circle"),
               " Please select at least 2 genes, or tick 'Heatmap: use ALL genes'."))
    }
    n_rows <- if (use_all) nrow(vals$data) else length(vals$ts_active_gene)
    # `ts_heatmap_height_override == 0` selects auto-fit based on the row
    # count; any non-zero value overrides the auto-fit so the user can pin
    # a height that fits their browser window.
    h <- if (isTRUE(input$ts_heatmap_height_override > 0))
      input$ts_heatmap_height_override else ts_heatmap_height_px(n_rows)
    note <- if (use_all) tagList(
      p(style = "color: #e67e22; font-size: 12px;",
        icon("exclamation-triangle"),
        sprintf(" All-genes mode: rendering up to %d genes (after dropping constant rows). This can take 10-30 s.", n_rows))
    ) else NULL
    tagList(
      note,
      withSpinner(plotOutput("ts_heatmap_plot", height = paste0(h, "px")),
                  type = 6, color = "#667eea")
    )
  })
  
  output$ts_heatmap_plot <- renderPlot({ res <- ts_heatmap_logic(); req(res); res$plot })
  
  output$ts_heatmap_download_ui <- renderUI({
    req(vals$ts_active_tps)
    use_all <- isTRUE(vals$ts_active_all_genes)
    if (!use_all) req(vals$ts_active_gene, length(vals$ts_active_gene) >= 2)
    downloadButton("ts_heatmap_download_plot", "Download Heatmap",
                   class = "btn-success")
  })
  
  output$ts_heatmap_download_plot <- downloadHandler(
    filename = function() {
      tag <- if (isTRUE(vals$ts_active_all_genes)) "_allgenes" else ""
      paste0("heatmap_", vals$ts_active_condition, tag, ".png")
    },
    content = function(file) {
      # I size the saved PNG by the same px-per-row scheme as the on-screen
      # plot, then convert to inches at 300 dpi. Capped so the file doesn't
      # blow past ggsave's safety limits in all-genes mode.
      n_rows <- if (isTRUE(vals$ts_active_all_genes)) nrow(vals$data)
      else length(vals$ts_active_gene)
      h_in <- max(6, min(ts_heatmap_height_px(n_rows) / 100, 60))
      hm <- ts_heatmap_logic(); req(hm)
      ggsave(file, plot = hm$plot,
             width = 10, height = h_in, dpi = 300, limitsize = FALSE)
    }
  )
  # --- Interesting-cluster callout (most dynamic gene cluster) ---
  output$ts_cluster_summary_ui <- renderUI({
    res <- ts_heatmap_logic()
    req(res, !is.null(res$cluster_stats), res$k > 1)
    top <- res$cluster_stats[1, ]
    div(class = "success-box", style = "margin-bottom: 12px;",
        h4(icon("star"), " Most dynamic cluster"),
        p(sprintf("%s has the largest mean-expression swing across timepoints (%d genes) and is the most temporally interesting set. Use the selector below the heatmap to download its genes for DEA.",
                  top$label, top$n_genes)))
  })
  
  # --- Per-cluster gene export (pick a cluster, download or send to DEA) ---
  # The CSV download is the portable path; the "Send to DEA" button is the
  # one-click bridge that loads the gene list straight into the DEA tab so
  # the volcano + table are restricted to that cluster's genes immediately.
  output$ts_cluster_export_ui <- renderUI({
    res <- ts_heatmap_logic()
    req(res, !is.null(res$cluster_stats), res$k > 1)
    cl_choices <- stats::setNames(res$cluster_stats$cluster, res$cluster_stats$label)
    div(style = "margin: 10px 0; padding: 12px; background-color: #f8f9fa; border-radius: 4px;",
        tags$label("Export a cluster's genes (download for later, or send straight to DEA):", style = "font-weight: bold;"),
        fluidRow(
          column(4, selectInput("ts_selected_cluster", NULL,
                                choices = cl_choices, selected = cl_choices[1])),
          column(4, downloadButton("ts_cluster_genes_download",
                                   "Download (CSV)", class = "btn-info",
                                   style = "width: 100%;")),
          column(4, actionButton("ts_send_cluster_to_dea",
                                 "Send to DEA",
                                 icon = icon("arrow-right"),
                                 class = "btn-success",
                                 style = "width: 100%;"))
        ),
        helpText("Send to DEA restricts the volcano plot and results table to this cluster's genes. Run DEA first if you have not already."))
  })
  
  output$ts_cluster_genes_download <- downloadHandler(
    filename = function() paste0("cluster_", input$ts_selected_cluster %||% "NA", "_genes.csv"),
    content  = function(file) {
      res <- ts_heatmap_logic(); req(res)
      cl    <- input$ts_selected_cluster
      genes <- names(res$clusters)[as.character(res$clusters) == as.character(cl)]
      utils::write.csv(data.frame(gene = genes), file, row.names = FALSE)
    }
  )
  
  # One-click cluster -> DEA bridge: copies the selected cluster's gene
  # symbols into vals$dea_custom_genes (the same slot the DEA upload writes
  # to) and switches to the DEA tab. No file round-trip required.
  observeEvent(input$ts_send_cluster_to_dea, {
    res <- ts_heatmap_logic(); req(res)
    cl    <- input$ts_selected_cluster
    genes <- names(res$clusters)[as.character(res$clusters) == as.character(cl)]
    if (length(genes) == 0) {
      showNotification("Selected cluster contains no genes.", type = "warning")
      return()
    }
    label <- tryCatch(
      res$cluster_stats$label[res$cluster_stats$cluster == cl][1],
      error = function(e) paste0("cluster ", cl)
    )
    vals$dea_custom_genes <- genes
    vals$dea_custom_genes_filename <- paste0("Time Series ", label %||% paste("cluster", cl))
    updateNavbarPage(session, "main_nav", selected = "DEA")
    showNotification(
      sprintf("Loaded %d genes from %s into DEA. Volcano + table now restricted to this list.",
              length(genes), label %||% paste("cluster", cl)),
      type = "message", duration = 6)
  })
  
  output$ts_plot_ui <- renderUI({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_gene)
    gene1 <- vals$ts_active_gene[1]
    
    if (is.null(vals$ts_active_tps) || length(vals$ts_active_tps) == 0) {
      return(p(style="color:#e67e22;", icon("info-circle"), " Please select at least one timepoint and click Show Plot."))
    }
    
    if(!gene1 %in% rownames(vals$data)) {
      return(p("Selected gene not found in dataset."))
    }
    
    withSpinner(plotOutput("ts_expression_plot",
                           height = paste0(input$ts_violin_height %||% 600, "px")),
                type = 6, color = "#667eea")
  })
  
  output$ts_download_ui <- renderUI({
    req(vals$ts_active_gene, vals$ts_active_tps)
    downloadButton("ts_download_plot", "Download Plot", class = "btn-success")
  })
  
  output$ts_download_plot <- downloadHandler(
    filename = function() { paste0("timeseries_", vals$ts_active_gene, "_", vals$ts_active_condition, ".png") },
    content = function(file) {
      ggsave(file, plot = ts_plot_logic(), width = 12, height = 8, dpi = 300)
    }
  )
  
  output$ts_expression_plot <- renderPlot({
    ts_plot_logic()
  })
  
  ts_plot_logic <- reactive({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_gene, vals$ts_active_tps)
    r <- vals$ts_active_roles %||% ts_roles()
    req(!is.null(r$timepoint), !is.null(r$condition), !is.null(r$celltype))
    
    # Always include the auto-detected control condition alongside the
    # user-selected condition so the treatment trajectory has a reference.
    meta <- vals$data@meta.data
    ds_mask <- ts_dataset_filter()
    ds_conds <- unique(as.character(meta[[r$condition]][ds_mask]))
    ctrl <- detect_control_conditions(ds_conds)
    conditions_use <- unique(c(vals$ts_active_condition, ctrl))
    mask <- ds_mask &
      meta[[r$condition]] %in% conditions_use &
      meta[[r$celltype]]  == vals$ts_active_celltype &
      meta[[r$timepoint]] %in% vals$ts_active_tps
    if (sum(mask) == 0) return(NULL)
    
    sub_obj <- vals$data[, mask]
    gene1 <- vals$ts_active_gene[1]
    expr  <- GetAssayData(sub_obj, layer = "data")[gene1, ]
    cond_vec <- as.character(sub_obj@meta.data[[r$condition]])
    df <- data.frame(
      Expression = as.numeric(expr),
      Timepoint  = factor(sub_obj@meta.data[[r$timepoint]], levels = vals$ts_active_tps),
      Condition  = factor(cond_vec, levels = conditions_use)
    )
    
    n_conds <- length(levels(df$Condition))
    cond_palette <- if (n_conds <= 2) c("#3498db", "#e74c3c")
    else RColorBrewer::brewer.pal(max(3, n_conds), "Set2")[seq_len(n_conds)]
    subtitle_str <- paste("Condition:", paste(levels(df$Condition), collapse = " + "),
                          "| Cell type:", vals$ts_active_celltype)
    
    ggplot(df, aes(x = Timepoint, y = Expression, fill = Condition)) +
      geom_violin(position = position_dodge(width = 0.85),
                  alpha = 0.7, scale = "width", trim = FALSE) +
      geom_boxplot(position = position_dodge(width = 0.85),
                   width = 0.15, fill = "white", alpha = 0.85, outlier.shape = NA) +
      scale_fill_manual(values = setNames(cond_palette, levels(df$Condition))) +
      labs(
        title    = paste("Expression of", gene1),
        subtitle = subtitle_str,
        y = "Log-normalised expression",
        x = paste("Timepoint",
                  if (!is.null(r$timepoint)) paste0("(", r$timepoint, ")") else ""),
        fill = "Condition"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position    = "right",
        axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
        panel.grid.major.x = element_blank()
      )
  })
  
  # === BOXPLOT (single-gene, per timepoint x condition) ===
  # Matches paper §3.3 ("boxplots and violin plots across time points"). The
  # Violin Plot subtab carries an inset boxplot, but this view exposes the
  # boxplot as the primary visualization for readers who prefer the IQR /
  # median summary over the kernel-density shape. Same data slice + palette
  # as ts_plot_logic so the two views are directly comparable.
  ts_boxplot_logic <- reactive({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype,
        vals$ts_active_gene, vals$ts_active_tps)
    r <- vals$ts_active_roles %||% ts_roles()
    req(!is.null(r$timepoint), !is.null(r$condition), !is.null(r$celltype))
    
    meta <- vals$data@meta.data
    ds_mask <- ts_dataset_filter()
    ds_conds <- unique(as.character(meta[[r$condition]][ds_mask]))
    ctrl <- detect_control_conditions(ds_conds)
    conditions_use <- unique(c(vals$ts_active_condition, ctrl))
    mask <- ds_mask &
      meta[[r$condition]] %in% conditions_use &
      meta[[r$celltype]]  == vals$ts_active_celltype &
      meta[[r$timepoint]] %in% vals$ts_active_tps
    if (sum(mask) == 0) return(NULL)
    
    sub_obj <- vals$data[, mask]
    gene1 <- vals$ts_active_gene[1]
    expr  <- as.numeric(GetAssayData(sub_obj, layer = "data")[gene1, ])
    cond_vec <- as.character(sub_obj@meta.data[[r$condition]])
    df <- data.frame(
      Expression = expr,
      Timepoint  = factor(sub_obj@meta.data[[r$timepoint]], levels = vals$ts_active_tps),
      Condition  = factor(cond_vec, levels = conditions_use)
    )
    
    n_conds <- length(levels(df$Condition))
    cond_palette <- if (n_conds <= 2) c("#3498db", "#e74c3c")
    else RColorBrewer::brewer.pal(max(3, n_conds), "Set2")[seq_len(n_conds)]
    subtitle_str <- paste("Condition:", paste(levels(df$Condition), collapse = " + "),
                          "| Cell type:", vals$ts_active_celltype)
    
    ggplot(df, aes(x = Timepoint, y = Expression, fill = Condition)) +
      geom_boxplot(position = position_dodge(width = 0.85),
                   width = 0.65, alpha = 0.85,
                   outlier.size = 0.6, outlier.alpha = 0.5) +
      scale_fill_manual(values = setNames(cond_palette, levels(df$Condition))) +
      labs(
        title    = paste("Expression of", gene1),
        subtitle = subtitle_str,
        y = "Log-normalised expression",
        x = paste("Timepoint",
                  if (!is.null(r$timepoint)) paste0("(", r$timepoint, ")") else ""),
        fill = "Condition"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position    = "right",
        axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
        panel.grid.major.x = element_blank()
      )
  })
  
  output$ts_boxplot_ui <- renderUI({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_gene)
    if (is.null(vals$ts_active_tps) || length(vals$ts_active_tps) == 0) {
      return(p(style = "color:#e67e22;", icon("info-circle"),
               " Please select at least one timepoint and click Show Plot."))
    }
    if (!vals$ts_active_gene[1] %in% rownames(vals$data)) {
      return(p("Selected gene not found in dataset."))
    }
    withSpinner(plotOutput("ts_boxplot",
                           height = paste0(input$ts_violin_height %||% 600, "px")),
                type = 6, color = "#667eea")
  })
  
  output$ts_boxplot <- renderPlot({ ts_boxplot_logic() })
  
  output$ts_boxplot_download_ui <- renderUI({
    req(vals$ts_active_gene, vals$ts_active_tps)
    downloadButton("ts_boxplot_download", "Download Plot", class = "btn-success")
  })
  
  output$ts_boxplot_download <- downloadHandler(
    filename = function() paste0("timeseries_boxplot_", vals$ts_active_gene, "_",
                                 vals$ts_active_condition, ".png"),
    content  = function(file) {
      ggsave(file, plot = ts_boxplot_logic(), width = 12, height = 8, dpi = 300)
    }
  )
  
  # === TEMPORAL TREND (LINE PLOT) ===
  # Mean log-normalised expression with standard-error bars at each
  # timepoint. The single-cell distribution is collapsed to one summary
  # point per timepoint so the cross-timepoint trajectory is the dominant
  # signal. A cell-count annotation is overlaid for sample-size context.
  ts_trend_plot_logic <- reactive({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_gene, vals$ts_active_tps)
    r <- vals$ts_active_roles %||% ts_roles()
    req(!is.null(r$timepoint), !is.null(r$condition), !is.null(r$celltype))
    
    meta <- vals$data@meta.data
    ds_mask <- ts_dataset_filter()
    ds_conds <- unique(as.character(meta[[r$condition]][ds_mask]))
    ctrl <- detect_control_conditions(ds_conds)
    conditions_use <- unique(c(vals$ts_active_condition, ctrl))
    mask <- ds_mask &
      meta[[r$condition]] %in% conditions_use &
      meta[[r$celltype]]  == vals$ts_active_celltype &
      meta[[r$timepoint]] %in% vals$ts_active_tps
    if (sum(mask) == 0) return(NULL)
    
    sub_obj <- vals$data[, mask]
    expr <- as.numeric(GetAssayData(sub_obj, layer = "data")[vals$ts_active_gene[1], ])
    cond_vec <- as.character(sub_obj@meta.data[[r$condition]])
    df <- data.frame(
      Expression = expr,
      Timepoint  = factor(sub_obj@meta.data[[r$timepoint]], levels = vals$ts_active_tps),
      Condition  = factor(cond_vec, levels = conditions_use)
    )
    
    summary_df <- df %>%
      group_by(Condition, Timepoint) %>%
      summarise(mean_expr = mean(Expression, na.rm = TRUE),
                se        = sd(Expression, na.rm = TRUE) / sqrt(n()),
                n_cells   = n(), .groups = "drop")
    
    n_conds <- length(levels(df$Condition))
    cond_palette <- if (n_conds <= 2) c("#3498db", "#e74c3c")
    else RColorBrewer::brewer.pal(max(3, n_conds), "Set2")[seq_len(n_conds)]
    subtitle_str <- paste("Condition:", paste(levels(df$Condition), collapse = " + "),
                          "| Cell type:", vals$ts_active_celltype)
    
    ggplot(summary_df,
           aes(x = Timepoint, y = mean_expr, group = Condition, color = Condition)) +
      geom_line(linewidth = 1.2) +
      geom_point(aes(fill = Condition), size = 4, shape = 21,
                 color = "white", stroke = 1.5) +
      geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se),
                    width = 0.2, alpha = 0.7) +
      geom_text(aes(label = paste0("n=", n_cells)),
                vjust = -1.5, size = 3, color = "#666", show.legend = FALSE) +
      scale_color_manual(values = setNames(cond_palette, levels(df$Condition))) +
      scale_fill_manual (values = setNames(cond_palette, levels(df$Condition))) +
      labs(
        title    = paste("Temporal Trend:", vals$ts_active_gene[1]),
        subtitle = subtitle_str,
        y = "Mean log-normalised expression (± SE)",
        x = "Timepoint",
        color = "Condition", fill = "Condition"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position    = "right",
        axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold")
      )
  })
  
  output$ts_trend_plot_ui <- renderUI({
    req(vals$data, vals$ts_active_condition, vals$ts_active_celltype, vals$ts_active_gene)
    if (is.null(vals$ts_active_tps) || length(vals$ts_active_tps) == 0) {
      return(p(style="color:#e67e22;", icon("info-circle"), " Please select at least one timepoint and click Show Plot."))
    }
    if(!vals$ts_active_gene[1] %in% rownames(vals$data)) {
      return(p("Selected gene not found in dataset."))
    }
    withSpinner(plotOutput("ts_trend_expression_plot",
                           height = paste0(input$ts_trend_height %||% 600, "px")),
                type = 6, color = "#667eea")
  })
  
  output$ts_trend_expression_plot <- renderPlot({
    ts_trend_plot_logic()
  })
  
  output$ts_trend_download_ui <- renderUI({
    req(vals$ts_active_gene, vals$ts_active_tps)
    downloadButton("ts_trend_download_plot", "Download Trend Plot", class = "btn-success")
  })
  
  output$ts_trend_download_plot <- downloadHandler(
    filename = function() { paste0("trend_", vals$ts_active_gene, "_", vals$ts_active_condition, ".png") },
    content = function(file) {
      ggsave(file, plot = ts_trend_plot_logic(), width = 12, height = 8, dpi = 300)
    }
  )
  # === ZOOM LOGIC ===
  explore_zoom_xlim <- reactiveVal(NULL)
  explore_zoom_ylim <- reactiveVal(NULL)
  
  observeEvent(input$umap_brush, {
    brush <- input$umap_brush
    if (!is.null(brush)) {
      explore_zoom_xlim(c(brush$xmin, brush$xmax))
      explore_zoom_ylim(c(brush$ymin, brush$ymax))
    }
  })
  
  observeEvent(input$umap_dblclick, {
    explore_zoom_xlim(NULL)
    explore_zoom_ylim(NULL)
  })
  
  observeEvent(input$reset_umap_zoom, {
    explore_zoom_xlim(NULL)
    explore_zoom_ylim(NULL)
  })
  
  # Snapshot the UMAP settings when the user clicks "Show Plot".
  # This isolates the plots from reactive filter changes until the button
  # is pressed again, matching the UX of the other subtabs.
  umap_snapshot <- eventReactive(input$show_umap_plot, {
    list(
      mask     = vals$explore_mask,
      meta_col = input$explore_filter_meta_col,
      gene     = input$explore_gene,
      subtitle = vals$explore_active_subtitle
    )
  }, ignoreNULL = FALSE)
  
  output$explore_umap_meta <- renderPlot({
    snap <- umap_snapshot()
    req(vals$global_plot_data, snap$meta_col)

    # Subset global DF (FAST)
    umap_df <- if (!is.null(snap$mask)) vals$global_plot_data[snap$mask, ] else vals$global_plot_data

    # Force explicit factor conversion
    group_data <- umap_df[[snap$meta_col]]
    umap_df$group <- as.factor(group_data)

    n_groups <- length(unique(umap_df$group))
    my_palette <- get_expanded_palette(n_groups)

    # The in-plot legend is suppressed here: with many cell-type levels it
    # used to overflow the image border or shrink the plotting panel to an
    # awkward non-square. The full categorical legend is rendered separately
    # below the UMAP row (output$explore_umap_meta_legend) as a wide CSS
    # grid of colour-chip / label pairs where every entry is fully visible.
    # The PNG export still includes the in-plot legend so the saved image
    # stays self-contained (see explore_umap_meta_download).
    p <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=group)) +
      fast_points() +
      scale_color_manual(values = my_palette) +
      theme_minimal(base_size = 14) +
      labs(title = paste("By", snap$meta_col), subtitle = snap$subtitle) +
      theme(legend.position = "none")

    if (!is.null(explore_zoom_xlim()) && !is.null(explore_zoom_ylim())) {
      p <- p + coord_cartesian(xlim = explore_zoom_xlim(), ylim = explore_zoom_ylim())
    }
    p
  }, res = 110, height = function() input$explore_umap_height %||% 1000)

  # Separate full-width metadata legend. Mirrors the metadata UMAP exactly:
  # same factor-level order (the order ggplot uses to map colours), same
  # palette from get_expanded_palette(). Rendered as a CSS grid so the
  # browser auto-flows the chips into as many columns as fit, then wraps
  # vertically; long labels are never truncated.
  output$explore_umap_meta_legend <- renderUI({
    snap <- umap_snapshot()
    req(vals$global_plot_data, snap$meta_col)
    umap_df <- if (!is.null(snap$mask)) vals$global_plot_data[snap$mask, ] else vals$global_plot_data
    if (is.null(umap_df[[snap$meta_col]])) return(NULL)
    # Use levels(as.factor(...)) so the chip order matches the plot's
    # internal factor levels byte-for-byte - same locale, same NA handling.
    group_factor <- as.factor(umap_df[[snap$meta_col]])
    groups <- levels(group_factor)
    n_groups <- length(groups)
    if (n_groups == 0) return(NULL)
    palette <- get_expanded_palette(n_groups)
    chips <- lapply(seq_along(groups), function(i) {
      div(style = "display: flex; align-items: center; gap: 10px; padding: 8px 12px; background: #ffffff; border: 1px solid #e1e5ea; border-radius: 6px; min-height: 36px;",
          div(style = sprintf(
            "width: 22px; height: 22px; flex-shrink: 0; border-radius: 4px; background: %s; border: 1px solid rgba(0,0,0,0.12); box-shadow: inset 0 0 0 1px rgba(255,255,255,0.25);",
            palette[i])),
          span(style = "font-size: 13px; color: #2c3e50; line-height: 1.3; word-break: break-word;",
               groups[i]))
    })
    wellPanel(style = "margin-top: 18px; padding: 18px; background: #f8f9fa;",
              div(style = "display: flex; align-items: baseline; justify-content: space-between; margin-bottom: 12px;",
                  h4(icon("palette"), paste(" Legend —", snap$meta_col),
                     style = "margin: 0; color: #2c3e50;"),
                  span(style = "color: #7f8c8d; font-size: 13px;",
                       sprintf("%d %s", n_groups,
                               if (n_groups == 1) "category" else "categories"))),
              div(style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 8px;",
                  chips))
  })

  output$explore_umap_expr <- renderPlot({
    snap <- umap_snapshot()
    req(vals$data, vals$global_plot_data)
    if(is.null(snap$gene) || snap$gene == "") {
      return(ggplot() + annotate("text", x=0.5, y=0.5, label="Select a gene to view expression", size=6, color="#bdc3c7") + theme_void())
    }
    
    # Check if gene exists in DATA (not filtered object)
    if(!snap$gene %in% rownames(vals$data)) return(NULL); 
    
    # Fetch expression data for ALL cells (Fast)
    expr_all <- GetAssayData(vals$data, layer="data")[snap$gene, ]
    
    # Subset if mask exists
    if (!is.null(snap$mask)) {
      umap_df <- vals$global_plot_data[snap$mask, ]
      expr <- expr_all[snap$mask]
    } else {
      umap_df <- vals$global_plot_data
      expr <- expr_all
    }
    
    umap_df$expr <- expr
    umap_df <- umap_df %>% arrange(expr)
    
    p <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=expr)) +
      fast_points() +
      scale_color_viridis_c(option="plasma", name="Expression") +
      theme_minimal(base_size = 14) +
      labs(title = paste("Gene:", snap$gene))
    
    if (!is.null(explore_zoom_xlim()) && !is.null(explore_zoom_ylim())) {
      p <- p + coord_cartesian(xlim = explore_zoom_xlim(), ylim = explore_zoom_ylim())
    }
    p
  }, res = 110, height = function() input$explore_umap_height %||% 1000)

  # === CLICK METADATA LOGIC ===
  last_explore_click <- reactiveVal(NULL)
  
  observeEvent(input$umap_meta_click, { last_explore_click(input$umap_meta_click) })
  observeEvent(input$umap_expr_click, { last_explore_click(input$umap_expr_click) })
  
  output$explore_click_info <- renderUI({
    click_event <- last_explore_click()
    
    if (is.null(click_event)) {
      return(p(style="color:#999;", icon("mouse-pointer"), " No cell selected. Click on a point in the UMAP plots above."))
    }
    
    req(vals$global_plot_data)
    
    # Target the currently filtered dataframe so coordinates align perfectly
    umap_df <- if (!is.null(vals$explore_mask)) vals$global_plot_data[vals$explore_mask, ] else vals$global_plot_data
    
    # Find the nearest point within a 15-pixel radius
    cell <- nearPoints(umap_df, click_event, xvar = "UMAP_1", yvar = "UMAP_2", maxpoints = 1, threshold = 15)
    
    if (nrow(cell) == 0) {
      return(p(style="color:#e67e22;", icon("exclamation-circle"), " No cell found near click. Try clicking closer to a point."))
    }
    
    # Fetch exact expression if a gene is selected
    expr_text <- "N/A"
    if(!is.null(input$explore_gene) && input$explore_gene != "" && input$explore_gene %in% rownames(vals$data)) {
      expr_val <- GetAssayData(vals$data, layer="data")[input$explore_gene, rownames(cell)]
      expr_text <- round(expr_val, 3)
    }
    
    # Dump every metadata column dynamically. UMAP coordinate columns and
    # the internal geometry column (used by some Seurat objects) are skipped.
    skip_cols <- c("UMAP_1", "UMAP_2", "geometry")
    meta_cols <- setdiff(names(cell), skip_cols)
    meta_rows <- vapply(meta_cols, function(col) {
      val <- cell[[col]]
      if (is.numeric(val)) val <- signif(val, 4)
      sprintf("<tr><th style='white-space:nowrap;'>%s</th><td>%s</td></tr>",
              htmltools::htmlEscape(col),
              htmltools::htmlEscape(as.character(val)))
    }, character(1))
    expr_row <- if (!is.null(input$explore_gene) && nzchar(input$explore_gene))
      sprintf("<tr><th>Expression (%s)</th><td>%s</td></tr>",
              htmltools::htmlEscape(input$explore_gene), expr_text)
    else ""
    HTML(paste0(
      "<table class='table table-condensed' style='margin-bottom: 0; font-size: 14px;'>",
      paste(meta_rows, collapse = ""),
      expr_row,
      "</table>"
    ))
  })
  
  # =========================================================================
  # === Dataset Exploration — Metadata UMAP download handlers ==============
  # =========================================================================
  output$explore_umap_meta_download <- downloadHandler(
    filename = function() paste0("umap_meta_", isolate(input$explore_filter_meta_col %||% "metadata"), ".png"),
    content  = function(file) {
      req(vals$global_plot_data, input$explore_filter_meta_col)
      umap_df <- if (!is.null(vals$explore_mask)) vals$global_plot_data[vals$explore_mask, ] else vals$global_plot_data
      umap_df$group <- as.factor(umap_df[[input$explore_filter_meta_col]])
      n_groups <- length(levels(umap_df$group))
      pal <- get_expanded_palette(n_groups)
      # Pick a legend column count so the legend tiles into rows-of-N rather
      # than one giant single column. Numbers tuned for the wide export
      # canvas below.
      legend_ncol <- if (n_groups <= 8) n_groups
                     else if (n_groups <= 30) 5L
                     else 7L
      n_rows <- as.integer(ceiling(n_groups / legend_ncol))
      p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = group)) +
        geom_point(size = 1.8, alpha = 0.9) +
        scale_color_manual(values = pal,
                           labels = function(x) stringr::str_wrap(x, width = 35)) +
        theme_minimal(base_size = 18) +
        labs(title = paste("By", input$explore_filter_meta_col)) +
        theme(legend.position = "bottom",
              legend.box      = "horizontal",
              plot.title      = element_text(face = "bold", size = 22),
              legend.text     = element_text(size = 11),
              legend.title    = element_text(size = 12),
              legend.key.size = unit(0.55, "cm")) +
        guides(color = guide_legend(ncol = legend_ncol, byrow = TRUE,
                                    override.aes = list(size = 4)))
      # Canvas: the PLOT PANEL is held at ~14 inches square regardless of
      # legend size, and extra vertical room (n_rows * 0.42 in) is added on
      # top for the legend rows. Width widens with legend.ncol so labels do
      # not clip horizontally. dpi = 400 + limitsize = FALSE give a "very
      # big" export - eg a 50-level cell-type column lands as ~6400 x 7200
      # pixels with the plot panel itself preserved at full square area.
      plot_in   <- 14
      legend_in <- 0.8 + n_rows * 0.42
      width_in  <- max(16, 2.4 * legend_ncol + 4)
      total_h   <- plot_in + legend_in
      ggsave(file, p,
             width = width_in, height = total_h,
             dpi = 400, limitsize = FALSE)
    }
  )
  output$explore_umap_expr_download <- downloadHandler(
    filename = function() paste0("umap_expr_", isolate(input$explore_gene %||% "gene"), ".png"),
    content  = function(file) {
      req(vals$data, vals$global_plot_data, input$explore_gene)
      if (!input$explore_gene %in% rownames(vals$data)) return(NULL)
      expr_all <- GetAssayData(vals$data, layer = "data")[input$explore_gene, ]
      if (!is.null(vals$explore_mask)) {
        umap_df <- vals$global_plot_data[vals$explore_mask, ]
        expr    <- expr_all[vals$explore_mask]
      } else {
        umap_df <- vals$global_plot_data
        expr    <- expr_all
      }
      # Expression UMAP has a continuous colour bar (slim), so the canvas
      # just needs to be big and square. dpi = 400 matches the metadata
      # export so the two PNGs look like a coherent pair.
      ggsave(file, build_expression_umap(umap_df, input$explore_gene, expr, raster = FALSE) +
               theme_minimal(base_size = 18) +
               theme(plot.title = element_text(face = "bold", size = 22)),
             width = 16, height = 14, dpi = 400, limitsize = FALSE)
    }
  )
  
  # =========================================================================
  # === Dataset Exploration — COEXPRESSION SUBTAB ==========================
  # =========================================================================
  # Just two side-by-side gene-expression UMAPs. No metadata-orientation
  # panel, no blend overlay. Plots wait for "Show Plot" and refresh when
  # either the button OR the shared filter "Update View" is clicked.
  output$coexp_g1_title <- renderText({
    if (is.null(input$coexp_gene_a) || input$coexp_gene_a == "")
      "Gene 1 (pick a gene)" else paste("Gene 1:", input$coexp_gene_a)
  })
  output$coexp_g2_title <- renderText({
    if (is.null(input$coexp_gene_b) || input$coexp_gene_b == "")
      "Gene 2 (pick a gene)" else paste("Gene 2:", input$coexp_gene_b)
  })
  
  # Shared colour-scale ceiling so the two coexpression UMAPs use the SAME
  # expression range: a given expression value maps to the same colour in
  # both panels, which is what makes them visually comparable. Reads whichever
  # genes are currently selected (respecting the shared cell filter) and
  # returns c(0, max) across both, or NULL when no valid gene is picked (in
  # which case each panel falls back to its own auto-range).
  coexp_shared_limits <- function() {
    if (is.null(vals$data)) return(NULL)
    genes <- c(input$coexp_gene_a, input$coexp_gene_b)
    genes <- genes[!vapply(genes, is.null, logical(1))]
    genes <- genes[nzchar(genes) & genes %in% rownames(vals$data)]
    if (length(genes) == 0) return(NULL)
    mask <- if (!is.null(vals$explore_mask)) vals$explore_mask else rep(TRUE, ncol(vals$data))
    m <- 0
    for (g in genes) {
      e <- GetAssayData(vals$data, layer = "data")[g, ][mask]
      if (length(e)) m <- max(m, max(e, na.rm = TRUE))
    }
    if (!is.finite(m) || m <= 0) NULL else c(0, m)
  }

  coexp_gene_plot <- function(gene, limits = NULL, raster = TRUE) {
    req(vals$data, vals$global_plot_data)
    if (is.null(gene) || gene == "") {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                 label = "Select a gene in the sidebar.",
                                 size = 6, color = "#bdc3c7") + theme_void())
    }
    if (!gene %in% rownames(vals$data)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                 label = paste("Gene", gene, "not in dataset."),
                                 size = 6, color = "#e67e22") + theme_void())
    }
    expr_all <- GetAssayData(vals$data, layer = "data")[gene, ]
    if (!is.null(vals$explore_mask)) {
      df   <- vals$global_plot_data[vals$explore_mask, ]
      expr <- expr_all[vals$explore_mask]
    } else {
      df   <- vals$global_plot_data
      expr <- expr_all
    }
    build_expression_umap(df, gene, expr, limits = limits, raster = raster)
  }

  output$coexp_umap_g1 <- renderPlot({
    req(input$show_coexp_plot, input$show_coexp_plot > 0)
    isolate(coexp_gene_plot(input$coexp_gene_a, limits = coexp_shared_limits()))
  }, res = 110, height = function() input$coexp_plot_height %||% 1000)
  output$coexp_umap_g2 <- renderPlot({
    req(input$show_coexp_plot, input$show_coexp_plot > 0)
    isolate(coexp_gene_plot(input$coexp_gene_b, limits = coexp_shared_limits()))
  }, res = 110, height = function() input$coexp_plot_height %||% 1000)

  output$coexp_g1_download <- downloadHandler(
    filename = function() paste0("coexp_", isolate(input$coexp_gene_a %||% "gene1"), ".png"),
    content  = function(file) ggsave(file, coexp_gene_plot(isolate(input$coexp_gene_a),
                                                           limits = isolate(coexp_shared_limits()),
                                                           raster = FALSE),
                                     width = 10, height = 9, dpi = 300))
  output$coexp_g2_download <- downloadHandler(
    filename = function() paste0("coexp_", isolate(input$coexp_gene_b %||% "gene2"), ".png"),
    content  = function(file) ggsave(file, coexp_gene_plot(isolate(input$coexp_gene_b),
                                                           limits = isolate(coexp_shared_limits()),
                                                           raster = FALSE),
                                     width = 10, height = 9, dpi = 300))
  
  # =========================================================================
  # === Dataset Exploration — MULTI-GENE DOT PLOT SUBTAB ===================
  # =========================================================================
  # Genes come from the selectize plus an optional uploaded gene-list .txt.
  # The shared filter (vals$explore_metadata_filter) is passed through so
  # Update View actually narrows the cells included in the plot. The plot
  # is button-triggered to avoid auto-rerunning on every keystroke.
  observeEvent(input$show_dotplot, {
    req(vals$data, input$dotplot_group_col)
    file_path <- if (!is.null(input$dotplot_genes_file)) input$dotplot_genes_file$datapath else NULL
    genes <- unique(c(input$dotplot_genes, parse_gene_input("", file_path)))
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0) {
      showNotification("Pick at least one gene or upload a .txt list.", type = "warning")
      return()
    }
    p <- create_multi_gene_dotplot(
      seurat_obj      = vals$data,
      genes           = genes,
      group_col       = input$dotplot_group_col,
      filter_list     = vals$explore_metadata_filter,
      cluster_genes   = isTRUE(input$dotplot_cluster_genes),
      scale_expression = isTRUE(input$dotplot_scale_expr)
    )
    vals$dotplot_active_genes <- genes
    vals$dotplot_plot_cached  <- p
    if (inherits(p, "ggplot")) vals$dotplot_data_cached <- p$data
  })
  
  output$dotplot_status_msg <- renderUI({
    if (is.null(vals$dotplot_active_genes))
      return(p(style = "color: #999; font-style: italic;", icon("info-circle"),
               " No genes plotted yet. Select genes (or upload a .txt) and press ",
               tags$b("Show Plot"), "."))
    div(class = "success-box",
        p(icon("check-circle"), strong(" Plotted "),
          length(vals$dotplot_active_genes), " gene(s) across the selected groups."))
  })
  
  output$explore_dotplot <- renderPlot({
    if (!is.null(vals$dotplot_plot_cached)) return(vals$dotplot_plot_cached)
    ggplot() + annotate("text", x = 0.5, y = 0.5,
                        label = "Press 'Plot dot plot' in the sidebar.",
                        size = 6, color = "#bdc3c7") + theme_void()
  }, height = function() input$dotplot_plot_height %||% 700)
  
  output$explore_dotplot_download_png <- downloadHandler(
    filename = function() "multigene_dotplot.png",
    content  = function(file) { req(vals$dotplot_plot_cached); ggsave(file, vals$dotplot_plot_cached, width = 12, height = 9, dpi = 300) })
  output$explore_dotplot_download_pdf <- downloadHandler(
    filename = function() "multigene_dotplot.pdf",
    content  = function(file) { req(vals$dotplot_plot_cached); ggsave(file, vals$dotplot_plot_cached, width = 12, height = 9, device = "pdf") })
  output$explore_dotplot_download_csv <- downloadHandler(
    filename = function() "multigene_dotplot_data.csv",
    content  = function(file) { req(vals$dotplot_data_cached); write.csv(vals$dotplot_data_cached, file, row.names = FALSE) })
  
  # =========================================================================
  # === Dataset Exploration — QC VIOLINS SUBTAB ============================
  # =========================================================================
  # Each violin renders when the user clicks "Show Plot". The current
  # values of the grouping column, log-y toggle, and shared filter are
  # snapshotted via isolate() so plots do not re-fire mid-edit.
  output$qc_status_msg <- renderUI({
    req(vals$data)
    needs <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    missing_cols <- setdiff(needs, colnames(vals$data@meta.data))
    if (length(missing_cols) == 0)
      div(class = "success-box",
          p(icon("check-circle"), " All three QC columns are present (computed at startup if missing)."))
    else
      div(class = "warning-box",
          p(icon("exclamation-triangle"), " Missing QC columns: ",
            paste(missing_cols, collapse = ", ")))
  })
  
  qc_plot_for <- function(qc_col) {
    create_qc_violin_plot(
      seurat_obj  = vals$data,
      qc_col      = qc_col,
      group_col   = if (is.null(input$qc_group_col) || !nzchar(input$qc_group_col)) NULL else input$qc_group_col,
      filter_list = vals$explore_metadata_filter,
      log_y       = isTRUE(input$qc_log_y)
    )
  }
  
  output$qc_violin_nFeature  <- renderPlot({
    req(input$show_qc_plot, input$show_qc_plot > 0)
    isolate(qc_plot_for("nFeature_RNA"))
  }, height = function() input$qc_plot_height %||% 550)
  output$qc_violin_nCount    <- renderPlot({
    req(input$show_qc_plot, input$show_qc_plot > 0)
    isolate(qc_plot_for("nCount_RNA"))
  }, height = function() input$qc_plot_height %||% 550)
  output$qc_violin_percentMt <- renderPlot({
    req(input$show_qc_plot, input$show_qc_plot > 0)
    isolate(qc_plot_for("percent.mt"))
  }, height = function() input$qc_plot_height %||% 550)
  
  output$qc_download_png <- downloadHandler(
    filename = function() "qc_violins.png",
    content  = function(file) {
      p1 <- qc_plot_for("nFeature_RNA"); p2 <- qc_plot_for("nCount_RNA"); p3 <- qc_plot_for("percent.mt")
      combined <- p1 / p2 / p3
      ggsave(file, combined, width = 12, height = 18, dpi = 300)
    })
  output$qc_download_pdf <- downloadHandler(
    filename = function() "qc_violins.pdf",
    content  = function(file) {
      p1 <- qc_plot_for("nFeature_RNA"); p2 <- qc_plot_for("nCount_RNA"); p3 <- qc_plot_for("percent.mt")
      combined <- p1 / p2 / p3
      ggsave(file, combined, width = 12, height = 18, device = "pdf")
    })
  
  output$filter_controls_ui <- renderUI({ req(vals$data); div(fluidRow(column(6, selectInput("filter_col_select", NULL, choices = get_valid_metadata_columns(vals$data, show_hidden = isTRUE(input$show_hidden_metadata)), width = "100%")), column(6, uiOutput("filter_vals_select_ui"))), actionButton("add_filter", "Add Filter", icon = icon("plus"), class = "btn-xs btn-info", style = "width: 100%; margin-top: 5px;")) })
  output$filter_vals_select_ui <- renderUI({ req(input$filter_col_select, vals$data); selectInput("filter_vals_select", NULL, choices = sort(unique(as.character(vals$data@meta.data[[input$filter_col_select]]))), multiple = TRUE, width = "100%") })
  observeEvent(input$add_filter, { col <- input$filter_col_select; val <- input$filter_vals_select; if(is.null(col) || is.null(val)) return(); curr <- vals$active_filters; exist <- which(vapply(curr, function(x) x$col == col, logical(1))); if(length(exist)>0) curr[[exist]]$vals <- val else curr[[length(curr)+1]] <- list(col=col, vals=val); vals$active_filters <- curr })
  observeEvent(input$remove_filter, { idx <- as.numeric(input$remove_filter); if(length(vals$active_filters)>=idx) vals$active_filters[[idx]] <- NULL })
  output$active_filters_ui <- renderUI({ if(length(vals$active_filters)==0) return(NULL); div(class = "active-filters", lapply(seq_along(vals$active_filters), function(i) { f <- vals$active_filters[[i]]; span(span(class="filter-badge", strong(f$col), ": ", paste(head(f$vals, 2), collapse=","), if(length(f$vals)>2)"..."), actionLink("remove_filter", "✕", onclick = paste0("Shiny.onInputChange('remove_filter', ", i, ");"), style="color: red; margin-left: 5px;")) })) })
  output$analysis_groups_ui <- renderUI({ req(vals$data, input$analysis_meta_col); obj <- vals$data; keep_cells <- rep(TRUE, ncol(obj)); for(f in vals$active_filters) { keep_cells <- keep_cells & (obj@meta.data[[f$col]] %in% f$vals) }; if (sum(keep_cells) > 0) { available_groups <- sort(unique(as.character(obj@meta.data[[input$analysis_meta_col]][keep_cells]))) } else { available_groups <- character(0) }; current_sel <- isolate(input$analysis_groups); selected_groups <- intersect(current_sel, available_groups); if (!is.null(vals$cocoa_history_pending)) { if (all(vals$cocoa_history_pending %in% available_groups)) { selected_groups <- vals$cocoa_history_pending; vals$cocoa_history_pending <- NULL } else if (length(intersect(vals$cocoa_history_pending, available_groups)) > 0) { selected_groups <- intersect(vals$cocoa_history_pending, available_groups); vals$cocoa_history_pending <- NULL } else { vals$cocoa_history_pending <- NULL } }; vals$analysis_current_choices <- available_groups; div(div(style = "margin-bottom: 5px; text-align: right;", actionLink("toggle_analysis_groups", "Select / Deselect All")), checkboxGroupInput("analysis_groups", NULL, choices = available_groups, selected = selected_groups)) })
  
  observeEvent(input$toggle_analysis_groups, {
    req(vals$analysis_current_choices)
    if (length(input$analysis_groups) == length(vals$analysis_current_choices)) {
      updateCheckboxGroupInput(session, "analysis_groups", selected = character(0))
    } else {
      updateCheckboxGroupInput(session, "analysis_groups", selected = vals$analysis_current_choices)
    }
  })
  output$analysis_cell_count_ui <- renderUI({ req(vals$data, input$analysis_meta_col); if(length(input$analysis_groups)==0) return(div(class = "cell-count-box", icon("calculator"), " Total Cells: 0")); obj <- vals$data; mask <- obj@meta.data[[input$analysis_meta_col]] %in% input$analysis_groups; for(f in vals$active_filters) mask <- mask & (obj@meta.data[[f$col]] %in% f$vals); div(class = "cell-count-box", icon("users"), paste(" Total Cells (Filtered):", format(sum(mask), big.mark = ","))) })
  output$run_btn_ui <- renderUI({
    # geneCOCOA matches genes against symbol-based pathway gene sets, so an
    # Ensembl-ID dataset must be converted to gene symbols before it can run.
    if (is_ensembl_object(vals$data)) {
      if (isTRUE(vals$cocoa_is_converting)) {
        return(actionButton("cocoa_convert_busy", "Converting gene IDs...",
                            class = "btn-warning btn-block btn-lg disabled",
                            icon = icon("spinner", class = "fa-spin")))
      }
      return(div(
        div(class = "error-box", style = "margin-bottom: 10px;",
            icon("exclamation-circle"),
            " This dataset uses Ensembl gene IDs. geneCOCOA needs gene symbols",
            " (its pathway gene sets are symbol-based) - convert them first."),
        actionButton("cocoa_convert_ensembl", "Convert gene IDs to symbols (biomaRt)",
                     class = "btn-info btn-block btn-lg", icon = icon("dna"),
                     style = "font-weight: bold;")
      ))
    }
    if (length(input$analysis_groups) == 0) return(actionButton("run_analysis_disabled", "Select Groups First", class = "btn-secondary btn-block btn-lg disabled", icon = icon("hand-pointer")))
    if (vals$is_analyzing) {
      tagList(
        actionButton("run_analysis_busy", "Processing...", class = "btn-warning btn-lg disabled", icon = icon("spinner", class="fa-spin"), style="width: 70%;"),
        actionButton("cancel_cocoa", "Stop", class = "btn-danger btn-lg", icon = icon("stop"), style="width: 28%; float: right;")
      )
    } else {
      actionButton("run_analysis", "Run geneCOCOA Analysis", class = "btn-primary btn-block btn-lg", icon = icon("play-circle"), style = "font-weight: bold; padding: 12px;")
    }
  })

  # Freeze the COCOA sidebar config in two situations: (1) the dataset uses
  # Ensembl gene IDs and must be converted to symbols first - the conversion
  # button lives outside this wrapper (in run_btn_ui), so the user is funnelled
  # to it; (2) an analysis is running, so the run cannot be reconfigured
  # mid-flight (the Stop button is also outside the wrapper).
  observe({
    if (is_ensembl_object(vals$data) || isTRUE(vals$is_analyzing)) {
      shinyjs::addClass(id = "cocoa_config_wrap", class = "panel-disabled")
    } else {
      shinyjs::removeClass(id = "cocoa_config_wrap", class = "panel-disabled")
    }
  })
  
  # geneCOCOA Ensembl -> symbol conversion (biomaRt), triggered from run_btn_ui.
  observeEvent(input$cocoa_convert_ensembl, {
    req(vals$data, !isTRUE(vals$cocoa_is_converting))
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      showNotification("biomaRt or its annotation database is not available in your active R library path. Please ensure your local environment is fully configured.",
                       type = "error", duration = 8)
      return()
    }
    obj      <- vals$data
    species  <- detect_species(obj)
    orig_ids <- rownames(obj)
    ens_keys <- sub("\\.[0-9]+$", "", orig_ids)
    vals$cocoa_is_converting <- TRUE
    showNotification("Querying Ensembl BioMart - this can take up to a minute...",
                     type = "message", duration = 6)
    
    future({
      run_biomart_query(unique(ens_keys), species)
    }, globals = list(run_biomart_query = run_biomart_query,
                      ens_keys = ens_keys, species = species),
    packages = "biomaRt", seed = TRUE) %...>% (function(bm) {
      if (is.null(bm) || nrow(bm) == 0) {
        showNotification("biomaRt returned no mappings; gene IDs were left unchanged.",
                         type = "warning", duration = 8)
        return()
      }
      sym       <- bm$external_gene_name[match(ens_keys, bm$ensembl_gene_id)]
      new_names <- ifelse(is.na(sym) | !nzchar(sym), orig_ids, sym)
      new_names <- make.unique(new_names)
      n_mapped  <- sum(new_names != orig_ids)
      new_obj   <- tryCatch(rename_genes_in_object(obj, new_names),
                            error = function(e) NULL)
      if (is.null(new_obj)) {
        showNotification("Could not rebuild the object after conversion.",
                         type = "error", duration = 8)
        return()
      }
      vals$data <- new_obj
      vals$dataset_hash <- paste0(vals$dataset_hash %||% "", "_sym")
      showNotification(sprintf("Converted %d of %d genes to symbols. geneCOCOA is ready.",
                               n_mapped, length(orig_ids)),
                       type = "message", duration = 8)
    }) %...!% (function(err) {
      showNotification(paste("biomaRt conversion failed:", conditionMessage(err)),
                       type = "error", duration = 10)
    }) %>% finally(function() {
      vals$cocoa_is_converting <- FALSE
    })
  })
  
  observeEvent(input$run_analysis, {
    req(input$analysis_gene, input$analysis_groups, !vals$is_analyzing)
    
    obj <- vals$data; meta_col <- input$analysis_meta_col; groups <- input$analysis_groups; gene <- input$analysis_gene; filters <- vals$active_filters
    global_mask <- rep(TRUE, ncol(obj))
    if (!is.null(filters) && length(filters) > 0) { masks <- lapply(filters, function(f) obj@meta.data[[f$col]] %in% f$vals); global_mask <- Reduce("&", masks) }
    
    if (sum(global_mask) == 0) { vals$error_msg <- "Filters resulted in 0 cells."; vals$is_analyzing <- FALSE; shinyjs::enable("run_analysis"); return() }
    
    # Estimate time roughly: ~0.1s per cell in parallel + overhead
    total_est_cells <- 0
    for(g in groups) {
      n_c <- sum(obj@meta.data[[meta_col]] == g & global_mask)
      total_est_cells <- total_est_cells + min(n_c, COCOA_MAX_CELLS)
    }
    
    
    shinyjs::show("cocoa_loading_overlay"); session$sendCustomMessage("start_timer", list(id = "cocoa_timer"))
    vals$analysis_res <- NULL; vals$error_msg <- NULL; vals$is_analyzing <- TRUE; vals$cancel_analysis <- FALSE
    
    # === OPTIMIZATION: Direct Sampling from Metadata (Avoids intermediate large subset) ===
    # Identify cells to keep based on metadata + global_mask without subsetting the object yet.
    
    final_cells_to_keep <- c()
    meta_df <- obj@meta.data
    
    for(g in groups) {
      # Filter metadata for this group AND global mask
      # We use rownames(meta_df) assuming they match colnames(obj) which is standard in Seurat
      in_group_mask <- meta_df[[meta_col]] == g & global_mask
      cells_in_group <- rownames(meta_df)[in_group_mask]
      
      # Downsample if needed
      if (length(cells_in_group) > COCOA_MAX_CELLS) {
        set.seed(42)
        
        # Stratified Sampling: If clusters exist, sample proportionally
        if ("seurat_clusters" %in% names(meta_df)) {
          # Use the already subsetted metadata for this group (virtual subset)
          # We only need the cluster column for these cells
          clusters <- meta_df[cells_in_group, "seurat_clusters"]
          
          # Calculate proportion of each cluster in this group
          cluster_props <- table(clusters) / length(clusters)
          
          # Determine target count for each cluster (minimum 1 cell if present)
          cluster_counts <- ceiling(cluster_props * COCOA_MAX_CELLS)
          
          sampled_cells <- c()
          for (cl in names(cluster_counts)) {
            # Identify cells in this specific cluster within the group
            cells_in_cluster <- cells_in_group[clusters == cl]
            
            if (length(cells_in_cluster) > 0) {
              size_to_pick <- min(length(cells_in_cluster), cluster_counts[[cl]])
              sampled_cells <- c(sampled_cells, sample(cells_in_cluster, size_to_pick))
            }
          }
          
          # Adjust if rounding errors
          if (length(sampled_cells) > COCOA_MAX_CELLS) {
            sampled_cells <- sample(sampled_cells, COCOA_MAX_CELLS)
          } else if (length(sampled_cells) < COCOA_MAX_CELLS) {
            remaining <- setdiff(cells_in_group, sampled_cells)
            needed <- COCOA_MAX_CELLS - length(sampled_cells)
            if(length(remaining) > 0) {
              sampled_cells <- c(sampled_cells, sample(remaining, min(length(remaining), needed)))
            }
          }
          cells_in_group <- sampled_cells
          
        } else {
          # Fallback to random if no clusters
          cells_in_group <- sample(cells_in_group, COCOA_MAX_CELLS)
        }
      }
      final_cells_to_keep <- c(final_cells_to_keep, cells_in_group)
    }
    
    if (length(final_cells_to_keep) == 0) { vals$error_msg <- "No cells remaining after sampling."; vals$is_analyzing <- FALSE; shinyjs::hide("cocoa_loading_overlay"); return() }
    
    # Subset the Seurat object to ONLY the sampled cells (~600), then extract
    # counts from the tiny object. This avoids touching the full 393k-cell
    # sparse matrix and matches the validated v160 approach.
    obj_subset_small <- obj[, final_cells_to_keep]
    DefaultAssay(obj_subset_small) <- "RNA"
    if (inherits(obj_subset_small[["RNA"]], "Assay5")) obj_subset_small[["RNA"]] <- JoinLayers(obj_subset_small[["RNA"]])
    counts_matrix <- GetAssayData(obj_subset_small, layer = "counts")
    meta_data_subset <- obj@meta.data[final_cells_to_keep, , drop = FALSE]
    # Log-transform paper-grade mode is currently disabled in the UI; the
    # worker is therefore always called with the fast configuration.
    log_transform_flag <- FALSE
    n_sims_cfg <- 100L
    # Species inference for msigdbr gene-set lookup.
    species_cfg <- detect_species(vals$data)
    
    future({
      # msigdbr runs INSIDE the future — with multicore (fork), this avoids
      # serializing the gene-set list over the future channel.
      suppressWarnings({ m_df <- msigdbr::msigdbr(species = species_cfg, category = "H") })
      gene_sets <- split(x = m_df$gene_symbol, f = m_df$gs_name)
      
      res_list <- list(); error_log <- c()
      for (g in groups) {
        cells_in_group <- rownames(meta_data_subset)[meta_data_subset[[meta_col]] == g]
        
        if (length(cells_in_group) >= 20) {
          grp_counts <- counts_matrix[, cells_in_group]
          res <- run_cocoa_worker(grp_counts, gene, gene_sets, 10, 100)
          if (!is.null(res)) res_list[[g]] <- res
        } else { error_log <- c(error_log, paste("Group", g, "skipped (<20 cells)")) }
      }
      return(list(res = res_list, logs = error_log))
    }, globals = list(run_cocoa_worker = run_cocoa_worker,
                      counts_matrix = counts_matrix, gene = gene, groups = groups,
                      meta_data_subset = meta_data_subset, meta_col = meta_col,
                      COCOA_MAX_CELLS = COCOA_MAX_CELLS,
                      species_cfg = species_cfg),
    packages = c("Matrix", "geneCOCOA", "dplyr", "Seurat", "msigdbr"), seed = TRUE) %...>% (function(result) {
      if (length(result$res) > 0) { vals$analysis_res <- result$res; vals$analysis_meta <- list(col = meta_col, groups = groups); vals$cocoa_gene_used <- gene; vals$cocoa_filters_used <- filters; save_cocoa_cache(vals$dataset_hash, result$res, filters, vals$analysis_meta, gene) } else { vals$error_msg <- if(length(result$logs)>0) paste("Failed:", paste(result$logs, collapse="; ")) else "No significant pathways." }
    }) %...!% (function(err) { vals$error_msg <- paste("Async Error:", err$message) }) %>% finally(function() {
      vals$is_analyzing <- FALSE; session$sendCustomMessage("stop_timer", list()); shinyjs::hide("cocoa_loading_overlay"); 
    })
  })
  
  
  output$analysis_status_msg <- renderUI({ if(vals$is_analyzing) return(NULL); if(!is.null(vals$error_msg)) return(div(class = "error-box", p(icon("exclamation-triangle"), strong("Error:"), vals$error_msg))); if(!is.null(vals$analysis_res)) return(div(class = "success-box", p(icon("check-circle"), strong("Complete!"), paste("Analyzed", length(vals$analysis_res), "group(s)")))); NULL })
  output$analysis_results_ui <- renderUI({ req(vals$analysis_res); wellPanel(style = "background-color: #f8f9fa; padding: 10px;", div(id = "plot_container", withSpinner(plotOutput("analysis_plot", height = paste0(input$cocoa_plot_height %||% 800, "px")), type = 6, color = "#3498db"))) })
  output$analysis_plot <- renderPlot({ req(vals$analysis_res); gene_to_show <- if(!is.null(vals$cocoa_gene_used)) vals$cocoa_gene_used else "Unknown"; filters_to_show <- if(!is.null(vals$cocoa_filters_used)) vals$cocoa_filters_used else list(); res <- vals$analysis_res; plot_df <- bind_rows(lapply(names(res), function(g) { d <- res[[g]]; if(is.null(d)) return(NULL); d <- as.data.frame(d); if(!("geneset" %in% colnames(d))) d$geneset <- rownames(d); tibble(Pathway = gsub("HALLMARK_", "", d$geneset), PVal = if("p.adj" %in% colnames(d)) d$p.adj else d$p, NLP = -log10(PVal), Group = g) })); if (nrow(plot_df) == 0) { plot.new(); text(0.5, 0.5, "No significant pathways.", cex = 1.5); return() }; plot_df_filtered <- plot_df %>% filter(PVal < input$p_val_thresh); if (nrow(plot_df_filtered) == 0) { plot.new(); text(0.5, 0.5, "No pathways below P-value threshold.", cex = 1.5); return() }; filter_str <- "Global Filters: None"; if (length(filters_to_show) > 0) { f_parts <- sapply(filters_to_show, function(f) { val_str <- paste(head(f$vals, 2), collapse=","); if(length(f$vals)>2) val_str <- paste0(val_str, "..."); paste0(f$col, "=(", val_str, ")") }); filter_str <- paste("Global Filters:", paste(f_parts, collapse="; ")) }; ggplot(plot_df_filtered, aes(x = NLP, y = reorder(Pathway, NLP), fill = Group)) + geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) + scale_fill_brewer(palette = "Set2") + labs(title = paste("Co-regulation:", gene_to_show), subtitle = paste0("Comparison: ", vals$analysis_meta$col, "\n", filter_str), x = "-log10(Adjusted P-Value)", y = NULL) + theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 18), legend.position = "bottom") }, height = 800)
  output$download_ui <- renderUI({
    req(vals$analysis_res)
    div(
      downloadButton("download_analysis", "Download Plot",
                     class = "btn-success btn-block",
                     style = "font-weight: bold;"),
      downloadButton("download_analysis_script", "Download R script",
                     class = "btn-default btn-block",
                     style = "margin-top: 6px;",
                     title = "Self-contained R script that reproduces this COCOA run from the same dataset.")
    )
  })
  output$download_analysis <- downloadHandler(filename = function() paste0("cocoa_", input$analysis_gene, ".png"), content = function(file) ggsave(file, width = 12, height = 9))
  output$download_analysis_script <- downloadHandler(
    filename = function() paste0("atlaslens_cocoa_", input$analysis_gene %||% "gene", "_reproduce.R"),
    content  = function(file) {
      writeLines(generate_cocoa_script(
        gene            = vals$cocoa_gene_used %||% input$analysis_gene %||% "",
        comparison_meta = vals$analysis_meta %||% list(),
        filter_list     = vals$cocoa_filters_used %||% list()
      ), file)
    }
  )
  
  observe({ 
    req(vals$data, input$dea_meta_col)
    obj <- vals$data
    keep_cells <- rep(TRUE, ncol(obj))
    for(f in vals$dea_filter_list) { keep_cells <- keep_cells & (obj@meta.data[[f$col]] %in% f$vals) }
    
    available_groups <- character(0)
    if (sum(keep_cells) > 0) { 
      available_groups <- sort(unique(as.character(obj@meta.data[[input$dea_meta_col]][keep_cells]))) 
    }
    
    # Handle History Restoration for Groups
    curr_g1 <- isolate(input$dea_group1)
    curr_g2 <- isolate(input$dea_group2)
    selected_g1 <- NULL
    selected_g2 <- NULL
    
    if (!is.null(vals$dea_history_pending)) {
      if (vals$dea_history_pending$g1 %in% available_groups && vals$dea_history_pending$g2 %in% available_groups) {
        selected_g1 <- vals$dea_history_pending$g1
        selected_g2 <- vals$dea_history_pending$g2
        vals$dea_history_pending <- NULL # Reset when fully applied
      } else {
        selected_g1 <- vals$dea_history_pending$g1
        selected_g2 <- vals$dea_history_pending$g2
      }
    } else {
      # Default behavior with preservation
      selected_g1 <- if(!is.null(curr_g1) && curr_g1 %in% available_groups) curr_g1 else NULL 
      selected_g2 <- if(!is.null(curr_g2) && curr_g2 %in% available_groups) curr_g2 else if(length(available_groups)>1) available_groups[2] else available_groups[1]
    }
    
    updateSelectInput(session, "dea_group1", choices = available_groups, selected = selected_g1)
    updateSelectInput(session, "dea_group2", choices = available_groups, selected = selected_g2)
  })
  
  output$dea_filter_controls_ui <- renderUI({ req(vals$data); div(fluidRow(column(6, selectInput("dea_filter_col_select", NULL, choices = get_valid_metadata_columns(vals$data, show_hidden = isTRUE(input$show_hidden_metadata)), width = "100%")), column(6, uiOutput("dea_filter_vals_select_ui"))), actionButton("dea_add_filter", "Add Filter", icon = icon("plus"), class = "btn-xs btn-info", style = "width: 100%; margin-top: 5px;")) })
  output$dea_filter_vals_select_ui <- renderUI({ req(input$dea_filter_col_select, vals$data); selectInput("dea_filter_vals_select", NULL, choices = sort(unique(as.character(vals$data@meta.data[[input$dea_filter_col_select]]))), multiple = TRUE, width = "100%") })
  observeEvent(input$dea_add_filter, { col <- input$dea_filter_col_select; val <- input$dea_filter_vals_select; if(is.null(col) || is.null(val)) return(); curr <- vals$dea_filter_list; exist <- which(vapply(curr, function(x) x$col == col, logical(1))); if(length(exist)>0) curr[[exist]]$vals <- val else curr[[length(curr)+1]] <- list(col=col, vals=val); vals$dea_filter_list <- curr })
  observeEvent(input$dea_remove_filter, { idx <- as.numeric(input$dea_remove_filter); if(length(vals$dea_filter_list)>=idx) vals$dea_filter_list[[idx]] <- NULL })
  output$dea_active_filters_ui <- renderUI({ if(length(vals$dea_filter_list)==0) return(NULL); div(class = "active-filters", lapply(seq_along(vals$dea_filter_list), function(i) { f <- vals$dea_filter_list[[i]]; span(span(class="filter-badge", strong(f$col), ": ", paste(head(f$vals, 2), collapse=","), if(length(f$vals)>2)"..."), actionLink("dea_remove_filter", "✕", onclick = paste0("Shiny.onInputChange('dea_remove_filter', ", i, ");"), style="color: red; margin-left: 5px;")) })) })
  output$dea_cell_count_ui <- renderUI({ req(vals$data, input$dea_meta_col, input$dea_group1, input$dea_group2); obj <- vals$data; mask1 <- obj@meta.data[[input$dea_meta_col]] == input$dea_group1; mask2 <- obj@meta.data[[input$dea_meta_col]] == input$dea_group2; for(f in vals$dea_filter_list) { mask1 <- mask1 & (obj@meta.data[[f$col]] %in% f$vals); mask2 <- mask2 & (obj@meta.data[[f$col]] %in% f$vals) }; div(class = "cell-count-box", p(icon("users"), " Group 1: ", format(sum(mask1), big.mark = ",")), p(icon("users"), " Group 2: ", format(sum(mask2), big.mark = ",")), p(strong("Total: "), format(sum(mask1)+sum(mask2), big.mark = ","))) })
  output$dea_run_btn_ui <- renderUI({ 
    if (is.null(input$dea_group1) || is.null(input$dea_group2) || input$dea_group1 == input$dea_group2) {
      actionButton("dea_run_disabled", "Complete Fields First", class = "btn-secondary btn-block btn-lg disabled", icon = icon("hand-pointer"))
    } else if(vals$dea_is_analyzing) {
      actionButton("dea_run_analysis_busy", "Processing...", class = "btn-warning btn-block btn-lg disabled", icon = icon("spinner", class="fa-spin"))
    } else {
      actionButton("dea_run_analysis", "Run Differential Analysis", class = "btn-primary btn-block btn-lg", icon = icon("play-circle"), style = "font-weight: bold;")
    }
  })

  # Freeze every DEA sidebar control (filters, group pickers, thresholds,
  # upload) while an analysis is running so the run cannot be reconfigured
  # mid-flight. The run button itself already shows a disabled "Processing..."
  # state, so it sits outside this wrapper.
  observe({
    if (isTRUE(vals$dea_is_analyzing)) {
      shinyjs::addClass(id = "dea_config_wrap", class = "panel-disabled")
    } else {
      shinyjs::removeClass(id = "dea_config_wrap", class = "panel-disabled")
    }
  })
  
  observeEvent(input$dea_run_analysis, {
    req(input$dea_group1, input$dea_group2, !vals$dea_is_analyzing)
    
    obj <- vals$data; filters <- vals$dea_filter_list; group1 <- input$dea_group1; group2 <- input$dea_group2; meta_col <- input$dea_meta_col; gene_in <- input$dea_gene; hash_in <- vals$dataset_hash
    cells_mask <- rep(TRUE, ncol(obj))
    if (!is.null(filters) && length(filters) > 0) { masks <- lapply(filters, function(f) obj@meta.data[[f$col]] %in% f$vals); cells_mask <- Reduce("&", masks) }
    
    n_total <- sum(cells_mask)
    
    
    shinyjs::show("dea_loading_overlay"); session$sendCustomMessage("start_timer", list(id = "dea_timer"))
    vals$dea_results <- NULL; vals$dea_error_msg <- NULL; vals$dea_is_analyzing <- TRUE
    
    # === v157 approach: extract sparse matrix + meta on main thread ===
    # Only the small subset gets serialized to the worker, NOT the full Seurat object.
    count_mat <- tryCatch(GetAssayData(obj, layer="counts"), error = function(e) GetAssayData(obj, layer="data"))
    subset_counts <- count_mat[, cells_mask, drop=FALSE]
    subset_meta <- obj@meta.data[cells_mask, , drop=FALSE]
    
    future({
      run_dea_worker(subset_counts, subset_meta, meta_col, group1, group2) 
    }, globals=list(run_dea_worker=run_dea_worker, subset_counts=subset_counts, subset_meta=subset_meta, group1=group1, group2=group2, meta_col=meta_col), packages = c("Seurat", "presto", "Matrix"), seed = TRUE) %...>% (function(res) {
      if (!is.null(res) && nrow(res) > 0) { vals$dea_results <- res; vals$dea_filters_used <- filters; vals$dea_meta_used <- list(col=meta_col, g1=group1, g2=group2); save_dea_cache(hash_in, res, filters, list(col=meta_col, groups=c(group1, group2)), gene_in) } else { vals$dea_error_msg <- "No results returned." }
    }) %...!% (function(err) { vals$dea_error_msg <- paste("Async Error:", err$message) }) %>% finally(function() {
      vals$dea_is_analyzing <- FALSE; session$sendCustomMessage("stop_timer", list()); shinyjs::hide("dea_loading_overlay")
    })
  })
  
  
  
  output$dea_status_msg <- renderUI({
    if(vals$dea_is_analyzing) return(NULL)
    if(!is.null(vals$dea_error_msg)) return(div(class = "error-box", h4(icon("exclamation-triangle"), " Analysis Failed"), p(vals$dea_error_msg)))
    if(!is.null(vals$dea_results)) {
      msg <- paste("Analyzed", nrow(vals$dea_results), "genes.")
      filt <- dea_filtered_results()
      if (!is.null(vals$dea_custom_genes) && !is.null(filt)) {
        msg <- paste0(msg, " Restricted to your uploaded gene list (",
                      nrow(filt), " / ", length(vals$dea_custom_genes),
                      " uploaded genes matched).")
      }
      return(div(class = "success-box", p(icon("check-circle"), strong("Complete! "), msg)))
    }
    NULL
  })
  
  # --- Custom gene list upload (parses CSV/TXT into vals$dea_custom_genes) ---
  # Accepts (a) CSV with a column named "gene" or "Gene", (b) any single-column
  # CSV/TSV, (c) a plain text file with one gene per line. Empty rows, NAs,
  # and the literal string "gene" header are dropped. The filename is kept
  # so the status UI can show it back to the user.
  observeEvent(input$dea_gene_upload, {
    f <- input$dea_gene_upload
    req(f)
    genes <- tryCatch({
      raw <- readLines(f$datapath, warn = FALSE)
      raw <- raw[nzchar(trimws(raw))]
      # If the file looks like a CSV/TSV, take the first column.
      if (any(grepl("[,\t]", raw))) {
        parts <- strsplit(raw, "[,\t]")
        col   <- vapply(parts, function(x) trimws(x[1]), character(1))
        # Drop a header row that is literally "gene"
        if (length(col) > 0 && tolower(col[1]) == "gene") col <- col[-1]
        col
      } else {
        trimws(raw)
      }
    }, error = function(e) {
      showNotification(paste("Could not parse uploaded file:", conditionMessage(e)), type = "error", duration = NULL)
      character(0)
    })
    # Strip quotes left over from CSV escaping, drop empties / NAs.
    genes <- gsub('^"|"$', "", genes)
    genes <- genes[nzchar(genes) & !is.na(genes)]
    genes <- unique(genes)
    if (length(genes) == 0) {
      showNotification("Uploaded file contained no usable gene names.", type = "warning")
      vals$dea_custom_genes <- NULL
      vals$dea_custom_genes_filename <- NULL
    } else {
      vals$dea_custom_genes <- genes
      vals$dea_custom_genes_filename <- f$name
      showNotification(paste("Loaded", length(genes), "genes from", f$name), type = "message", duration = 4)
    }
  })
  
  # Clear button handler — also resets the fileInput so the same file can be
  # re-selected (Shiny does not fire input$dea_gene_upload twice for the same
  # path without a reset).
  observeEvent(input$dea_clear_custom_genes, {
    vals$dea_custom_genes <- NULL
    vals$dea_custom_genes_filename <- NULL
    shinyjs::reset("dea_gene_upload")
  })
  
  output$dea_custom_genes_status_ui <- renderUI({
    if (is.null(vals$dea_custom_genes)) return(NULL)
    matched <- if (!is.null(vals$dea_results)) {
      sum(vals$dea_custom_genes %in% vals$dea_results$gene)
    } else NA_integer_
    div(style = "margin-top: 8px; padding: 8px; background-color: #e8f5e9; border-left: 3px solid #27ae60; border-radius: 3px;",
        div(icon("filter"), strong(" Restricting outputs to uploaded list")),
        div(style = "font-size: 0.9em; color: #2c3e50;",
            sprintf("%d genes from %s", length(vals$dea_custom_genes), vals$dea_custom_genes_filename %||% "uploaded file")),
        if (!is.na(matched)) {
          div(style = "font-size: 0.9em; color: #2c3e50;",
              sprintf("%d of %d match the DEA results.", matched, length(vals$dea_custom_genes)))
        })
  })
  
  output$dea_custom_genes_clear_ui <- renderUI({
    if (is.null(vals$dea_custom_genes)) return(NULL)
    actionButton("dea_clear_custom_genes", "Clear uploaded list",
                 icon = icon("times"), class = "btn-warning btn-sm",
                 style = "margin-top: 6px;")
  })
  
  # Volcano + table + downloads all read from this reactive. When no upload
  # is active it is a passthrough; when an upload is active it filters by
  # the uploaded gene symbols.
  dea_filtered_results <- reactive({
    res <- vals$dea_results
    if (is.null(res)) return(NULL)
    if (is.null(vals$dea_custom_genes)) return(res)
    res[res$gene %in% vals$dea_custom_genes, , drop = FALSE]
  })
  output$dea_results_ui <- renderUI({
    req(vals$dea_results)
    tagList(
      # Always-visible reproducibility bar so the R-script export is reachable
      # from any subtab (e.g. right after restoring a run from History, which
      # lands on the Volcano Plot subtab).
      div(style = "margin-bottom: 12px; padding: 10px 14px; background-color: #eef5fb; border-left: 4px solid #3498db; border-radius: 3px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
          span(icon("file-code"), strong(" Reproducibility:"),
               " download a self-contained R script that re-runs this exact comparison."),
          downloadButton("dea_download_script", "Download R script", class = "btn-primary btn-sm")),
      tabsetPanel(
        tabPanel("Volcano Plot", icon = icon("mountain"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("Genome-wide Differential Expression"),
                           withSpinner(plotOutput("dea_volcano_plot", height = "600px"), type = 6, color = "#3498db")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("dea_download_plot", "Download Volcano (PNG)", class = "btn-primary")
                 )
        ),
        tabPanel("Single Gene Expression", icon = icon("vial"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4(textOutput("dea_violin_title")),
                           p("This plot shows the actual expression distribution of your highlighted gene in the two conditions."),
                           withSpinner(plotOutput("dea_violin_plot", height = "500px"), type = 6, color = "#9b59b6")
                 )
        ),
        tabPanel("Results Table", icon = icon("table"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("Detailed Results (All Genes)"),
                           withSpinner(DT::dataTableOutput("dea_results_table"), type = 6, color = "#3498db")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("dea_download_table", "Download Table (CSV)", class = "btn-info")
                 )
        )
      )
    )
  })
  output$dea_volcano_plot <- renderPlot({
    fr <- dea_filtered_results()
    req(fr)
    validate(need(nrow(fr) > 0,
                  "None of the uploaded genes match the DEA results. Clear the upload or run DEA on a comparison that includes these genes."))
    create_dea_volcano_plot(fr, input$dea_gene, input$dea_p_threshold,
                            input$dea_logfc_threshold, vals$dea_filters_used,
                            input$dea_show_significant)
  }, height = 600)
  output$dea_violin_title <- renderText({ if(is.null(input$dea_gene) || input$dea_gene == "") "Select a gene to view expression" else paste("Expression of", input$dea_gene) })
  output$dea_violin_plot <- renderPlot({ req(vals$dea_results, input$dea_gene, vals$data, vals$dea_meta_used); create_gene_violin_plot(vals$data, input$dea_gene, vals$dea_meta_used$col, vals$dea_meta_used$g1, vals$dea_meta_used$g2, vals$dea_filters_used) }, height = function() input$dea_violin_height %||% 500)
  output$dea_results_table <- DT::renderDataTable({
    fr <- dea_filtered_results()
    req(fr)
    validate(need(nrow(fr) > 0,
                  "None of the uploaded genes match the DEA results."))
    p_thr  <- input$dea_p_threshold     %||% 0.05
    fc_thr <- input$dea_logfc_threshold %||% 0.585
    raw <- fr
    if (isTRUE(input$dea_show_significant)) {
      sig_mask <- !is.na(raw$p_val_adj) & raw$p_val_adj < p_thr &
        !is.na(raw$avg_log2FC) & abs(raw$avg_log2FC) >= fc_thr
      raw <- raw[sig_mask, , drop = FALSE]
    }
    raw <- raw[order(raw$p_val_adj, na.last = TRUE), , drop = FALSE]
    table_df <- data.frame(
      gene       = raw$gene,
      avg_log2FC = round(raw$avg_log2FC, 4),
      p_val      = formatC(raw$p_val,     format = "e", digits = 2),
      p_val_adj  = formatC(raw$p_val_adj, format = "e", digits = 2),
      pct.1      = round(raw$pct.1, 3),
      pct.2      = round(raw$pct.2, 3),
      Significant = ifelse(!is.na(raw$p_val_adj) & raw$p_val_adj < p_thr &
                             !is.na(raw$avg_log2FC) & abs(raw$avg_log2FC) >= fc_thr,
                           "Yes", ""),
      stringsAsFactors = FALSE
    )
    dt <- DT::datatable(
      table_df,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      caption  = tags$caption(
        style = "caption-side: top; text-align: left; color: #2c3e50; font-weight: bold;",
        sprintf("Rows marked 'Yes' in Significant: p_val_adj < %.3g and |log2FC| >= %.3g.",
                p_thr, fc_thr))
    )
    DT::formatStyle(dt, "Significant",
                    backgroundColor = DT::styleEqual("Yes", "#e8f5e9"),
                    fontWeight = DT::styleEqual("Yes", "bold"))
  }, server = TRUE)
  output$dea_download_plot <- downloadHandler(
    filename = function() {
      suffix <- if (!is.null(vals$dea_custom_genes)) "_custom_list" else ""
      paste0("volcano_", input$dea_gene, suffix, ".png")
    },
    content = function(file) {
      fr <- dea_filtered_results()
      ggsave(file,
             create_dea_volcano_plot(fr, input$dea_gene, input$dea_p_threshold,
                                     input$dea_logfc_threshold,
                                     show_significant = input$dea_show_significant),
             width = 12, height = 7)
    }
  )
  output$dea_download_table <- downloadHandler(
    filename = function() {
      if (!is.null(vals$dea_custom_genes)) "dea_results_custom_list.csv" else "dea_results.csv"
    },
    content = function(file) write.csv(dea_filtered_results(), file, row.names = FALSE)
  )
  
  # Reproducibility export: hand the reviewer a self-contained R script that
  # re-runs the exact comparison from the same source dataset.
  output$dea_download_script <- downloadHandler(
    filename = function() "atlaslens_dea_reproduce.R",
    content  = function(file) {
      writeLines(generate_dea_script(
        comparison_meta = vals$dea_meta_used %||% list(),
        filter_list     = vals$dea_filters_used %||% list(),
        highlight_gene  = input$dea_gene,
        p_thr  = input$dea_p_threshold %||% 0.05,
        fc_thr = input$dea_logfc_threshold %||% 0.585
      ), file)
    }
  )
  
  # ===========================================================================
  # GO ENRICHMENT SERVER LOGIC
  # ===========================================================================
  
  # Status indicator showing whether DEA results are available
  output$go_dea_status_ui <- renderUI({
    if (!is.null(vals$dea_results) && nrow(vals$dea_results) > 0) {
      n_sig <- sum(vals$dea_results$p_val_adj < 0.05 & abs(vals$dea_results$avg_log2FC) >= 0.585, na.rm = TRUE)
      div(style = "color: #27ae60; font-weight: bold;", icon("check-circle"),
          paste0(" DEA results ready (", nrow(vals$dea_results), " genes, ", n_sig, " significant)"))
    } else {
      div(style = "color: #e67e22;", icon("exclamation-circle"), " Run DEA first, or upload a gene list.")
    }
  })
  
  # Handle the uploaded gene list. Accepts a plain one-gene-per-line file or a
  # delimited table; if a header naming a gene column is present it is used,
  # and an optional fold-change column enables the Up-vs-Down comparison.
  observeEvent(input$go_upload_file, {
    req(input$go_upload_file)
    tryCatch({
      raw <- readLines(input$go_upload_file$datapath, warn = FALSE)
      raw <- trimws(raw)
      raw <- raw[nzchar(raw)]
      if (length(raw) == 0) {
        showNotification("Uploaded file is empty.", type = "error", duration = NULL)
        vals$go_upload_data <- NULL; return()
      }
      cells <- strsplit(raw, "[[:space:],]+")
      first <- tolower(cells[[1]])
      gene_names <- c("gene", "genes", "symbol", "gene_symbol", "gene_id")
      fc_names   <- c("avg_log2fc", "log2fc", "logfc", "avg_logfc", "log2foldchange")
      has_header <- any(first %in% gene_names)
      if (has_header) {
        body <- cells[-1]
        gcol <- which(first %in% gene_names)[1]
        fcol <- which(first %in% fc_names)[1]
      } else {
        body <- cells; gcol <- 1L; fcol <- NA_integer_
      }
      getcol <- function(i) vapply(body, function(x)
        if (!is.na(i) && length(x) >= i) x[i] else NA_character_, character(1))
      df <- data.frame(gene = getcol(gcol), stringsAsFactors = FALSE)
      if (!is.na(fcol)) df$avg_log2FC <- suppressWarnings(as.numeric(getcol(fcol)))
      df <- df[!is.na(df$gene) & nzchar(trimws(df$gene)), , drop = FALSE]
      if (nrow(df) == 0) {
        showNotification("No genes found in the file.", type = "error", duration = NULL)
        vals$go_upload_data <- NULL; return()
      }
      vals$go_upload_data <- df
      msg <- paste("Uploaded", nrow(df), "genes.")
      if (is.null(df$avg_log2FC))
        msg <- paste(msg, "No fold-change column - only the combined analysis is available.")
      showNotification(msg, type = "message")
    }, error = function(e) {
      showNotification(paste("File read error:", conditionMessage(e)), type = "error", duration = NULL)
      vals$go_upload_data <- NULL
    })
  })
  
  # Run button
  output$go_run_btn_ui <- renderUI({
    has_data <- if (input$go_data_source == "dea")
      (!is.null(vals$dea_results) && nrow(vals$dea_results) > 0)
    else
      (!is.null(vals$go_upload_data) && nrow(vals$go_upload_data) > 0)
    if (!has_data) {
      actionButton("go_run_disabled", "No Data Available", class = "btn-secondary btn-block btn-lg disabled", icon = icon("hand-pointer"))
    } else if (vals$go_is_analyzing) {
      actionButton("go_run_busy", "Processing...", class = "btn-warning btn-block btn-lg disabled", icon = icon("spinner", class = "fa-spin"))
    } else {
      actionButton("go_run_analysis", "Run GO Enrichment", class = "btn-success btn-block btn-lg", icon = icon("play-circle"), style = "font-weight: bold;")
    }
  })

  # Freeze the GO sidebar config (data source, mode, thresholds, ontology)
  # while an enrichment run is in progress so it cannot be reconfigured
  # mid-flight. The run button sits outside this wrapper and shows its own
  # disabled "Processing..." state.
  observe({
    if (isTRUE(vals$go_is_analyzing)) {
      shinyjs::addClass(id = "go_config_wrap", class = "panel-disabled")
    } else {
      shinyjs::removeClass(id = "go_config_wrap", class = "panel-disabled")
    }
  })

  # Main analysis trigger
  observeEvent(input$go_run_analysis, {
    req(!vals$go_is_analyzing)
    
    from_dea <- input$go_data_source == "dea"
    src_df   <- if (from_dea) vals$dea_results else vals$go_upload_data
    req(src_df, nrow(src_df) > 0)
    
    p_cut  <- input$go_p_cutoff     %||% 0.05
    fc_cut <- input$go_logfc_cutoff %||% 0.585
    mode   <- input$go_mode
    has_fc <- "avg_log2FC" %in% colnames(src_df)
    
    # Select significant genes. DEA results are filtered by adjusted p-value
    # and |log2FC|; an uploaded list is already significant, so it is used
    # as-is (its fold-change column, if present, only sets gene direction).
    if (from_dea) {
      sig_df <- src_df[!is.na(src_df$p_val_adj) & src_df$p_val_adj < p_cut &
                         !is.na(src_df$avg_log2FC) & abs(src_df$avg_log2FC) >= fc_cut, , drop = FALSE]
    } else {
      sig_df <- src_df
    }
    if (nrow(sig_df) == 0) {
      showNotification("No significant genes pass the thresholds. Try relaxing the cutoffs.", type = "warning")
      return()
    }
    if (mode == "compare" && !has_fc) {
      showNotification("Compare mode needs a fold-change column. Upload a list with avg_log2FC, or choose 'All DE genes (combined)'.", type = "warning")
      return()
    }
    
    # Build the gene list(s) for the selected mode.
    if (mode == "compare") {
      up_genes   <- unique(sig_df$gene[!is.na(sig_df$avg_log2FC) & sig_df$avg_log2FC > 0])
      down_genes <- unique(sig_df$gene[!is.na(sig_df$avg_log2FC) & sig_df$avg_log2FC < 0])
      gene_list  <- list(Upregulated = up_genes, Downregulated = down_genes)
      gene_list  <- gene_list[vapply(gene_list, length, integer(1)) > 0]
      if (length(gene_list) == 0) {
        showNotification("No up- or down-regulated genes to compare.", type = "warning")
        return()
      }
    } else {
      gene_list <- list(genes = unique(sig_df$gene))
    }
    
    # Auto-detect species and gene-ID type from the genes themselves.
    all_genes <- unique(unlist(gene_list, use.names = FALSE))
    species <- if (!is.null(vals$data)) detect_species(vals$data) else {
      if (mean(grepl("^ENSMUS", all_genes)) > 0.5) "Mus musculus"
      else if (mean(grepl("^ENSG[0-9]", all_genes)) > 0.5) "Homo sapiens"
      else {
        alpha <- all_genes[grepl("^[A-Za-z]", all_genes)]
        if (length(alpha) > 0 && mean(alpha == toupper(alpha)) > 0.7) "Homo sapiens" else "Mus musculus"
      }
    }
    key_type <- if (length(all_genes) > 0 &&
                    mean(grepl("^ENS", all_genes)) > 0.5) "ENSEMBL" else "SYMBOL"
    
    ontology <- input$go_ontology
    q_cutoff <- input$go_q_cutoff %||% 0.1
    
    # Settings captured for the History cache so a saved run can be restored.
    hash_in <- vals$dataset_hash
    go_settings <- list(data_source = input$go_data_source, mode = mode,
                        ontology = ontology, p_cut = p_cut,
                        q_cut = q_cutoff, logfc_cut = fc_cut)
    
    shinyjs::show("go_loading_overlay")
    session$sendCustomMessage("start_timer", list(id = "go_timer"))
    vals$go_results <- NULL
    vals$go_error_msg <- NULL
    vals$go_is_analyzing <- TRUE
    
    future({
      run_go_worker(gene_list, species, ontology, p_cut, q_cutoff, mode, key_type)
    }, globals = list(
      run_go_worker = run_go_worker, gene_list = gene_list, species = species,
      ontology = ontology, p_cut = p_cut, q_cutoff = q_cutoff,
      mode = mode, key_type = key_type
    ), packages = c("clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
                    "AnnotationDbi", "rrvgo", "GOSemSim"),
    seed = TRUE) %...>% (function(res) {
      if (!is.null(res)) {
        vals$go_results <- res
        tryCatch(save_go_cache(hash_in, res, go_settings), error = function(e) NULL)
      } else {
        vals$go_error_msg <- "No results returned from GO enrichment."
      }
    }) %...!% (function(err) {
      vals$go_error_msg <- paste("GO Enrichment Error:", err$message)
    }) %>% finally(function() {
      vals$go_is_analyzing <- FALSE
      session$sendCustomMessage("stop_timer", list())
      shinyjs::hide("go_loading_overlay")
    })
  })
  
  # Total enriched terms across whatever data frames the result holds.
  go_total_terms <- function(res) {
    if (is.null(res)) return(0)
    sum(vapply(list(res$all_df, res$up_df, res$down_df),
               function(d) if (is.null(d)) 0L else nrow(d), integer(1)))
  }
  
  # Status message
  output$go_status_msg <- renderUI({
    if (vals$go_is_analyzing) return(NULL)
    if (!is.null(vals$go_error_msg)) {
      return(div(class = "error-box", h4(icon("exclamation-triangle"), " Analysis Failed"), p(vals$go_error_msg)))
    }
    if (!is.null(vals$go_results)) {
      return(div(class = "success-box",
                 p(icon("check-circle"), strong("Complete!"),
                   paste0("Found ", go_total_terms(vals$go_results), " enriched GO terms."))))
    }
    NULL
  })
  
  # Named list of enrichment data frames the dot plot / table consume.
  go_df_list <- reactive({
    res <- vals$go_results
    req(res)
    l <- if (identical(res$mode, "compare"))
      list(Upregulated = res$up_df, Downregulated = res$down_df)
    else
      list(`All DE genes` = res$all_df)
    l[!vapply(l, is.null, logical(1))]
  })
  
  # Results UI: Dot Plot + Results Table.
  output$go_results_ui <- renderUI({
    req(vals$go_results)
    if (go_total_terms(vals$go_results) == 0) {
      return(div(class = "error-box",
                 h4(icon("info-circle"), " No Enriched GO Terms"),
                 p("No GO terms passed the enrichment cutoffs for this gene set. ",
                   "Try relaxing the P-value / Q-value cutoffs, or switching the ",
                   "ontology (BP / MF / CC).")))
    }
    tagList(
      # Always-visible reproducibility bar (mirrors the DEA tab) so the
      # R-script export is reachable from any subtab, including right after a
      # History restore which lands on the Dot Plot subtab.
      div(style = "margin-bottom: 12px; padding: 10px 14px; background-color: #f4ecf7; border-left: 4px solid #8e44ad; border-radius: 3px; display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
          span(icon("file-code"), strong(" Reproducibility:"),
               " download a self-contained R script that re-runs this GO enrichment from a DEA results table."),
          downloadButton("go_download_script", "Download R script", class = "btn-primary btn-sm")),
      tabsetPanel(
        tabPanel("Dot Plot", icon = icon("braille"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("GO Enrichment Dot Plot"),
                           p("Top enriched GO terms. GeneRatio on the x-axis; dot size = gene count, colour = adjusted p-value."),
                           withSpinner(plotOutput("go_dotplot", height = "650px"), type = 6, color = "#8e44ad")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("go_download_dotplot", "Download Dot Plot (PNG)", class = "btn-primary"))
        ),
        tabPanel("Treemap", icon = icon("th-large"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("Semantic Cluster Treemap"),
                           p("Redundant GO terms collapsed via rrvgo. Each tile is a parent term; tile area scales with -log10(p.adjust)."),
                           withSpinner(plotOutput("go_treemap", height = "650px"), type = 6, color = "#8e44ad")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("go_download_treemap", "Download Treemap (PNG)", class = "btn-primary"))
        ),
        tabPanel("Scatter", icon = icon("circle-nodes"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("Semantic Similarity Scatter"),
                           p("2D MDS projection of the GO term similarity matrix. Points are coloured by parent term; the spread reflects how semantically distant the GO terms are."),
                           withSpinner(plotOutput("go_scatter", height = "650px"), type = 6, color = "#8e44ad")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("go_download_scatter", "Download Scatter (PNG)", class = "btn-primary"))
        ),
        tabPanel("Results Table", icon = icon("table"), br(),
                 wellPanel(style = "background-color: #f8f9fa; padding: 15px;",
                           h4("Full GO Enrichment Results"),
                           withSpinner(DT::dataTableOutput("go_full_table"), type = 6, color = "#3498db")
                 ),
                 div(style = "margin-top: 15px;",
                     downloadButton("go_download_full_table", "Download Results (CSV)", class = "btn-info"))
        )
      )
    )
  })
  
  # Dot plot
  output$go_dotplot <- renderPlot({
    dl <- go_df_list(); req(length(dl) > 0)
    create_go_dotplot(dl)
  }, height = 650)
  
  # Named list of rrvgo reductions (one per panel) for the treemap + scatter.
  # Mirrors go_df_list() so the panel labels stay in sync between subtabs.
  go_rrvgo_list <- reactive({
    res <- vals$go_results
    req(res)
    l <- if (identical(res$mode, "compare"))
      list(Upregulated = res$up_rrvgo, Downregulated = res$down_rrvgo)
    else
      list(`All DE genes` = res$all_rrvgo)
    l[!vapply(l, is.null, logical(1))]
  })
  
  # Treemap. In compare mode there are two panels (Up + Down). rrvgo's
  # treemap draws via *grid*, so graphics::par(mfrow) has no effect on it -
  # both treemaps would land on the full device and the second (Down) would
  # overwrite the first (Up). Instead we push a 1xN grid layout and hand each
  # treemap its own column viewport via create_go_treemap(vp = ...).
  output$go_treemap <- renderPlot({
    rl <- go_rrvgo_list()
    # Keep only panels that actually produced a reduction so the layout width
    # matches the number of treemaps we can draw.
    rl <- rl[!vapply(rl, function(x) is.null(x) || is.null(x$reduced), logical(1))]
    validate(need(length(rl) > 0,
                  "No semantic-similarity reduction available. This usually means too few enriched terms or the ontology does not support rrvgo reduction (BP/MF/CC only)."))
    if (length(rl) == 1) {
      create_go_treemap(rl[[1]], panel_title = names(rl)[1])
    } else {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(rl))))
      for (i in seq_along(rl)) {
        create_go_treemap(
          rl[[i]], panel_title = names(rl)[i],
          vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
      }
    }
  }, height = 650)
  
  # Scatter. ggplot output, so faceting with patchwork is trivial.
  output$go_scatter <- renderPlot({
    rl <- go_rrvgo_list()
    validate(need(length(rl) > 0,
                  "No semantic-similarity reduction available. This usually means too few enriched terms or the ontology does not support rrvgo reduction (BP/MF/CC only)."))
    plots <- lapply(names(rl), function(nm) create_go_scatter(rl[[nm]], panel_title = nm))
    plots <- plots[!vapply(plots, is.null, logical(1))]
    validate(need(length(plots) > 0, "Could not render the scatter plot for this reduction."))
    if (length(plots) == 1) {
      plots[[1]]
    } else if (requireNamespace("patchwork", quietly = TRUE)) {
      Reduce(`+`, plots) + patchwork::plot_layout(ncol = length(plots))
    } else {
      plots[[1]]
    }
  }, height = 650)
  
  # Combined results table (a Direction column flags up vs down)
  go_results_table <- reactive({
    dl <- go_df_list(); req(length(dl) > 0)
    do.call(rbind, lapply(names(dl), function(nm)
      data.frame(Direction = nm, dl[[nm]], stringsAsFactors = FALSE, check.names = FALSE)))
  })
  
  output$go_full_table <- DT::renderDataTable({
    DT::datatable(go_results_table(), options = list(pageLength = 25, scrollX = TRUE),
                  rownames = FALSE)
  }, server = TRUE)
  
  # --- DOWNLOAD HANDLERS ---
  output$go_download_dotplot <- downloadHandler(
    filename = function() "go_dotplot.png",
    content  = function(file) {
      dl <- go_df_list()
      ggplot2::ggsave(file, create_go_dotplot(dl), width = 12, height = 8, dpi = 300)
    }
  )
  output$go_download_treemap <- downloadHandler(
    filename = function() "go_treemap.png",
    content  = function(file) {
      rl <- go_rrvgo_list()
      rl <- rl[!vapply(rl, function(x) is.null(x) || is.null(x$reduced), logical(1))]
      if (length(rl) == 0) return()
      # treemap draws via grid; render to a PNG device directly. Widen the
      # canvas for the side-by-side compare-mode layout.
      grDevices::png(file, width = 8 * length(rl), height = 8, units = "in", res = 300)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (length(rl) == 1) {
        create_go_treemap(rl[[1]], panel_title = names(rl)[1])
      } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, length(rl))))
        for (i in seq_along(rl)) {
          create_go_treemap(
            rl[[i]], panel_title = names(rl)[i],
            vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
        }
      }
    }
  )
  output$go_download_scatter <- downloadHandler(
    filename = function() "go_scatter.png",
    content  = function(file) {
      rl <- go_rrvgo_list()
      plots <- lapply(names(rl), function(nm) create_go_scatter(rl[[nm]], panel_title = nm))
      plots <- plots[!vapply(plots, is.null, logical(1))]
      if (length(plots) == 0) return()
      if (length(plots) == 1) {
        ggplot2::ggsave(file, plots[[1]], width = 10, height = 8, dpi = 300)
      } else if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- Reduce(`+`, plots) + patchwork::plot_layout(ncol = length(plots))
        ggplot2::ggsave(file, combined, width = 10 * length(plots), height = 8, dpi = 300)
      } else {
        ggplot2::ggsave(file, plots[[1]], width = 10, height = 8, dpi = 300)
      }
    }
  )
  output$go_download_full_table <- downloadHandler(
    filename = function() "go_enrichment_results.csv",
    content  = function(file) utils::write.csv(go_results_table(), file, row.names = FALSE)
  )
  
  # Reproducibility export: re-run this GO enrichment from a DEA CSV with the
  # same mode / ontology / cutoffs. Cell-filter is inherited from the DEA
  # run that produced the source CSV, so we pass list() here.
  output$go_download_script <- downloadHandler(
    filename = function() "atlaslens_go_reproduce.R",
    content  = function(file) {
      writeLines(generate_go_script(
        go_settings = list(
          data_source = input$go_data_source,
          mode        = input$go_mode,
          ontology    = input$go_ontology,
          p_cut       = input$go_p_cutoff,
          q_cut       = input$go_q_cutoff,
          logfc_cut   = input$go_logfc_cutoff
        ),
        filter_list = list()
      ), file)
    }
  )
  
  # ===========================================================================
  # END GO ENRICHMENT
  # ===========================================================================
  
  observeEvent(input$trigger_load, {
    file_path <- file.path(CACHE_DIR, input$trigger_load)
    if (file.exists(file_path)) {
      tryCatch({
        cache_data <- if(grepl("\\.qs$", file_path)) qs::qread(file_path) else readRDS(file_path)
        if (!is.null(cache_data$type) && cache_data$type == "COCOA") {
          updateSelectizeInput(session, "analysis_gene", selected = cache_data$gene)
          if (!is.null(cache_data$comparison_meta$groups)) { vals$cocoa_history_pending <- cache_data$comparison_meta$groups }
          updateSelectInput(session, "analysis_meta_col", selected = cache_data$comparison_meta$col)
          vals$active_filters <- cache_data$filter_list; vals$analysis_res <- cache_data$results; vals$analysis_meta <- cache_data$comparison_meta
          vals$cocoa_gene_used <- cache_data$gene; vals$cocoa_filters_used <- cache_data$filter_list
          updateNavbarPage(session, "main_nav", selected = "Analysis")
        } else if (!is.null(cache_data$type) && cache_data$type == "GO") {
          s <- cache_data$go_settings
          if (!is.null(s)) {
            updateRadioButtons(session, "go_data_source",  selected = s$data_source)
            updateRadioButtons(session, "go_mode",         selected = s$mode)
            updateSelectInput(session,  "go_ontology",     selected = s$ontology)
            updateNumericInput(session, "go_p_cutoff",     value = s$p_cut)
            updateNumericInput(session, "go_q_cutoff",     value = s$q_cut)
            updateNumericInput(session, "go_logfc_cutoff", value = s$logfc_cut)
          }
          vals$go_results <- cache_data$results
          updateNavbarPage(session, "main_nav", selected = "GO")
        } else {
          if(is.null(cache_data$gene) || cache_data$gene == "All Genes" || cache_data$gene == "Unknown") { updateSelectizeInput(session, "dea_gene", selected = "") } else { updateSelectizeInput(session, "dea_gene", selected = cache_data$gene) }
          
          # Store group selections in pending state for the observer
          if (!is.null(cache_data$comparison_meta$groups) && length(cache_data$comparison_meta$groups) >= 2) {
            vals$dea_history_pending <- list(
              g1 = cache_data$comparison_meta$groups[1],
              g2 = cache_data$comparison_meta$groups[2]
            )
          }
          
          updateSelectInput(session, "dea_meta_col", selected = cache_data$comparison_meta$col)
          vals$dea_filter_list <- cache_data$filter_list; vals$dea_filters_used <- cache_data$filter_list; vals$dea_meta_used <- list(col=cache_data$comparison_meta$col, g1=cache_data$comparison_meta$groups[1], g2=cache_data$comparison_meta$groups[2]); vals$dea_results <- cache_data$results
          updateNavbarPage(session, "main_nav", selected = "DEA")
        }
        showNotification("Restored from history!", type = "message")
      }, error = function(e) showNotification("Error loading history.", type = "error", duration = NULL))
    }
  })
  output$history_list <- renderUI({ input$refresh_history; hist <- list_cached_analyses(); if (is.null(hist)) return(p(style = "text-align: center; color: #999;", icon("info-circle"), " No history.")); tagList(lapply(1:nrow(hist), function(i) { h <- hist[i, ]; type_badge <- if(h$type == "COCOA") span(class="type-badge COCOA", "COCOA") else if(h$type == "GO") span(class="type-badge GO", "GO") else span(class="type-badge DEA", "DEA"); div(class = "history-item-row", onclick = paste0("Shiny.setInputValue('trigger_load', '", h$cache_file, "', {priority: 'event'});"), fluidRow(column(8, type_badge, strong(h$gene), br(), span(style="font-size:12px;", h$desc)), column(4, style="text-align:right;", span(style="font-size:10px;", format(h$timestamp, "%m-%d %H:%M"))))) })) })
  
  session$onSessionEnded(function() {
    gc()
  })
}

shinyApp(ui, server)