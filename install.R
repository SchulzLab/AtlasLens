#!/usr/bin/env Rscript
# =============================================================================
# AtlasLens - R dependency installer
# =============================================================================
# Installs every R package required to run AtlasLens (app.R).
#
#   Validated environment :  R 4.3.3  +  Bioconductor 3.18
#   Usage                 :  Rscript install.R
#
# For a reproducible, pinned environment we recommend the conda route instead (environment.yml).
# This script is the portable path for any machine with a working R 4.3.x.
# =============================================================================

message("=====================================================")
message(" AtlasLens dependency installer")
message("=====================================================")

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Bioconductor release that matches the validated R 4.3.x environment.
# Change this only if you deliberately run a different R version.
bioc_version <- "3.18"

# --- R version sanity check --------------------------------------------------
r_ver <- getRversion()
if (r_ver < "4.3" || r_ver >= "4.4") {
  message("")
  message("NOTE: AtlasLens was validated on R 4.3.3 with Bioconductor ", bioc_version, ".")
  message("      You are running R ", r_ver, ". If the Bioconductor step fails,")
  message("      set 'bioc_version' above to the release matching your R version")
  message("      (see https://bioconductor.org/about/release-announcements/).")
}

# --- helper: install only what is missing ------------------------------------
install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0L) {
    message("  already installed: ", paste(pkgs, collapse = ", "))
    return(invisible(TRUE))
  }
  message("  installing: ", paste(missing, collapse = ", "))
  installer(missing)
  invisible(TRUE)
}

# --- 1. bootstrap tools ------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")

# --- 2. CRAN packages --------------------------------------------------------
cran_pkgs <- c(
  # Shiny UI / interactivity
  "shiny", "shinyjs", "shinycssloaders", "DT",
  # single-cell core
  "Seurat", "Matrix",
  # plotting
  "ggplot2", "ggrepel", "plotly", "viridis", "RColorBrewer", "patchwork",
  # fast rasterised point layer for the big UMAPs (geom_point_rast)
  "ggrastr",
  # data wrangling
  "dplyr", "tidyr",
  # async / caching / misc
  "future", "promises", "qs", "digest",
  # landing-page config loader (reads landing_config.json)
  "jsonlite",
  # gene sets
  "msigdbr",
  # rendering back-ends used by rrvgo's treemap / heatmap plots
  "treemap", "pheatmap"
)
message("")
message("[1/3] CRAN packages")
install_if_missing(cran_pkgs, function(p) install.packages(p))

# --- 3. Bioconductor packages (GO enrichment) --------------------------------
bioc_pkgs <- c(
  "clusterProfiler",   # GO over-representation analysis (enrichGO / compareCluster)
  "rrvgo",             # semantic similarity reduction of GO terms
  "GOSemSim",          # GO semantic similarity back-end for rrvgo
  "enrichplot",        # dotplot() for enrichResult / compareClusterResult
  "DOSE",              # defines the enrichResult S4 class
  "AnnotationDbi",     # OrgDb access layer
  "biomaRt",           # Ensembl ID -> gene symbol conversion (geneCOCOA tab)
  "GO.db",             # Gene Ontology term database
  "org.Hs.eg.db",      # human gene annotation
  "org.Mm.eg.db"       # mouse gene annotation
)
message("")
message("[2/3] Bioconductor packages (Bioconductor ", bioc_version, ")")
install_if_missing(bioc_pkgs, function(p) {
  BiocManager::install(p, version = bioc_version, update = FALSE, ask = FALSE)
})

# --- 4. GitHub-only packages -------------------------------------------------
# Neither package is on CRAN, Bioconductor, or conda - they must come from git.
message("")
message("[3/3] GitHub packages")
if (!requireNamespace("presto", quietly = TRUE)) {
  message("  installing: presto (immunogenomics/presto)")
  remotes::install_github("immunogenomics/presto", upgrade = "never")
}
if (!requireNamespace("geneCOCOA", quietly = TRUE)) {
  message("  installing: geneCOCOA (si-ze/geneCOCOA)")
  remotes::install_github("si-ze/geneCOCOA", upgrade = "never")
}

# --- 5. verify ---------------------------------------------------------------
message("")
message("=====================================================")
message(" Verifying installation")
message("=====================================================")
all_pkgs <- c(cran_pkgs, bioc_pkgs, "presto", "geneCOCOA")
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
for (p in all_pkgs) message(sprintf("  %-18s %s", p, if (ok[[p]]) "OK" else "MISSING"))

if (all(ok)) {
  message("")
  message("All ", length(all_pkgs), " packages installed. AtlasLens is ready to run:")
  message("  R -e 'shiny::runApp(\"app.R\", launch.browser = TRUE)'")
} else {
  message("")
  message("WARNING: ", sum(!ok), " package(s) failed to install - see the log above.")
  if (!interactive()) quit(status = 1L)
}
