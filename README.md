# AtlasLens

**The source  code behind the Shiny application for reproducible exploration of single-cell RNA-seq atlases.**

---

## Table of contents

- [Overview](#overview)
- [Features](#features)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Dataset and configuration](#dataset-and-configuration)
- [Running AtlasLens](#running-atlaslens)
- [Usage guide](#usage-guide)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Overview

AtlasLens is a single-file R/Shiny application (`app.R`). It loads one single-cell RNA-seq atlas or dataset (a Seurat object) and exposes a suite
of analysis tools through a tabbed interface. Heavy computations (differential expression, geneCOCOA, GO enrichment) run **asynchronously** in
background workers (via `future` / `promises`), so the interface stays
responsive, and results are cached on disk (`qs`) for instant retrieval.

The application auto-detects metadata role columns (timepoint, condition, cell type, dataset) so it can
be pointed at any compatible atlas.

The reference dataset is the Tabula Muris atlas; and a time-resolved whole-lung single-cell atlas of bleomycin-induced lung injury
and fibrosis is used to demonstrate the time-series and
condition-specific modules.

## Features

| Tab | Purpose |
|-----|---------|
| **Introduction** | Landing page; dataset load status and app overview. |
| **Dataset Exploration** | Interactive UMAP projections with a shared multi-dimensional cell filter. Sub-tabs: metadata-coloured UMAP, gene co-expression overlay, multi-gene dot plot, and QC violin plots. |
| **Gene Function (COCOA)** | Runs [geneCOCOA](https://github.com/si-ze/geneCOCOA) to identify biological pathways co-regulated with a gene of interest across user-defined groups. Results are displayed as comparative bar plots. |
| **Differential Expression** | Genome-wide two-group comparison using Seurat `FindMarkers` (Wilcoxon test, accelerated by `presto`). Volcano plot, per-gene expression violins, and a downloadable results table. Optional restriction of outputs to an uploaded gene list (e.g. a Time Series cluster export). |
| **GO Enrichment** | Gene Ontology over-representation analysis on DEA results or an uploaded gene table. clusterProfiler `enrichGO` / `compareCluster` plus rrvgo semantic clustering. See [below](#go-enrichment). |
| **Time Series Information** | Temporal expression analysis across timepoints, k-means clustering of genes by their temporal patterns, mean-expression heatmap, mean-trajectory line plot, per-timepoint boxplot, and per-timepoint violin plot. |
| **History** | One-click restore of previously cached COCOA / DEA / GO analyses. Restoring a saved run also re-enables its **Download R script** export. |

### Multi-dimensional metadata filtering

Cells can be subset iteratively by any combination of metadata columns
(e.g. *Condition = AMI AND Sex = Female AND Day = 7*). The filter applies
across the Dataset Exploration, DEA, COCOA and GO tabs, so the same cohort
of cells is analysed consistently throughout a session. Non-informative
columns (per-cell barcode, raw QC counters, columns with a single value
across the whole object, and any non-atomic list-columns) are hidden by
default; the Dataset Overview card lists them explicitly and exposes a
toggle to reveal them.

### Cluster → DEA bridge

The Time Series heatmap clusters genes by temporal pattern. Each cluster is
labelled with the timepoint at which it peaks. A one-click **Send to DEA**
button copies the selected cluster's genes into the DEA tab and restricts
the volcano plot and results table to that gene set; a parallel
**Download (CSV)** button produces a portable file that can be re-uploaded
later.

### Reproducible R-script export

Every DEA, COCOA and GO results panel exposes a **Download R script** button
that emits a self-contained R script reproducing the exact analysis from
the same source dataset. Each script:

- loads the Seurat object from `DATASET_PATH`,
- re-applies the multi-dimensional cell filter that was active at run time,
- re-runs the analysis with the same comparison, thresholds and mode,
- writes the results to disk.

The History tab restores any cached analysis to its tab; the script export
button on that tab then reflects the restored settings, so previously saved
runs are also recoverable as code.

### GO Enrichment

The GO Enrichment tab uses either the **Differential Expression results
table** generated in-app or a **user-uploaded CSV** (DEA-style: `gene`,
`avg_log2FC`, `p_val_adj`). It offers:

- **Two analysis modes**: a *single gene list* (up-, down-, or both-regulated)
  via `clusterProfiler::enrichGO`, or a *comparison* of up- vs down-regulated
  genes via `clusterProfiler::compareCluster`.
- **Standard visualisation**: a clusterProfiler enrichment dot plot.
- **Semantic clustering** with [rrvgo](https://www.bioconductor.org/packages/rrvgo/):
  treemap and scatter plot that collapse redundant GO terms into
  interpretable parent groups.
- **Results table**: the full enrichGO output, downloadable as CSV.

## Repository structure

```
AtlasLens/
├── app.R               # The AtlasLens Shiny application
├── Dockerfile         
├── environment.yml     
├── install.R          
├── README.md
└── .gitignore
```

## Requirements

- **R 4.3.x** (validated on R 4.3.3) with **Bioconductor 3.18**.
- **Operating system:** Linux or macOS recommended. 
- **Memory:** depends on Atlas size.
- **RAM**: depends on Atlas size.
- **System libraries** (Linux): `libcurl`, `libxml2`, `libssl`, `libpng`,
  `libtiff`, `libhdf5`, `libglpk`, `libgsl`, `libbz2`, and a C / C++ / Fortran
  toolchain. These are provided automatically by the Docker image and by
  the conda environment.

## Installation

Choose **one** of the three routes below.

### Option A: Docker (recommended)

```bash
# From the repository root
docker build -t atlaslens .
docker run --rm -p 3838:3838 \
    -v /absolute/path/to/your_dataset.rds:/data/dataset.rds:ro \
    -e DATASET_PATH=/data/dataset.rds \
    atlaslens
```

Then open `http://localhost:3838`. The `Dockerfile` pins R 4.3.1 +
Bioconductor 3.18 and installs every dependency, so the container is
self-contained.

### Option B: YAML 


```bash
conda env create -f environment.yml
conda activate atlaslens_env
```


### Option C: portable R installer

For any machine with a working R 4.3.x toolchain:

```bash
Rscript install.R
```

## Dataset and configuration

AtlasLens expects **one integrated Seurat object** saved as an `.rds` file.
Gene symbols (or Ensembl IDs  converted on-the-fly via biomaRt) are the row
names; cell-level metadata columns are used for grouping and filtering. The
app supports both Seurat v3 / v4 `Assay` and Seurat v5 `Assay5` objects;
multi-layer v5 objects are joined automatically at startup.

The dataset path is resolved at startup from the `DATASET_PATH` environment
variable:

```bash
export DATASET_PATH=/path/to/your_integrated_object.rds
```

A cache directory for analysis results is created at `~/tmp` (configurable
via `CACHE_DIR`).

## Running AtlasLens

```bash
# Activate the environment first (conda activate atlaslens_env, if using conda)
export DATASET_PATH=/path/to/your_integrated_object.rds

R -e 'shiny::runApp("app.R", host = "0.0.0.0", port = 3838, launch.browser = FALSE)'
```

Then open `http://localhost:3838` (or the server address) in a browser.

## Example workflow

1. **Load**: open the app; the Introduction tab confirms the dataset is
   ready and the Dataset Overview card lists the available metadata columns.
2. **Explore**: in *Dataset Exploration*, colour the UMAP by any metadata
   column, overlay gene expression, and apply a shared multi-dimensional
   cell filter that propagates to every other tab.
3. **Time-series**: in *Time Series Information*, pick a
   condition, cell type and timepoints, then inspect the mean-expression
   heatmap. Genes are k-means-clustered by temporal pattern; each cluster is
   labelled with its peak timepoint. Use **Send to DEA** to carry a cluster
   directly into differential expression.
4. **Differential expression**: in *Differential Expression*, choose a
   metadata column and two groups, then run the comparison. Optionally
   restrict the volcano plot and results table to an uploaded gene list.
5. **GO enrichment**: in *GO Enrichment*, use the DEA results directly (or
   upload a CSV), choose single-list or comparison mode, and run. Review the
   dot plot and the rrvgo treemap / scatter.
6. **Analyse gene function**: in *Gene Function (COCOA)*, pick a gene,
   define the groups to compare, and run the co-regulation analysis.
7. **History**: revisit any previous COCOA / DEA / GO run from the
   *History* tab.

Long-running steps display a progress overlay; results are cached, so
repeating an analysis with the same settings is instant.


## Troubleshooting

- **`DATASET_PATH` not found.** Confirm the `.rds` file exists and the
  variable is exported in the same shell that launches R (or passed to
  `docker run` via `-e DATASET_PATH=...`).
- **Analyses run slowly on Windows.** The `multicore` plan is Unix-only; use
  Linux / macOS, or expect sequential fallback.
- **Ensembl-to-symbol conversion stalls.** biomaRt queries Ensembl over the
  network. Confirm the host has outbound HTTPS access; on a restricted
  server, pre-convert IDs offline and reload the object.

## Citation


[Seurat](https://satijalab.org/seurat/),
[geneCOCOA](https://github.com/si-ze/geneCOCOA),
[clusterProfiler](https://bioconductor.org/packages/clusterProfiler/),
[rrvgo](https://bioconductor.org/packages/rrvgo/),
[presto](https://github.com/immunogenomics/presto).


## License

MIT
