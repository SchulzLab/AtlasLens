# Getting started

AtlasLens is a single-file R/Shiny application (`app.R`) for reproducible
exploration and analysis of single-cell RNA-seq atlases. Heavy computations (differential
expression, geneCOCOA, GO enrichment) run **asynchronously** in background
workers and results are cached on disk for instant retrieval, so the interface
stays responsive. The app auto-detects metadata role columns (timepoint,
condition, cell type, dataset), so it can be pointed at any compatible atlas.

This page covers everything you need to install it, prepare a dataset, and launch
the app.

## Requirements

- **R 4.3.x** (validated on R 4.3.3) with **Bioconductor 3.18**.
- **Operating system:** Linux or macOS recommended.
- **Memory / RAM:** depends on the size of your atlas.
- **System libraries** (Linux): `libcurl`, `libxml2`, `libssl`, `libpng`,
  `libtiff`, `libhdf5`, `libglpk`, `libgsl`, `libbz2`, and a C / C++ / Fortran
  toolchain. These are provided automatically by the Docker image and by the
  conda environment.

## Repository structure

```
AtlasLens/
├── app.R                        # The AtlasLens Shiny application
├── anndata_to_seurat.R          # Convert an AnnData (.h5ad) object to a Seurat .rds
├── build_seurat_from_files.R    # Build a Seurat .rds from raw matrix / metadata files
├── Dockerfile
├── environment.yml
├── install.R
├── README.md
└── .gitignore
```

## Installation

Choose **one** of the three routes below.

=== "Docker (recommended)"

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

=== "Conda (environment.yml)"

    ```bash
    conda env create -f environment.yml
    conda activate atlaslens_env
    ```

=== "Portable R installer"

    For any machine with a working R 4.3.x toolchain:

    ```bash
    Rscript install.R
    ```

## Preparing a Seurat object

AtlasLens expects **one integrated [Seurat](https://satijalab.org/seurat/) object** saved as an `.rds` file.
If you do not already have one, two helper scripts in the repository can build
it for you:

- **From an AnnData object** (a `.h5ad` file, e.g. exported from Scanpy / Python)
  — use
  [`anndata_to_seurat.R`](https://github.com/SchulzLab/AtlasLens/blob/main/anndata_to_seurat.R).
  It reads the file with [anndataR](https://github.com/scverse/anndataR)'s native
  R reader (no Python required) and writes a Seurat `.rds` with a **classic (v3)
  assay**, which loads under any Seurat version:

    ```bash
    # In your own R (requires R >= 4.5 and BiocManager::install(c("anndataR","rhdf5")))
    Rscript anndata_to_seurat.R  input.h5ad  output.rds
    ```

    The script copies `obs` to cell metadata, uses `layers['counts']` (or `X`) as
    the expression matrix, and carries over `obsm` embeddings (`X_umap`, `X_pca`);
    it prints the resulting assay class — confirm it says `Assay` (v3).

- **From raw / matrix files** (10x, CSV, …) — use
  [`build_seurat_from_files.R`](https://github.com/SchulzLab/AtlasLens/blob/main/build_seurat_from_files.R)
  to assemble a Seurat object from raw input files, then load the `.rds` the same
  way.

Point `DATASET_PATH` (below) at the resulting `.rds` file.

!!! info "About the object"
    Gene symbols (or Ensembl IDs, converted on-the-fly using the offline
    Bioconductor annotation databases `org.Hs.eg.db` / `org.Mm.eg.db` /
    `org.Dr.eg.db`, with Ensembl [biomaRt](https://bioconductor.org/packages/biomaRt/)
    as a network fallback) are the row names; cell-level metadata columns are used
    for grouping and filtering. Both Seurat v3 / v4 `Assay` and Seurat v5 `Assay5`
    objects are supported; multi-layer v5 objects are joined automatically at
    startup.

## Dataset configuration

The dataset path is resolved at startup from the `DATASET_PATH` environment
variable. Point it at the `.rds` file you prepared above:

```bash
export DATASET_PATH=/path/to/your_integrated_object.rds
```

A cache directory for analysis results is created at `~/tmp` (configurable via
`CACHE_DIR`).

## Configuration file (`landing_config.json`)

Alongside the dataset, AtlasLens reads an optional **`landing_config.json`** file.
Use it to personalise the landing page (an introduction and dataset description
shown on the **Introduction** tab) and, if needed, to tell the **Time Series** tab
which metadata columns to use.

### Time Series column mapping (optional)

The Time Series tab auto-detects the timepoint, condition, cell-type and dataset
columns from your metadata. If auto-detection picks the wrong column (or your
timepoint column has an unusual name), declare it explicitly in
`landing_config.json` — a declared column always wins:

```json
"timeseries_column": "Day",
"condition_column": "",
"celltype_column": "celltype",
"dataset_column": ""
```

Leave a field as `""` to keep auto-detecting it. Names are **case-sensitive** and
must match a column in the object's `meta.data`.

!!! tip
    If the Time Series tab shows no timepoints or picks the wrong column, set
    `timeseries_column` here to the exact name of your timepoint column.

## Running AtlasLens

```bash
# Activate the environment first (conda activate atlaslens_env, if using conda)
export DATASET_PATH=/path/to/your_integrated_object.rds

R -e 'shiny::runApp("app.R", host = "0.0.0.0", port = 3838, launch.browser = FALSE)'
```

Then open `http://localhost:3838` (or the server address) in a browser.

!!! success "Ready"
    Once the app is running, the **Introduction** tab confirms the dataset is
    loaded. Head to the **[Features](features.md)** overview, then follow the
    **[User guide](step-by-step.md)**.

## Troubleshooting

- **`DATASET_PATH` not found.** Confirm the `.rds` file exists and the
  variable is exported in the same shell that launches R (or passed to
  `docker run` via `-e DATASET_PATH=...`).
- **Analyses run slowly on Windows.** The `multicore` plan is Unix-only; use
  Linux / macOS, or expect sequential fallback.
- **Background analyses crash** (e.g. on some macOS setups, a
  `MultisessionFuture ... non-exportable reference (externalptr)` error). Set
  `ATLASLENS_FUTURE_PLAN=sequential` to run analyses on the main thread instead.
  This is the most portable mode and works on every OS; the only trade-off is
  that the interface is blocked while a long analysis runs. Accepted values are
  `sequential`, `multisession`, and `multicore`.
- **Ensembl-to-symbol conversion stalls.** biomaRt queries Ensembl over the
  network. The app uses the offline OrgDb databases first and only falls back to
  biomaRt; if it still stalls, confirm the host has outbound HTTPS access, or
  pre-convert IDs offline and reload the object.
- **`invalid class "Assay5" object: All layers must have a record in the cells
  map`.** The object is a Seurat v5 object written under a different Seurat
  version. Re-save it with a **classic v3 assay** — `anndata_to_seurat.R` does
  this automatically, or in R run `obj[["RNA"]] <- as(obj[["RNA"]], "Assay")`
  and `saveRDS()` again.
- **Time Series tab shows no timepoints, or picks the wrong column.** Declare the
  timepoint column explicitly via `timeseries_column` in `landing_config.json`
  (see [Time Series column mapping](#time-series-column-mapping-optional)).
