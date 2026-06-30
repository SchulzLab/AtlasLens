# FAQ

## Do I need to register or create an account to use AtlasLens?

No. AtlasLens runs locally and requires no registration.

## Can I try AtlasLens without installing it?

Yes. An online demo using the Tabula Muris dataset is available at
[atlaslens.uni-frankfurt.de](https://atlaslens.uni-frankfurt.de/). To explore your
own data, install AtlasLens locally (see [Getting started](getting-started.md)).

## Is AtlasLens available offline? Is my data uploaded anywhere?

AtlasLens is designed to run **locally** — your dataset stays on your own machine
and is never uploaded to a remote server, which preserves data privacy. The online
demo above is provided only for testing with the public Tabula Muris dataset.

## What input data does AtlasLens need?

One integrated **Seurat object** saved as an `.rds` file. If your data is in
another format, convert it first: use `anndata_to_seurat.R` for an AnnData
(`.h5ad`) file or `build_seurat_from_files.R` for raw matrix/metadata files. See
[Preparing a Seurat object](getting-started.md#preparing-a-seurat-object).

## My dataset uses Ensembl IDs — can I still run geneCOCOA?

Yes. geneCOCOA needs gene symbols, so AtlasLens converts the Ensembl IDs for you
the first time you open the **Gene Function (COCOA)** tab. After conversion you can
switch between symbols and Ensembl IDs from the Introduction page. See the
[Gene Function (COCOA)](step-by-step.md#gene-function-cocoa) guide.

## The Time Series tab shows no timepoints, or picks the wrong column. What do I do?

Declare the timepoint column explicitly in `landing_config.json` via
`timeseries_column` — a declared column always wins over auto-detection. See
[Time Series column mapping](getting-started.md#time-series-column-mapping-optional).

## How do I make my analyses reproducible?

Every DEA, GO and geneCOCOA results panel has a **Download R script** button that
exports a self-contained script reproducing that exact analysis. The **History**
tab also logs every run during a session so you can reopen it later.

## An analysis runs slowly or crashes — what can I do?

See the [Troubleshooting](getting-started.md#troubleshooting) section. On some
systems, setting `ATLASLENS_FUTURE_PLAN=sequential` runs analyses on the main
thread, which is the most portable mode.

## How do I cite AtlasLens?

See the [References](references.md) page for citation details, the datasets used,
and the underlying tools.

## How do I contact support?

Email [ashrafiyan@med.uni-frankfurt.de](mailto:ashrafiyan@med.uni-frankfurt.de).

!!! question "Didn't find your answer?"
    Reach out at [ashrafiyan@med.uni-frankfurt.de](mailto:ashrafiyan@med.uni-frankfurt.de).
