# References

## How to cite AtlasLens

AtlasLens is not yet formally published — a preprint will be available on
bioRxiv soon. In the meantime, please reference the source code and the archived
release:

- Source code: [github.com/SchulzLab/AtlasLens](https://github.com/SchulzLab/AtlasLens)
- Archived release & processed datasets (Zenodo):
  [doi.org/10.5281/zenodo.20666624](https://doi.org/10.5281/zenodo.20666625)

## Datasets

AtlasLens is demonstrated using two publicly available single-cell datasets. The
processed Seurat objects for both are deposited on Zenodo
([zenodo.org/records/20666624](https://zenodo.org/records/20666624)).

- **Tabula Muris** — a mouse multi-organ single-cell atlas (23 tissues, 21,025
  genes, 110,824 cells). Tabula Muris Consortium. Single-cell transcriptomics of
  20 mouse organs creates a Tabula Muris. *Nature*, 2018;562:367–372.
  [doi:10.1038/s41586-018-0590-4](https://doi.org/10.1038/s41586-018-0590-4)
- **Whole-lung bleomycin injury time course** — a time-resolved whole-lung
  single-cell atlas of bleomycin-induced lung injury and fibrosis (23,400 genes,
  29,297 cells across six timepoints). Strunz M, *et al.* Alveolar regeneration
  through a Krt8+ transitional stem cell state that persists in human lung
  fibrosis. *Nature Communications*, 2020;11:3559.
  [doi:10.1038/s41467-020-17358-3](https://doi.org/10.1038/s41467-020-17358-3)

## Tools and packages

AtlasLens builds on the following open-source tools — please cite them as
appropriate when reporting analyses produced with AtlasLens:

- **[Seurat](https://satijalab.org/seurat/)** — single-cell data structure,
  visualisation, and `FindMarkers` differential expression. Stuart T, *et al.*
  Comprehensive integration of single-cell data. *Cell*, 2019.
  [doi:10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
- **[presto](https://github.com/immunogenomics/presto)** — fast Wilcoxon rank-sum
  test used to accelerate differential expression.
- **[geneCOCOA](https://github.com/si-ze/geneCOCOA)** — context-dependent gene
  function analysis. Zehr S, *et al.* GeneCOCOA: detecting context-specific
  functions of individual genes using co-expression data. *PLOS Computational
  Biology*, 2025;21(3):e1012278.
  [doi:10.1371/journal.pcbi.1012278](https://doi.org/10.1371/journal.pcbi.1012278)
- **[clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)** — GO
  over-representation analysis (`enrichGO` / `compareCluster`). Yu G, *et al.*
  clusterProfiler: an R package for comparing biological themes among gene
  clusters. *OMICS*, 2012;16(5):284–287.
  [doi:10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118)
- **[rrvgo](https://www.bioconductor.org/packages/rrvgo/)** — redundancy reduction
  and summarisation of GO terms (treemap / scatter). Pérez-Silva JG, *et al.*
  rrvgo: a Bioconductor package for interpreting lists of Gene Ontology terms.
  *Bioinformatics*, 2021;37(19):3166–3167.
  [doi:10.1093/bioinformatics/btab123](https://doi.org/10.1093/bioinformatics/btab123)
- **[biomaRt](https://bioconductor.org/packages/biomaRt/)** — Ensembl ID ↔ gene
  symbol conversion (via [Ensembl BioMart](https://www.ensembl.org/info/data/biomart/index.html)).
- **[anndataR](https://github.com/scverse/anndataR)** — native R reader for
  AnnData (`.h5ad`) files, used by `anndata_to_seurat.R`.
