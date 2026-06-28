# Use the rocker/shiny-verse image (matches R 4.3.x requirement)
FROM rocker/shiny-verse:4.3.1

# Install system dependencies required for Seurat, devtools, and other bio-packages
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libglpk-dev \
    libxml2-dev \
    zlib1g-dev \
    cmake \
    libgsl-dev \
    libbz2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    git \
    && rm -rf /var/lib/apt/lists/*


# Install remotes + BiocManager first so the exact versions below can be pinned.
RUN R -e "install.packages(c('remotes', 'devtools', 'BiocManager'))"

# Pin the Seurat stack to the versions AtlasLens was validated on (these match
# the conda environment.yml used on the cluster: Seurat 5.3.0 / SeuratObject
# 5.2.0 / Matrix 1.6-5). The rocker base image's frozen CRAN snapshot otherwise
# installs Seurat v4, which - loaded next to SeuratObject v5 - silently enters a
# broken "v3/v4 compatibility mode" and mangles large objects on load
# ("replacement has 21025 rows, data has 110824"). Matrix is pinned first
# because SeuratObject and Seurat compile against its ABI; Matrix 1.6-5 is also
# the newest release that supports R 4.3.x (1.7.x needs R >= 4.4). cloud.r-project.org
# is used so these exact versions resolve from CRAN (incl. the archive).
RUN R -e "remotes::install_version('Matrix', version = '1.6-5', repos = 'https://cloud.r-project.org', upgrade = 'never')"
RUN R -e "remotes::install_version('SeuratObject', version = '5.2.0', repos = 'https://cloud.r-project.org', upgrade = 'never')"
RUN R -e "remotes::install_version('Seurat', version = '5.3.0', repos = 'https://cloud.r-project.org', upgrade = 'never')"

# Sanity check: fail the build loudly if the wrong major version slipped in.
RUN R -e "stopifnot(packageVersion('Seurat') == '5.3.0', packageVersion('SeuratObject') == '5.2.0'); cat('Seurat stack OK:', as.character(packageVersion('Seurat')), '/', as.character(packageVersion('SeuratObject')), '\n')"

# Remaining CRAN packages (Seurat is pinned above).
RUN R -e "install.packages(c('msigdbr', 'shinycssloaders', 'viridis', 'ggrepel', 'digest', 'future', 'promises', 'DT', 'shinyjs', 'qs', 'RColorBrewer', 'plotly', 'patchwork', 'treemap', 'pheatmap', 'ggrastr'))"

RUN R -e "BiocManager::install('gemma.R')"

# GO Enrichment + biomaRt stack (Bioconductor 3.18 matches R 4.3.x).
# clusterProfiler / rrvgo / GOSemSim / enrichplot / DOSE / AnnotationDbi : GO over-representation + semantic-similarity reduction
# biomaRt                                                                 : Ensembl ID -> gene symbol conversion (geneCOCOA + GO tabs)
# GO.db / org.Hs.eg.db / org.Mm.eg.db / org.Dr.eg.db                      : GO term + human / mouse / zebrafish gene annotation
RUN R -e "BiocManager::install(c('clusterProfiler', 'rrvgo', 'GOSemSim', 'enrichplot', 'DOSE', 'AnnotationDbi', 'biomaRt', 'GO.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dr.eg.db'), version = '3.18', update = FALSE, ask = FALSE)"

RUN R -e "remotes::install_github('si-ze/geneCOCOA', upgrade = 'never')"

# Install Presto for fast Wilcoxon tests
RUN R -e "remotes::install_github('immunogenomics/presto')"




# Re-pin the future/promises ecosystem to the versions AtlasLens was validated
# with (conda environment.yml). The frozen base-image CRAN snapshot installs a
# 2023-era `future` that SeuratObject 5.2.0 cannot load ("object
# 'FutureInterruptError' is not exported by 'namespace:future'"). Done here,
# AFTER the heavy installs, so it overrides the stale copies WITHOUT invalidating
# their build cache; dependencies precede the packages that need them.
RUN R -e "for (pv in list(c('globals','0.18.0'), c('parallelly','1.45.1'), c('listenv','0.9.1'), c('later','1.4.4'), c('future','1.67.0'), c('future.apply','1.20.0'), c('promises','1.3.3'))) remotes::install_version(pv[1], pv[2], repos='https://cloud.r-project.org', upgrade='never')"

# Pin Shiny + ggplot2 to the mutually-compatible pair used in environment.yml.
# The base image / dependency resolution otherwise leaves a too-old Shiny next to
# a modern ggplot2. When that happens Shiny cannot read ggplot2's panel ranges,
# so the click-and-drag brush on the UMAPs returns NORMALISED [0,1] coordinates
# instead of data coordinates - zoom-to-brush then selects the wrong region
# (see rstudio/shiny#1420). Done AFTER the heavy installs so it overrides the
# stale copies without invalidating their build cache.
RUN R -e "for (pv in list(c('ggplot2','3.5.2'), c('shiny','1.11.1'))) remotes::install_version(pv[1], pv[2], repos='https://cloud.r-project.org', upgrade='never')"
RUN R -e "stopifnot(packageVersion('shiny') == '1.11.1', packageVersion('ggplot2') == '3.5.2'); cat('Shiny/ggplot2 pinned OK\n')"

# Create app directory
RUN mkdir -p /srv/shiny-server/app

# Copy the app source code and its landing-page / species / Time Series column
# configuration. landing_config.json is optional at runtime, but baking it in
# lets the deployer ship intro text, an explicit species, and Time Series column
# overrides without a separate mount.
COPY app.R /srv/shiny-server/app/app.R
COPY landing_config.json /srv/shiny-server/app/landing_config.json

# Ensure the cache directory structure exists and is writable by the 'shiny' user
# The app uses ~/tmp by default, but having this ensures permissions are correct if
# the internal cache fallback is triggered.
RUN mkdir -p /srv/shiny-server/app/cache && chown -R shiny:shiny /srv/shiny-server/app/cache

# Expose the Shiny port
EXPOSE 3838

# Start the application
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app/app.R', host = '0.0.0.0', port = 3838)"]
