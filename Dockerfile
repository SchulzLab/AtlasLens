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
    git \
    && rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('remotes', 'devtools', 'BiocManager', 'Seurat', 'msigdbr', 'shinycssloaders', 'viridis', 'ggrepel', 'digest', 'future', 'promises', 'DT', 'shinyjs', 'qs', 'RColorBrewer', 'plotly', 'patchwork', 'treemap', 'pheatmap'))"

RUN R -e "BiocManager::install('gemma.R')"

# GO Enrichment + biomaRt stack (Bioconductor 3.18 matches R 4.3.x).
# clusterProfiler / rrvgo / GOSemSim / enrichplot / DOSE / AnnotationDbi : GO over-representation + semantic-similarity reduction
# biomaRt                                                                 : Ensembl ID -> gene symbol conversion (geneCOCOA + GO tabs)
# GO.db / org.Hs.eg.db / org.Mm.eg.db                                     : GO term + human / mouse gene annotation
RUN R -e "BiocManager::install(c('clusterProfiler', 'rrvgo', 'GOSemSim', 'enrichplot', 'DOSE', 'AnnotationDbi', 'biomaRt', 'GO.db', 'org.Hs.eg.db', 'org.Mm.eg.db'), version = '3.18', update = FALSE, ask = FALSE)"

RUN R -e "remotes::install_github('si-ze/geneCOCOA', upgrade = 'never')"

# Install Presto for fast Wilcoxon tests
RUN R -e "remotes::install_github('immunogenomics/presto')"




# Create app directory
RUN mkdir -p /srv/shiny-server/app

# Copy the app source code
COPY app.R /srv/shiny-server/app/app.R

# Ensure the cache directory structure exists and is writable by the 'shiny' user
# The app uses ~/tmp by default, but having this ensures permissions are correct if
# the internal cache fallback is triggered.
RUN mkdir -p /srv/shiny-server/app/cache && chown -R shiny:shiny /srv/shiny-server/app/cache

# Expose the Shiny port
EXPOSE 3838

# Start the application
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app/app.R', host = '0.0.0.0', port = 3838)"]
