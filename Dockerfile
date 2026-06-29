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

# Redefine R's "C++11" compiler flag to gnu++17. fgsea (a dependency of
# clusterProfiler via DOSE/enrichplot) declares SystemRequirements: C++11, so R
# compiles it with -std=gnu++11. The BH (Boost) headers in this image now require
# C++14+, so that compile fails with "Boost.Math requires C++14" /
# "'is_final' has not been declared in 'std'". Because BiocManager::install()
# does NOT return a non-zero exit code on a failed package, the build would
# otherwise continue and ship an image silently missing clusterProfiler. Setting
# CXX11STD here makes any C++11-requesting package compile under C++17.
#
# Also raise R's download timeout from its 60s default to 600s, for ALL R
# invocations from here on (Rprofile.site is read by every `R -e` below). The
# Bioconductor annotation packages (esp. org.Hs.eg.db) and some source tarballs
# are large enough to exceed 60s on a normal connection; when a download times
# out, BiocManager::install() silently skips that package and the build ships
# without it (org.Hs.eg.db / enrichplot / clusterProfiler went missing this way).
RUN echo 'CXX11STD = -std=gnu++17' >> /usr/local/lib/R/etc/Makevars.site \
    && echo 'options(timeout = 600)' >> /usr/local/lib/R/etc/Rprofile.site

# Install ggplot2 3.5.2 BEFORE the GO stack. ggtree 3.10.1 (a dependency of
# enrichplot -> clusterProfiler) calls ggplot2::check_linewidth, which only
# exists in ggplot2 >= 3.5.0; the frozen 2023-10-30 CRAN snapshot ships ggplot2
# 3.4.4, so ggtree fails to byte-compile ("object 'check_linewidth' not found")
# and takes enrichplot + clusterProfiler down with it. 3.5.2 is the version
# AtlasLens is validated on (and is re-asserted with shiny lower down, where it
# is now a no-op). update=FALSE in the GO step below preserves this version.
RUN R -e "remotes::install_version('ggplot2', '3.5.2', repos='https://cloud.r-project.org', upgrade='never')"

# GO Enrichment + biomaRt stack (Bioconductor 3.18 matches R 4.3.x).
# clusterProfiler / rrvgo / GOSemSim / enrichplot / DOSE / AnnotationDbi : GO over-representation + semantic-similarity reduction
# biomaRt                                                                 : Ensembl ID -> gene symbol conversion (geneCOCOA + GO tabs)
# GO.db / org.Hs.eg.db / org.Mm.eg.db / org.Dr.eg.db                      : GO term + human / mouse / zebrafish gene annotation
#
# This is done as a self-healing + self-diagnosing R script rather than a bare
# BiocManager::install() because that call returns exit code 0 even when a
# package fails, so the build would otherwise ship an image whose GO Enrichment
# tab dies at runtime with "there is no package called 'clusterProfiler'". The
# script:
#   1. installs the stack with update=FALSE (protects the pinned Seurat/Matrix);
#   2. LOAD-tests each package and prints the REAL error (e.g. a failed shared
#      object or a "namespace ... is required" version clash) - not just "missing";
#   3. retries the failures with update=TRUE. enrichplot/clusterProfiler do NOT
#      depend on Seurat, and Matrix 1.6-5 is already the newest build for R 4.3,
#      so allowing dependency updates here cannot disturb the pinned stack. This
#      fixes the common case where a too-old pre-installed dependency blocked
#      enrichplot under update=FALSE;
#   4. stops the build (non-zero exit) if anything still won't load.
RUN R --no-save <<'EOF'
pkgs <- c('clusterProfiler','rrvgo','GOSemSim','enrichplot','DOSE','AnnotationDbi',
          'biomaRt','GO.db','org.Hs.eg.db','org.Mm.eg.db','org.Dr.eg.db')
BiocManager::install(pkgs, version = '3.18', update = FALSE, ask = FALSE)

# Load-test each package; print the actual error for any that fail to load.
loads <- function(p) tryCatch({ loadNamespace(p); TRUE },
  error = function(e) { message('LOAD FAIL [', p, ']: ', conditionMessage(e)); FALSE })

bad <- pkgs[!vapply(pkgs, loads, logical(1))]
if (length(bad)) {
  message('Retrying with dependency updates allowed: ', paste(bad, collapse = ', '))
  BiocManager::install(bad, version = '3.18', update = TRUE, ask = FALSE)
  bad <- pkgs[!vapply(pkgs, loads, logical(1))]
}
if (length(bad))
  stop('GO/enrichment packages still failing after retry: ', paste(bad, collapse = ', '),
       ' (see the LOAD FAIL lines above for the underlying cause)')
cat('GO/enrichment stack OK\n')
EOF

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

# Pin Shiny + patchwork to the versions mutually-compatible with ggplot2 3.5.2
# (installed above, before the GO stack); these match the validated conda
# environment.yml (shiny 1.11.1, patchwork 1.3.2). The base image's frozen
# 2023-10-30 snapshot otherwise leaves them too old next to the modern ggplot2:
#  - too-old Shiny can't read ggplot2's panel ranges, so the click-and-drag brush
#    on the UMAPs returns NORMALISED [0,1] coords and zoom-to-brush selects the
#    wrong region (rstudio/shiny#1420);
#  - patchwork 1.1.3 predates ggplot2 3.5 support, so combining plots (the GO
#    Enrichment scatter, which uses patchwork::plot_layout) dies in add_guides
#    with "'==' only defined for equally-sized data frames".
# The assertion re-confirms all three versions made it into the final image.
RUN R -e "for (pv in list(c('shiny','1.11.1'), c('patchwork','1.3.2'))) remotes::install_version(pv[1], pv[2], repos='https://cloud.r-project.org', upgrade='never')"
RUN R -e "stopifnot(packageVersion('shiny') == '1.11.1', packageVersion('ggplot2') == '3.5.2', packageVersion('patchwork') == '1.3.2'); cat('Shiny/ggplot2/patchwork pinned OK\n')"

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
