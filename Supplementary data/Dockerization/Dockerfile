FROM bioconductor/bioconductor_docker:RELEASE_3_22

RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyjs', 'rintrojs', 'DT', 'waiter', 'RColorBrewer', 'ggplot2', 'ggvenn', 'ggrepel', 'remotes', 'future', 'httpuv', 'devtools', 'dplyr', 'purrr', 'readxl', 'stringr', 'tidyverse', 'cowplot', 'gridExtra', 'igraph', 'plyr', 'reshape2', 'scales', 'wrapr'), repos='http://cran.rstudio.com')"

RUN R -e "BiocManager::install(c('iSEE', 'iSEEu', 'iSEEde', 'SummarizedExperiment', 'SingleCellExperiment', 'DeeDeeExperiment', 'edgeR', 'karyoploteR', 'circlize', 'ComplexHeatmap', 'GO.db', 'maftools', 'AnnotationDbi', 'biomaRt', 'Biostrings', 'clusterProfiler', 'DOSE', 'enrichplot', 'fgsea', 'ggtree', 'GOSemSim', 'limma', 'org.Hs.eg.db', 'org.Mm.eg.db', 'ReactomePA', 'treeio'), ask=FALSE)"

# Instalar paquetes VHIOinformatics
RUN R -e "remotes::install_github('VHIOinformatics/BasicPlots')"
RUN R -e "remotes::install_github('VHIOinformatics/BasicFunctions')"
RUN R -e "remotes::install_github('VHIOinformatics/fromparse2onco@maria')"

# Copiar apps
COPY VVO_rnaseq /app/VVO_rnaseq/
COPY VVO_genomics /app/VVO_genomics/
COPY landing /app/landing/

EXPOSE 8180 8181 8000

COPY start.sh /app/start.sh
RUN chmod +x /app/start.sh

CMD ["/app/start.sh"]