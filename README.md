# VVO – VHIO's Visual Omics

Shiny platform for interactive omics data visualization, powered by iSEE (Bioconductor) and extended with custom S4 panel classes.

## Modules

- **VVO** – Bulk RNA-seq exploration. Requires a `DeeDeeExperiment` `.rds` file.
- **VVO_genomics** – Mutation profiling and oncoplot visualization.

Both applications are accessed through a shared landing page. The Docker image is currently intended for internal use at VHIO; the `Dockerfile` and `start.sh` script are provided in this repository as reference.

## Supplementary Data

Test datasets and application testing documentation used for this Master's Thesis are available in this repository.

## Author

Master's Thesis (TFM) – Maria Suau Oliver
Bioinformatics Unit, Vall d'Hebron Institute of Oncology (VHIO)
