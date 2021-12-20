# Analysis steps for running multi-omic analyses

## Before running

* Ensure that all [required software](#software-requirements) is installed.
* Ensure that a `data` folder exists at the root level of the project directory structure, containing the [required data files](#required-data-files) (see below).

Note: the expected locations of the `results` and `data` folder may be changed by editing the respective definitions of the `results_folder` and `data_folder` variables in `make.R`.

## Running

This analysis subpipeline can be invoked by running

```
make.sh
```

(in this working directory).

This simply invokes `targets::tar_make()` in R.
Consult the [documentation of the `targets` package](https://books.ropensci.org/targets) for finer control over execution.

## Output produced

* `results/ksp_outliers.csv` Table of molecule candidates for knockout-specific phenotypes sorted by distance
* `results/tsne/dbscan-selected.html` A visualization of the tSNE embedding with HDBSCAN clusters and their extensions outlined.
* `results/tsne/qardm-tsne.html` A visualization of the tSNE embedding (details in the Methods section)
* `results/tsne/tsne_neighbors.csv` A table listing all neighbors of MXPs on the tSNE plot that fall within one unit in all tSNE runs.

## Software requirements

    * R version 4.0.2
    * R packages:
        * dbscan 1.1.8.1
        * ggplot2 3.3.2
        * here 1.0.1
        * jsonlite 1.6.1
        * magrittr 2.0.1
        * patchwork 1.1.1
        * qvalue 2.20.0
        * readxl 1.3.1
        * Rtsne 0.15
        * targets 0.8.1
        * tidyverse 1.3.0

## Required data files
Paths are relative to project root

* `data/metadata/sample_control_matches.txt` - download from MassIVE
* `data/metadata/mxp_genes.csv` - download from MassIVE
* `data/metadata/protein_metadata.csv` - download from MassIVE
* `data/data tables/Fold changes, p-values, SDs/20210819_Lipidomics_FCsPvaluesSDs.csv`
* `data/data tables/Fold changes, p-values, SDs/20210819_Metabolomics_FCsPvaluesSDs.csv`
* `data/data tables/Fold changes, p-values, SDs/20210819_Proteomics_FCsPvaluesSDs.csv`

