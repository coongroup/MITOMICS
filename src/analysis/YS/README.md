# Analysis steps to compute log2-fold-changes, q-values, and run tSNE

## Before running

* Ensure that all [required software](#software-requirements) is installed.
* Ensure that a `results` folder exists at the root level of the project directory structure.
* Ensure that a `data` folder exists at the root level of the project directory structure, containing the [required data files](#required-data-files) (see below).

Note: the expected locations of the `results` and `data` folder may be changed by editing the respective definitions of the `results_folder` and `data_folder` variables in `make.R`.

## Running

This analysis subpipeline can be invoked by running

```
Rscript make.R
```

Specific targets in the analysis subpipeline can be run by passing their names as subsequent arguments to the same command. `plan.R` defines the target names. Consult the [documentation of the `drake` package](https://books.ropensci.org/drake/) for finer control over execution.

## Output produced

* `results/molecule_data.xlsx` An excel workbook containing mean log-fold changes, p-values, and three sets of q-values (in-condition, single-ome, multi-ome).
* `results/tsne_figure.html` A visualization of the tSNE embedding (details in the Methods section)
* `results/tsne_neighbors.csv` A table listing all neighbors of MXPs on the tSNE plot that fall within one unit.

## Software requirements

    * R version 4.0.2
    * R packages:
        * drake 7.12.5
        * here 0.1
        * jsonlite 1.6.1
        * magrittr 1.5
        * qvalue 2.20.0
        * tidyverse 1.3.0
        * tsne 0.1.3
        * writexl 1.3.1

## Required data files
Paths are relative to project root

* `data/metadata/sample_control_matches.txt` - download from MassIVE
* `data/metadata/mxp_genes.csv` - download from MassIVE
* `data/metadata/protein_metadata.csv` - download from MassIVE
* `data/combined_LFQ_filtered_imputed_combat.tsv` - output from proteomics processing
* `data/metabolomics_data_w_batch.csv` - output from metabolomics processing
* `data/combined_lipidomics_data_filtered.tsv` - output from lipidomics processing

