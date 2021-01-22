# Analysis steps to compute log2-fold-changes, q-values, and run tSNE

## Before running

* Ensure that all [required software](Software requirements) is installed.
* Ensure that a `results` folder exists at the root level of the project directory structure.
* Ensure that a `data` folder exists at the root level of the project directory structure, containing the following files:
    * `combined_LFQ_filtered_imputed_combat.tsv` - output from proteomics processing
    * `metabolomics_data_w_batch.csv` - output from metabolomics processing
    * `combined_lipidomics_data_filtered.tsv` - output from lipidomics processing
    * `20200514_H3K_proteomics_avg_log2_FC_w_metadata.csv`
    * `molecule_data_201001_JWR.xlsx`
    * `protein_metadata.csv`

Note: the expected locations of the `results` and `data` folder may be changed by editing the respective definitions of the `results_folder` and `data_folder` variables in `make.R`.

## Running

This analysis subpipeline can be invoked by running

```
Rscript make.R
```

Specific targets in the analysis subpipeline can be run by passing their names as subsequent arguments to the same command. `plan.R` defines the target names. Consult the [documentation of the `drake` package](https://books.ropensci.org/drake/) for finer control over execution.

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

## Required data files
Paths are relative to project root

* `src/analysis/YS/metadata/sample_control_matches.txt` - present in git repository
* `data/combined_LFQ_filtered_imputed_combat.tsv` - output from proteomics processing
* `data/metabolomics_data_w_batch.csv` - output from metabolomics processing
* `data/combined_lipidomics_data_filtered.tsv` - output from lipidomics processing
* `data/20200514_H3K_proteomics_avg_log2_FC_w_metadata.csv`
* `data/molecule_data_201001_JWR.xlsx`
* `data/protein_metadata.csv`

