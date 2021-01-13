## Script, data, and processing pipeline details

### 1.0-proteomics_processing.py
- description: Code to merge searches and QC filter by percentile rank of pWT measurements.
- input:
    - proteinGroups_1st_search.txt (1.6 Gb)
    - proteinGroups_2nd_search.txt (1.6 Gb)
    - sample_batch_info.tsv
- output:
    - combined_pWT_LFQ_df_filtered.csv
    - combined_LFQ_df_filtered.csv
- notes: output serves as input for imputation code in /src/processing/DRB.


### 2.0-proteomics_processing.py
- description: Code to combine pWT KNN and KO left-censored imputed values.
- input:
    - combined_pWT_LFQ_filtered_imputed.csv
    - combined_LFQ_filtered_imputed.csv
- output:
    - combined_LFQ_filtered_imputed_for_combat.tsv
- notes: Input from DRB imputation pipeline.


### 3.0-proteomics_processing.R
- description: Code for batch correction via ComBat.
- input:
    - combined_LFQ_filtered_imputed_for_combat.tsv
    - sample_w_pWT_batch_info.tsv
- output:
    - combined_LFQ_filtered_imputed_combat.tsv
- notes: 


### 4.0-proteomics_processing.R
- description: Code to standardize KO names and drop unused growth WT samples.
- input:
    - combined_LFQ_filtered_imputed_combat.tsv
    - H3K_Project_Master_Lists.xlsx
- output:
    - proteomics_tier2.csv
- notes: 


## Additional info on software version, dependencies, and required data files

### Dependencies for Python 3.7.6 running under: macOS  10.15.7
    - pandas==1.0.1
    - re (builtin)
    - numpy==1.18.1
    - collections (builtin)
    - tqdm==4.42.1
    - scipy==1.4.1

### Dependencies for R version 3.5.1 (2018-07-02) -- "Feather Spray" running under: macOS  10.15.7
    - sva_3.30.1
    - data.table_1.12.8

## Required data files to rerun processing
    - proteinGroups_1st_search.txt (1.6 Gb) - download via MassIVE
    - proteinGroups_2nd_search.txt (1.6 Gb) - download via MassIVE
    - sample_batch_info.tsv
    - sample_w_pWT_batch_info.tsv
    - combined_LFQ_filtered_imputed.csv (output from imputation pipeline)
    - combined_pWT_LFQ_filtered_imputed.csv (output from imputation pipeline)
    - H3K_Project_Master_Lists.xlsx
    