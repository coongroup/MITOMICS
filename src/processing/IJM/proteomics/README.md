## Script, data, and pipeline details

### 1.0-proteomics_processing.py
- description: Code to merge searches and QC filter by percentile rank of pWT measurements.
- input:
    - proteinGroups_1st_search.txt (1.6 Gb)
    - proteinGroups_2nd_search.txt (1.6 Gb)
    - sample_batch_info.tsv
- output:
    - combined_pWT_LFQ_df_filtered.csv
    - combined_LFQ_df_filtered.csv
- notes: output serves as input for imputation code in /src/processing/DRB