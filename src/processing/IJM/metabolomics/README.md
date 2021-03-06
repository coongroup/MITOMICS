## Script, data, and processing pipeline details

### 1.0-metabolomics_processing.py
- description: Code to standardize sample names.
- input:
    - QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20191101.csv
    - ../proteomics/H3K_Project_Master_Lists.xlsx
- output:
    - metabolomics_data_w_batch.csv
- notes: This data was input for Yuriy's FC (tier3) calculations.

### 2.0-metabolomics_processing.py
- description: Code to standardize sample names.
- input:
    - metabolomics_data_w_batch.csv
    - ../proteomics/H3K_Project_Master_Lists.xlsx
- output:
    - metabolomics_tier2.csv
- notes:


## Additional info on software version, dependencies, and required data files

### Dependencies for Python 3.7.6 running under: macOS  10.15.7
    - pandas==1.0.1
    - re (builtin)
    - numpy==1.18.1
    - collections (builtin)

## Required data files to rerun processing
    - QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20191101.csv
    - H3K_Project_Master_Lists.xlsx
    