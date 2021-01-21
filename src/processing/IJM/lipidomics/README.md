## Script, data, and processing pipeline details

# combine lipidomics data
# collapse lipidomcis data
# filter CoQ species data
# standardize names

### 1.0-lipidomics_processing.py
- description: Code to combine lipidomics batches.
- input:
    - Batch*txt
- output:
    - combined_lipidomics_data.tsv
- notes: 

### 2.0-lipidomics_processing.py
- description: Code to collapse complementary cardiolipin species.
- input:
    - collapse_data.tsv
    - combined_lipidomics_data.tsv
- output:
    - combined_lipidomics_data_collapsed.tsv
- notes: 

### 3.0-lipidomics_processing.py
- description: Code to filter out redundant CoQ species.
- input:
    - combined_lipidomics_data_collapsed.tsv
    - coq_species_to_drop.xlsx
- output:
    - combined_lipidomics_data_filtered.tsv
- notes:

### 4.0-lipidomics_processing.py
- description: Code to standardize KO sample names and drop growth WTs
- input:
    - combined_lipidomics_data_filtered.tsv
    - ../proteomics/H3K_Project_Master_Lists.xlsx
- output:
    - lipidomics_tier2.csv
- notes: 


## Additional info on software version, dependencies, and required data files

### Dependencies for Python 3.7.6 running under: macOS  10.15.7
    - pandas==1.0.1
    - re (builtin)
    - numpy==1.18.1
    - collections (builtin)
    - tqdm==4.42.1
    - scipy==1.4.1
    - glob2==0.7
    
## Required data files to rerun processing
    - Batch 1 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 2 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 3 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 4 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 5 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 6 - Lipidomics (Branch_ Lipidomics).txt
    - Batch 7 - Lipidomics (Branch_ Lipidomics).txt
    - collapse_data.tsv
    - coq_species_to_drop.xlsx
    - H3K_Project_Master_Lists.xlsx