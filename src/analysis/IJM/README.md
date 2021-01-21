## Script, data, and processing pipeline details

### dynamic_range_boxplots.R
- description: Code to plot multi-omics dynamic range box plots
- input:
    - tier3_combined.xlsx
- output:
    - multiomics dynamic range boxplots
    
### multi_omics_heatmap.R
- description: Code to plot multi-omics heatmap
- input:
    - tier3_combined.xlsx
    - mitomics_target_data.xlsx
- output:
    - multiomics heatmap
    
## Additional info on software version, dependencies, and required data files

### Dependencies for R version 3.5.1 (2018-07-02) -- "Feather Spray" running under: macOS  10.15.7
    - sva_3.30.1
    - data.table_1.12.8
    - RColorBrewer_1.1-2
    - readxl_1.3.1
    - dplyr_0.8.5
    - pheatmap_1.0.12
    - ggplot2_3.3.1

## Required data files to rerun processing
    - tier3_combined.xlsx
    - mitomics_target_data.xlsx