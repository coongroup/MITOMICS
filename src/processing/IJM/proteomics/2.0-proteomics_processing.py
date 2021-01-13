import pandas as pd
import numpy as np
import re

### Combine KNN imputations from pWT with left-tail imputation for knockouts ###

# get KNN data for pWTs, note these are log2 values
pWT_KNN_path = ("combined_pWT_LFQ_filtered_imputed.csv")
pWT_KNN_df = pd.read_csv(pWT_KNN_path, index_col="Majority protein ID")

# convert back to LFQ (from log2 values)
def power2(log2_val):
    return 2**log2_val

# transpose and get LFQ values (from log2 values)
pWT_KNN_df = pWT_KNN_df.apply(power2)

# load low-end imputed data for knockouts
combined_LFQ_path = "combined_LFQ_filtered_imputed.csv"
combined_LFQ_df = pd.read_csv(combined_LFQ_path, index_col="Majority protein ID")
# filter by protein list used for pWT KNN imputation
combined_LFQ_df = combined_LFQ_df.loc[pWT_KNN_df.index]
print(combined_LFQ_df.shape)
combined_LFQ_df.head() # shape = (7409, 772)

# update LFQ values for pWTs
for column in combined_LFQ_df.columns:
    if "pWT" in column:
        combined_LFQ_df[column] = pWT_KNN_df[column]
        
# sort columns alphabetically, so that their batch info can be more easily matched in sample_w_pWT_batch_info.tsv
combined_LFQ_df.sort_index(axis=1, inplace=True)

# write out for combat correction
output_path = "combined_LFQ_filtered_imputed_for_combat.tsv"
combined_LFQ_df.to_csv(output_path, sep="\t")