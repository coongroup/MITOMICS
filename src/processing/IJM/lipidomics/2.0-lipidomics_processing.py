
# collapse inconsistently quantified lipid classes (e.g., cardiolipinds)
# add new class assignments from Evgenia

import pandas as pd
import numpy as np

# load new class info from Evgenia
lipid_class_info_path = "collapse_data.tsv"
lipid_class_info_df = pd.read_csv(lipid_class_info_path, sep="\t")

# load combined lipidomics dataset (downloaded from data viewer)
combined_lipidomics_path = "combined_lipidomics_data.tsv"
combined_lipidomics_df = pd.read_csv(combined_lipidomics_path, sep="\t",
                                        index_col=0)
combined_lipidomics_df.index.name = "replicate"

# updated lipid class list
updated_lipid_class_list = []

for index, row in lipid_class_info_df.iterrows():
    new_lipid_class = row['new_lipid_class']
    lipid_class = row['lipid_class']

    if pd.notnull(new_lipid_class):
        lipid_class = new_lipid_class

    updated_lipid_class_list.append(lipid_class)

lipid_class_info_df['updated_lipid_class'] = updated_lipid_class_list
# write the lipid class out to a reference file
new_lipid_class_reference_path = "reference/lipid_class.csv"
use_cols = ['lipid', 'lipid_class', 'updated_lipid_class']
# # NOTE: may to make list non-redunant for easy mapping in R
#lipid_class_info_df[use_cols].to_csv(new_lipid_class_reference_path)

# get a dictionary with cardiolipins, etc. to collapse
collapse_info_dict = {}

for index, row in lipid_class_info_df.iterrows():
    new_lipid = row['new_lipid']
    lipid = row['lipid']

    if pd.notnull(new_lipid):
        if not new_lipid in collapse_info_dict:
            collapse_info_dict[new_lipid] = [lipid]
        else:
            collapse_info_dict[new_lipid].append(lipid)

# now iterate through combined_lipidomics_df and collapse relevant lipids
collpased_lipidomics_dict = {}
# record which columns (individual lipids to drop)
drop_cols = []

for new_lipid, old_lipid_list in collapse_info_dict.items():
    new_lipid_sum_list = []

    for index, row in combined_lipidomics_df[old_lipid_list].iterrows():
        new_lipid_sum = np.log2(sum([2**row[old_lipid] for old_lipid in old_lipid_list]))
        new_lipid_sum_list.append(new_lipid_sum)

        # add the old lipid names to be dropped
        drop_cols.extend(old_lipid_list)

    combined_lipidomics_df[new_lipid] = new_lipid_sum_list

combined_lipidomics_df.drop(columns=drop_cols, inplace=True)
# combined_lipidomics_df.shape -> (775, 3579) (-1 for 'batch' = 3578)

# reorder columns alphabetically
combined_lipidomics_df = combined_lipidomics_df[combined_lipidomics_df.columns.sort_values()]

# output the lipidomics data
combined_lipidomics_path = "combined_lipidomics_data_collapsed.tsv"
combined_lipidomics_df = combined_lipidomics_df.to_csv(combined_lipidomics_path, sep="\t")
