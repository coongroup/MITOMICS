
# load modules
import pandas as pd
import re
import numpy as np
from collections import Counter

# load proteomics data
proteomics_path = "combined_LFQ_filtered_imputed_combat.tsv"
proteomics_df = pd.read_csv(proteomics_path, sep="\t", index_col=0)

# load in proteomics key from master list
proteomics_key_path = "H3K_Project_Master_Lists.xlsx"
proteomics_key_df = pd.read_excel(proteomics_key_path, sheet_name='Protein Samples')

# first need to create a list of cell lines with duplicate KOs (i.e., with "*-2")
duplicate_ko_list = []

for ko in list(set(proteomics_key_df['H3K Cell Line Name'])):
    
    # check for "-2"
    if ko.split("-")[1] == "2":
        
        # add the prefix to the list
        duplicate_ko_list.append(ko.split("-")[0])
        
# populate a dict to map sample id to cell line name (in proper format: ADCK2-KO1 vs ADCK2-KO2)
id_to_cell_line_dict = {}

# iterate through proteomics key
for index, row in proteomics_key_df.iterrows():
    
    # get cell line
    cell_line = row['H3K Cell Line Name']
    
    # get sample id
    sample_id = row['Sample ID']
    
    # pull out sample as int / string (to align with proteomics format of)
    sample_pattern = "[A-Z](\d+)P"
    sample_id = re.search(sample_pattern, sample_id).group(1)
    sample_id = str(int(sample_id))
    
    # rename if not WT
    if not re.match("^WT-", cell_line):
    
        # standardize name (i.e., 'ADCK2-KO2' -> 'ADCK2-KO2')
        cell_line_split = cell_line.split("-")
        if cell_line_split[0] in duplicate_ko_list:

            cell_line = cell_line_split[0] + "-KO" + cell_line_split[1] # ADCK2 + -KO + 2

        # else just add "-KO"
        else:

            cell_line = cell_line_split[0] + "-KO"
            
    else:
        cell_line = "WT"

    # populate dictionary
    # check that it's consistent if repeated
    if not sample_id in id_to_cell_line_dict:
        id_to_cell_line_dict[sample_id] = cell_line
    else:
        # will drop all WT samples anyway, so preserving precise mapping should not matter here
        if not cell_line == "WT":
            assert(id_to_cell_line_dict[sample_id] == cell_line)
        
    # print the cell line name with the new name
    #print("{}\t{}".format(row['H3K Cell Line Name'], cell_line))
    
# get pWT cols
pWT_cols = [col for col in proteomics_df.columns if "pWT" in col]

# create subset without pWT cols
proteomics_subset_df = proteomics_df.drop(pWT_cols, axis=1)

# create an update colname list
updated_colnames = []

# iterate through colnames
for column in proteomics_subset_df.columns:
    
    # skip pWT samples (will average those)
    if not "pWT" in column:
    
        # parse out sample id (might have to handle pWT separately)
        sample_id_pattern = "(\d+)_"
        sample_id = re.search(sample_id_pattern, column).group(1)

        # parse out replicate id
        replicate_pattern = "_rep([ABCD])"
        replicate = re.search(replicate_pattern, column).group(1)
        
        # get updated cell line
        cell_line = id_to_cell_line_dict[sample_id]
        
        # add new colnames
        updated_colname = cell_line + "-" + replicate
        
        # add to updated_colnames list
        updated_colnames.append(updated_colname)
        
# rename columns
proteomics_subset_df.columns = updated_colnames

# now drop growth WT cols
gWT_cols = [col for col in proteomics_subset_df.columns if "WT-" in col]
proteomics_subset_df = proteomics_subset_df.drop(gWT_cols, axis=1)

# now get average values for pWTs
pWT_dict = {}

# iterate through dataframe
for protein, row in proteomics_df[pWT_cols].iterrows():
    
    # add protein to dict
    pWT_dict[protein] = {}
    
    # iterate through pWT samples
    for column in proteomics_df[pWT_cols].columns:

        # get replicate from name
        # parse out replicate id
        replicate_pattern = "_rep([ABCD])"
        replicate = re.search(replicate_pattern, column).group(1)

        # get LFQ value
        LFQ_value = int(round(2**row[column]))

        # add to dictionary
        if not replicate in pWT_dict[protein]:
            pWT_dict[protein][replicate] = [LFQ_value]
            
        else:
            pWT_dict[protein][replicate].append(LFQ_value)


# now calculate averages for each replicate
pWT_replicate_mean_dict = {}

for protein, replicate_dict in pWT_dict.items():
    
    # set upper level
    pWT_replicate_mean_dict[protein] = {}
    
    # get average for each replicate
    for replicate, LFQ_list in replicate_dict.items():
        
        # set replicate to WT-Replicate
        replicate = "WT-" + replicate
        
        # get log2 of mean intensity value
        pWT_replicate_mean_dict[protein][replicate] = np.log2(np.mean(LFQ_list))
        
# convert to df
pWT_replicate_mean_df = pd.DataFrame(pWT_replicate_mean_dict).T

# check the alignment of indices
Counter(pWT_replicate_mean_df.index == proteomics_subset_df.index)

# rejoin df by index
proteomics_subset_df = proteomics_subset_df.join(pWT_replicate_mean_df)

# write out data
outpath = "proteomics_tier2.csv"
proteomics_df = proteomics_subset_df.to_csv(outpath, sep=",")