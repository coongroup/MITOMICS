# load modules
import pandas as pd
import re
import numpy as np
from collections import Counter

# load proteomics data
lipidomics_path = "combined_lipidomics_data_filtered.tsv"
lipidomics_df = pd.read_csv(lipidomics_path, sep="\t", index_col=0)

# drop batch column and transpose
lipidomics_df = lipidomics_df.drop(['batch'], axis=1).T

# load in lipidomics key from master list
lipidomics_key_path = "../proteomics/H3K_Project_Master_Lists.xlsx"
lipidomics_key_df = pd.read_excel(lipidomics_key_path, sheet_name='Lipid Samples')

# first need to create a list of cell lines with duplicate KOs (i.e., with "*-2")
duplicate_ko_list = []

for ko in list(set(lipidomics_key_df['H3K Cell Line Name'])):
    
    # skip BLANK
    if not ko == "BLANK":
    
        # check for "-2"
        if ko.split("-")[1] == "2":

            # add the prefix to the list
            duplicate_ko_list.append(ko.split("-")[0])
            
# populate a dict to map sample id to cell line name (in proper format: ADCK2-KO1 vs ADCK2-KO2)
id_to_cell_line_dict = {}

# iterate through proteomics key
for index, row in lipidomics_key_df.iterrows():
    
    # get cell line
    cell_line = row['H3K Cell Line Name']
    
    # skip BLANK
    if not cell_line == "BLANK":
    
        # get sample id
        sample_id = row['Sample ID']

        # pull out sample as int / string 
        sample_pattern = "[A-Z](\d+)L"
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
                
# pWT samples
pWT_cols = [col for col in lipidomics_df.columns if "PWT" in col]

# create an update colname list
updated_colnames = []

# iterate through colnames
for column in lipidomics_df.columns:
    
    # skip pWT samples (will average those)
    if not "PWT" in column:
    
        # parse out sample id (might have to handle pWT separately)
        sample_id_pattern = "[A-Z](\d+)L" # "A136L Quant Values"
        sample_id = str(int(re.search(sample_id_pattern, column).group(1)))

        # parse out replicate id
        replicate_pattern = "([ABCDE])\d+L" # "A136L Quant Values"
        replicate = re.search(replicate_pattern, column).group(1)
        
        # get updated cell line
        cell_line = id_to_cell_line_dict[sample_id]
        
        # add new colnames
        updated_colname = cell_line + "-" + replicate
        
        # add to updated_colnames list
        updated_colnames.append(updated_colname)
        
    else:
        
        # get batch number
        batch_pattern = "PWT-(\d+)-\d+"
        batch = re.search(batch_pattern, column).group(1)
        
        # get replicate number
        replicate_pattern = "PWT-\d+-(\d+)"
        replicate = re.search(replicate_pattern, column).group(1)
        
        # rename
        updated_column = "PWT-BATCH{}-{}".format(batch, replicate)
        updated_colnames.append(updated_column) 
        
# create an update colname list
updated_colnames = []

# iterate through colnames
for column in lipidomics_df.columns:
    
    # skip pWT samples (will average those)
    if not "PWT" in column:
    
        # parse out sample id (might have to handle pWT separately)
        sample_id_pattern = "[A-Z](\d+)L" # "A136L Quant Values"
        sample_id = str(int(re.search(sample_id_pattern, column).group(1)))

        # parse out replicate id
        replicate_pattern = "([ABCDE])\d+L" # "A136L Quant Values"
        replicate = re.search(replicate_pattern, column).group(1)
        
        # get updated cell line
        cell_line = id_to_cell_line_dict[sample_id]
        
        # add new colnames
        updated_colname = cell_line + "-" + replicate
        
        # add to updated_colnames list
        updated_colnames.append(updated_colname)
        
    else:
        
        # get batch number
        batch_pattern = "PWT-(\d+)-\d+"
        batch = re.search(batch_pattern, column).group(1)
        
        # get replicate number
        replicate_pattern = "PWT-\d+-(\d+)"
        replicate = re.search(replicate_pattern, column).group(1)
        
# rename
updated_column = "PWT-BATCH{}-{}".format(batch, replicate)
updated_colnames.append(updated_column)

# now drop growth WT cols
gWT_cols = [col for col in lipidomics_df.columns if re.match("^WT-", col)]
lipidomics_subset_df = lipidomics_df.drop(gWT_cols, axis=1)

# write out data
outpath = "lipidomics_tier2.csv"
lipidomics_subset_df.to_csv(outpath, sep=",")



                