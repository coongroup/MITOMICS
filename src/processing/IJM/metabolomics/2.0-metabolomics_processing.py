# load modules
import pandas as pd
import re
import numpy as np
from collections import Counter

# load metabolomics data
metabolomics_path = "metabolomics_data_w_batch.csv"
metabolomics_df = pd.read_csv(metabolomics_path, index_col=0)

# get column list
column_list = metabolomics_df.columns.tolist()

# replace one column name
column_list[7] = "Adenosine-5-monophosphate 3 TMS" # currently "denosine-5-monophosphate 3 TMS_I would call this a sugar phosphate too_TR"

# update columns
metabolomics_df.columns = column_list

# drop batch and cell line columns and transpose
metabolomics_df = metabolomics_df.drop(['batch', 'H3K_cell_line'], axis=1).T

# load in lipidomics key from master list
metabolomics_key_path = "../proteomics/H3K_Project_Master_Lists.xlsx"
metabolomics_key_df = pd.read_excel(metabolomics_key_path, sheet_name='Metabolite Samples')

# first need to create a list of cell lines with duplicate KOs (i.e., with "*-2")
duplicate_ko_list = []

for ko in list(set(metabolomics_key_df['H3K Cell Line Name'])):
    
    # skip BLANK
    if not ko == "BLANK":
    
        # check for "-2"
        if ko.split("-")[1] == "2":

            # add the prefix to the list
            duplicate_ko_list.append(ko.split("-")[0])
            
# populate a dict to map sample id to cell line name (in proper format: ADCK2-KO1 vs ADCK2-KO2)
id_to_cell_line_dict = {}

# iterate through proteomics key
for index, row in metabolomics_key_df.iterrows():
    
    # get cell line
    cell_line = row['H3K Cell Line Name']
    
    # skip BLANK
    if not cell_line == "BLANK":
    
        # get sample id
        sample_id = row['Sample ID']

        # pull out sample as int / string 
        sample_pattern = "[A-Z](\d+)M" 
        sample_id = re.search(sample_pattern, sample_id).group(1)
        sample_id = str(int(sample_id))

        # rename if not WT
        if not re.match("^WT-", cell_line):

            # standardize name (i.e., 'ADCK2-KO2' : 'ADCK2-KO2')
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
                
# pPT samples
pWT_cols = [col for col in metabolomics_df.columns if "PWT" in col]

# create an update colname list
updated_colnames = []

# translate replicate name
replicate_dict = {
    "1":"A",
    "2":"B",
    "3":"C"
}

# translate batch names
batch_dict = {
    "Batch1and2":"BATCH1",
    "Batch3and4":"BATCH2",
    "Batch5and6":"BATCH3",
    "Batch7and8":"BATCH4",
    "Batch9and10":"BATCH5",
    "Batch11and12":"BATCH6",
    "Batch14and15":"BATCH7"
}

# iterate through colnames
for column in metabolomics_df.columns:
    
    # skip pWT samples (will average those)
    if not "PWT" in column:
    
        # parse out sample id (might have to handle pWT separately)
        sample_id_pattern = "_Day[1-3]_(\d+)M" # "Batch11and12_Day1_168M"
        sample_id = str(int(re.search(sample_id_pattern, column).group(1)))

        # parse out replicate id
        replicate_pattern = "_Day([1-3])_\d+M" # "Batch11and12_Day1_168M"
        replicate = re.search(replicate_pattern, column).group(1)
        
        # translate replicate days to letter
        replicate = replicate_dict[replicate]
        
        # get updated cell line
        cell_line = id_to_cell_line_dict[sample_id]
        
        # add new colnames
        updated_colname = cell_line + "-" + replicate
        
        # add to updated_colnames list
        updated_colnames.append(updated_colname)
        
    else:
        
        # get batch number
        batch_pattern = "(Batch\d+and\d+)_Day[1-3]_PWT" # 'Batch11and12_Day2_PWT'
        batch = re.search(batch_pattern, column).group(1)
        
        # replace batch name
        batch = batch_dict[batch]
        
        # get replicate number
        replicate_pattern = "_Day([1-3])_PWT" # 'Batch11and12_Day2_PWT'
        replicate = re.search(replicate_pattern, column).group(1)
        
        # translate replicate days to letter
        replicate = replicate_dict[replicate]
        
        # rename
        updated_column = "PWT-{}-{}".format(batch, replicate)
        updated_colnames.append(updated_column) 
        
# rename columns
metabolomics_df.columns = updated_colnames

# now drop growth WT cols
gWT_cols = [col for col in metabolomics_df.columns if re.match("^WT-", col)]
metabolomics_subset_df = metabolomics_df.drop(gWT_cols, axis=1)

# write out data
outpath = "metabolomics_tier2.csv"
metabolomics_subset_df.to_csv(outpath, sep=",")
