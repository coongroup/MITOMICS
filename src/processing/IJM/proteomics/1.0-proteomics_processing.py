import pandas as pd
import numpy as np
import re
from collections import Counter
from scipy import stats
from tqdm import tqdm

#### 1. Merge search data ####

# load data from proteinGroups files
first_search_path = "proteinGroups_1st_search.txt"
first_search_df = pd.read_csv(first_search_path, sep="\t")

# start with all columns of interest except LFQ for second search
second_search_path = "proteinGroups_2nd_search.txt"
second_search_df = pd.read_csv(second_search_path, sep="\t")

# filter data
def filter_df(proteinGroup_df):
    proteinGroup_df = proteinGroup_df[proteinGroup_df['Reverse'] != '+']
    proteinGroup_df = proteinGroup_df[proteinGroup_df['Only identified by site'] != '+']
    proteinGroup_df = proteinGroup_df[proteinGroup_df['Potential contaminant'] != '+']

    ## TODO: Filter by sparsity
    # get list of LFQ value columns
    LFQ_columns = []
    for column in proteinGroup_df.columns:
        if "LFQ" in column:
            LFQ_columns.append(column)

    return proteinGroup_df

first_search_df = filter_df(first_search_df)
second_search_df = filter_df(second_search_df)

# populate LFQ dict, key by majority protein ID
LFQ_dict = {}

for index,row in first_search_df.iterrows():
    first_search_id = row['id']
    majority_protein_group = row['Majority protein IDs']
    majority_protein = majority_protein_group.split(";")[0]
    gene_names = row['Gene names']

    # pull out LFQ values
    for column in first_search_df.columns:
        if "LFQ" in column:

            LFQ_value = row[column]
            if LFQ_value > 0:

                if not majority_protein in LFQ_dict:
                    LFQ_dict[majority_protein] = {}
                    LFQ_dict[majority_protein]['first_search_id'] = first_search_id
                    LFQ_dict[majority_protein]['first_search_gene_names'] = gene_names
                    LFQ_dict[majority_protein]['first_search_majority_proteins'] = majority_protein_group
                    LFQ_dict[majority_protein][column] = LFQ_value
                    # keep track of non-zero values
                    LFQ_dict[majority_protein]['count_in_first_search'] = 1

                else:
                    LFQ_dict[majority_protein][column] = LFQ_value
                    LFQ_dict[majority_protein]['count_in_first_search'] += 1

# now for second search
for index,row in second_search_df.iterrows():
    second_search_id = row['id']
    majority_protein_group = row['Majority protein IDs']
    majority_protein = majority_protein_group.split(";")[0]
    gene_names = row['Gene names']
    if not pd.isnull(gene_names):
        first_gene_name = gene_names.split(";")[0]

    # pull out LFQ values
    for column in second_search_df.columns:
        if "LFQ" in column:

            LFQ_value = row[column]
            if LFQ_value > 0:

                if not majority_protein in LFQ_dict:
                    LFQ_dict[majority_protein] = {}
                    LFQ_dict[majority_protein]['second_search_id'] = second_search_id
                    LFQ_dict[majority_protein]['second_search_gene_names'] = gene_names
                    LFQ_dict[majority_protein]['second_search_majority_proteins'] = majority_protein_group
                    LFQ_dict[majority_protein][column] = LFQ_value
                    # keep track of non-zero values
                    LFQ_dict[majority_protein]['count_in_second_search'] = 1

                # if the majority protein is in the LFQ_dict, but this is the
                # first sample from the second search
                elif 'second_search_id' not in LFQ_dict[majority_protein]:
                    LFQ_dict[majority_protein]['merged_by'] = 'majority_protein'
                    LFQ_dict[majority_protein]['second_search_id'] = second_search_id
                    LFQ_dict[majority_protein]['second_search_gene_names'] = gene_names
                    LFQ_dict[majority_protein]['second_search_majority_proteins'] = majority_protein_group
                    LFQ_dict[majority_protein][column] = LFQ_value

                    if not 'count_in_second_search' in LFQ_dict[majority_protein]:
                        LFQ_dict[majority_protein]['count_in_second_search'] = 1
                    else:
                        LFQ_dict[majority_protein]['count_in_second_search'] += 1

                else:
                    LFQ_dict[majority_protein][column] = LFQ_value
                    LFQ_dict[majority_protein]['count_in_second_search'] += 1

combined_LFQ_df = pd.DataFrame(LFQ_dict).T
combined_LFQ_df.index.name = "Majority protein ID"
combined_LFQ_df['count_in_first_search'].fillna(0, inplace=True)
combined_LFQ_df['count_in_second_search'].fillna(0, inplace=True)
combined_LFQ_df['total_count'] = combined_LFQ_df['count_in_first_search'] + \
    combined_LFQ_df['count_in_second_search']
combined_LFQ_df.sort_values(ascending=False,inplace=True,
    by=['total_count', 'count_in_first_search', 'count_in_second_search'])

# get gene names
first_search_gene_names_dict = Counter(combined_LFQ_df['first_search_gene_names'])
second_search_gene_names_dict = Counter(combined_LFQ_df['second_search_gene_names'])

# filter for 0 count in first search and >= 200 in second search
combined_LFQ_df_first_subset = combined_LFQ_df[combined_LFQ_df['count_in_second_search'] == 0]
combined_LFQ_df_second_subset = combined_LFQ_df[combined_LFQ_df['count_in_first_search'] == 0]

# iterate through combined_LFQ_df and look for merge opportunities by shared gene name
merge_dict = {}

for majority_protein_i,row_i in tqdm(combined_LFQ_df_first_subset.iterrows(),total=combined_LFQ_df_first_subset.shape[0]):

    first_search_gene_names = row_i['first_search_gene_names']
    if not pd.isnull(first_search_gene_names):

        for majority_protein_j,row_j in combined_LFQ_df_second_subset.iterrows():
            second_search_gene_names = row_j['second_search_gene_names']

            if not pd.isnull(second_search_gene_names):

                if second_search_gene_names.split(";")[0] == \
                    first_search_gene_names.split(";")[0]:

                    # Make sure there is only one protein with this gene name
                    if first_search_gene_names_dict[first_search_gene_names] == 1 \
                        and second_search_gene_names_dict[second_search_gene_names]:

                        merge_dict[majority_protein_i] = majority_protein_j

# repeat for second subset
for majority_protein_i,row_i in tqdm(combined_LFQ_df_second_subset.iterrows(),total=combined_LFQ_df_second_subset.shape[0]):
    second_search_gene_names = row_i['second_search_gene_names']

    if not pd.isnull(second_search_gene_names):

        for majority_protein_j,row_j in combined_LFQ_df_second_subset.iterrows():

            first_search_gene_names = row_j['first_search_gene_names']

            if not pd.isnull(first_search_gene_names):

                if second_search_gene_names.split(";")[0] == \
                    first_search_gene_names.split(";")[0]:

                    # Make sure there is only one protein with this gene name
                    if first_search_gene_names_dict[first_search_gene_names] == 1 \
                        and second_search_gene_names_dict[second_search_gene_names]:

                        merge_dict[majority_protein_i] = majority_protein_j

# drop any rows that values in the merge dictionary, then repopulate
# rows that are keys, using LFQ_dict and merge_dict mapping
combined_LFQ_df.drop(merge_dict.values(), inplace=True)

# in order to overwrite exist data via .loc, need to turned off chained assignment
# https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
pd.options.mode.chained_assignment = None

for majority_protein in merge_dict:
    second_search_majority_protein = merge_dict[majority_protein]
    for key in LFQ_dict[second_search_majority_protein]:
        combined_LFQ_df.loc[majority_protein, key] = LFQ_dict[second_search_majority_protein][key]
    combined_LFQ_df.loc[majority_protein, 'merged_by'] = 'gene_name'

combined_LFQ_df['total_count'] = combined_LFQ_df['count_in_first_search'] + \
    combined_LFQ_df['count_in_second_search']
combined_LFQ_df.sort_values(ascending=False,inplace=True,
    by=['total_count', 'count_in_first_search', 'count_in_second_search',
        'merged_by'])

# get LFQ cols
combined_LFQ_cols = [col for col in combined_LFQ_df.columns if "LFQ" in col]

# filter based on coverage across combined dataset. And record the protein groups that are dropped.
cutoff = round(len(combined_LFQ_cols)/2)

# now populate a dictionary with info regarding the proteins that were dropped
dropped_protein_group_dict = {}
for majority_protein_id, row in combined_LFQ_df[combined_LFQ_df['total_count'] < cutoff].iterrows():
    dropped_protein_group_dict[majority_protein_id] = {}
    
    count_in_first_search = row['count_in_first_search']
    count_in_second_search = row['count_in_second_search']
    total_count = row['total_count']
    merged_by = row['merged_by']
    dropped_protein_group_dict[majority_protein_id]['count_in_first_search'] = count_in_first_search
    dropped_protein_group_dict[majority_protein_id]['count_in_second_search'] = count_in_second_search
    dropped_protein_group_dict[majority_protein_id]['total_count'] = total_count
    dropped_protein_group_dict[majority_protein_id]['merged_by'] = merged_by
    dropped_protein_group_dict[majority_protein_id]['reason_dropped'] = "total_count < {}".format(cutoff)
    
# drop any proteins that don't appear at all in the first search
for majority_protein_id, row in combined_LFQ_df[combined_LFQ_df['count_in_first_search']==0].iterrows():
    
    dropped_protein_group_dict[majority_protein_id] = {}
    
    count_in_first_search = row['count_in_first_search']
    count_in_second_search = row['count_in_second_search']
    total_count = row['total_count']
    merged_by = row['merged_by']
    dropped_protein_group_dict[majority_protein_id]['count_in_first_search'] = count_in_first_search
    dropped_protein_group_dict[majority_protein_id]['count_in_second_search'] = count_in_second_search
    dropped_protein_group_dict[majority_protein_id]['total_count'] = total_count
    dropped_protein_group_dict[majority_protein_id]['merged_by'] = merged_by
    dropped_protein_group_dict[majority_protein_id]['reason_dropped'] = "count_in_first_search == 0".format(cutoff)

### 2. Percentile rank filtering ###

# Record list and metadata of dropped measurements
# For pWT samples, calculate percentile rank of average LFQ value by batch.

# need a map of pWT sample name and batch

batch_info_path = "sample_batch_info.tsv"
batch_info_df = pd.read_csv(batch_info_path, sep="\t")

# get list of pWT sample names
pWT_sample_list = [sample_name for sample_name in batch_info_df['sample'] if "pWT" in sample_name]

# get dictionary with pWT batch info
pWT_sample_dict = {}
for index, row in batch_info_df.iterrows():
    sample_name = row['sample']
    maxquant_batch = row['maxquant_batch']
    batch = row['batch']
    
    if 'pWT' in sample_name:
        pWT_sample_dict[sample_name] = {'maxquant_batch':maxquant_batch,
                          'batch':batch}

# get pWT batch LFQ data
pWT_batch_LFQ_dict = {}

# filter out measurements in less than half of data
majority_filtered_combined_LFQ_df = combined_LFQ_df[combined_LFQ_df['total_count'] > cutoff]

# replace NAs with 0
for majority_protein_id, row in majority_filtered_combined_LFQ_df[pWT_sample_list].replace(np.nan, 0).iterrows():
    
    for pWT_sample, batch_info_dict in pWT_sample_dict.items():
        
        quant_value = row[pWT_sample]
        
        #batch = pWT_sample_dict[pWT_sample]
        # actually batch is baked into pWT sample name
        batch = pWT_sample.split()[-1].split("_")[0]
        
        # collect LFQ data for pWTs by batch
        if not batch in pWT_batch_LFQ_dict:
            pWT_batch_LFQ_dict[batch] = {}
            pWT_batch_LFQ_dict[batch][majority_protein_id] = [quant_value]
            
        elif majority_protein_id not in pWT_batch_LFQ_dict[batch]:
            pWT_batch_LFQ_dict[batch][majority_protein_id] = [quant_value]
            
        else:
            pWT_batch_LFQ_dict[batch][majority_protein_id].append(quant_value)
            
# for each batch, for each protein group, get pWT percentile rank
percentile_rank_dict = {}

for batch, protein_group_quant_dict in tqdm(pWT_batch_LFQ_dict.items(), total=len(pWT_batch_LFQ_dict), desc='outer loop'):
    
    # calculate mean for each for protein
    protein_group_mean_dict = {}
    for protein_group, quant_value_list in protein_group_quant_dict.items():
        
        # filter out zeros from quant_value_list
        quant_value_list = [i for i in quant_value_list if i > 0]
    
        if len(quant_value_list) == 0:
            protein_group_mean = np.nan
        
        else:
            protein_group_mean = np.mean(quant_value_list)
            protein_group_mean_dict[protein_group] = protein_group_mean
        
    # calculate percentile rank
    percentile_rank_dict[batch] = {}
    protein_group_mean_list = list(protein_group_mean_dict.values())
    for protein_group, mean in protein_group_mean_dict.items():
        percentile_rank = round(stats.percentileofscore(protein_group_mean_list, mean), 1)
        percentile_rank_dict[batch][protein_group] = percentile_rank

# convert to df
percentile_rank_df = pd.DataFrame(percentile_rank_dict)
        
# now add range and merge by data
percentile_range_list = []
merged_by_list = []

# create a protein group merge by dict, to prevent nested for loop
merge_by_dict = {}
for majority_protein_id, row in majority_filtered_combined_LFQ_df.iterrows():
    merged_by = row['merged_by']
    merge_by_dict[majority_protein_id] = merged_by
    
for majority_protein_id, row in percentile_rank_df.iterrows():
    percentile_max = max(row.dropna())
    percentile_min = min(row.dropna())
    percentile_rank_range = percentile_max - percentile_min
    percentile_range_list.append(percentile_rank_range)
    
    merged_by = merge_by_dict[majority_protein_id]
    merged_by_list.append(merged_by)
    
# append range and merged by info
percentile_rank_df['range'] = percentile_range_list
percentile_rank_df['merged_by'] = merged_by_list

# populate dictionary from combined_LFQ_df
combined_LFQ_dict = {}
cols_of_interst = ['merged_by', 'count_in_first_search', 'count_in_second_search', 'total_count']

for majority_protein_id, row in combined_LFQ_df.iterrows():
    combined_LFQ_dict[majority_protein_id] = {}
    for col in cols_of_interst:
        combined_LFQ_dict[majority_protein_id][col] = row[col]

# populate dictionary with info regarding the proteins that were dropped
for majority_protein_id_i in percentile_rank_df[percentile_rank_df['merged_by'].isnull()].index:
    first_search_count = combined_LFQ_dict[majority_protein_id_i]['count_in_first_search']
    second_search_count = combined_LFQ_dict[majority_protein_id_i]['count_in_second_search']
    total_count = combined_LFQ_dict[majority_protein_id_i]['total_count']
    merged_by = combined_LFQ_dict[majority_protein_id_i]['merged_by']
    
    if not majority_protein_id_i in dropped_protein_group_dict:
        dropped_protein_group_dict[majority_protein_id_i] = {}
        dropped_protein_group_dict[majority_protein_id_i]['count_in_first_search'] = first_search_count
        dropped_protein_group_dict[majority_protein_id_i]['count_in_second_search'] = second_search_count
        dropped_protein_group_dict[majority_protein_id_i]['total_count'] = total_count
        dropped_protein_group_dict[majority_protein_id_i]['merged_by'] = merged_by
        dropped_protein_group_dict[majority_protein_id_i]['reason_dropped'] = "Missing in first search"
        
    else:
        dropped_protein_group_dict[majority_protein_id_i]['reason_dropped'] += " & Missing in first search"

    # print(majority_protein_id_i, first_search_count, second_search_count)
    
# now filter out any proteins with percentile rank ranges > 50
range_rank_cutoff = 50

rank_range_filter_total = percentile_rank_df[percentile_rank_df['range'] > range_rank_cutoff].shape[0]

for majority_protein_id_i, row_i in tqdm(percentile_rank_df[percentile_rank_df['range'] > range_rank_cutoff].iterrows(), total=rank_range_filter_total):
    percentile_rank_range = round(row_i['range'], 1)
    first_search_count = combined_LFQ_dict[majority_protein_id_i]['count_in_first_search']
    second_search_count = combined_LFQ_dict[majority_protein_id_i]['count_in_second_search']
    total_count = combined_LFQ_dict[majority_protein_id_i]['total_count']
    merged_by = combined_LFQ_dict[majority_protein_id_i]['merged_by']

    if not majority_protein_id_i in dropped_protein_group_dict:
        dropped_protein_group_dict[majority_protein_id_i] = {}
        dropped_protein_group_dict[majority_protein_id_i]['count_in_first_search'] = first_search_count
        dropped_protein_group_dict[majority_protein_id_i]['count_in_second_search'] = second_search_count
        dropped_protein_group_dict[majority_protein_id_i]['total_count'] = total_count
        dropped_protein_group_dict[majority_protein_id_i]['merged_by'] = merged_by
        dropped_protein_group_dict[majority_protein_id_i]['reason_dropped'] = \
            "Percentile rank range > {}".format(range_rank_cutoff)
    else:
        dropped_protein_group_dict[majority_protein_id_i]['reason_dropped'] += \
            " & Percentile rank range > {}".format(range_rank_cutoff)

# covert to df
dropped_protein_group_df = pd.DataFrame(dropped_protein_group_dict).T
        
# write out qc filtering info
#dropped_protein_group_path = "protein_groups_to_drop.csv"
#dropped_protein_group_df.index.name = "majority_protein_id"
#dropped_protein_group_df.to_csv(dropped_protein_group_path)
        
# write out data for imputation

# writeout pWT data (for KNN)
LFQ_cols = [col for col in combined_LFQ_df.columns if "LFQ" in col]
non_LFQ_cols = [col for col in combined_LFQ_df.columns if not "LFQ" in col]
pWT_LFQ_cols = [col for col in combined_LFQ_df.columns if "LFQ" in col and "pWT" in col]
combined_LFQ_filtered_outpath = "combined_pWT_LFQ_df_filtered.csv"
combined_LFQ_df.drop(dropped_protein_group_df.index)[pWT_LFQ_cols].to_csv(combined_LFQ_filtered_outpath)

# write out full data
LFQ_cols = [col for col in combined_LFQ_df.columns if "LFQ" in col]
non_LFQ_cols = [col for col in combined_LFQ_df.columns if not "LFQ" in col]
pWT_LFQ_cols = [col for col in combined_LFQ_df.columns if "LFQ" in col and "pWT" in col]
combined_LFQ_filtered_outpath = "combined_LFQ_df_filtered.csv"
combined_LFQ_df.drop(dropped_protein_group_df.index)[LFQ_cols].to_csv(combined_LFQ_filtered_outpath)