
# 1. Get lipid batch info
# 2. Get knockout info metadata

# load modules
import pandas as pd

# load metabolomics data
metabolomics_path = "QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20191101.csv"
metabolomics_df = pd.read_csv(metabolomics_path, index_col='Feature_ID')
# drop the first 4 columns (RT, Pubchem, HMDB, Class) and transpose
metabolomics_df = metabolomics_df[metabolomics_df.columns[4:]].T

# get batch info from sample name
batch_list = []
for sample, row in metabolomics_df.iterrows():
    sample_id = sample.split("_")[-1]
    batch = "_".join(sample.split("_")[:-1])
    batch_list.append(batch)
    print(sample_id,batch)
# add batch info to the dataframe
metabolomics_df['batch'] = batch_list

# load sample info
metabolite_samples_path = "../proteomics/H3K_Project_Master_Lists.xlsx"
metabolite_samples_df = pd.read_excel(metabolite_samples_path, sheet_name="Metabolite Samples")
metabolite_samples_df.head()

# generate dict to translate sample id to cell line name
cell_line_dict = {}
for index, row_j in metabolite_samples_df.iterrows():
    # e.g., "A193M"
    jarred_sample_id = row_j['Sample ID']
    h3k_cell_line = row_j['H3K Cell Line Name']

    cell_line_dict[jarred_sample_id[1:]] = h3k_cell_line

# for each sample translate to cell line name
cell_line_list = []
for sample, row in metabolomics_df.iterrows():
    # e.g., "193M"
    sample_id = sample.split("_")[-1]

    if not "PWT" in sample:
        cell_line = cell_line_dict[sample_id]
    else:
        # e.g., "Batch9and10_Day3_PWT"
        cell_line = "_".join(sample.split("_")[-2:])

    cell_line_list.append(cell_line)
    print(cell_line)

# add batch info to the dataframe
metabolomics_df['H3K_cell_line'] = cell_line_list

# write out the file for ComBat normalization
metabolomics_df.to_csv("metabolomics_data_w_batch.csv")
