
import pandas as pd
from glob import glob
import re

lipidomics_data_paths = "Batch*txt"

batch_list = []
df_list = []

for path in glob(lipidomics_data_paths):
    batch = re.search('Batch\s(\d{1,2})', path).group(1)
    df = pd.read_csv(path, sep="\t", index_col="Identification")
    df_list.append(df.T)
    num_samples = df.shape[1]
    batch_list.extend([batch for i in range(num_samples)])
    print(path, df.shape)

combined_df = pd.concat(df_list, ignore_index=False)
#combined_df = combined_df.T # combined_df.shape -> (782, 3631)
combined_df['batch'] = batch_list
# drop rows with NAs (should be just the last one)
combined_df.dropna(inplace=True) # combined_df.shape -> (775, 3631)

# write out the combined file
combined_output_path = "combined_lipidomics_data.tsv"
combined_df.to_csv(combined_output_path, sep="\t")
