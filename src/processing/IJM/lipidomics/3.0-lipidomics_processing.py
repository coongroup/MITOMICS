
# load modules
import pandas as pd

# load data
combined_collapsed_path = "combined_lipidomics_data_collapsed.tsv"
combined_collapsed_df = pd.read_csv(combined_collapsed_path, sep="\t", index_col="replicate")

# load redundant CoQ details
Coq10_species_path = "coq_species_to_drop.xlsx"
Coq10_species_df = pd.read_excel(Coq10_species_path)

# get list of identifiers
lipid_drop_list = Coq10_species_df['Identifier'].tolist()

# write out data for Yuriy 6/4/20
combined_collapsed_df.drop(lipid_drop_list, axis=1).to_csv("combined_lipidomics_data_filtered.tsv", sep="\t")