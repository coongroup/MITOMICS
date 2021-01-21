
# load modules
library(ggplot2)
library(readxl)
library(dplyr)

# set working directory
setwd("~/Desktop/UW_2020/H3K/MITOMICS/")

# set data path
data_path = "src/analysis/IJM/tier3_combined.xlsx"

# load proteomics data
proteomics_df <- read_excel(data_path, sheet = "Protein L2FC")

# subset quant values
proteomics_meta_cols <- names(proteomics_df)[2:13]
proteomics_quant_df <- proteomics_df[,-which(names(proteomics_df) %in% proteomics_meta_cols)]

# add molecule type
proteomics_quant_df$type <- "Protein"

# load lipidomics data (doesn't appear to be any metadata)
lipidomics_df <- read_excel(data_path, sheet = "Lipids L2FC")

# add molecule type
lipidomics_df$type <- "Lipid"

# load metabolomics data
metabolomics_df <- read_excel(data_path, sheet = "Metabolites L2FC")

# add molecule type
metabolomics_df$type <- "Metabolite"

# "MICOS/MIB" missing from metabolomics
setdiff(names(proteomics_quant_df), names(metabolomics_df))

# and lipidomics
setdiff(names(proteomics_quant_df), names(lipidomics_df))

# drop "MICOS/MIB" from proteomics_quant_df
proteomics_quant_df <- proteomics_quant_df[,-which(names(proteomics_quant_df) %in% c("MICOS/MIB"))]

# rbind the dataframes together
combined_df <- rbind(proteomics_quant_df, lipidomics_df, metabolomics_df)

# define dynamic range function
dynamic_range <- function(row) {
  row_max <- max(row)
  row_min <- min(row)
  dynamic_range <- row_max - row_min
  return(dynamic_range)
}

# test out with example row

# drop the first (molecule name) and last columns (molecule type)
example_row <- combined_df[1,2:204]
max(example_row)
min(example_row)
max(example_row) - min(example_row)
dynamic_range(example_row)

# create new column
combined_df$dynamic_range <- apply(combined_df[,2:204], 1, dynamic_range)

# make boxplot
ggplot(combined_df, aes(x=type, y=dynamic_range, fill=type)) +
  geom_boxplot(show.legend = F) +
  theme_light() +
  labs(x = "Measurement Type", y = "log2 FC Dynamic Range") +
  scale_fill_manual(values = c("#2c966a", "#e3ca4a", "#003947")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face="bold"))

# get values for boxplot
combined_df %>%
  select(dynamic_range, type) %>%
  group_by(type) %>%
  summarize(
    min = min(dynamic_range),
    max = max(dynamic_range),
    median = median(dynamic_range),
    q1 = quantile(dynamic_range, 0.25),
    q3 = quantile(dynamic_range, 0.75),
    IQR = IQR(dynamic_range)
  )
