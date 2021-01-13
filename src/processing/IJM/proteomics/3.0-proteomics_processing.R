
# load libraries
library(sva)
library(data.table)

# set wd
setwd("~/Desktop/UW_2020/H3K/MITOMICS/src/processing/IJM/proteomics/")

# load data
LFQ_w_imputed_data_df <- read.table("combined_LFQ_filtered_imputed_for_combat.tsv", sep="\t", header = T)
LFQ_batch_df <- read.table("sample_w_pWT_batch_info.tsv", sep="\t", header=T)

# drop protein group names, then transpose 
LFQ_w_imputed_data_df_t <- as.data.frame(t(LFQ_w_imputed_data_df[,2:ncol(LFQ_w_imputed_data_df)]))
colnames(LFQ_w_imputed_data_df_t) <- LFQ_w_imputed_data_df$Majority.protein.ID

# order by batch
ordered_batch <- c('1', '2', '3', '4', '6', '7', '8', '9',
                   '10', '11', '12', '13', '14', '15', 'pWT')
LFQ_batch_df$batch <- factor(LFQ_batch_df$batch, levels=ordered_batch)

#### comBAT #### 

# set batch
batch = LFQ_batch_df$ES_batch

# covert to matrix
dat_matrix = as.matrix(log2(LFQ_w_imputed_data_df_t))

# set rownames
rownames(dat_matrix) = rownames(LFQ_w_imputed_data_df_t)

# transpose and run ComBat
combat_dat_matrix = ComBat(dat=t(dat_matrix), batch=batch)

# tranpose back
combat_dat_matrix_t = t(combat_dat_matrix)

#### write out ComBat corrected values ####

# check that protein IDs still line up
table(rownames(combat_dat_matrix) == LFQ_w_imputed_data_df$Majority.protein.ID)

# convert to data table
combat_dat_table <- as.data.table(combat_dat_matrix)

# assign rownames
rownames(combat_dat_table) <-rownames(combat_dat_matrix)

# write out file
fwrite(file="combined_LFQ_filtered_imputed_combat.tsv", x=combat_dat_table, sep="\t", row.names = T)
