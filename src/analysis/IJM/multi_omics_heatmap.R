
# load modules
library(pheatmap)
library(dplyr)
library(readxl)
library(data.table) # for transpose of dfs
library(RColorBrewer)

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
       
#### Plot Heatmap ####

# load MXP data from Jarred's file
mxp_df <- read_excel("src/analysis/IJM/mitomics_target_data.xlsx", sheet = "HAP1 KO cell lines")

# check overlap of KO names
setdiff(names(combined_df), mxp_df$`H3K Cell Line Name (updated 200929)`)

# get df by KO name and MXP status
mxp_ko_df <- mxp_df %>%
  select("H3K Cell Line Name (updated 200929)", "MXP")

# get mat_col
mat_col <- data.frame(Status = mxp_ko_df$MXP)
rownames(mat_col) <- mxp_ko_df$`H3K Cell Line Name (updated 200929)`

# drop molecule id and type columns
#> colnames(combined_df[,ncol(combined_df)])
#[1] "type"

# > colnames(combined_df[,1])
# [1] "Molecule ID"

# create subset for testing
n = nrow(combined_df)
#n = 100
#combined_df_sample <- sample_n(combined_df[,2:204], n)
combined_df_sample <- sample_n(combined_df, n)

# set up colors for rows
mat_row <- data.frame(Molecule = combined_df_sample$type)
row.names(mat_row) <- row.names(combined_df_sample)

# Protein: R0 G57 B71 -> #003947
# lipid: R44 G150 B106 -> #2c966a
# Metabolite: R227 G202 B74 -> #e3ca4a
# Sentinel: R35 G31 B32 -> #231f20
# MXP: R241 G45 B49 -> #f12d30

# set up color map
mat_colors <- list(Molecule = c("#003947", "#2c966a", "#e3ca4a"),
                   Status = c("#231f20", "#f12d30"))
names(mat_colors$Molecule) <- c("Protein", "Lipid", "Metabolite")
names(mat_colors$Status) <- c("Sentinal", "MXP")

# set input to matrix
mat <- as.matrix(combined_df_sample[,2:204])
row.names(mat) <- row.names(combined_df_sample)

# create heatmap
out <- pheatmap(
  mat               = mat,
  annotation_row    = mat_row,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  border_color      = NA,
  fontsize          = 14,
  show_colnames     = F,
  show_rownames     = F,
  main              = "",
  color             = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                        "RdBu")))(101),
  breaks            = seq(-3, 3, length.out = 100)
)

# get dendrogram order
col_order_df <- as.data.frame(colnames(mat[,out$tree_col[["order"]]]))
colnames(col_order_df) <- c("H3K Cell Line")
#write.table(col_order_df, "heatmap_col_order.tsv", sep="\t")