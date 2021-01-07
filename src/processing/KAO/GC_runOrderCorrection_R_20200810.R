library(MASS)
library(sva)

## function for RSD 
rsd <- function(featureRow, index = 1:length(featureRow)){
  sd(2^(featureRow[index]), na.rm = T)/mean(2^(featureRow[index]), na.rm = T) * 100
}

############## Run time correction for H3K ##################

## use time stamp of file "mtime" for doing run-order correction

file_info <- read.csv("P:/SRP_012019_H3K/KAO_GC_quant/file_time_stamp.csv") ## read in

file_info <- file_info[grep("D:/SRP/H3K/AllFiles_1_2019/",file_info$X),] ## keep only releveant files

file_info$X <- sub("D:/SRP/H3K/AllFiles_1_2019/","", file_info$X) ## remove early folder path information

file_info <- file_info[-grep("Hexane|MSTFA|Alkane",file_info$X), ] ## remove hexanes, alkanes, MSFTA files  

fileNames <- sub(".raw", "_LFQ", file_info$X, fixed = T) ## to match gc-quant output, substitute .raw with _LFQ

## fixing a couple file naming issues to match the gcquant output
fileNames[grep("Batch9and10/Day2/PWT_5", fileNames)] <- paste( fileNames[grep("Batch9and10/Day2/PWT_5", fileNames)], ".1", sep = "")
fileNames[grep("Batch14and15/Day2/PWT_7", fileNames)] <- paste( fileNames[grep("Batch14and15/Day2/PWT_7", fileNames)], ".1", sep = "")


file_info <- data.frame(mtime = as.POSIXct(file_info$mtime), as.data.frame(do.call(rbind, strsplit(fileNames, "/")), stringsAsFactors = F)) ## pull out relevante info from file path into dataframe
names(file_info)  <- c("mtime", "Batch", "Day", "Sample")

file_info$BatchDay <- as.factor(paste(file_info$Batch, file_info$Day, sep ="_"))
file_info$group <- sub("A|B|C", "", sapply(strsplit(file_info$Sample,"_"), function(x) x[[1]]))


## write cleaned up file 
#write.csv(file_info, "P:/SRP_012019_H3K/KAO_GC_quant/file_time_stamp_v2.csv")

################### Read in sample quant data ################

df <- read.csv("P:/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB.csv", stringsAsFactors = F)


sampleNames <- names(df)[grep("_LFQ", names(df))]

## Confirm matches
file_info <- file_info[file_info$Sample %in% sampleNames, ]
sampleNames[!sampleNames %in% file_info$Sample]
file_info$Sample[match(sampleNames,file_info$Sample)] == sampleNames
file_info <- file_info[match(sampleNames,file_info$Sample),]
dim(file_info)

write.csv(file_info, "P:/SRP_012019_H3K/KAO_GC_quant/file_time_stamp_v2.csv", row.names = F)

# Create a file where imputed values are NA
df_na <- df[grep("_LFQ", names(df))]
df_tier <- df[grep("_tier", names(df))]
df_na[df_tier == 5] <- NA
metabolite_quality <- rowSums(is.na(df_na))
sample_quality <- colSums(is.na(df_na))

index <- rev(order(metabolite_quality))
index2 <- rev(order(sample_quality))

PercTrue <- table(is.na(df_na))[1]/sum(table(is.na(df_na)))
df_na_trim <- df_na
i=1
j=1
while(PercTrue < 0.95){
  
  df_na_trim_i <- df_na_trim[-index[i],]
  df_na_trim_j <- df_na_trim[,-index2[j]]
  PercTrue_i <- table(is.na(df_na_trim_i))[1]/sum(table(is.na(df_na_trim_i)))
  PercTrue_j <- table(is.na(df_na_trim_j))[1]/sum(table(is.na(df_na_trim_j)))
  if(PercTrue_i > PercTrue_j){
    PercTrue <- PercTrue_i
    df_na_trim <- df_na_trim_i
    i = i + 1
  } else {
    PercTrue <- PercTrue_j
    df_na_trim <- df_na_trim_j
    j = j+1
  }
}
dim(df_na_trim)
metabolite_quality[index[i]] # 338
sample_quality[index2[j]] # 40

dim(df_na) # rows = metabolites (223), # columns = samples (770) 

## Tranform DF to be samples x metabolties 
df_t_na <- as.data.frame(t(df_na))
row.names(df_t_na) == file_info$Sample

## merge runtime and LFQ by sample names 
df_merge_na <- cbind(file_info, df_t_na)


## version of merge w/out NA
df_t <- as.data.frame(t(df[grep("_LFQ", names(df))]))
row.names(df_t) == file_info$Sample
df_merge <- cbind(file_info, df_t)

############### Plotting data by run order #############
library(pheatmap)

pheatmap(t(df_merge[order(df_merge$mtime),-c(1:6) ]), 
         cluster_cols = F, scale = "row",
         breaks = c(-15, -5, seq(-3,3,length.out = 16),5, 15),
         color = colorRampPalette(c(4,"white",2))(20))

## example plots
plot(V18 ~ mtime , data = df_merge)
plot(V3 ~ mtime , data = df_merge_na[df_merge_na$group == "PWT",])

plot(V70 ~ BatchDay, data = df_merge, las = 2)
plot(V3 ~ BatchDay, data = df_merge_na, las=2)
df[3,1]

########### Calculating PWT RSDs prior to correation #########
# PWT = pooled wild type cells (control QC samples)

as.numeric(df_merge$group == "PWT") + 1
# split by Batch-Day, imputed values are NA
df_split_na <- split(df_merge_na,  df_merge_na$BatchDay)
df_split <- split(df_merge, df_merge$BatchDay)

# calculate RSDs by Batch-Day, imputed values are NA
PWT_RSDs_pre_na <- cbind(sapply(df_split_na, function(x) apply(x[x$group == "PWT", -c(1:6)],2,rsd)),apply(df_merge_na[df_merge_na$group == "PWT",-c(1:6)], 2, rsd))
boxplot(PWT_RSDs_pre_na)

PWT_RSDs_pre <- cbind(sapply(df_split, function(x) apply(x[x$group == "PWT", -c(1:6)],2, rsd)), apply(df_merge[df_merge$group == "PWT", -c(1:6)], 2,rsd))
boxplot(PWT_RSDs_pre)


######### Correcting based on run time ###########
df_split_corrected <- df_split #creating new data frame to put corrected data into

names(df_split) 
#day = 21 i = 229 df_split[[9]][,10], test
# normalizing each day for within day linear changes to intensity values and mean of PWT values, each compound at a time.  
for(day in 1:length(df_split)){
    
    ## find numeric columns
    for(i in which(lapply(df_split_na[[day]], class) == "numeric" )){ ## take only numeric columns 
      df_split_na[[day]][["mtime"]] <- log2(rank(df_split_na[[day]][["mtime"]]))
      #df_split_na[[day]][,i] <- 2^df_split_na[[day]][,i]
      
      model <- formula(paste(names(df_split_na[[day]])[i], "mtime", sep = " ~ "))
      
      ## fit linear model to data table w/ NAs to prevent skewing based on imputation
      fit <- tryCatch(rlm(model, data = df_split_na[[day]], weight = as.numeric(group == "PWT")+1),warning = function(w) fit =NULL, error = function(e) fit = NULL) 
      if(!is.null(fit) ){
        intercept <- summary(fit)$coef[1,1]
        slope <- summary(fit)$coef[2,1]
      } else {
        intercept <- mean(df_split[[day]][,i], na.rm = T)
        slope <- 0
      }
        
        # normalize by run-order
        # sampleValue - ((slope * runOrder) + intercept) + mean
      
        df_split_corrected[[day]][,i] <- (df_split[[day]][,i] - ((slope * df_split_na[[day]][["mtime"]]) + intercept)) + mean(df_merge_na[df_merge_na$group == "PWT",i], na.rm = T)
    
    }
}

######## RSDs post correction ############# 
PWT_RSDs_post <- cbind(sapply(df_split_corrected, function(x) apply(x[x$group == "PWT", -c(1:6)],2, rsd)), apply(df_merge[df_merge$group == "PWT", -c(1:6)], 2,rsd))
names(df_split_corrected)

pdf("P:/SRP_012019_H3K/KAO_GC_quant/Pre_post_runOrderCorrection_20191030.pdf")
boxplot(cbind(PWT_RSDs_pre, PWT_RSDs_post), col = c(rep('#F93A38',ncol(PWT_RSDs_pre)),rep( '#AB298D', ncol(PWT_RSDs_post))),
        xaxt = "n", las = 1, ylab = "% RSD by day")

legend("topright", c("Pre run order correction","Post run order correction"), bty = "n", col = c('#F93A38', "#AB298D"), pch = 15)
dev.off()

boxplot(PWT_RSDs_post, las = 2)
PWT_RSDs_post

## Unsplit 
df_merge_correct <- unsplit(df_split_corrected, df_merge$BatchDay)
table(is.na(df_merge_correct))
table(is.na(df_merge))

table(df_merge_correct$Sample == df_merge$Sample)# check match

## example plots
par(mfrow=c(1,2))

pheatmap(t(df_merge_correct[order(df_merge_correct$mtime),-c(1:6) ]), 
         cluster_cols = F, scale = "row",
         breaks = c(-15, -5, seq(-3,3,length.out = 16),5, 15),
         color = colorRampPalette(c(4,"white",2))(20))


plot(V4 ~ mtime , data = df_merge)
plot(V4 ~ mtime , data = df_merge_correct[df_merge_correct$group == "PWT",])


plot(V36 ~ BatchDay, data = df_merge)
plot(V36 ~ BatchDay, data = df_merge_correct)
par(mfrow = c(1,2))
plot(V70 ~ mtime, data = df_split[[10]])
plot(V70 ~ mtime, data = df_split_corrected[[10]])
dev.off()

df[36,1]
names(df_split_corrected)
boxplot(t(df_merge_correct[,-c(1:6)])[,order(df_merge_correct$Batch)], outline = F, border = as.factor(df_merge_correct$Batch)[order(df_merge_correct$Batch)])

########## remove poor features ############

rowSums(df_tier)
hist(rowSums(df_tier[df_merge$group == "PWT"]))
#questionable_metabolites <- which(rowSums(is.na(df_na[df_merge$group == "PWT"])) > 10)
questionable_metabolites <- which(rowSums(df_tier[df_merge$group == "PWT"]) > 150)

table(df_merge$group  == "PWT")
df[questionable_metabolites, 1]

poorQuant_metabolites <- which(rowMeans(PWT_RSDs_post) > 30)
df[poorQuant_metabolites,1]

#41 = lactic acid
#70 = tryptophan
plot(df_merge_correct[df_merge$group == "PWT",70], pch = rep(c(1,19), each = 4), 
     col = as.factor(df_merge$Batch[which(df_merge$group == "PWT")] ))
legend("topright",pch = 19,  col = 1:7, levels(as.factor(df_merge$Batch)))

dev.off()

plot(df_merge[df_merge$group == "PWT",70 ], pch = rep(c(1,19), each = 4), col = as.factor(df_merge$Batch[which(df_merge$group == "PWT")] ))
plot(file_info$mtime[order(file_info$mtime)][1:50], col = as.factor( file_info$group == "PWT")[order(file_info$mtime)][1:50])

metabolites_to_remove <- unique(c(poorQuant_metabolites, questionable_metabolites))
df[metabolites_to_remove,1]

########## Write csv of run order Corrected 
t_df_merge_corrected <- t(df_merge_correct[,-c(1:6)])
colnames(t_df_merge_corrected) <- paste(df_merge_correct$Sample, "_RunOrderCorrected", sep= "")
df_corrected <- cbind(df, t_df_merge_corrected)
df_corrected <- df_corrected[-metabolites_to_remove, ]

write.csv(df_corrected, "P:/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20200810.csv")

########## PCA ############
help(prcomp)
pca_pre <- prcomp((df_merge[,-c(1:6)]))
plot(pca_pre$x, col = as.factor(df_merge$group == "PWT"), pch =19)

pca_post <- prcomp((df_merge_correct[,-c(1:6)][,-metabolites_to_remove]))
plot(pca_post$x, col = as.factor(df_merge$group == "PWT"), pch =19)
dev.off()
pca_post_noRemove <-prcomp((df_merge_correct[,-c(1:6)]))
plot(pca_post_noRemove$x, col = as.factor(df_merge$group == "PWT"), pch = 19)

plot(pca_pre)
plot(pca_post)
pca_post_PWT <- prcomp((df_merge_correct[df_merge_correct$group == "PWT",-c(1:6)][,-metabolites_to_remove]))
plot(pca_post_PWT$x, col = as.factor(df_merge_correct$Batch[df_merge_correct$group == "PWT"]), pch =19)


plot(pca_post$rotation[,1:2])
df_correct[which(pca_post$rotation[,1] > 0.15),1]

pdf("P:/SRP_012019_H3K/KAO_GC_quant/PCA_Pre_post_runOrderCorrection_20200810.pdf")
par(las = 1, mfrow = c(2,2))

plot(pca_pre$x, col = as.factor(df_merge$Batch), pch =19)
legend("topleft", levels(as.factor(df_merge$Batch)), col = 1:7, pch = 19, bty = "n")

plot(pca_post$x, col = as.factor(df_merge_correct$Batch), pch =19)
legend("topleft", levels(as.factor(df_merge_correct$Batch)), col = 1:7, pch = 19, bty = "n")

plot(pca_post_PWT$x, col = as.factor(df_merge_correct$Batch[df_merge_correct$group == "PWT"]), pch =19)
legend("topright", levels(as.factor(df_merge_correct$Batch)), col = 1:7, pch = 19, bty = "n")


plot(pca_post$rotation[,1:2] , pch = 19)
text(pca_post$rotation[,1:2][which(pca_post$rotation[,1] > 0.15),], df_corrected[which(pca_post$rotation[,1] > 0.15),1])
dev.off()

############# Filter duplicates and remove batch 3/4 day 1 #############

file_info  <- read.csv("P:/SRP_012019_H3K/KAO_GC_quant/file_time_stamp_v2.csv")

df_corrected <- read.csv("P:/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20200810.csv", stringsAsFactors = F)

dim(file_info)
dim(df_corrected)

df_values <- df_corrected[,grep("_RunOrderCorrected", names(df_corrected))]
dim(df_values)
df_tier <- df_corrected[,grep("_tier", names(df_corrected))]

sample_duplicates <- which(duplicated(paste(file_info$group,file_info$BatchDay)) & file_info$group != "PWT")
sample_remove <- unique(c(sample_duplicates, which(file_info$BatchDay == "Batch3and4_Day1")))

df_values_unique <- df_values[,-sample_remove]
file_info_unique <- file_info[-sample_remove,]
dim(df_values_unique)

table(file_info_unique$Sample == sub("_RunOrderCorrected","",names(df_values_unique)))

t_df_values <- t(df_values_unique)

values_batchDay_PWT_mean <- aggregate(t_df_values, by = list(file_info_unique$group, file_info_unique$BatchDay), mean)
dim(values_batchDay_PWT_mean)

plot(values_batchDay_PWT_mean[, 9])
df_values_unique_PWT_mean <- t(values_batchDay_PWT_mean[,-c(1:2)])
dim(df_values_unique_PWT_mean)

colnames(df_values_unique_PWT_mean) <- paste(values_batchDay_PWT_mean$Group.2, values_batchDay_PWT_mean$Group.1,  sep = "_")

#final <- cbind(df_corrected[,1:6], df_values_unique_PWT_mean)
#write.csv(final, "P:/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20191101.csv", row.names = F)

final_2 <- cbind(df_corrected[,1:6], df_values_unique)
write.csv(final_2, "P:/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20200810_2.csv", row.names = F)

######### PCA after filtering ##########
pca_pre <- prcomp(t(df_values))
plot(pca_pre$x, col = as.factor(file_info$BatchDay), pch =19)
dev.off()

pca_post <- prcomp(t(final_2[,-c(1:6)]))
plot(pca_post$x, col = as.factor(values_batchDay_PWT_mean$Group.2), pch =19)
legend("topleft", levels(as.factor(values_batchDay_PWT_mean$Group.2)), col = 1:7, pch = 19, bty = "n")
plot(pca_post)
summary(pca_post)
final_2[,2][pca_post$rotation[,1] > -0.05]
final_2[,2][pca_post$rotation[,2] > 0.05]
range(pca_post$rotation[,2])
batch_final <-  as.factor(sub("_Day[c(1,2,3)]", "", values_batchDay_PWT_mean$Group.2))
combat_final <- ComBat(as.matrix(final_2[,-c(1:6)]), batch = as.numeric(batch_final))


pca_combat <- prcomp(t(combat_final))
plot(pca_combat$x, col = batch_final, pch=19, main = "ComBat")
plot(pca_post$x, col = batch_final, pch =19)
