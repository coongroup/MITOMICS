## comparing GC vs. LC ## 

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

#### Load in GC and LC data ####

setwd("P:/ES_20170625_H3K/SRP_012019_H3K/")
df_gc <- read.csv("P:/ES_20170625_H3K/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20200810_2.csv", stringsAsFactors = F)

tf <- read.csv("P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_TIC_GROWTH.csv", stringsAsFactors = F)

## metabolite info table ## 

metabolite_info <- sub("X", "", names(tf)[-c(1:4)])

## expression matrix ##
exprs <- tf[,-c(1:4)]

## sample info table ##
names(tf)
sample_info_tf <- data.frame(names = tf$File_Name, file_names = tf$File_Name, label = tf$Label, TIC = tf$TIC, growth = tf$Growth)
sample_info_tf$date <- apply(sample_info_tf, 1, function(x) strsplit(x[2], '_')[[1]][[1]])
sample_info_tf$sampleNumber <- apply(sample_info_tf, 1, function(x) strsplit(x[2], '_')[[1]][2])
sample_info_tf$sampleID <- apply(sample_info_tf, 1, function(x) strsplit(x[2], '_')[[1]][3])
sample_info_tf$group <- apply(sample_info_tf, 1, function(x) strsplit(x[8], "-", fixed = T)[[1]][[1]])
sample_info_tf$group[sample_info_tf$group == 1] <- "Control"
sample_info_tf$group[sample_info_tf$group == 2] <- "Control"
sample_info_tf$is_sample <- grepl("[0-9]{3}M",sample_info_tf$group)
sample_info_tf$gc_name <- apply(sample_info_tf, 1, function(x) paste(strsplit(x[8], "-", fixed = T)[[1]][2],strsplit(x[8], "-", fixed = T)[[1]][[1]], sep = ""))
sample_info_tf$gc_name <- sub("NAControl", "PWT", sample_info_tf$gc_name)
names(sample_info_tf)

##
gc_name <- sub("_LFQ_RunOrderCorrected","", names(df_gc), fixed = T)
gc_name <- sub("_1901........", "", gc_name)
gc_name <- sub("_._._.", "", gc_name)
gc_name <- sub("_LFQ.1_RunOrderCorrected","" , gc_name)


gc_exprs <- aggregate(t(df_gc[,-c(1:6)]), by = list(gc_name[-c(1:6)]), mean)
gc_exprs$Group.1 
names(gc_exprs)[-1] <- df_gc$ï..Feature.ID

## impute 0 values from bottom 5 percent within compound ##
# 
# imputeZeros <- function(exprs, percent = 0.05) {
#         newExprs <- exprs  
#         for(i in names(exprs)){
#                 lowPerc <- quantile(exprs[[i]][exprs[[i]] != 0], prob = percent)
#                 index <- which(exprs[[i]] == 0) 
#                 newExprs[[i]][index] <- runif(length(exprs[[i]][exprs[[i]] == 0]), min(exprs[[i]]), lowPerc)
#         }
#         newExprs
# }
# set.seed(123)
# exprs_impute <- imputeZeros(exprs)

### 
boxplot(log2(exprs))

png("E:/MITOMICS/plots/No_batch_normalization.png")
boxplot(t(log2(exprs)), border = as.numeric(as.factor(sample_info_tf$date)),
        main = "No batch correction")
dev.off()
x <- boxplot(t(exprs), col = as.numeric(sample_info_tf$date))

df_exprs <- aggregate(log2(exprs), by = list(sample_info_tf$gc_name), mean)
df_exprs$Group.1
names(df_exprs)[-1] <- metabolite_info

## Function for strain RSD calculation 
strainRSD <- function(matrix, strain_info){
  # matrix is samples in columns, features in rows 
  # strain info is vector with length == to nrow(matrix)
  strain_sd <- aggregate(matrix, by = list(strain_info), sd)
  strain_mean <- aggregate(matrix, by = list(strain_info), mean)
  rsd <- data.frame(names = strain_mean[,1], strain_sd[,-1]/strain_mean[,-1] *100)
  rsd
}

noNorm_rsd <- strainRSD(exprs[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])

boxplot(rowMeans(noNorm_rsd[,-1]))

names(noNorm_rsd[,-1])[noNorm_rsd[noNorm_rsd$names == "Control",-1] > 30 & colMeans(noNorm_rsd[,-1]) > 30]

#### Batch Normaliztion #### 

## Option 1: Combat normalization ## 

library(sva)
exprs_combat <- ComBat(t(log2(exprs)), batch = as.factor(sample_info_tf$date))
df_combat <- aggregate(t(exprs_combat), by = list(sample_info_tf$gc_name), mean)
names(df_combat)[-1] <- metabolite_info

png("E:/MITOMICS/plots/Combat_normalization.png")
boxplot(exprs_combat, border = as.numeric(as.factor(sample_info_tf$date)), 
        main = "Combat batch normalization")
dev.off()

Combat_rsd <- strainRSD(t(2^exprs_combat)[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
boxplot(rowMeans(Combat_rsd[,-1]))

## Option 2: Batch Median normalization ## 

batch_median <- aggregate(t(x$stats), by = list(sample_info_tf$date), mean)
ratio <- batch_median[,3]/mean(batch_median[,3])
sample_info_tf$batchMedian <- ratio[as.numeric(as.factor(sample_info_tf$date))]

exprs_medianNorm <- exprs/sample_info_tf$batchMedian
df_medianNorm <- aggregate(log2(exprs_medianNorm), by = list(sample_info_tf$gc_name), mean)
names(df_medianNorm)[-1] <- metabolite_info

png("E:/MITOMICS/plots/Batch_median_normalization.png")
boxplot(t(log2(exprs_medianNorm)), border = as.numeric(as.factor(sample_info_tf$date)),
        main = "Batch median batch normalization" )
dev.off()

medianNorm_rsd <- strainRSD(exprs_medianNorm[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
boxplot(rowMeans(medianNorm_rsd[,-1]))

#### Normalizing for sample extraction amount #### 
# Option 0.1: TIC normalizization 
TIC_ratio <- sample_info_tf$TIC/mean(sample_info_tf$TIC[sample_info_tf$sampleID == "Control"])
exprs_TIC <- exprs/TIC_ratio

exprs_TIC_rsd <- strainRSD(exprs_TIC[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
x <- boxplot(rowMeans(exprs_TIC_rsd[,-1]))
text(1.2, 40, paste("median", round(x$stats[3,], 2)))

# Option 1.1: Combat + TIC normalization
TIC_ratio <- sample_info_tf$TIC/mean(sample_info_tf$TIC[sample_info_tf$sampleID == "Control"])
exprs_combat_TIC <- t(2^exprs_combat)/TIC_ratio

combat_TIC_rsd <- strainRSD(exprs_combat_TIC[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
x <- boxplot(rowMeans(combat_TIC_rsd[,-1]))
text(1.2, 40, paste("median", round(x$stats[3,], 2)))
plot(exprs_combat_growth[,combat_TIC_rsd[combat_TIC_rsd$names == "Control",-1] > 40 & colMeans(combat_TIC_rsd[,-1]) > 40])

# Option 1.2: Combat + strain growth normalization
growth_ratio <- sample_info_tf$growth/mean(sample_info_tf$growth[sample_info_tf$sampleID == "Control"])
exprs_combat_growth <- t(2^exprs_combat)/growth_ratio

combat_growth_rsd <- strainRSD(exprs_combat_growth[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
x <- boxplot(rowMeans(combat_growth_rsd[,-1]))
text(1.2, 40, paste("median", round(x$stats[3,], 2)))

df_combat_growth <- aggregate(exprs_combat_growth, by = list(sample_info_tf$gc_name), mean)
names(combat_growth_rsd[,-1])[combat_growth_rsd[combat_growth_rsd$names == "Control",-1] > 40 & colMeans(combat_growth_rsd[,-1]) > 40]
boxplot(combat_growth_rsd[,-1], las =2)

df_combat_growth_filter <- df_combat_growth[,!(combat_growth_rsd[combat_growth_rsd$names == "Control",-1] > 40 & colMeans(combat_growth_rsd[,-1]) > 40)]

# Option 2.1: MedianNorm + TIC normalization
exprs_medianNorm_TIC <- exprs_medianNorm/TIC_ratio

medianNorm_TIC_rsd <- strainRSD(exprs_medianNorm_TIC[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
x <- boxplot(rowMeans(medianNorm_TIC_rsd[,-1]))
text(1.2, 40, paste("median", round(x$stats[3,], 2)))
names(medianNorm_TIC_rsd[,-1])[medianNorm_TIC_rsd[medianNorm_TIC_rsd$names == "Control",-1] > 30 & colMeans(medianNorm_TIC_rsd[,-1]) > 30]

# Option 2.2: MedianNorm + strain grown normalization
exprs_medianNorm_growth <- exprs_medianNorm/growth_ratio

medianNorm_growth_rsd <- strainRSD(exprs_medianNorm_growth[grep("A|B|C", sample_info_tf$sampleID),], sample_info_tf$group[grep("A|B|C", sample_info_tf$sampleID)])
x <- boxplot(rowMeans(medianNorm_growth_rsd[,-1]))
text(1.2, 40, paste("median", round(x$stats[3,], 2)))
names(medianNorm_growth_rsd[,-1])[medianNorm_growth_rsd[medianNorm_growth_rsd$names == "Control",-1] > 30 & colMeans(medianNorm_growth_rsd[,-1]) > 30]

#### Write out normalized tables #### 

write.csv(cbind(sample_info_tf, exprs_combat_TIC), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_CombatBatchNorm_TICsampleNorm.csv")
write.csv(cbind(sample_info_tf, exprs_combat_growth), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_CombatBatchNorm_growthsampleNorm.csv")
write.csv(cbind(sample_info_tf, exprs_medianNorm_TIC), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_medianBatchNorm_TICsampleNorm.csv")
write.csv(cbind(sample_info_tf, exprs_medianNorm_growth), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_medianBatchNorm_growthsampleNorm.csv")
write.csv(cbind(sample_info_tf, exprs_TIC), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_TICsampleNorm.csv")
write.csv(cbind(sample_info_tf, t(exprs_combat)), file = "P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_CombatBatchNorm.csv")
write.csv(cbind(sample_info_tf, exprs_))
#### Growth norm for GC #### 
gc_sample_info <- data.frame(gc_name = gc_name)
gc_sample_info$growth <- sample_info_tf$growth[match(gc_name, sample_info_tf$gc_name)]
gc_sample_info <- gc_sample_info[-c(1:6),]

gc_growthNorm <- 2^t(df_gc[,-c(1:6)])/(gc_sample_info$growth/gc_sample_info$growth[gc_sample_info$gc_name == "PWT"][1])
df_gc_growthNorm <- cbind(df_gc[,1:6], t(log2(gc_growthNorm)))
names(df_gc_growthNorm) <- sub("RunOrderCorrected","GrowthNormalized", names(df_gc_growthNorm))
write.csv(df_gc_growthNorm, file = "P:/ES_20170625_H3K/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_growthNormalized_20210614.csv")

## alternative
merge_gc_growth <- merge(gc_exprs, sample_info_tf, by.x = "Group.1", by.y = "gc_name" )
names(merge_gc_growth)

gc_growth_ratio <- merge_gc_growth$growth/mean(merge_gc_growth$growth[merge_gc_growth$sampleID == "Control"])
gc_exprs_growth <- 2^merge_gc_growth[,2:165]/gc_growth_ratio 
df_gc_growth <- aggregate(gc_exprs_growth, by = list(merge_gc_growth[,1]), mean)


#### Merging GC and LC data to make comparisons easier #### 
merge_exprs <- merge(gc_exprs, aggregate(t(exprs_combat), by = list(sample_info_tf$gc_name), mean), by = "Group.1")
merge_exprs_growth <- merge(df_gc_growth, df_combat_growth, by = "Group.1")
merge_exprs_growth[,-1] <- log2(merge_exprs_growth[,-1])

dim(merge_exprs)
dim(merge_exprs_growth)

group <- as.factor(sub("A|B|C", "", merge_exprs$Group.1))

### Plot proline ### 

pdf("Proline_GCMS_LCMS_replicates.pdf")
plot(merge_exprs$`Proline 2 TMS`, merge_exprs$Proline,
     main = "H3K replicates",
     ylab = "LCMS Proline (Tracefinder)",
     xlab = "GCMS Proline", col = as.numeric(group), pch =19)
linearReg <- lm(merge_exprs$Proline ~ merge_exprs$`Proline 2 TMS`, )
abline(linearReg, col = "gray30", lty = 2)
text(23, 30.5, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
dev.off()

pdf("Proline_GCMS_LCMS_replicates_growthnorm.pdf")
plot(merge_exprs_growth$`Proline 2 TMS`, merge_exprs_growth$Proline,
     main = "H3K replicates",
     ylab = "LCMS Proline (COMBAT norm)",
     xlab = "GCMS Proline", col = as.numeric(group), pch =19)
linearReg <- lm(merge_exprs_growth$Proline ~ merge_exprs_growth$`Proline 2 TMS`, )
abline(linearReg, col = "gray30", lty = 2)
text(23, 30.5, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
dev.off()

merge_exprs_mean <- aggregate(merge_exprs, by = list(group), mean)
merge_exprs_growth_mean <- aggregate(merge_exprs_growth, by =list(group), mean)


pdf("E:/MITOMICS/plots/Proline_GCMS_LCMS_strainsFC.pdf")
plot(merge_exprs_mean$`Proline 2 TMS`-merge_exprs_mean$`Proline 2 TMS`[204],
     merge_exprs_mean$Proline-merge_exprs_mean$Proline[204], 
     ylab = "log2(Strain/WT) LCMS Proline (Tracefinder)",
     xlab = "log2(Strain/WT) GCMS Proline", main = "H3K strains",
     col = as.factor(mean_group$Group.1), pch=19, cex=1.5)
linearReg <- lm(FC_LCMS_Proline ~ FC_GCMS_Proline )
abline(linearReg, col = "gray30", lty = 2)
text(-2, 1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))

dev.off()


pdf("E:/MITOMICS/plots/Proline_GCMS_LCMS_strainsFC_COMBAT.pdf")
plot(merge_exprs_growth_mean$`Proline 2 TMS`-merge_exprs_growth_mean$`Proline 2 TMS`[204],
     merge_exprs_growth_mean$Proline-merge_exprs_growth_mean$Proline[204], 
     ylab = "log2(Strain/WT) LCMS Proline (COMBAT)",
     xlab = "log2(Strain/WT) GCMS Proline", main = "H3K strains",
     col = as.factor(mean_group$Group.1), pch=19, cex=1.5)
linearReg <- lm(FC_LCMS_ProlineCombat ~ FC_GCMS_Proline )
abline(linearReg, col = "gray30", lty = 2)
text(-2, 1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
dev.off()

# List of matching GC and LC 
matchingNames <- rbind(
c("2-Aminoadipic acid 3 TMS", "X2_Aminoadipic_acid"),
c("2-Ketoglutaric acid 3 TMS" , "Alpha_ketoglutaric_acid"),
c("2-Hydroxyglutaric acid 3 TMS", "Alpha_Hydroxyglutaric_acid"),
c("Alanine 2 TMS", "Alanine"),
c("Aspartic acid 3 TMS", "Aspartic_acid"),
c("Citric acid 4 TMS ", "Citric_acid"),
c("Dihydroorotate 3TMS", "Dihydroorotate"),
c("Fumaric acid 2 TMS", "Fumaric_acid"),
c("Glutamic acid 3 TMS", "Glutamic_Acid"),
c("Glutamine 3 TMS", "Glutamine"),
c("Glutathione reduced", "Glutathione_reduced"),
c("Glycerol 3-Phosphate 4 TMS", "Glycerol_2_or_3_phosphate"),
c("Isoleucine TMS", "Isoleucine"),
c("Lysine 4 TMS", "Lysine"),
c("Malic acid 3 TMS", "Malic_acid"),
c("Pantothenic acid 3 TMS", "Pantothenic_acid"),
c("Phenylalanine 2 TMS", "Phenylalanine"),
c("Proline 2 TMS", "Proline"),
c("Putrescine 4 TMS", "Putrescine"),
c("Pyruvic acid 2 TMS", "Pyruvic_acid"),
c("Serine 3 TMS", "Serine"),
c("Spermine 6TMS", "Spermine"),
c("Succinic acid 2 TMS", "Methylmalonic_acid_or_Succinic_acid"),
c("Threonine 3 TMS", "Threonine"),
c("Tryptophan 3TMS", "Tryptophan"),
c("Tyrosine 3 TMS", "Tyrosine"),
c("Uric acid 4 TMS", "Uric_acid"))

dim(matchingNames) ## 27

merge_exprs_mean <- aggregate(merge_exprs, by = list(group), mean)
merge_exprs_growth_mean <- aggregate(merge_exprs_growth, by =list(group), mean)

plotLC_vs_GC <- function(nameGC, nameLC) {
  par(mfrow = c(2,2))
  
  plot(merge_exprs[[nameGC]], merge_exprs[[nameLC]],
       main = "H3K replicates",
       ylab = paste("LCMS",nameLC),
       xlab = paste("GCMS", nameGC), 
       col = as.numeric(group), pch =19)
  linearReg <- lm(merge_exprs[[nameLC]] ~ merge_exprs[[nameGC]] )
  abline(linearReg, col = "gray30", lty = 2)
  text(min(merge_exprs[[nameGC]])+1, max(merge_exprs[[nameLC]])-1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
  
  plot(merge_exprs_growth[[nameGC]], merge_exprs_growth[[nameLC]],
       main = "H3K replicates - growth normalized",
       ylab = paste("LCMS", nameLC),
       xlab = paste("GCMS", nameGC), col = as.numeric(group), pch =19)
  linearReg <- lm(merge_exprs_growth[[nameLC]] ~ merge_exprs_growth[[nameGC]])
  abline(linearReg, col = "gray30", lty = 2)
  text(min(merge_exprs_growth[[nameGC]])+1, max(merge_exprs_growth[[nameLC]])-1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
  
  FC_GC <- merge_exprs_mean[[nameGC]]-merge_exprs_mean[[nameGC]][204]
  FC_LC <- merge_exprs_mean[[nameLC]]-merge_exprs_mean[[nameLC]][204]
  plot(FC_GC,
       FC_LC, 
       ylab = paste("log2(Strain/WT) LCMS", nameLC),
       xlab = paste("log2(Strain/WT) GCMS", nameGC), main = "H3K strains",
       col = as.factor(mean_group$Group.1), pch=19, cex=1.5)
  linearReg <- lm(FC_LC ~ FC_GC )
  abline(linearReg, col = "gray30", lty = 2)
  text(-2, 1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
 
  FC_GC <- merge_exprs_growth_mean[[nameGC]]-merge_exprs_growth_mean[[nameGC]][204]
  FC_LC <- merge_exprs_growth_mean[[nameLC]]-merge_exprs_growth_mean[[nameLC]][204]
  
  plot(FC_GC,
       FC_LC, 
       ylab = paste("log2(Strain/WT) LCMS", nameLC),
       xlab = paste("log2(Strain/WT) GCMS", nameGC), main = "H3K strains - growth normalized",
       col = as.factor(mean_group$Group.1), pch=19, cex=1.5)
  linearReg <- lm(FC_LC ~ FC_GC )
  abline(linearReg, col = "gray30", lty = 2)
  text(-2, 1, paste("r2 = ", round(summary(linearReg)$r.squared, 2)))
  
}


i = 25
plotLC_vs_GC(nameGC = matchingNames[i,1], nameLC = matchingNames[i,2])

for(i in 1:nrow(matchingNames)){
  fileName <- paste("E:/MITOMICS/plots/GCMS_LCMS_", matchingNames[i,2],".pdf", sep = "")
  pdf(fileName)
  plotLC_vs_GC(matchingNames[i,1], matchingNames[i,2])
  dev.off()
}


correl <- cor(df_exprs[,-c(1)])

#correl_2 <- two_omic_cor_fast(t(merge_exprs[,2:165]), t(merge_exprs[,-c(1:165)]))

library(pheatmap)
pdf("E:/MITOMICS/plots/LCMS_possibleDuplicates.pdf", height = 12, width = 9)
pheatmap(correl, show_rownames = F)
dev.off()

#pheatmap(correl_2[grepl(".",metabolite_info$Name),!grepl("unknown", colnames(correl_2))],show_rownames = T, show_colnames = T,
#         labels_col = strtrim(df_gc$i..Feature.ID[!grepl("unknown", colnames(correl_2))], 10),
#         labels_row = strtrim(metabolite_info$Name[grepl(".",metabolite_info$Name)], 10))

