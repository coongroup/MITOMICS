###### Combinging GC-MS data and LC-MS data with 3 PWT samples (for website) ######

## Function for strain RSD calculation 
strainRSD <- function(matrix, strain_info){
  # matrix is samples in columns, features in rows 
  # strain info is vector with length == to nrow(matrix)
  strain_sd <- aggregate(matrix, by = list(strain_info), sd)
  strain_mean <- aggregate(matrix, by = list(strain_info), mean)
  rsd <- data.frame(names = strain_mean[,1], strain_sd[,-1]/strain_mean[,-1] *100)
  rsd
}

##### Read in GC-MS Data ######

## GCMS file info and sample info
## This gc table still has replicate samples (reinjected samples and batch 3 and 4 day 1)
file_info  <- read.csv("P:/ES_20170625_H3K/SRP_012019_H3K/KAO_GC_quant/file_time_stamp_v2.csv")
gc <- read.csv("P:/ES_20170625_H3K/SRP_012019_H3K/KAO_GC_quant/QuantResults_20191030_KAO_H3K_HMDB_runOrderCorrected_20200810.csv", stringsAsFactors = F)


dim(file_info)
dim(gc)
names(gc)

gc_values <- gc[,grep("_RunOrderCorrected", names(gc))]
row.names(gc_values) <- gc$ï..Feature.ID
dim(gc_values)

file_info$Sample[match(sub("_RunOrderCorrected", "", names(gc_values)),file_info$Sample)] == sub("_RunOrderCorrected", "", names(gc_values))
file_info <- file_info[match(sub("_RunOrderCorrected", "", names(gc_values)),file_info$Sample), ]
file_info$Sample == sub("_RunOrderCorrected", "", names(gc_values))


sample_duplicates <- which(duplicated(paste(file_info$group,file_info$BatchDay)) & file_info$group != "PWT")
sample_remove <- unique(c(sample_duplicates))

gc_values_unique <- gc_values[,-sample_remove]

sampleNames <- names(gc_values_unique)[grep("_LFQ", names(gc_values_unique))]
sampleNames <- sub("_RunOrderCorrected", "", sampleNames)

## Confirm matches
file_info <- file_info[file_info$Sample %in% sampleNames, ]
sampleNames[!sampleNames %in% file_info$Sample]
file_info$Sample[match(sampleNames,file_info$Sample)] == sampleNames
dim(file_info)

table(file_info$Sample == sub("_RunOrderCorrected","",names(gc_values_unique))) ## All True

## transpose value table
t_gc_values <- t(gc_values_unique)

## Calculated PWT mean for each Day to yield 3 averaged PWT; exclude outlier Batch3and4_Day1
values_Day_mean <- aggregate(t_gc_values[file_info$BatchDay != "Batch3and4_Day1",], by = list(file_info$group[file_info$BatchDay != "Batch3and4_Day1"], file_info$Day[file_info$BatchDay != "Batch3and4_Day1"]), mean)

gc_Day_PWT_mean <- t(values_Day_mean[values_Day_mean[,1] == "PWT",-c(1:2)])
names(gc_Day_PWT_mean) <- c("A_PWT","B_PWT","C_PWT")


file_info$ID <- sub("_LFQ","", file_info$Sample, fixed = T)
file_info$ID <- sub("_1901........", "", file_info$ID)
file_info <- rbind(file_info, c(NA,NA,NA,"Day1","A_PWT",NA,"PWT","A_PWTM"))
file_info <- rbind(file_info, c(NA,NA,NA,"Day2","B_PWT",NA,"PWT","B_PWTM"))
file_info <- rbind(file_info, c(NA,NA,NA,"Day3","C_PWT",NA,"PWT","C_PWTM"))
dim(file_info) # 739 8

gc_values_unique <- cbind(gc_values_unique, gc_Day_PWT_mean)
dim(gc_values_unique) # 164 739

## clean up data by removing extra PWT (only keep averaged PWT)
gc_values_triplicates <- gc_values_unique[,-grep('PWT_._._.', names(gc_values_unique))]
file_info_triplicates <- file_info[-grep("PWT_._._", file_info$Sample),]

dim(gc_values_triplicates) # 164 660
dim(file_info_triplicates) # 660 8

###### Load in LC-MS data ######

lc <- read.csv("P:/AZS_20212604_H3K_Metabolites/TraceFinder/Peak Areas/For Katie_H3K_Growth_TIC/20210603_H3K_Overall_Master_File_Imputed_Perseus_2p1v20210604AZS_CombatBatchNorm_growthsampleNorm.csv", stringsAsFactors = )
# compounds to remove due to high variability
metabolites_to_remove <- c("X4_Guanidinobutyric_acid","Dihydroorotate", "Pimelic_acid", "Suberic_acid",             
  "Sedoheptulose_7_phosphate", "Pyruvic_acid")

lc <- lc[,-match(metabolites_to_remove, names(lc))]

controls1 <- grep("Control", lc$sampleID)[c(FALSE, FALSE, TRUE)] #sample every 3rd control sample
controls2 <- grep("Control", lc$sampleID)[c(FALSE, TRUE, FALSE)]
controls3 <- grep("Control", lc$sampleID)[c(TRUE, FALSE, FALSE)]

labels1 <- c(NA, "Control_A","Control_A", "Control", NA, 182357.0, NA, NA, "A_PWTM", "PWT",TRUE, "A_PWTM", NA)
labels2 <- c(NA, "Control_B","Control_B", "Control", NA, 182357.0, NA, NA, "B_PWTM", "PWT",TRUE, "B_PWTM", NA)
labels3 <- c(NA, "Control_C","Control_C", "Control", NA, 182357.0, NA, NA, "C_PWTM", "PWT",TRUE, "C_PWTM", NA)

lc_PWT_means <- rbind(colMeans(lc[controls1,-c(1:13)]), colMeans(lc[controls2,-c(1:13)]),colMeans(lc[controls3,-c(1:13)]))

#replace first 3 controls with average control information 
control_index <- grep("Control", lc$group)[1:3]
lc_PWT <- lc
lc_PWT[control_index[1],1:13] <- labels1
lc_PWT[control_index[1],-c(1:13)] <- colMeans(lc[controls1,-c(1:13)])
lc_PWT[control_index[2],1:13] <- labels2
lc_PWT[control_index[2],-c(1:13)] <- colMeans(lc[controls2,-c(1:13)])
lc_PWT[control_index[3],1:13] <- labels3
lc_PWT[control_index[3],-c(1:13)] <- colMeans(lc[controls3,-c(1:13)])


lc_PWT <- lc_PWT[lc_PWT$is_sample == T,]

dim(lc_PWT)

### Growth normalize the GC-MS data ####
file_info_triplicates$growth <- NA
file_info_triplicates$growth <- as.numeric(lc_PWT$growth[match(file_info_triplicates$ID,lc_PWT$gc_name)])

gc_values_triplicates_growthNorm <- 2^t(gc_values_triplicates)/(file_info_triplicates$growth/file_info_triplicates$growth[file_info_triplicates$group == "PWT"][1])

gc_growthNorm <- cbind(file_info_triplicates, gc_values_triplicates_growthNorm)

## Combine GC and LC data ###########

gc_growthNorm$ID[is.na(match(file_info_triplicates$ID, lc_PWT$gc_name))]
lc_PWT$gc_name[is.na(match(lc_PWT$gc_name, file_info_triplicates$ID))]
gc_lc_merge <- merge(gc_growthNorm, lc_PWT, by.x = "ID", by.y = "gc_name")

write.csv(gc_lc_merge, "E:/Mitomics/data/gc_lc_merge_20210706_KAO.csv")
names(gc_lc_merge) #10:173 = GC, 186:265 = LC

rsd_gc <- strainRSD(gc_lc_merge[,10:173], gc_lc_merge$group.x)
rsd_lc <- strainRSD(gc_lc_merge[,186:265], gc_lc_merge$group.y)
rsd_merge<- merge(rsd_gc, rsd_lc, by = "names")
write.csv(rsd_merge, "E:/Mitomics/data/gc_lc_merge_rsds_20210706_KAO.csv")
