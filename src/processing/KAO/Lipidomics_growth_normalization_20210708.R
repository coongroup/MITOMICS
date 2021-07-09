###### Growth normalize Lipidomics Data ###### 


files <- c("E:/MITOMICS/data/Lipidomics Batch 1 (Branch_ Lipidomics).txt" ,
  "E:/MITOMICS/data/Lipidomics Batch 2 (Branch_ Lipidomics).txt" ,
  "E:/MITOMICS/data/Lipidomics Batch 3 (Branch_ Lipidomics).txt" ,
  "E:/MITOMICS/data/Lipidomics Batch 4 (Branch_ Lipidomics).txt" ,
  "E:/MITOMICS/data/Lipidomics Batch 5 (Branch_ Lipidomics).txt" ,
  "E:/MITOMICS/data/Lipidomics Batch 6 (Branch_ Lipidomics).txt", 
  "E:/MITOMICS/data/Lipidomics Batch 7 (Branch_ Lipidomics).txt")

# lipidNames <- unique(c(names(lipid_batch1), names(lipid_batch2), names(lipid_batch3), names(lipid_batch4), names(lipid_batch5), names(lipid_batch6), names(lipid_batch7)))
# lipidGroups <- sub(".A","", lipidNames, fixed =T)
# lipidGroups <- sub(".B","", lipidGroups, fixed =T)
# lipidGroups <- sub(".C","", lipidGroups, fixed =T)
# lipidGroups <- sub(".D","", lipidGroups, fixed =T)
# 
# write.csv(unique(lipidGroups), "E:/MITOMICS/data/namesLipidSamples.csv")

mapping <- read.csv("E:/MITOMICS/data/LipidName_growth.csv") 
#write.csv(mapping, "E:/MITOMICS/data/LipidName_growth.csv", row.names = F)

### Growth normalize for each Batch #####

fileName <- files[1]
growthNormalize <- function(fileName) {
  df <- read.delim(fileName)
  df_groups<- substr(names(df),1,nchar(names(df))-2)
  
  df_growth_ratio <- mapping$growth[match(df_groups, mapping$LipidName)]/mapping$growth[mapping$Sample.ID == "PWT"][1]
  growthNorm <- 2^t(df[,-c(1:2)])/df_growth_ratio[-c(1:2)]

  df_growth <- cbind(df[,1:2], log2(t(growthNorm)))
  write.csv(df_growth, paste(sub(".txt","",fileName),"_growthNorm.csv", sep =""), row.names = F)
}

for(i in files){
  growthNormalize(i)
}
