details <- file.info(list.files("D:/SRP/H3K" ,pattern = ".raw", full.names = T, recursive = T))
write.csv(details, "D:/SRP/H3K/file_time_stamp.csv")
