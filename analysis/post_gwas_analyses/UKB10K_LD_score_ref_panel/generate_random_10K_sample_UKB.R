#Read in library
library(dplyr) 


#Read in sample file 
data_chr1_22_sample <- read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/id_mapping/data.chr1-22.sample", header = T)

#Extract 10,000 random participant IDs and create new output
keep_ids <- data_chr1_22_sample %>% filter(ID_1 %in% sample(unique(ID_1),10000))

write.table(keep_ids, "/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/keep_ids.txt", sep = " ", col.names=T,row.names=F,quote=F)
