### Create lipid and blood pressure instruments for MVMR ###

# Read in relevant libraries
library(dplyr)

setwd('/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr')

# Read in datasets 
ldl <- readRDS("ldlc_pop_extract_common_snps__mvmr_clumped.rds")
hdl <- readRDS("hdl_pop_extract_common_snps_mvmr_clumped.rds")
tg <- readRDS("tg_pop_extract_common_snps_mvmr_clumped.rds")

# Format exposure column 
ldl$trait <- "LDLC"
hdl$trait <- "HDL"
tg$trait <- "TG"

# Combine datasets
lipid_com <- rbind(ldl, hdl,tg) 
lipid_com <- lipid_com %>% filter(pop != "CROSS_ANCESTRY")

# Save output 
saveRDS(lipid_com,"lipid_pop_extract_common_snps_mvmr_clumped.rds")
