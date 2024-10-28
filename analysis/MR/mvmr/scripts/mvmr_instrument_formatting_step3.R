### Create lipid and blood pressure instruments for MVMR ###

# Read in relevant libraries
library(dplyr)

setwd('/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr')

# Read in relevant data 
lipid <- readRDS("lipid_pop_extract_common_snps_mvmr_clumped.rds")

# Format data 
lipid <- lipid %>% rename(
	SNP = rsid, 
	exposure = trait, 
	id.exposure = trait, 
	effect_allele.exposure = ea, 
	other_allele.exposure = oa, 
	eaf.exposure = eaf, 
	beta.exposure = beta, 
	se.exposure = se, 
	pval.exposure = pval, 
	pop = pop) 

lipid <- lipid %>% select(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure) 

lipid <- as.data.frame(lipid) 

# Save output 
saveRDS(lipid, "lipid_pop_extract_common_snps_clumped_for_mvmr.rds")

	
	
	
