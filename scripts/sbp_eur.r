library(gwasvcf)
library(here)

a <- read.table(here("data", "gill_sbp.txt"), header=T)

chrpos <- paste0(a$Chr, ":", a$Pos)

set_bcftools()
sbp <- query_gwas("/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/panukbb_gwas_import/processed/ukb-e-SBP_p2_EUR/ukb-e-SBP_p2_EUR.vcf.gz", chrompos=chrpos)

sbp <- vcf_to_tibble(sbp)
save(sbp, file=here("data", "ukb-e-SBP_p2_EUR.rdata"))

