library(gwasvcf)
set_bcftools()
library(here)

a <- read.table(here("data", "gill_sbp.txt"), header=T)
chrpos1 <- paste0(a$Chr, ":", a$Pos)
load(here("data", "sbp_tophits.rdata"))
chrpos2 <- paste0(sbp_tophits$chr, ":", sbp_tophits$position)
wojcik <- c("2:26932031", "4:81174592", "7:27243221", "4:81169912", "10:104616663", "4:81184341", "2:127016740", "6:151004770", "8:142389954", "4:81164723", "10:104846178")

load(here("data", "eqtlgen.rdata"))
chrpos3 <- paste0(eqtlgen$chr, ":", eqtlgen$position)

table(chrpos1 %in% chrpos2)
table(wojcik %in% chrpos2)
table(wojcik %in% chrpos1)
table(chrpos3 %in% wojcik)
chrpos <- unique(c(chrpos1, chrpos2, wojcik, chrpos3))

sbp <- query_gwas("/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/panukbb_gwas_import/processed/ukb-e-SBP_p2_EUR/ukb-e-SBP_p2_EUR.vcf.gz", chrompos=chrpos)
sbp <- vcf_to_tibble(sbp)
save(sbp, file=here("data", "ukb-e-SBP_p2_EUR.rdata"))

sbp_p3 <- query_gwas("/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/panukbb_gwas_import/processed/ukb-e-SBP_p3_EUR/ukb-e-SBP_p3_EUR.vcf.gz", chrompos=chrpos)
sbp_p3 <- vcf_to_tibble(sbp_p3)
sbp_p3
save(sbp_p3, file=here("data", "ukb-e-SBP_p3_EUR.rdata"))


