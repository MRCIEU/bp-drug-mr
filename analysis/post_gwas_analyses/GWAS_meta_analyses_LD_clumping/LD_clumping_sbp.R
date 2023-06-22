library(data.table)
library(TwoSampleMR)
library(ieugwasr)

#Read in directory 
dir <- c("/user/work/ac14629/MRC_network_project/results/METAL") 
diro <- c("/user/work/ac14629/MRC_network_project/results/METAL")

#Read in GWAS meta-analysis 
dat <- fread(paste0(dir, "/METAANALYSIS_SBP_UKB_UGR_gwsigthres.TBL"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)

#LD clumping 
dat <- ld_clump(
    dplyr::tibble(rsid=dat$MarkerName, pval=dat$P.value),
    plink_bin = "/mnt/storage/software/languages/R/4.2.1/lib64/R/library/plinkbinr/bin/plink_Linux",
    bfile = "/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink/merged_chr1_22/chr1-22_merged"
)

#Write output 
write.table(dat, "METAANALYSIS_SBP_UKB_UGR_clumped", sep = "\t", col.names=T,row.names=F,quote=F)
