# Outline
# This script extracts the LDL-lowering  drug target regions from each population and harmonises them

# 1. Use multi-ancestry meta-analysis to identify ancestry-agnostic top-hits (need to add chr:pos from dbSNP VCF)
# 2. Use ldl-lowering drug target regions to extract best SNP from each potential region
# 3. Use UKBB cross-ancestry panel to clump
# 4. Standardise and harmonise across ancestries using chr:pos_a1_a2 (a1,a2 alphabetical)
# 5. Save extracted data


library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(ieugwasr)
library(gwasvcf)
library(parallel)

gwasvcf::set_bcftools("/mnt/storage/software/apps/BCFTOOLS/bcftools/bcftools")
bcftools <- options()[["tools_bcftools"]]

# Read in dataset
dat <- fread("/user/work/ac14629/MRC_network_project/data/GLGC/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1")

# Need to get chr:pos
# Use dbsnp vcf file to lookup all rsids and then extract chr, pos, ID via bcftools
rsid <- dat$rsID
vcffile <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.vcf.gz"
tmp <- tempfile()
write.table(unique(rsid), file=paste0(tmp, ".snplist"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cmd <- glue("{bcftools} query -f '%CHROM %POS %ID\\n' -i'ID=@{tmp}.snplist' {vcffile} > {tmp}.txt")
system(cmd)

map <- fread(glue("{tmp}.txt"))
dat <- left_join(dat, map, by=c("rsID" = "V3"))

# Sort out non rsid IDs
dat$pvalue_GC <- as.numeric(dat$pvalue_GC)
dat$V2 <- as.numeric(dat$V2)

table(is.na(dat$V2))
table(is.na(dat$V1))
subset(dat, is.na(V1))

# Extract ldl-lowering drug target gene region by chr:pos

ldl <- fread("/user/work/ac14629/MRC_network_project/scripts/ldl-analysis/ldl_lowering_drug_target_regions.csv") %>%
  as_tibble() %>%
  tidyr::separate(`Position (hg19)`, sep = "-", into = c("start", "end"))

dat2 <- subset(dat, pvalue_GC < 1e-4)
dim(dat2)

extr <- mclapply(1:nrow(ldl), \(i) {
  a <- subset(dat2, V1 == ldl$Chromosome[i] & V2 >= ldl$start[i] & V2 <= ldl$end[i])
  a$`Target gene` <- ldl$`Target gene`[i]
  a$`Function` <- ldl$`Function`[i]
  a$`Drug` <- ldl$`Drug`[i]
  a
}, mc.cores = 10) %>%
  bind_rows()

# Clump using European reference panel
extr$rsid <- extr$rsID
extr$pval <- extr$`pvalue_GC`
bfile <- "/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink/merged_chr1_22/chr1-22_merged"
plink_bin <- "plink"
extr_clumped <- ld_clump(extr, plink_bin = plink_bin, bfile = bfile) %>%
  filter(!duplicated(rsid))

dim(extr_clumped)
table(extr_clumped$Drug)
table(ldl$Drug)

# Write
write.csv(extr_clumped, file = "/user/work/ac14629/MRC_network_project/results/ldl-analysis/meta_ldl_extract.csv")

# Need to extract relevant SNPs from each of the ancestry-specific GWASs

clumped <- fread("/user/work/ac14629/MRC_network_project/results/ldl-analysis/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1_revised.clumped")

sel_rsid <- unique(c(clumped$rsid, extr_clumped$rsid))
tmp <- tempfile()
write.table(sel_rsid, file = tmp, row = FALSE, col = FALSE, quote = FALSE)

pops <- c("AFR", "EAS", "EUR", "HIS", "SAS")

extract <- lapply(pops, function(i) {
  cmd <- glue("grep -wf {tmp} /user/work/ac14629/MRC_network_project/data/GLGC/LDL_INV_{i}_HRC_1KGP3_others_ALL.meta.singlevar.results > {tmp}_{i}")
  system(cmd)
  fread(glue("{tmp}_{i}"))
})

ex <- list()

a <- scan("/user/work/ac14629/MRC_network_project/data/GLGC/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results", what = "character", nlines = 1)
names(extract[[1]]) <- a
ex$afr <- extract[[1]] %>%
  as_tibble() %>%
  dplyr::select(rsid = rsID, chr = CHROM, pos = POS_b37, eaf = POOLED_ALT_AF, ea = ALT, oa = REF, pval = pvalue_GC, beta = EFFECT_SIZE, se = SE, n = N) %>%
  mutate(pop = "AFR")

names(extract[[2]]) <- a
ex$eas <- extract[[2]] %>%
  as_tibble() %>%
  dplyr::select(rsid = rsID, chr = CHROM, pos = POS_b37, eaf = POOLED_ALT_AF, ea = ALT, oa = REF, pval = pvalue_GC, beta = EFFECT_SIZE, se = SE, n = N) %>%
  mutate(pop = "EAS")

names(extract[[3]]) <- a
ex$eur <- extract[[3]] %>%
  as_tibble() %>%
  dplyr::select(rsid = rsID, chr = CHROM, pos = POS_b37, eaf = POOLED_ALT_AF, ea = ALT, oa = REF, pval = pvalue_GC, beta = EFFECT_SIZE, se = SE, n = N) %>%
  mutate(pop = "EUR")

names(extract[[4]]) <- a
ex$his <- extract[[4]] %>%
  as_tibble() %>%
  dplyr::select(rsid = rsID, chr = CHROM, pos = POS_b37, eaf = POOLED_ALT_AF, ea = ALT, oa = REF, pval = pvalue_GC, beta = EFFECT_SIZE, se = SE, n = N) %>%
  mutate(pop = "HIS")

names(extract[[5]]) <- a
ex$sas <- extract[[5]] %>%
  as_tibble() %>%
  dplyr::select(rsid = rsID, chr = CHROM, pos = POS_b37, eaf = POOLED_ALT_AF, ea = ALT, oa = REF, pval = pvalue_GC, beta = EFFECT_SIZE, se = SE, n = N) %>%
  mutate(pop = "SAS")

standardise <- function(d) {
  toflip <- d$ea > d$oa
  d$eaf[toflip] <- 1 - d$eaf[toflip]
  d$beta[toflip] <- d$beta[toflip] * -1
  temp <- d$oa[toflip]
  d$oa[toflip] <- d$ea[toflip]
  d$ea[toflip] <- temp
  d$snpid <- paste0(d$chr, ":", d$pos, "_", d$ea, "_", d$oa)
  d %>% arrange(chr, pos)
}

ex <- lapply(ex, standardise) %>%
  bind_rows()

ex <- left_join(ex, extr_clumped %>% dplyr::select(rsid, Drug, `Target gene`, Function), by = "rsid")

saveRDS(ex, file = "/user/work/ac14629/MRC_network_project/results/ldl-analysis/pop_extract_clumped.rds")
