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

# Read in extracted GLGC LDL GWAS meta-analysis
extr_clumped <- fread("/user/work/ac14629/MRC_network_project/results/ldl-analysis/meta_ldl_extract.csv")

# Need to extract relevant SNPs from each of the ancestry-specific GWASs
clumped <- fread("/user/work/ac14629/MRC_network_project/data/GLGC/GLGC_cross_ancestry_ldlc_gwas_metanalysis_clumped_dt.txt")
sel_rsid <- unique(c(clumped$rsid, extr_clumped$rsid))
tmp <- tempfile()
write.table(sel_rsid, file = tmp, row = FALSE, col = FALSE, quote = FALSE)

pops <- c("afr", "eas", "eur", "sas", "ugr")

extract <- lapply(pops, function(i) {
  cmd <- glue("grep -wf {tmp} /user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis/hp/GWAS_assoc_{i}_hypertension_out.txt > {tmp}_{i}")
  system(cmd)
  fread(glue("{tmp}_{i}"))
})

ex <- list()

a <- scan("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis/hp/GWAS_assoc_afr_hypertension_out.txt", what = "character", nlines = 1)
names(extract[[1]]) <- a
ex$afr <- extract[[1]] %>%
  as_tibble() %>%
  dplyr::select(rsid = SNP, chr = chr, pos = pos, eaf = freq, ea = effect.allele, oa = other.allele, pval = SPA.pval, beta = logOR, se = logOR_SE, n = n.obs) %>%
  mutate(pop = "AFR")

names(extract[[2]]) <- a
ex$eas <- extract[[2]] %>%
  as_tibble() %>%
  dplyr::select(rsid = SNP, chr = chr, pos = pos, eaf = freq, ea = effect.allele, oa = other.allele, pval = SPA.pval, beta = logOR, se = logOR_SE, n = n.obs) %>%
  mutate(pop = "EAS")

names(extract[[4]]) <- a
ex$sas <- extract[[4]] %>%
  as_tibble() %>%
  dplyr::select(rsid = SNP, chr = chr, pos = pos, eaf = freq, ea = effect.allele, oa = other.allele, pval = SPA.pval, beta = logOR, se = logOR_SE, n = n.obs) %>%
  mutate(pop = "SAS")

names(extract[[5]]) <- a
ex$sas <- extract[[5]] %>%
  as_tibble() %>%
  dplyr::select(rsid = SNP, chr = chr, pos = pos, eaf = freq, ea = effect.allele, oa = other.allele, pval = SPA.pval, beta = logOR, se = logOR_SE, n = n.obs) %>%
  mutate(pop = "UGR")

a <- scan("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis/hp/GWAS_assoc_eur_hypertension_out.txt", what="character", nlines=1)
names(extract[[3]]) <- a
ex$eur <- extract[[3]] %>%
  as_tibble() %>%
  dplyr::select(rsid = SNP, chr = CHR, pos = BP, eaf = A1FREQ, ea = ALLELE1, oa = ALLELE0, pval = P_BOLT_LMM_INF, beta = BETA, se = SE, n = N) %>%
  mutate(pop = "EUR")

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

saveRDS(ex, file = "/user/work/ac14629/MRC_network_project/results/ldl-analysis/pop_extract_ldlc_hp_clumped.rds")
