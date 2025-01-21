###Multivariable cis-MR using PCs###

# Load relevant libraries
library(simulateGP)
library(dplyr)
library(data.table)
library(tidyr)
library(MendelianRandomization)

# Get LD matrix
load(url("https://github.com/explodecomputer/simulateGP/raw/master/data/ldetect.rdata"))
head(ldetect)
a <- subset(ldetect, pop == "EUR" & chr == "chr5") %>%
     mutate(len = stop - start) %>%
     arrange(len) %>%
     slice(100)

afrbfile <- "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/AFR"
amrbfile <- "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/AMR"
easbfile <- "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/EAS"
eurbfile <- "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/EUR"
sasbfile <- "/user/work/ac14629/MRC_network_project/data/UKB/LD_ref_dataset/SAS"

afrld <- get_ld(a$chr, from=72465174, to=77176458, afrbfile)
amrld <- get_ld(a$chr, from=72465174, to=77176458, amrbfile)
easld <- get_ld(a$chr, from=72465174, to=77176458, easbfile)
eurld <- get_ld(a$chr, from=72465174, to=77176458, eurbfile)
sasld <- get_ld(a$chr, from=72465174, to=77176458, sasbfile)

## Read in summary statistics for lipid traits
afr_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_afr_ldlc_common_snps_formatted.txt")
amr_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/mcps_ldlc_common_snps_formatted.txt")
eas_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ckb_ldlc_common_snps_formatted.txt")
eur_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_eur_ldlc_common_snps_formatted.txt")
sas_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ukb_sas_ldlc_common_snps_formatted.txt")
ugr_ldlc <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/ugr_ldlc_common_snps_formatted.txt")

afr_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/ukb_afr_hdl_common_snps_formatted.txt")
amr_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/mcps_hdl_common_snps_formatted.txt")
eas_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/ckb_hdl_common_snps_formatted.txt")
eur_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/ukb_eur_hdl_common_snps_formatted.txt")
sas_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/ukb_sas_hdl_common_snps_formatted.txt")
ugr_hdl <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/ugr_hdl_common_snps_formatted.txt")

afr_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/ukb_afr_tg_common_snps_formatted.txt")
amr_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/mcps_tg_common_snps_formatted.txt")
eas_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/ckb_tg_common_snps_formatted.txt")
eur_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/ukb_eur_tg_common_snps_formatted.txt")
sas_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/ukb_sas_tg_common_snps_formatted.txt")
ugr_tg <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/ugr_tg_common_snps_formatted.txt")

# Read in disease outcome summary statistics
afr_chd <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/chd/ukb_afr_chd_formatted_rescaled.txt")
amr_chd <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/chd/mcps_chd_formatted_rescaled.txt")
eas_chd <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/chd/ckb_chd_formatted_rescaled.txt")
eur_chd <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/chd/ukb_eur_chd_formatted.txt")
sas_chd <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/chd/ukb_sas_chd_formatted_rescaled.txt")

afr_stroke <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/stroke/ukb_afr_stroke_formatted_rescaled.txt")
amr_stroke <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/stroke/mcps_stroke_formatted_rescaled.txt")
eas_stroke <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/stroke/ckb_stroke_formatted_rescaled.txt")
eur_stroke <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/stroke/ukb_eur_stroke_formatted.txt")
sas_stroke <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/stroke/ukb_sas_stroke_formatted_rescaled.txt")

afr_hp <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hp/ukb_afr_hp_formatted_rescaled.txt")
amr_hp <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hp/mcps_hp_formatted_rescaled.txt")
eas_hp <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hp/ckb_hp_formatted_rescaled.txt")
eur_hp <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hp/ukb_eur_hp_formatted.txt")
sas_hp <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hp/ukb_sas_hp_formatted_rescaled.txt")

## Select drug target gene region 
hmgcr_ldlc_afr <- afr_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_ldlc_amr <- amr_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_ldlc_eas <- eas_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_ldlc_eur <- eur_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_ldlc_sas <- sas_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_ldlc_ugr <- amr_ldlc %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)

hmgcr_hdl_afr <- afr_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_hdl_amr <- amr_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_hdl_eas <- eas_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_hdl_eur <- eur_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_hdl_sas <- sas_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_hdl_ugr <- amr_hdl %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)

hmgcr_tg_afr <- afr_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_tg_amr <- amr_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_tg_eas <- eas_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_tg_eur <- eur_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_tg_sas <- sas_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)
hmgcr_tg_ugr <- amr_tg %>% filter(chr == 5, pos >= 74632154 - 250000, pos <= 74632154 + 250000)

## Get common snp ids
common_snp_ids_afr <- Reduce(intersect, list(hmgcr_ldlc_afr$snp.id, 
                                         hmgcr_hdl_afr$snp.id, 
                                         hmgcr_tg_afr$snp.id))


common_snp_ids_amr <- Reduce(intersect, list(hmgcr_ldlc_amr$snp.id,
                                         hmgcr_hdl_amr$snp.id,
                                         hmgcr_tg_amr$snp.id))

common_snp_ids_eas <- Reduce(intersect, list(hmgcr_ldlc_eas$snp.id,
                                         hmgcr_hdl_eas$snp.id,
                                         hmgcr_tg_eas$snp.id))

common_snp_ids_eur <- Reduce(intersect, list(hmgcr_ldlc_eur$snp.id,
                                         hmgcr_hdl_eur$snp.id,
                                         hmgcr_tg_eur$snp.id))

common_snp_ids_sas <- Reduce(intersect, list(hmgcr_ldlc_sas$snp.id,
                                         hmgcr_hdl_sas$snp.id,
                                         hmgcr_tg_sas$snp.id))

common_snp_ids_ugr <- Reduce(intersect, list(hmgcr_ldlc_ugr$snp.id,
                                         hmgcr_hdl_ugr$snp.id,
                                         hmgcr_tg_ugr$snp.id))


## Subset by common snp ids 
hmgcr_ldlc_afr_common <- hmgcr_ldlc_afr[hmgcr_ldlc_afr$snp.id %in% common_snp_ids_afr, ]
hmgcr_hdl_afr_common <- hmgcr_hdl_afr[hmgcr_hdl_afr$snp.id %in% common_snp_ids_afr, ]
hmgcr_tg_afr_common <- hmgcr_tg_afr[hmgcr_tg_afr$snp.id %in% common_snp_ids_afr, ]
afr_chd_common <- afr_chd[afr_chd$snp %in% common_snp_ids_afr, ]
afr_stroke_common <- afr_stroke[afr_stroke$snp %in% common_snp_ids_afr, ]
afr_hp_common <- afr_hp[afr_hp$snp %in% common_snp_ids_afr, ]

hmgcr_ldlc_amr_common <- hmgcr_ldlc_amr[hmgcr_ldlc_amr$snp.id %in% common_snp_ids_amr, ]
hmgcr_hdl_amr_common <- hmgcr_hdl_amr[hmgcr_hdl_amr$snp.id %in% common_snp_ids_amr, ]
hmgcr_tg_amr_common <- hmgcr_tg_amr[hmgcr_tg_amr$snp.id %in% common_snp_ids_amr, ]
amr_chd_common <- amr_chd[amr_chd$variant.id %in% common_snp_ids_afr, ]
amr_stroke_common <- amr_stroke[amr_stroke$variant.id %in% common_snp_ids_afr, ]
amr_hp_common <- amr_hp[amr_hp$variant.id %in% common_snp_ids_afr, ]

hmgcr_ldlc_eas_common <- hmgcr_ldlc_eas[hmgcr_ldlc_eas$snp.id %in% common_snp_ids_eas, ]
hmgcr_hdl_eas_common <- hmgcr_hdl_eas[hmgcr_hdl_eas$snp.id %in% common_snp_ids_eas, ]
hmgcr_tg_eas_common <- hmgcr_tg_eas[hmgcr_tg_eas$snp.id %in% common_snp_ids_eas, ]
eas_chd_common <- eas_chd[eas_chd$variant.id %in% common_snp_ids_eas, ]
eas_stroke_common <- eas_stroke[eas_stroke$variant.id %in% common_snp_ids_eas, ]
eas_hp_common <- eas_hp[eas_hp$variant.id %in% common_snp_ids_eas, ]

hmgcr_ldlc_eur_common <- hmgcr_ldlc_eur[hmgcr_ldlc_eur$snp.id %in% common_snp_ids_eur, ]
hmgcr_hdl_eur_common <- hmgcr_hdl_eur[hmgcr_hdl_eur$snp.id %in% common_snp_ids_eur, ]
hmgcr_tg_eur_common <- hmgcr_tg_eur[hmgcr_tg_eur$snp.id %in% common_snp_ids_eur, ]
eur_chd_common <- eur_chd[eur_chd$variant.id %in% common_snp_ids_eur, ]
eur_stroke_common <- eur_stroke[eur_stroke$variant.id %in% common_snp_ids_eur, ]
eur_hp_common <- eur_hp[eur_hp$variant.id %in% common_snp_ids_eur, ]

hmgcr_ldlc_sas_common <- hmgcr_ldlc_sas[hmgcr_ldlc_sas$snp.id %in% common_snp_ids_sas, ]
hmgcr_hdl_sas_common <- hmgcr_hdl_sas[hmgcr_hdl_sas$snp.id %in% common_snp_ids_sas, ]
hmgcr_tg_sas_common <- hmgcr_tg_sas[hmgcr_tg_sas$snp.id %in% common_snp_ids_sas, ]
sas_chd_common <- sas_chd[sas_chd$snp %in% common_snp_ids_sas, ]
sas_stroke_common <- sas_stroke[sas_stroke$snp %in% common_snp_ids_sas, ]
sas_hp_common <- sas_hp[sas_hp$snp %in% common_snp_ids_sas, ]

hmgcr_ldlc_ugr_common <- hmgcr_ldlc_ugr[hmgcr_ldlc_ugr$snp.id %in% common_snp_ids_ugr, ]
hmgcr_hdl_ugr_common <- hmgcr_hdl_ugr[hmgcr_hdl_ugr$snp.id %in% common_snp_ids_ugr, ]
hmgcr_tg_ugr_common <- hmgcr_tg_ugr[hmgcr_tg_ugr$snp.id %in% common_snp_ids_ugr, ]

# Subset matrix by matching to GWAS common SNPs 
afr_map_subset <- afrld$map %>%
  filter(snp %in% common_snp_ids_afr)
afr_selected_indices <- which(afrld$map$snp %in% common_snp_ids_afr)
afr_ld_subset <- afrld$ld[selected_indices, selected_indices]

amr_map_subset <- amrld$map %>%
  filter(snp %in% common_snp_ids_amr)
amr_selected_indices <- which(amrld$map$snp %in% common_snp_ids_amr)
amr_ld_subset <- amrld$ld[selected_indices, selected_indices]

eas_map_subset <- easld$map %>%
  filter(snp %in% common_snp_ids_eas)
eas_selected_indices <- which(easld$map$snp %in% common_snp_ids_eas)
eas_ld_subset <- easld$ld[selected_indices, selected_indices]

eur_map_subset <- eurld$map %>%
  filter(snp %in% common_snp_ids_eur)
eur_selected_indices <- which(eurld$map$snp %in% common_snp_ids_eur)
eur_ld_subset <- eurld$ld[selected_indices, selected_indices]

sas_map_subset <- sasld$map %>%
  filter(snp %in% common_snp_ids_sas)
sas_selected_indices <- which(sasld$map$snp %in% common_snp_ids_sas)
sas_ld_subset <- sasld$ld[selected_indices, selected_indices]

ugr_map_subset <- afrld$map %>%
  filter(snp %in% common_snp_ids_ugr)
ugr_selected_indices <- which(afrld$map$snp %in% common_snp_ids_ugr)
ugr_ld_subset <- afrld$ld[selected_indices, selected_indices]

# Ensure that GWASs and matrix have the same dimensions 
hmgcr_ldlc_afr_common <- hmgcr_ldlc_afr_common %>%
  filter(variant.id %in% afr_map_subset$snp)
hmgcr_hdl_afr_common <- hmgcr_hdl_afr_common %>%
  filter(variant.id %in% afr_map_subset$snp)
hmgcr_tg_afr_common <- hmgcr_tg_afr_common %>%
  filter(variant.id %in% afr_map_subset$snp)
afr_chd_common <- afr_chd_common %>%
  filter(variant.id %in% afr_map_subset$snp)
afr_stroke_common <- afr_stroke_common %>%
  filter(variant.id %in% afr_map_subset$snp)
afr_hp_common <- afr_hp_common %>%
  filter(variant.id %in% afr_map_subset$snp)

hmgcr_ldlc_amr_common <- hmgcr_ldlc_amr_common %>%
  filter(variant.id %in% amr_map_subset$snp)
hmgcr_hdl_amr_common <- hmgcr_hdl_amr_common %>%
  filter(variant.id %in% amr_map_subset$snp)
hmgcr_tg_amr_common <- hmgcr_tg_amr_common %>%
  filter(variant.id %in% amr_map_subset$snp)
amr_chd_common <- amr_chd_common %>%
  filter(variant.id %in% amr_map_subset$snp)
amr_stroke_common <- amr_stroke_common %>%
  filter(variant.id %in% amr_map_subset$snp)
amr_hp_common <- amr_hp_common %>%
  filter(variant.id %in% amr_map_subset$snp)

hmgcr_ldlc_eas_common <- hmgcr_ldlc_eas_common %>%
  filter(variant.id %in% eas_map_subset$snp)
hmgcr_hdl_eas_common <- hmgcr_hdl_eas_common %>%
  filter(variant.id %in% eas_map_subset$snp)
hmgcr_tg_eas_common <- hmgcr_tg_eas_common %>%
  filter(variant.id %in% eas_map_subset$snp)
eas_chd_common <- eas_chd_common %>%
  filter(variant.id %in% eas_map_subset$snp)
eas_stroke_common <- eas_stroke_common %>%
  filter(variant.id %in% eas_map_subset$snp)
eas_hp_common <- eas_hp_common %>%
  filter(variant.id %in% eas_map_subset$snp)

hmgcr_ldlc_eur_common <- hmgcr_ldlc_eur_common %>%
  filter(variant.id %in% eur_map_subset$snp)
hmgcr_hdl_eur_common <- hmgcr_hdl_eur_common %>%
  filter(variant.id %in% eur_map_subset$snp)
hmgcr_tg_eur_common <- hmgcr_tg_eur_common %>%
  filter(variant.id %in% eur_map_subset$snp)
eur_chd_common <- eur_chd_common %>%
  filter(variant.id %in% eur_map_subset$snp)
eur_stroke_common <- eur_stroke_common %>%
  filter(variant.id %in% eur_map_subset$snp)
eur_hp_common <- eur_hp_common %>%
  filter(variant.id %in% afr_map_subset$snp)

hmgcr_ldlc_sas_common <- hmgcr_ldlc_sas_common %>%
  filter(variant.id %in% sas_map_subset$snp)
hmgcr_hdl_sas_common <- hmgcr_hdl_sas_common %>%
  filter(variant.id %in% sas_map_subset$snp)
hmgcr_tg_sas_common <- hmgcr_tg_sas_common %>%
  filter(variant.id %in% sas_map_subset$snp)
sas_chd_common <- sas_chd_common %>%
  filter(variant.id %in% sas_map_subset$snp)
sas_stroke_common <- sas_stroke_common %>%
  filter(variant.id %in% sas_map_subset$snp)
sas_hp_common <- sas_hp_common %>%
  filter(variant.id %in% afr_map_subset$snp)

hmgcr_ldlc_ugr_common <- hmgcr_ldlc_ugr_common %>%
  filter(variant.id %in% ugr_map_subset$snp)
hmgcr_hdl_ugr_common <- hmgcr_hdl_ugr_common %>%
  filter(variant.id %in% ugr_map_subset$snp)
hmgcr_tg_ugr_common <- hmgcr_tg_ugr_common %>%
  filter(variant.id %in% ugr_map_subset$snp)
ugr_chd_common <- afr_chd_common %>%
  filter(variant.id %in% ugr_map_subset$snp)
ugr_stroke_common <- afr_stroke_common %>%
  filter(variant.id %in% ugr_map_subset$snp)
ugr_hp_common <- afr_hp_common %>%
  filter(variant.id %in% ugr_map_subset$snp)

# Harmonise matrix with data 
flip <- ifelse(afr_map_subset == hmgcr_ldlc_afr_common$effect.allele, 1, 0)
m <- flip %*% t(flip)
afrld <- afr_ld_subset * m

flip <- ifelse(amr_map_subset == hmgcr_ldlc_amr_common$A1, 1, 0)
m <- flip %*% t(flip)
amrld <- amr_ld_subset * m

flip <- ifelse(eas_map_subset == hmgcr_ldlc_eas_common$EA, 1, 0)
m <- flip %*% t(flip)
easld <- eas_ld_subset * m

flip <- ifelse(eur_map_subset == hmgcr_ldlc_eur_common$ALLELE1, 1, 0)
m <- flip %*% t(flip)
eurld <- eur_ld_subset * m

flip <- ifelse(sas_map_subset == hmgcr_ldlc_sas_common$effect.allele, 1, 0)
m <- flip %*% t(flip)
sasld <- sas_ld_subset * m

## MV-PCA method (Batool et al. 2022) 

Psi_afr <- ((abs(hmgcr_ldlc_afr_common$Est)+abs(hmgcr_hdl_afr_common$Est)+abs(hmgcr_tg_afr_common$Est))/afr_chd_common$Est.SE_rescaled)%o%
((abs(hmgcr_ldlc_afr_common$Est)+abs(hmgcr_hdl_afr_common$Est)+abs(hmgcr_tg_afr_common$Est))/afr_chd_common$Est.SE_rescaled)*afrld_subset

Psi_amr <- ((abs(hmgcr_ldlc_amr_common$BETA)+abs(hmgcr_hdl_amr_common$BETA)+abs(hmgcr_tg_amr_common$BETA))/amr_chd_common$SE_rescaled)%o%
((abs(hmgcr_ldlc_amr_common$BETA)+abs(hmgcr_hdl_amr_common$BETA)+abs(hmgcr_tg_amr_common$BETA))/amr_chd_common$SE_rescaled)*amrld_subset

Psi_eas <- ((abs(hmgcr_ldlc_eas_common$BETA)+abs(hmgcr_hdl_eas_common$BETA)+abs(hmgcr_tg_eas_common$BETA))/eas_chd_common$SE_rescaled)%o%
((abs(hmgcr_ldlc_eas_common$BETA)+abs(hmgcr_hdl_eas_common$BETA)+abs(hmgcr_tg_eas_common$BETA))/eas_chd_common$SE_rescaled)*easld_subset

Psi_eur <- ((abs(hmgcr_ldlc_eur_common$BETA)+abs(hmgcr_hdl_eur_common$BETA)+abs(hmgcr_tg_eur_common$BETA))/eur_chd_common$logOR.SE)%o%
((abs(hmgcr_ldlc_eur_common$BETA)+abs(hmgcr_hdl_eur_common$BETA)+abs(hmgcr_tg_eur_common$BETA))/eur_chd_common$logOR.SE)*eurld_subset

Psi_sas <- ((abs(hmgcr_ldlc_sas_common$Est)+abs(hmgcr_hdl_sas_common$Est)+abs(hmgcr_tg_sas_common$Est))/sas_chd_common$Est.SE_rescaled)%o%
((abs(hmgcr_ldlc_sas_common$Est)+abs(hmgcr_hdl_sas_common$Est)+abs(hmgcr_tg_sas_common$Est))/sas_chd_common$Est.SE_rescaled)*sasld_subset

Psi_ugr <- ((abs(hmgcr_ldlc_ugr_common$Est)+abs(hmgcr_hdl_ugr_common$Est)+abs(hmgcr_tg_ugr_common$Est))/afr_chd_common$Est.SE_rescaled)%o%
((abs(hmgcr_ldlc_ugr_common$Est)+abs(hmgcr_hdl_ugr_common$Est)+abs(hmgcr_tg_ugr_common$Est))/afr_chd_common$Est.SE_rescaled)*afrld_subset

pcs_afr <- eigen(Psi_afr)
pcs_amr <- eigen(Psi_amr)
pcs_eas <- eigen(Psi_eas)
pcs_eur <- eigen(Psi_eur)
pcs_sas <- eigen(Psi_sas)
pcs_ugr <- eigen(Psi_ugr)

K_afr <- which(cumsum(pcs_afr$values^2) / sum(pcs_afr$values^2) > 0.99)[1]
K_afr

K_amr <- which(cumsum(pcs_amr$values^2) / sum(pcs_amr$values^2) > 0.99)[1]
K_amr

K_eas <- which(cumsum(pcs_eas$values^2) / sum(pcs_eas$values^2) > 0.99)[1]
K_eas

K_eur <- which(cumsum(pcs_eur$values^2) / sum(pcs_eur$values^2) > 0.99)[1]
K_eur

K_sas <- which(cumsum(pcs_sas$values^2) / sum(pcs_sas$values^2) > 0.99)[1]
K_sas

K_ugr <- which(cumsum(pcs_ugr$values^2) / sum(pcs_ugr$values^2) > 0.99)[1]
K_ugr

afrx1comp <- hmgcr_ldlc_afr_common$Est %*% pcs_afr$vectors[,1:K_afr] %>% drop() 
afrx2comp <- hmgcr_hdl_afr_common$Est %*% pcs_afr$vectors[,1:K_afr] %>% drop()  
afrx3comp <- hmgcr_tg_afr_common$Est %*% pcs_afr$vectors[,1:K_afr] %>% drop()
afrycomp <- afr_chd_common$Est_rescaled %*% pcs_afr$vectors[,1:K_afr] %>% drop()

amrx1comp <- hmgcr_ldlc_amr_common$BETA %*% pcs_amr$vectors[,1:K_amr] %>% drop()
amrx2comp <- hmgcr_hdl_amr_common$BETA %*% pcs_amr$vectors[,1:K_amr] %>% drop()
amrx3comp <- hmgcr_tg_amr_common$BETA %*% pcs_amr$vectors[,1:K_amr] %>% drop()
amrycomp <- amr_chd_common$BETA_rescaled %*% pcs_amr$vectors[,1:K_amr] %>% drop()

easx1comp <- hmgcr_ldlc_eas_common$BETA %*% pcs_eas$vectors[,1:K_eas] %>% drop()
easx2comp <- hmgcr_hdl_eas_common$BETA %*% pcs_eas$vectors[,1:K_eas] %>% drop()
easx3comp <- hmgcr_tg_eas_common$BETA %*% pcs_eas$vectors[,1:K_eas] %>% drop()
easycomp <- eas_chd_common$BETA_rescaled %*% pcs_eas$vectors[,1:K_eas] %>% drop()

eurx1comp <- hmgcr_ldlc_eur_common$BETA %*% pcs_eur$vectors[,1:K_eur] %>% drop()
eurx2comp <- hmgcr_hdl_eur_common$BETA %*% pcs_eur$vectors[,1:K_eur] %>% drop()
eurx3comp <- hmgcr_tg_eur_common$BETA %*% pcs_eur$vectors[,1:K_eur] %>% drop()
eurycomp <- eur_chd_common$logOR %*% pcs_eur$vectors[,1:K_eur] %>% drop()

sasx1comp <- hmgcr_ldlc_sas_common$Est %*% pcs_sas$vectors[,1:K_sas] %>% drop()
sasx2comp <- hmgcr_hdl_sas_common$Est %*% pcs_sas$vectors[,1:K_sas] %>% drop()
sasx3comp <- hmgcr_tg_sas_common$Est %*% pcs_sas$vectors[,1:K_sas] %>% drop()
sasycomp <- sas_chd_common$Est_rescaled %*% pcs_sas$vectors[,1:K_sas] %>% drop()

ugrx1comp <- (hmgcr_ldlc_ugr_common$Est %*% pcs_ugr$vectors[,1:K_ugr] %>% drop()
ugrx2comp <- (hmgcr_hdl_ugr_common$Est %*% pcs_ugr$vectors[,1:K_ugr] %>% drop()
ugrx3comp <- (hmgcr_tg_ugr_common$Est %*% pcs_ugr$vectors[,1:K_ugr] %>% drop()
ugrycomp <- (afr_chd_common$Est_rescaled %*% pcs_ugr$vectors[,1:K_ugr] %>% drop()

## Multivariable MR (unweighted)
afr_mvmr_unw <- summary(lm(afrycomp ~ 0 + afrx1comp + afrx2comp + afrx3comp))
amr_mvmr_unw <- summary(lm(amrycomp ~ 0 + amrx1comp + amrx2comp + amrx3comp))
eas_mvmr_unw <- summary(lm(easycomp ~ 0 + easx1comp + easx2comp + easx3comp))
eur_mvmr_unw <- summary(lm(afrycomp ~ 0 + afrx1comp + afrx2comp + afrx3comp))
sas_mvmr_unw <- summary(lm(sasycomp ~ 0 + sasx1comp + sasx2comp + sasx3comp))
ugr_mvmr_unw <- summary(lm(afrycomp ~ 0 + ugrx1comp + ugrx2comp + ugrx3comp))

#saveRDS(afr_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_afr.rds")
#saveRDS(amr_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_amr.rds")
#saveRDS(eas_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_eas.rds")
#saveRDS(eur_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_eur.rds")
#saveRDS(sas_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_sas.rds")
#saveRDS(ugr_mvmr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_cis_mvmr_unweighted_ugr.rds")

## Univariable MR (unweighted)
#afr_univariable_mr_unw <- summary(lm(afrycomp ~ 0 + afrx1comp))
#amr_univariable_mr_unw <- summary(lm(amrycomp ~ 0 + amrx1comp))
#eas_univariable_mr_unw <- summary(lm(easycomp ~ 0 + easx1comp))
#eur_univariable_mr_unw <- summary(lm(eurycomp ~ 0 + eurx1comp))
#sas_univariable_mr_unw <- summary(lm(sasycomp ~ 0 + sasx1comp))
#ugr_univariable_mr_unw <- summary(lm(afrycomp ~ 0 + ugrx1comp))
#
#saveRDS(afr_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_afr.rds")
#saveRDS(amr_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_amr.rds")
#saveRDS(eas_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_eas.rds")
#saveRDS(eur_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_eur.rds")
#saveRDS(sas_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_sas.rds")
#saveRDS(ugr_univariable_mr_unw, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_unweighted_ugr.rds")

## To do weights
Omega_afr <- (afr_chd_common$Est.SE_rescaled %o% afr_chd_common$Est.SE_rescaled) * afrld_subset
pcOmega_afr <- t(prcomp(Psi_afr, scale=FALSE)$rotation[,1:K_afr])%*%Omega_afr%*%
prcomp(Psi_afr, scale=FALSE)$rotation[,1:K_afr]

Omega_amr <- (amr_chd_common$SE_rescaled %o% amr_chd_common$SE_rescaled) * amrld_subset
pcOmega_amr <- t(prcomp(Psi_amr, scale=FALSE)$rotation[,1:K_amr])%*%Omega_amr%*%
prcomp(Psi_amr, scale=FALSE)$rotation[,1:K_amr]

Omega_eas <- (eas_chd_common$SE_rescaled %o% eas_chd_common$SE_rescaled) * easld_subset
pcOmega_eas <- t(prcomp(Psi_eas, scale=FALSE)$rotation[,1:K_eas])%*%Omega_eas%*%
prcomp(Psi_eas, scale=FALSE)$rotation[,1:K_eas]

Omega_eur <- (eur_chd_common$logOR.SE %o% eur_chd_common$logOR.SE) * eurld_subset
pcOmega_eur <- t(prcomp(Psi_eur, scale=FALSE)$rotation[,1:K_eur])%*%Omega_eur%*%
prcomp(Psi_eur, scale=FALSE)$rotation[,1:K_eur]

Omega_sas <- (sas_chd_common$Est.SE_rescaled %o% sas_chd_common$Est.SE_rescaled) * sasld_subset
pcOmega_sas <- t(prcomp(Psi_sas, scale=FALSE)$rotation[,1:K_sas])%*%Omega_sas%*%
prcomp(Psi_sas, scale=FALSE)$rotation[,1:K_sas]

Omega_ugr <- (afr_chd_common$Est.SE_rescaled %o% afr_chd_common$Est.SE_rescaled) * afrld_subset
pcOmega_ugr <- t(prcomp(Psi_ugr, scale=FALSE)$rotation[,1:K_ugr])%*%Omega_ugr%*%
prcomp(Psi_ugr, scale=FALSE)$rotation[,1:K_ugr]

# Multivariable MR (weighted)
#afr_mvmr_weighted <- summary(lm(afrycomp ~ 0 + afrx1comp + afrx2comp + afrx3comp, weight = 1/diag(pcOmega_afr)^2))
#amr_mvmr_weighted <- summary(lm(amrycomp ~ 0 + amrx1comp + amrx2comp + amrx3comp, weight = 1/diag(pcOmega_amr)^2))
#eas_mvmr_weighted <- summary(lm(easycomp ~ 0 + easx1comp + easx2comp + easx3comp, weight = 1/diag(pcOmega_eas)^2))
#eur_mvmr_weighted <- summary(lm(eurycomp ~ 0 + eurx1comp + eurx2comp + eurx3comp, weight = 1/diag(pcOmega_eur)^2))
#sas_mvmr_weighted <- summary(lm(sasycomp ~ 0 + sasx1comp + sasx2comp + sasx3comp, weight = 1/diag(pcOmega_sas)^2))
#ugr_mvmr_weighted <- summary(lm(afrycomp ~ 0 + ugrx1comp + ugrx2comp + ugrx3comp, weight = 1/diag(pcOmega_ugr)^2))

#saveRDS(afr_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_afr.rds")
#saveRDS(amr_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_amr.rds")
#saveRDS(eas_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_eas.rds")
#saveRDS(eur_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_eur.rds")
#saveRDS(sas_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_sas.rds")
#saveRDS(ugr_mvmr_weighted, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mvmr_weighted_ugr.rds")

# Univariable MR (weighted)
afr_univariable_mr_w <- summary(lm(afrycomp ~ 0 + afrx1comp, weight = 1/diag(pcOmega_afr)^2))
amr_univariable_mr_w <- summary(lm(amrycomp ~ 0 + amrx1comp, weight = 1/diag(pcOmega_amr)^2))
eas_univariable_mr_w <- summary(lm(easycomp ~ 0 + easx1comp, weight = 1/diag(pcOmega_eas)^2))
eur_univariable_mr_w <- summary(lm(eurycomp ~ 0 + eurx1comp, weight = 1/diag(pcOmega_eur)^2))
sas_univariable_mr_w <- summary(lm(sasycomp ~ 0 + sasx1comp, weight = 1/diag(pcOmega_sas)^2))
ugr_univariable_mr_w <- summary(lm(afrycomp ~ 0 + ugrx1comp, weight = 1/diag(pcOmega_ugr)^2))

saveRDS(afr_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.afr.rds")
saveRDS(amr_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.amr.rds")
saveRDS(eas_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.eas.rds")
saveRDS(eur_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.eur.rds")
saveRDS(sas_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.sas.rds")
saveRDS(ugr_univariable_mr_w, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_univ_mr_weighted.ugr.rds")

## MV-PCA (weighted)
mvpca_afr <- mr_mvivw(mr_mvinput(cbind(afrx1comp, afrx2comp, afrx3comp),
cbind(rep(1, length(afrx1comp)), rep(1, length(afrx1comp)), rep(1, length(afrx1comp))),
afrycomp, rep(1, length(afrx1comp)), corr=pcOmega_afr), model="fixed")
mvpca_afr$Estimate; mvpca_afr$StdError

mvpca_est_afr = solve(rbind(afrx1comp, afrx2comp, afrx3comp)%*%solve(pcOmega_afr)%*%
cbind(afrx1comp, afrx2comp, afrx3comp))%*%
rbind(afrx1comp, afrx2comp, afrx3comp)%*%solve(pcOmega_afr)%*%afrycomp
mvpca_se_afr = sqrt(diag(solve(rbind(afrx1comp, afrx2comp, afrx3comp)%*%solve(pcOmega_afr)%*%
cbind(afrx1comp, afrx2comp, afrx3comp))))

z_scores_afr <- mvpca_est_afr / mvpca_se_afr
p_values_afr <- 2 * pnorm(abs(z_scores_afr), lower.tail = FALSE)
p_values_afr

mvpca_est_afr <- matrix(mvpca_est_afr, ncol = 1)
rownames(mvpca_est_afr) <- c("afrx1comp", "afrx2comp", "afrx3comp")
colnames(mvpca_est_afr) <- "beta"

mvpca_se_afr <- matrix(mvpca_se_afr, ncol = 1)
rownames(mvpca_se_afr) <- names(mvpca_se_afr)
colnames(mvpca_se_afr) <- c("se")

mvpca_pval_afr <- matrix(p_values_afr, ncol = 1)
rownames(mvpca_pval_afr) <- c("afrx1comp", "afrx2comp", "afrx3comp")
colnames(mvpca_pval_afr) <- "pval"

res_afr <- cbind(mvpca_est_afr, mvpca_se_afr, mvpca_pval_afr)

write.table(res_afr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_afr.txt", sep = "\t", quote = F, row.names = F)

mvpca_amr <- mr_mvivw(mr_mvinput(cbind(amrx1comp, amrx2comp, amrx3comp),
cbind(rep(1, length(amrx1comp)), rep(1, length(amrx1comp)), rep(1, length(amrx1comp))),
amrycomp, rep(1, length(amrx1comp)), corr=pcOmega_amr), model="fixed")
mvpca_amr$Estimate; mvpca_amr$StdError

mvpca_est_amr = solve(rbind(amrx1comp, amrx2comp, amrx3comp)%*%solve(pcOmega_amr)%*%
cbind(amrx1comp, amrx2comp, amrx3comp))%*%
rbind(amrx1comp, amrx2comp, amrx3comp)%*%solve(pcOmega_amr)%*%amrycomp
mvpca_se_amr = sqrt(diag(solve(rbind(amrx1comp, amrx2comp, amrx3comp)%*%solve(pcOmega_amr)%*%
cbind(amrx1comp, amrx2comp, amrx3comp))))

z_scores_amr <- mvpca_est_amr / mvpca_se_amr
p_values_amr <- 2 * pnorm(abs(z_scores_amr), lower.tail = FALSE)
p_values_amr

mvpca_est_amr <- matrix(mvpca_est_amr, ncol = 1)
rownames(mvpca_est_amr) <- c("amrx1comp", "amrx2comp", "amrx3comp")
colnames(mvpca_est_amr) <- "beta"

mvpca_se_amr <- matrix(mvpca_se_amr, ncol = 1)
rownames(mvpca_se_amr) <- names(mvpca_se_amr)
colnames(mvpca_se_amr) <- c("se")

mvpca_pval_amr <- matrix(p_values_amr, ncol = 1)
rownames(mvpca_pval_amr) <- c("amrx1comp", "amrx2comp", "amrx3comp")
colnames(mvpca_pval_amr) <- "pval"

res_amr <- cbind(mvpca_est_amr, mvpca_se_amr, mvpca_pval_amr)

write.table(res_amr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_amr.txt", sep = "\t", quote = F, row.names = F)

mvpca_eas <- mr_mvivw(mr_mvinput(cbind(easx1comp, easx2comp, easx3comp),
cbind(rep(1, length(easx1comp)), rep(1, length(easx1comp)), rep(1, length(easx1comp))),
easycomp, rep(1, length(easx1comp)), corr=pcOmega_eas), model="fixed")
mvpca_eas$Estimate; mvpca_eas$StdError

mvpca_est_eas = solve(rbind(easx1comp, easx2comp, easx3comp)%*%solve(pcOmega_eas)%*%
cbind(easx1comp, easx2comp, easx3comp))%*%
rbind(easx1comp, easx2comp, easx3comp)%*%solve(pcOmega_eas)%*%easycomp
mvpca_se_eas = sqrt(diag(solve(rbind(easx1comp, easx2comp, easx3comp)%*%solve(pcOmega_eas)%*%
cbind(easx1comp, easx2comp, easx3comp))))

z_scores_eas <- mvpca_est_eas / mvpca_se_eas
p_values_eas <- 2 * pnorm(abs(z_scores_eas), lower.tail = FALSE)
p_values_eas

mvpca_est_eas <- matrix(mvpca_est_eas, ncol = 1)
rownames(mvpca_est_eas) <- c("easx1comp", "easx2comp", "easx3comp")
colnames(mvpca_est_eas) <- "beta"

mvpca_se_eas <- matrix(mvpca_se_eas, ncol = 1)
rownames(mvpca_se_eas) <- names(mvpca_se_eas)
colnames(mvpca_se_eas) <- c("se")

mvpca_pval_eas <- matrix(p_values_eas, ncol = 1)
rownames(mvpca_pval_eas) <- c("easx1comp", "easx2comp", "easx3comp")
colnames(mvpca_pval_eas) <- "pval"

res_eas <- cbind(mvpca_est_eas, mvpca_se_eas, mvpca_pval_eas)

write.table(res_eas, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_eas.txt", sep = "\t", quote = F, row.names = F)

mvpca_eur <- mr_mvivw(mr_mvinput(cbind(eurx1comp, eurx2comp, eurx3comp),
cbind(rep(1, length(eurx1comp)), rep(1, length(eurx1comp)), rep(1, length(eurx1comp))),
eurycomp, rep(1, length(eurx1comp)), corr=pcOmega_eur), model="fixed")
mvpca_eur$Estimate; mvpca_eur$StdError

mvpca_est_eur = solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega_eur)%*%
cbind(eurx1comp, eurx2comp, eurx3comp))%*%
rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega_eur)%*%eurycomp
mvpca_se_eur = sqrt(diag(solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega_eur)%*%
cbind(eurx1comp, eurx2comp, eurx3comp))))

z_scores_eur <- mvpca_est_eur / mvpca_se_eur
p_values_eur <- 2 * pnorm(abs(z_scores_eur), lower.tail = FALSE)
p_values_eur

mvpca_est_eur <- matrix(mvpca_est_eur, ncol = 1)
rownames(mvpca_est_eur) <- c("eurx1comp", "eurx2comp", "eurx3comp")
colnames(mvpca_est_eur) <- "beta"

mvpca_se_eur <- matrix(mvpca_se_eur, ncol = 1)
rownames(mvpca_se_eur) <- names(mvpca_se_eur)
colnames(mvpca_se_eur) <- c("se")

mvpca_pval_eur <- matrix(p_values_eur, ncol = 1)
rownames(mvpca_pval_eur) <- c("eurx1comp", "eurx2comp", "eurx3comp")
colnames(mvpca_pval_eur) <- "pval"

res_eur <- cbind(mvpca_est_eur, mvpca_se_eur, mvpca_pval_eur)

write.table(res_eur, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_eur.txt", sep = "\t", quote = F, row.names = F)

mvpca_sas <- mr_mvivw(mr_mvinput(cbind(sasx1comp, sasx2comp, sasx3comp),
cbind(rep(1, length(sasx1comp)), rep(1, length(sasx1comp)), rep(1, length(sasx1comp))),
sasycomp, rep(1, length(sasx1comp)), corr=pcOmega_sas), model="fixed")
mvpca_sas$Estimate; mvpca_sas$StdError

mvpca_est_sas = solve(rbind(sasx1comp, sasx2comp, sasx3comp)%*%solve(pcOmega_sas)%*%
cbind(sasx1comp, sasx2comp, sasx3comp))%*%
rbind(sasx1comp, sasx2comp, sasx3comp)%*%solve(pcOmega_sas)%*%sasycomp
mvpca_se_sas = sqrt(diag(solve(rbind(sasx1comp, sasx2comp, sasx3comp)%*%solve(pcOmega_sas)%*%
cbind(sasx1comp, sasx2comp, sasx3comp))))

z_scores_sas <- mvpca_est_sas / mvpca_se_sas
p_values_sas <- 2 * pnorm(abs(z_scores_sas), lower.tail = FALSE)
p_values_sas

mvpca_est_sas <- matrix(mvpca_est_sas, ncol = 1)
rownames(mvpca_est_sas) <- c("sasx1comp", "sasx2comp", "sasx3comp")
colnames(mvpca_est_sas) <- "beta"

mvpca_se_sas <- matrix(mvpca_se_sas, ncol = 1)
rownames(mvpca_se_sas) <- names(mvpca_se_sas)
colnames(mvpca_se_sas) <- c("se")

mvpca_pval_sas <- matrix(p_values_sas, ncol = 1)
rownames(mvpca_pval_sas) <- c("sasx1comp", "sasx2comp", "sasx3comp")
colnames(mvpca_pval_sas) <- "pval"

res_sas <- cbind(mvpca_est_sas, mvpca_se_sas, mvpca_pval_sas)

write.table(res_sas, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_sas.txt", sep = "\t", quote = F, row.names = F)

mvpca_ugr <- mr_mvivw(mr_mvinput(cbind(ugrx1comp, ugrx2comp, ugrx3comp),
cbind(rep(1, length(ugrx1comp)), rep(1, length(ugrx1comp)), rep(1, length(ugrx1comp))),
ugrycomp, rep(1, length(ugrx1comp)), corr=pcOmega_ugr), model="fixed")
mvpca_ugr$Estimate; mvpca_ugr$StdError

mvpca_est_ugr = solve(rbind(ugrx1comp, ugrx2comp, ugrx3comp)%*%solve(pcOmega_ugr)%*%
cbind(ugrx1comp, ugrx2comp, ugrx3comp))%*%
rbind(ugrx1comp, ugrx2comp, ugrx3comp)%*%solve(pcOmega_ugr)%*%ugrycomp
mvpca_se_ugr = sqrt(diag(solve(rbind(ugrx1comp, ugrx2comp, ugrx3comp)%*%solve(pcOmega_ugr)%*%
cbind(ugrx1comp, ugrx2comp, ugrx3comp))))

z_scores_ugr <- mvpca_est_ugr / mvpca_se_ugr
p_values_ugr <- 2 * pnorm(abs(z_scores_ugr), lower.tail = FALSE)
p_values_ugr

mvpca_est_ugr <- matrix(mvpca_est_ugr, ncol = 1)
rownames(mvpca_est_ugr) <- c("ugrx1comp", "ugrx2comp", "ugrx3comp")
colnames(mvpca_est_ugr) <- "beta"

mvpca_se_ugr <- matrix(mvpca_se_ugr, ncol = 1)
rownames(mvpca_se_ugr) <- names(mvpca_se_ugr)
colnames(mvpca_se_ugr) <- c("se")

mvpca_pval_ugr <- matrix(p_values_ugr, ncol = 1)
rownames(mvpca_pval_ugr) <- c("ugrx1comp", "ugrx2comp", "ugrx3comp")
colnames(mvpca_pval_ugr) <- "pval"

res_ugr <- cbind(mvpca_est_ugr, mvpca_se_ugr, mvpca_pval_ugr)

write.table(res_ugr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/hmgcr_chd_mvmr_pca_ugr.txt", sep = "\t", quote = F, row.names = F)

## Compare using wald ratio for a single snp in univariable MR
#library(TwoSampleMR)
#i <- which.min(hmgcr_ldlc_afr_common$Score.pval)
#afr_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_afr$Est[i], afr_chd_common$Est_rescaled[i], hmgcr_ldlc_afr$Est.SE[i], afr_chd_common$Est.SE_rescaled[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(afr_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_afr.rds")

#i <- which.min(hmgcr_ldlc_amr_common$P)
#amr_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_amr$BETA[i], amr_chd_common$BETA_rescaled[i], hmgcr_ldlc_amr$BETA[i], amr_chd_common$SE_rescaled[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(amr_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_amr.rds")

#i <- which.min(hmgcr_ldlc_eas_common$P)
#eas_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_eas$BETA[i], eas_chd_common$BETA_rescaled[i], hmgcr_ldlc_eas$BETA[i], eas_chd_common$SE_rescaled[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(eas_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_eas.rds")

#i <- which.min(hmgcr_ldlc_eur_common$P_BOLT_LMM_INF)
#eur_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_eur$BETA[i], eur_chd_common$logOR[i], hmgcr_ldlc_eur$BETA[i], eur_chd_common$logOR.SE[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(eur_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_eur.rds")

#i <- which.min(hmgcr_ldlc_sas_common$Score.pval)
#sas_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_sas$Est[i], sas_chd_common$Est_rescaled[i], hmgcr_ldlc_sas$Est.SE[i], sas_chd_common$Est.SE_rescaled[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(sas_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_sas.rds")

#i <- which.min(hmgcr_ldlc_afr_common$Score.pval)
#ugr_wald_ratio_univ_mr <- mr_wald_ratio(hmgcr_ldlc_ugr$Est[i], afr_chd_common$Est_rescaled[i], hmgcr_ldlc_ugr$Est.SE[i], afr_chd_common$Est.SE_rescaled[i]) %>% as_tibble %>% mutate(t=b/se)

#saveRDS(ugr_wald_ratio_univ_mr, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/pc_based_cis_mvmr/hmgcr_chd_mr_wald_ratio_ugr.rds")
