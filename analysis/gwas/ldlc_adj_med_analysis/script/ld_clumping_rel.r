#Outline
#This script will LD clump the GWAS outputs for relative medication adjustment and create a table hightlighting the number of clumped SNPs for each GWAS output

#Read in relevant packages
library(TwoSampleMR)
library(data.table)

#Set working directory
setwd('/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis')

#Read GWAS outputs
ldlc_unadj <- fread("GENESIS_assoc_ldlc_unadjusted_afr_out.txt", header = T)
ldlc_10 <- fread("GENESIS_assoc_ldlc_10_afr_out.txt", header = T)
ldlc_20 <- fread("GENESIS_assoc_ldlc_20_afr_out.txt", header = T)
ldlc_30 <- fread("GENESIS_assoc_ldlc_30_afr_out.txt", header = T)
ldlc_40 <- fread("GENESIS_assoc_ldlc_40_afr_out.txt", header = T)
ldlc_50 <- fread("GENESIS_assoc_ldlc_50_afr_out.txt", header = T)
ldlc_60 <- fread("GENESIS_assoc_ldlc_60_afr_out.txt", header = T)
ldlc_70 <- fread("GENESIS_assoc_ldlc_70_afr_out.txt", header = T)
ldlc_80 <- fread("GENESIS_assoc_ldlc_80_afr_out.txt", header = T)

# Filter GWAS outputs by P-value 
ldlc_unadj <- subset(ldlc_unadj, Score.pval < 5e-8)
ldlc_10 <- subset(ldlc_10, Score.pval < 5e-8)
ldlc_20 <- subset(ldlc_20, Score.pval < 5e-8)
ldlc_30 <- subset(ldlc_30, Score.pval < 5e-8)
ldlc_40 <- subset(ldlc_40, Score.pval < 5e-8)
ldlc_50 <- subset(ldlc_50, Score.pval < 5e-8)
ldlc_60 <- subset(ldlc_60, Score.pval < 5e-8)
ldlc_70 <- subset(ldlc_70, Score.pval < 5e-8)
ldlc_80 <- subset(ldlc_80, Score.pval < 5e-8)
                     
#Reformat GWAS outputs
ldlc_unadj <- format_data(ldlc_unadj, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_10 <- format_data(ldlc_10, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_20 <- format_data(ldlc_20, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_30 <- format_data(ldlc_30, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_40 <- format_data(ldlc_40, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_50 <- format_data(ldlc_50, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_60 <- format_data(ldlc_60, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_70 <- format_data(ldlc_70, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_80 <- format_data(ldlc_80, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")

# Perform LD clumping
ldlc_unadj <- clump_data(ldlc_unadj, pop = "AFR")
ldlc_10 <- clump_data(ldlc_10, pop = "AFR")
ldlc_20 <- clump_data(ldlc_20, pop = "AFR")
ldlc_30 <- clump_data(ldlc_30, pop = "AFR")
ldlc_40 <- clump_data(ldlc_40, pop = "AFR")
ldlc_50 <- clump_data(ldlc_50, pop = "AFR")
ldlc_60 <- clump_data(ldlc_60, pop = "AFR")
ldlc_70 <- clump_data(ldlc_70, pop = "AFR")
ldlc_80 <- clump_data(ldlc_80, pop = "AFR")

# Get number of clumped SNPs 
num_clumped_snps_ldlc_unadj <- print(nrow(ldlc_unadj))
num_clumped_snps_ldlc_10 <- print(nrow(ldlc_10))
num_clumped_snps_ldlc_20 <- print(nrow(ldlc_20))
num_clumped_snps_ldlc_30 <- print(nrow(ldlc_30))
num_clumped_snps_ldlc_40 <- print(nrow(ldlc_40))
num_clumped_snps_ldlc_50 <- print(nrow(ldlc_50))
num_clumped_snps_ldlc_60 <- print(nrow(ldlc_60))
num_clumped_snps_ldlc_70 <- print(nrow(ldlc_70))
num_clumped_snps_ldlc_80 <- print(nrow(ldlc_80))

# Organise results 
results <- data.frame(
  Model = c("Unadjusted", "10%", "20%", "30%", "40", "50", "60", "70", "80"),
  Number_LD_clumped_SNPs = c(num_clumped_snps_ldlc_unadj, num_clumped_snps_ldlc_10, num_clumped_snps_ldlc_20, num_clumped_snps_ldlc_30, num_clumped_snps_ldlc_40, num_clumped_snps_ldlc_50, num_clumped_snps_ldlc_60, num_clumped_snps_ldlc_70, num_clumped_snps_ldlc_80)
  )

  write.table(results, "LD_clumped_res_rel_afr.txt", sep = "\t", col.names=T,row.names=F,quote=F)
