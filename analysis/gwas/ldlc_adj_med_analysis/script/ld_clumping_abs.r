#Outline
#This script will LD clump the GWAS outputs for absolute medication adjustment and create a table hightlighting the number of clumped SNPs for each GWAS output

#Read in relevant packages
library(TwoSampleMR)
library(data.table)

#Set working directory
setwd('/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis')

#Read GWAS outputs
ldlc_0.3 <- fread("GENESIS_assoc_ldlc_0.3_afr_out.txt", header = T)
ldlc_0.5 <- fread("GENESIS_assoc_ldlc_0.5_afr_out.txt", header = T)
ldlc_0.7 <- fread("GENESIS_assoc_ldlc_0.7_afr_out.txt", header = T)
ldlc_0.9 <- fread("GENESIS_assoc_ldlc_0.9_afr_out.txt", header = T)
ldlc_1.0 <- fread("GENESIS_assoc_ldlc_1.0_afr_out.txt", header = T)
ldlc_1.2<- fread("GENESIS_assoc_ldlc_1.2_afr_out.txt", header = T)
ldlc_1.4 <- fread("GENESIS_assoc_ldlc_1.4_afr_out.txt", header = T)
ldlc_1.6 <- fread("GENESIS_assoc_ldlc_1.6_afr_out.txt", header = T)
ldlc_1.8 <- fread("GENESIS_assoc_ldlc_1.8_afr_out.txt", header = T)
ldlc_2.0 <- fread("GENESIS_assoc_ldlc_2.0_afr_out.txt", header = T)

# Filter GWAS outputs by P-value 
ldlc_0.3 <- subset(ldlc_0.3, Score.pval < 5e-8)
ldlc_0.5 <- subset(ldlc_0.5, Score.pval < 5e-8)
ldlc_0.7 <- subset(ldlc_0.7, Score.pval < 5e-8)
ldlc_0.9 <- subset(ldlc_0.9, Score.pval < 5e-8)
ldlc_1.0 <- subset(ldlc_1.0, Score.pval < 5e-8)
ldlc_1.2 <- subset(ldlc_1.2, Score.pval < 5e-8)
ldlc_1.4 <- subset(ldlc_1.4, Score.pval < 5e-8)
ldlc_1.6 <- subset(ldlc_1.6, Score.pval < 5e-8)
ldlc_1.8 <- subset(ldlc_1.8, Score.pval < 5e-8)
ldlc_2.0 <- subset(ldlc_2.0, Score.pval < 5e-8)                     

#Reformat GWAS outputs
ldlc_0.3 <- format_data(ldlc_0.3, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_0.5 <- format_data(ldlc_0.5, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_0.7 <- format_data(ldlc_0.7, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_0.9 <- format_data(ldlc_0.9, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_1.0 <- format_data(ldlc_1.0, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_1.2 <- format_data(ldlc_1.2, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_1.4 <- format_data(ldlc_1.4, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_1.6 <- format_data(ldlc_1.6, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_1.8 <- format_data(ldlc_1.8, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")
ldlc_2.0 <- format_data(ldlc_2.0, type = "exposure", header = T, snp_col = "variant.id", chr_col = "chr", pos_col = "pos", samplesize_col = "n.obs", eaf_col = "freq", beta_col = "Est", se_col = "Est.SE", effect_allele_col = "effect.allele", other_allele_col = "other.allele", pval_col = "Score.pval")

# Perform LD clumping
ldlc_0.3 <- clump_data(ldlc_0.3, pop = "AFR")
ldlc_0.5 <- clump_data(ldlc_0.5, pop = "AFR")
ldlc_0.7 <- clump_data(ldlc_0.7, pop = "AFR")
ldlc_0.9 <- clump_data(ldlc_0.9, pop = "AFR")
ldlc_1.0 <- clump_data(ldlc_1.0, pop = "AFR")
ldlc_1.2 <- clump_data(ldlc_1.2, pop = "AFR")
ldlc_1.4 <- clump_data(ldlc_1.4, pop = "AFR")
ldlc_1.6 <- clump_data(ldlc_1.6, pop = "AFR")
ldlc_1.8 <- clump_data(ldlc_1.8, pop = "AFR")
ldlc_2.0 <- clump_data(ldlc_2.0, pop = "AFR")

# Get number of clumped SNPs 
num_clumped_snps_ldlc_0.3 <- print(nrow(ldlc_0.3))
num_clumped_snps_ldlc_0.5 <- print(nrow(ldlc_0.5))
num_clumped_snps_ldlc_0.7 <- print(nrow(ldlc_0.7))
num_clumped_snps_ldlc_0.9 <- print(nrow(ldlc_0.9))
num_clumped_snps_ldlc_1.0 <- print(nrow(ldlc_1.0))
num_clumped_snps_ldlc_1.2 <- print(nrow(ldlc_1.2))
num_clumped_snps_ldlc_1.4 <- print(nrow(ldlc_1.4))
num_clumped_snps_ldlc_1.6 <- print(nrow(ldlc_1.6))
num_clumped_snps_ldlc_1.8 <- print(nrow(ldlc_1.8))
num_clumped_snps_ldlc_2.0 <- print(nrow(ldlc_2.0))

# Organise results 
results <- data.frame(
  Model = c("0.3", "0.5", "0.7", "0.9", "1.0", "1.2", "1.4", "1.6", "1.8","2.0"),
  Number_LD_clumped_SNPs = c(num_clumped_snps_ldlc_0.3, num_clumped_snps_ldlc_0.5, num_clumped_snps_ldlc_0.7, num_clumped_snps_ldlc_0.9, num_clumped_snps_ldlc_1.0, num_clumped_snps_ldlc_1.2, num_clumped_snps_ldlc_1.4, num_clumped_snps_ldlc_1.6, num_clumped_snps_ldlc_2.0)
)

  write.table(results, "LD_clumped_res_abs_afr.txt", sep = "\t", col.names=T,row.names=F,quote=F)
