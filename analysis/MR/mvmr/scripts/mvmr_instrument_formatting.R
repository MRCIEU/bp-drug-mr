### Identifying unique instruments from each lipid fraction and blood pressure traits 

# Read in relevant packages 
library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(ieugwasr)
library(gwasvcf)
library(parallel)

# Read in data for lipid fractions 
ldlc_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/ldlc/METAL_LDLC_clumped_common_snps_formatted.txt")
hdl_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/hdl/METAL_HDL_clumped_common_snps_formatted.txt")
tg_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/tg/METAL_TG_clumped_common_snps_formatted.txt")

# Read in data for blood pressure traits 
sbp_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/sbp/METAL_SBP_clumped_common_snps_formatted.txt")
dbp_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/dbp/METAL_DBP_clumped_common_snps_formatted.txt")
pp_clumped <- fread("/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/pp/METAL_PP_clumped_common_snps_formatted.txt")

# Combine data for lipids and BP with SNPs present in all three datasets 
ldlc_clumped$trait <- "LDLC"
hdl_clumped$trait <- "HDL"
tg_clumped$trait <- "TG"
sbp_clumped$trait <- "SBP"
dbp_clumped$trait <- "DBP"
pp_clumped$trait <- "PP"

lipid_com <- rbind(ldlc_clumped, hdl_clumped, tg_clumped)
bp_com <- rbind(sbp_clumped, dbp_clumped, pp_clumped)

lipid_com <- lipid_com %>% arrange(Chromosome)
bp_com <- bp_com %>% arrange(Chromosome)

lipid_com <- lipid_com %>%
  distinct(Chromosome, Position, .keep_all = TRUE)

bp_com <- bp_com %>%
  distinct(Chromosome, Position, .keep_all = TRUE)

# Match chr and pos to rsid
gwasvcf::set_bcftools("/mnt/storage/software/apps/BCFTOOLS/bcftools/bcftools")
bcftools <- options()[["tools_bcftools"]]

vcffile <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.vcf.gz"

chromosome <- lipid_com$Chromosome
position <- lipid_com$Position

tmp_positions <- tempfile()
write.table(data.frame(chromosome, position), file = tmp_positions, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

cmd <- glue("{bcftools} query -f '%ID %CHROM %POS\\n' -R {tmp_positions} {vcffile} > {tmp_positions}.rsids")
system(cmd)

extracted_data_lipid <- read.table(glue("{tmp_positions}.rsids"), header = FALSE, col.names = c("rsid", "Chromosome", "Position"), stringsAsFactors = FALSE)
extr_lipid <- merge(lipid_com, extracted_data_lipid, by = c("Chromosome", "Position"))

chromosome <- bp_com$Chromosome
position <- bp_com$Position

tmp_positions <- tempfile()
write.table(data.frame(chromosome, position), file = tmp_positions, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

cmd <- glue("{bcftools} query -f '%ID %CHROM %POS\\n' -R {tmp_positions} {vcffile} > {tmp_positions}.rsids")
system(cmd)

extracted_data_bp <- read.table(glue("{tmp_positions}.rsids"), header = FALSE, col.names = c("rsid", "Chromosome", "Position"), stringsAsFactors = FALSE)
extr_bp <- merge(bp_com, extracted_data_bp, by = c("Chromosome", "Position"))

#LD clumping
lipid_com_c <- ld_clump(
    dplyr::tibble(rsid=extr_lipid$rsid, pval=extr_lipid$`P-value`),
    plink_bin = "/mnt/storage/software/languages/R/4.2.1/lib64/R/library/plinkbinr/bin/plink_Linux",
    bfile = "/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink/merged_chr1_22/chr1-22_merged"
)

bp_com_c <- ld_clump(
    dplyr::tibble(rsid=extr_bp$rsid, pval=extr_bp$`P-value`),
    plink_bin = "/mnt/storage/software/languages/R/4.2.1/lib64/R/library/plinkbinr/bin/plink_Linux",
    bfile = "/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink/merged_chr1_22/chr1-22_merged"
)

# Join and extract relevant columns
lipid_j <- left_join(lipid_com_c, extr_lipid, by = "rsid")
bp_j <- left_join(bp_com_c, extr_bp, by = "rsid")

#Save output 
write.table(lipid_j, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/mvmr_instruments_lipids.txt", sep = " ", quote = F, row.names = F) 
write.table(bp_j, "/user/work/ac14629/MRC_network_project/results/MAIN_ANALYSIS/mvmr/mvmr_instruments_bp.txt", sep = " ", quote = F, row.names = F)


