#!/bin/bash
#SBATCH --job-name=afr_generate_merged_genotype_files
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28 
#SBATCH--time=20:00:00
#SBATCH --account=smed001801

module load apps/plink/1.90

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/plink/AFR/exclusion_multiallelic_snps
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/plink/AFR

plink --bed $data_dir/ukb_imp_filtered_excl_multiallelic_snps_chr1_v3_AFR.bed --bim $data_dir/ukb_imp_filtered_excl_multiallelic_snps_chr1_v3_AFR.bim --fam $data_dir/ukb_imp_filtered_excl_multiallelic_snps_chr1_v3_AFR.fam --merge-list $data_dir/list_chr2-22_AFR.txt --make-bed --out $out_dir/merged_chr1_22/chr1-22_merged
