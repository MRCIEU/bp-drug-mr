#!/bin/bash
#SBATCH --job-name=ukb10k_generate_merged_genotype_files
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28 
#SBATCH--time=20:00:00
#SBATCH --account=smed001801

module load apps/plink/1.90

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink

plink --bed $data_dir/ukb10k_imp_filtered_chr1_v3.bed --bim $data_dir/ukb10k_imp_filtered_chr1_v3.bim --fam $data_dir/ukb10k_imp_filtered_chr1_v3.fam --merge-list $data_dir/list_chr2-22_ukb10k.txt --make-bed --out $out_dir/merged_chr1_22/chr1-22_merged
