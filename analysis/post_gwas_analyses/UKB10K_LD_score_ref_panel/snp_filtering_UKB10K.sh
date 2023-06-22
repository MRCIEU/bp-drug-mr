#!/bin/bash
#SBATCH --job-name=ukb10k_filter_bgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/qctool/2.2.0

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen
export snp_list_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen/filtered/snp_list
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen

echo chr="${SLURM_ARRAY_TASK_ID}"

qctool \
-s $data_dir/ukb10k_imp.sample \
-g $data_dir/ukb10k_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen \
-og $out_dir/filtered/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen \
-os $out_dir/filtered/ukb_imp_filtered.sample \
-incl-rsids $snp_list_dir/snp_list_chr_"${SLURM_ARRAY_TASK_ID}".txt 
