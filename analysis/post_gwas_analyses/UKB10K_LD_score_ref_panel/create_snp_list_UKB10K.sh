#!/bin/bash
#SBATCH --job-name=ukb10k-snp-list
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28 
#SBATCH --time=10:00:00
#SBATCH --array=1-22
#SBATCH --account=smed001801

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen

echo chr="${SLURM_ARRAY_TASK_ID}"

awk '{if ($14>0.01 && $17>0.5){ print $2;}}' $data_dir/snp_stats/snp_stats_chr_"${SLURM_ARRAY_TASK_ID}".txt > $out_dir/filtered/snp_list/snp_list_chr_"${SLURM_ARRAY_TASK_ID}".txt
