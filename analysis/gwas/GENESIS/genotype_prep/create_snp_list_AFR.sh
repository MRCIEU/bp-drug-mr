#!/bin/bash
#SBATCH --job-name=afr-snp-list
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28 
#SBATCH --time=100:00:00
#SBATCH --array=1-22
#SBATCH --account=smed001801

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB

echo chr="${SLURM_ARRAY_TASK_ID}"

awk '{if ($14>0.01 && $17>0.5){ print $2;}}' $data_dir/AFR/snp_stats/snp_stats_chr"${SLURM_ARRAY_TASK_ID}"_AFR.txt > $out_dir/AFR/filtered/snp_list/snp_list_chr_"${SLURM_ARRAY_TASK_ID}"_AFR.txt
