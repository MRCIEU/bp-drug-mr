#!/bin/bash
#SBATCH --job-name=eas-snp-list
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28 
#SBATCH --time=30:00:00
#SBATCH --array=1-22
#SBATCH --account=smed001801

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB

echo chr="${SLURM_ARRAY_TASK_ID}"

awk '{if ($14>0.01 && $17>0.5){ print $2;}}' $data_dir/EAS/snp_stats/snp_stats_chr_"${SLURM_ARRAY_TASK_ID}"_EAS.txt > $out_dir/EAS/filtered/snp_list/snp_list_chr_"${SLURM_ARRAY_TASK_ID}"_EAS.txt
