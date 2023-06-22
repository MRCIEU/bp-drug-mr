#!/bin/bash
#SBATCH --job-name=ukb10k-bgen-plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/plink/2.00

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen/filtered
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink

echo chr="${SLURM_ARRAY_TASK_ID}"

plink2 --bgen $data_dir/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen --sample $data_dir/ukb_imp_filtered.sample --make-bed --out $out_dir/ukb10k_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3
