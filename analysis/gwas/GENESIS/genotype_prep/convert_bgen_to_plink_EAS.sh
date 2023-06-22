#!/bin/bash
#SBATCH --job-name=eas_plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/plink/2.3

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/plink

echo chr="${SLURM_ARRAY_TASK_ID}"

plink2 --bgen $data_dir/EAS/filtered/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3_EAS.bgen --sample $data_dir/EAS/filtered/ukb_imp_filtered_EAS.sample --make-bed --out $out_dir/EAS/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3_EAS
