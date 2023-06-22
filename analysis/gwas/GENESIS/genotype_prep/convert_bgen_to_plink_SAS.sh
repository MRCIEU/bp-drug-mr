#!/bin/bash
#SBATCH --job-name=sas_plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/plink/2.00

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/plink

echo chr="${SLURM_ARRAY_TASK_ID}"

plink2 --bgen $data_dir/SAS/filtered/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3_SAS.bgen --sample $data_dir/SAS/filtered/ukb_imp_filtered_SAS.sample --make-bed --out $out_dir/SAS/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3_SAS
