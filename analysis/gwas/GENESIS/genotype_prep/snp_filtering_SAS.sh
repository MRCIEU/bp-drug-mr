#!/bin/bash
#SBATCH --job-name=sas_filter_bgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/qctool/2.2.0

export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB
export snp_list_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB/SAS/filtered/snp_list
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB

echo chr="${SLURM_ARRAY_TASK_ID}"

qctool \
-s $data_dir/SAS/ukb_imp_SAS.sample \
-g $data_dir/SAS/ukb_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3_SAS.bgen \
-og $out_dir/SAS/filtered/ukb_imp_filtered_chr"${SLURM_ARRAY_TASK_ID}"_v3_SAS.bgen \
-os $out_dir/SAS/filtered/ukb_imp_filtered_SAS.sample \
-incl-rsids $snp_list_dir/snp_list_chr_"${SLURM_ARRAY_TASK_ID}"_SAS.txt 
