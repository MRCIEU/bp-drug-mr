#!/bin/bash
#SBATCH --job-name=qctool_sas_snp_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=330:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --account=smed001801 
#SBATCH --array=1-22

module load apps/qctool/2.2.0

export bgen_dir=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/raw_downloaded/bgen
export sample_dir=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/id_mapping
export dat_dir=/user/work/ac14629/MRC_network_project/scripts/GENESIS
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/dosage_bgen/UKB

echo chr="${SLURM_ARRAY_TASK_ID}"

qctool \
-s $sample_dir/data.chr1-22.sample \
-g $bgen_dir/ukb_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen \
-og $out_dir/SAS/ukb_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3_SAS.bgen \
-os $out_dir/SAS/ukb_imp_SAS.sample \
-incl-samples sample_ID_keep_SAS.txt \
-snp-stats \
-osnp $out_dir/SAS/snp_stats/snp_stats_chr_"${SLURM_ARRAY_TASK_ID}"_SAS.txt
