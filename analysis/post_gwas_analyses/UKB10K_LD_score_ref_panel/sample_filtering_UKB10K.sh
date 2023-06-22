#!/bin/bash
#SBATCH --job-name=ukb10k_snp_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-22
#SBATCH --account=smed001801

module load apps/qctool/2.2.0

export bgen_dir=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/raw_downloaded/bgen
export sample_dir=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/id_mapping
export data_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K
export out_dir=/user/work/ac14629/MRC_network_project/data/UKB/UKB10K/dosage_bgen

echo chr="${SLURM_ARRAY_TASK_ID}"

qctool \
-s $sample_dir/data.chr1-22.sample \
-g $bgen_dir/ukb_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen \
-og $out_dir/ukb10k_imp_chr"${SLURM_ARRAY_TASK_ID}"_v3.bgen \
-os $out_dir/ukb10k_imp.sample \
-incl-samples $data_dir/keep_ids.txt \
-snp-stats \
-osnp $out_dir/snp_stats/snp_stats_chr_"${SLURM_ARRAY_TASK_ID}".txt
