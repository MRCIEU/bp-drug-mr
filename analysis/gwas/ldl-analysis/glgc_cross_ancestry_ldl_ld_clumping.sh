#!/bin/bash
#SBATCH --job-name=plink-ldl-clumping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=05:00:0
#SBATCH --account=smed001801

module load apps/plink/1.90

plink --bfile /user/work/ac14629/MRC_network_project/data/UKB/UKB10K/plink/merged_chr1_22/chr1-22_merged --clump-p1 1 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 --clump /user/work/ac14629/MRC_network_project/data/GLGC/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1  --out /user/work/ac14629/MRC_network_project/results/ldl-analysis/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1_clumped
