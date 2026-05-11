#!/bin/bash
#SBATCH --job-name=plink_vcor_afr
#SBATCH --array=1-$(wc -l < snps_clean.txt)
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --account=smed001801

module load apps/plink2/2.00a68LM

mkdir -p vcor_out/AFR logs

# SNP for this task
SNP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" snps_clean.txt)

# Look up chromosome + position
read CHR BP <<< $(awk -v s="$SNP" '$1==s {print $2, $3}' AFR_snp_pos.txt)

# Safety check (skip missing SNPs)
if [ -z "$CHR" ]; then
  echo "SNP not found in bim: $SNP"
  exit 0
fi

# Compute ±1Mb window
START=$((BP - 1000000))
END=$((BP + 1000000))
if [ "$START" -lt 0 ]; then START=0; fi

RANGE=$(mktemp)

echo "${CHR}:${START}-${END}" > "$RANGE"

# Run PLINK2 LD correlation 
plink2 \
  --bfile AFR_modified \
  --extract range "$RANGE" \
  --r2-phased \
  --ld-window 99999
  --ld-window-kb 1000
  --ld-window-r2 0
  --out vcor_out/AFR/AFR_${SNP}

rm "$RANGE"
