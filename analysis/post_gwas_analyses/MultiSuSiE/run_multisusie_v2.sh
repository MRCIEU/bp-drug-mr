#!/bin/bash
#SBATCH --job-name=msusie
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:0
#SBATCH --output=logs/%A_%a.out
#SBATCH --no-requeue
#SBATCH --account=smed001801

module add languages/python/3.14.3

set -eo pipefail

source ~/.bashrc
conda activate MultiSuSiE

SNP=$(awk "NR==${SLURM_ARRAY_TASK_ID}" snps_clean.txt)

echo "Task: $SLURM_ARRAY_TASK_ID"
echo "SNP: '$SNP'"

if [ -z "$SNP" ]; then
    echo "ERROR: empty SNP"
    exit 1
fi

python run_multisusie_v2.py \
    --snp "$SNP" \
    --window 1000000 \
    --max_snps 5000 \
    --outdir results
