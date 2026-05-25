#!/bin/bash
#SBATCH --job-name=guide_align_qc
#SBATCH --partition=engreitz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err

set -euo pipefail

mkdir -p logs
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gem3-map

python scripts/run_guide_alignment_qc.py \
  --guides-tsv data/guides.tsv \
  --reference-fasta data/IGVFFI0653VCGH.fasta \
  --chromsizes data/hg38.chrom.sizes \
  --outdir results/guide_alignment_qc \
  --threads 8 \
  --pam NGG \
  --add-leading-g
