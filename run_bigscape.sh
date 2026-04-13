#!/bin/bash
#SBATCH --job-name=bigscape
#SBATCH --output=/data/habib/EcoSurfBGC/logs/bigscape_%j.out
#SBATCH --error=/data/habib/EcoSurfBGC/logs/bigscape_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=debug
#SBATCH --account=habib
#SBATCH --qos=user_habib

source /home/habib/miniconda3/etc/profile.d/conda.sh
conda activate bigscape_env

BASE="/data/habib/EcoSurfBGC"
PFAM_DIR="/home/habib/miniconda3/envs/antismash/lib/python3.9/site-packages/antismash/databases/pfam/35.0"
INPUT="${BASE}/bigscape_input"
OUTPUT="${BASE}/bigscape_output"

mkdir -p "${OUTPUT}"

echo "Starting BiG-SCAPE at $(date)"
echo "Input GBK files: $(ls ${INPUT} | wc -l)"
echo "Pfam directory: ${PFAM_DIR}"

bigscape \
  -i "${INPUT}" \
  -o "${OUTPUT}" \
  --pfam_dir "${PFAM_DIR}" \
  -c 16 \
  --mibig \
  --mix \
  --include_singletons \
  --cutoffs 0.3 0.5 0.7 \
  --mode auto \
  -l "EcoSurfBGC_2057genomes" \
  -v \
  2>&1 | tee "${BASE}/logs/bigscape_run.log"

echo "BiG-SCAPE finished at $(date)"
echo "Exit code: $?"
