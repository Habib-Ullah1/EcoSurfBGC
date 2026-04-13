#!/bin/bash
#SBATCH --job-name=ecosurf_as
#SBATCH --output=/data/habib/EcoSurfBGC/logs/antismash_logs/as_%A_%a.out
#SBATCH --error=/data/habib/EcoSurfBGC/logs/antismash_logs/as_%A_%a.err
#SBATCH --array=1-5%5
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=debug
#SBATCH --account=habib
#SBATCH --qos=user_habib

# ── Activate antiSMASH environment ───────────────────────
source /home/habib/miniconda3/etc/profile.d/conda.sh
conda activate antismash

# ── Paths ─────────────────────────────────────────────────
BASE="/data/habib/EcoSurfBGC"
ID_LIST="${BASE}/genome_lists/bacteria_fasta_ids.txt"
FASTA_DIR="${BASE}/genomes/fasta"
OUT_DIR="${BASE}/antismash_output"
LOGFILE="${BASE}/logs/antismash_logs/antismash_summary.log"

mkdir -p "${OUT_DIR}"
mkdir -p "${BASE}/logs/antismash_logs"

# ── Get genome ID for this array task ────────────────────
GENOME_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ID_LIST}")

if [ -z "${GENOME_ID}" ]; then
  echo "No genome ID for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

FASTA="${FASTA_DIR}/${GENOME_ID}.fasta"
GENOME_OUTDIR="${OUT_DIR}/${GENOME_ID}"

# ── Skip if already completed successfully ────────────────
if [ -d "${GENOME_OUTDIR}" ] && [ -f "${GENOME_OUTDIR}/index.html" ]; then
  echo "ALREADY_DONE: ${GENOME_ID}" >> "${LOGFILE}"
  exit 0
fi

# ── Skip if FASTA missing ─────────────────────────────────
if [ ! -f "${FASTA}" ] || [ ! -s "${FASTA}" ]; then
  echo "NO_FASTA: ${GENOME_ID}" >> "${LOGFILE}"
  exit 1
fi

# ── Run antiSMASH ─────────────────────────────────────────
echo "START: ${GENOME_ID} task ${SLURM_ARRAY_TASK_ID} at $(date)" >> "${LOGFILE}"

antismash \
  --taxon bacteria \
  --genefinding-tool prodigal \
  --cb-knownclusters \
  --asf \
  --pfam2go \
  --cpus 4 \
  --output-dir "${GENOME_OUTDIR}" \
  "${FASTA}" \
  2>> "${BASE}/logs/antismash_logs/as_${GENOME_ID}.err"

# ── Check result ──────────────────────────────────────────
EXIT_CODE=$?
if [ ${EXIT_CODE} -eq 0 ] && [ -f "${GENOME_OUTDIR}/index.html" ]; then
  echo "SUCCESS: ${GENOME_ID} at $(date)" >> "${LOGFILE}"
else
  echo "FAILED: ${GENOME_ID} exit=${EXIT_CODE} at $(date)" >> "${LOGFILE}"
  exit 1
fi
