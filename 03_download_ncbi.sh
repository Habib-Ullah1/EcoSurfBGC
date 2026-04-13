#!/bin/bash
#SBATCH --job-name=ecosurf_dl
#SBATCH --output=/data/habib/EcoSurfBGC/logs/download_logs/dl_%A_%a.out
#SBATCH --error=/data/habib/EcoSurfBGC/logs/download_logs/dl_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

FASTA_DIR=/data/habib/EcoSurfBGC/genomes/fasta
CSV=/data/habib/EcoSurfBGC/genome_lists/ncbi_accessions.csv
LOG=/data/habib/EcoSurfBGC/logs/download_logs/ncbi_summary.log

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ${CSV})
IMG_ID=$(echo $LINE | cut -d',' -f1)
ACC=$(echo $LINE | cut -d',' -f3)
OUTFILE=${FASTA_DIR}/${IMG_ID}.fasta

if [ -f "${OUTFILE}" ] && [ -s "${OUTFILE}" ]; then
  echo "SKIP: ${IMG_ID}" >> ${LOG}
  exit 0
fi

FTP_PATH=$(efetch -db assembly -id ${ACC} -format fasta -location ftp 2>/dev/null | head -1)

if [[ "${FTP_PATH}" == *"ftp"* ]]; then
  BASE=$(echo ${FTP_PATH} | tr -d '[:space:]')
  NAME=$(basename ${BASE})
  URL="${BASE}/${NAME}_genomic.fna.gz"
  wget -q -O ${OUTFILE}.gz "${URL}" && gunzip -f ${OUTFILE}.gz
fi

if [ -f "${OUTFILE}" ] && [ -s "${OUTFILE}" ]; then
  FIRST=$(head -c 1 ${OUTFILE})
  if [ "${FIRST}" == ">" ]; then
    echo "SUCCESS: ${IMG_ID} ${ACC}" >> ${LOG}
    exit 0
  fi
fi

echo "FAILED: ${IMG_ID} ${ACC}" >> ${LOG}
exit 1
EOF
