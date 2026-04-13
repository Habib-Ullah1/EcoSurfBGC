#!/bin/bash
# Submits antiSMASH batch 2 (JGI genomes) in chunks of 5

BASE="/data/habib/EcoSurfBGC"
ID_LIST="${BASE}/genome_lists/jgi_antismash_ids.txt"
TEMPLATE="${BASE}/scripts/04_run_antismash.sh"
TOTAL=$(wc -l < "${ID_LIST}")
CHUNK=5
START=1

echo "Batch 2: ${TOTAL} JGI genomes"
echo "Submitting chunks of ${CHUNK} starting from task ${START}"

while [ ${START} -le ${TOTAL} ]; do
  END=$((START + CHUNK - 1))
  if [ ${END} -gt ${TOTAL} ]; then
    END=${TOTAL}
  fi

  # Wait until 4 or fewer jobs in queue
  while true; do
    CURRENT_JOBS=$(squeue -u habib -h | wc -l)
    if [ ${CURRENT_JOBS} -le 4 ]; then
      break
    fi
    echo "$(date): ${CURRENT_JOBS} jobs in queue, waiting 60s..."
    sleep 60
  done

  # Create chunk script - use jgi ID list and update array range
  sed "s|#SBATCH --array=.*|#SBATCH --array=${START}-${END}%5|" \
    "${TEMPLATE}" | \
  sed "s|bacteria_fasta_ids.txt|jgi_antismash_ids.txt|" \
    > /tmp/antismash_batch2_chunk.sh

  # Verify
  ARRAY_LINE=$(grep "^#SBATCH --array" /tmp/antismash_batch2_chunk.sh)
  ID_LINE=$(grep "ID_LIST" /tmp/antismash_batch2_chunk.sh | head -1)
  echo "$(date): Submitting ${START}-${END} | ${ARRAY_LINE}"

  sbatch /tmp/antismash_batch2_chunk.sh
  SUBMIT_STATUS=$?

  if [ ${SUBMIT_STATUS} -ne 0 ]; then
    echo "ERROR: Submission failed for tasks ${START}-${END}. Retrying in 120s..."
    sleep 120
    continue
  fi

  echo "Submitted tasks ${START}-${END}"
  START=$((END + 1))
  sleep 10
done

echo "=== Batch 2 complete at $(date) ==="
