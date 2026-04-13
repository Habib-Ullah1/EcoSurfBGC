#!/bin/bash
BASE="/data/habib/EcoSurfBGC"
ID_LIST="${BASE}/genome_lists/final_missing_ids.txt"
TEMPLATE="${BASE}/scripts/04_run_antismash.sh"
TOTAL=$(wc -l < "${ID_LIST}")
CHUNK=5
START=1

echo "Remaining genomes: ${TOTAL}"
echo "Starting at $(date)"

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

  # Create chunk using final_missing_ids.txt and fresh array indices 1-302
  sed "s|#SBATCH --array=.*|#SBATCH --array=${START}-${END}%5|" \
    "${TEMPLATE}" | \
  sed "s|bacteria_fasta_ids.txt|final_missing_ids.txt|" \
    > /tmp/antismash_remaining_chunk.sh

  ARRAY_LINE=$(grep "^#SBATCH --array" /tmp/antismash_remaining_chunk.sh)
  echo "$(date): Submitting ${START}-${END} | ${ARRAY_LINE}"

  sbatch /tmp/antismash_remaining_chunk.sh
  STATUS=$?

  if [ ${STATUS} -ne 0 ]; then
    echo "ERROR: Failed tasks ${START}-${END}. Retrying in 60s..."
    sleep 60
    continue
  fi

  echo "Submitted tasks ${START}-${END}"
  START=$((END + 1))
  sleep 10
done

echo "=== All remaining genomes submitted at $(date) ==="
