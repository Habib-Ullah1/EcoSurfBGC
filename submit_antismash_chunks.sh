#!/bin/bash
# Submits antiSMASH jobs in chunks of 5, respecting 10-job QOS limit

BASE="/data/habib/EcoSurfBGC"
ID_LIST="${BASE}/genome_lists/bacteria_fasta_ids.txt"
TEMPLATE="${BASE}/scripts/04_run_antismash.sh"
TOTAL=$(wc -l < "${ID_LIST}")
CHUNK=5
START=6  # Tasks 1-5 already submitted

echo "Total genomes: ${TOTAL}"
echo "Submitting chunks of ${CHUNK} starting from task ${START}"

while [ ${START} -le ${TOTAL} ]; do
  END=$((START + CHUNK - 1))
  if [ ${END} -gt ${TOTAL} ]; then
    END=${TOTAL}
  fi

  # Wait until 5 or fewer jobs in queue (leaves room for our new chunk)
  while true; do
    CURRENT_JOBS=$(squeue -u habib -h | wc -l)
    if [ ${CURRENT_JOBS} -le 4 ]; then
      break
    fi
    echo "$(date): ${CURRENT_JOBS} jobs in queue, waiting 60s..."
    sleep 60
  done

  # Create chunk script by replacing array line
  sed "s|#SBATCH --array=.*|#SBATCH --array=${START}-${END}%5|" \
    "${TEMPLATE}" > /tmp/antismash_chunk.sh

  # Verify the array line was set correctly
  ARRAY_LINE=$(grep "^#SBATCH --array" /tmp/antismash_chunk.sh)
  echo "$(date): Submitting tasks ${START}-${END} | ${ARRAY_LINE}"

  sbatch /tmp/antismash_chunk.sh
  SUBMIT_STATUS=$?

  if [ ${SUBMIT_STATUS} -ne 0 ]; then
    echo "ERROR: Submission failed for tasks ${START}-${END}. Waiting 120s and retrying..."
    sleep 120
    continue  # Retry same chunk
  fi

  echo "Successfully submitted tasks ${START}-${END}"
  START=$((END + 1))
  sleep 10
done

echo "=== All ${TOTAL} genomes submitted at $(date) ==="
