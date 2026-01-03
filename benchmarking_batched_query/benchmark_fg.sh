#!/usr/bin/env bash

INPUT_FILE="functional_groups.txt"
SCRIPT="single_cohort_API.py"
RUNS=100
OUT="benchmark_results.csv"

echo "protein_id,implementation,run,time_seconds" > "$OUT"

tail -n +2 "$INPUT_FILE" | while read -r protein_id fg_query; do
  for impl in new old; do
    for run in $(seq 1 $RUNS); do

      if [[ "$impl" == "old" ]]; then
        FLAG="--use-old"
      else
        FLAG=""
      fi

      start=$(date +%s.%N)

      python "$SCRIPT" \
        $FLAG \
        --queryFG "$fg_query" \
        > /dev/null

      end=$(date +%s.%N)
      elapsed=$(awk "BEGIN {print $end - $start}")

      echo "$protein_id,$impl,$run,$elapsed" >> "$OUT"
    done
  done
done