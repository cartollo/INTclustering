#!/bin/bash

# FOLDER="results/subsample_10_1000_16S_fam"

# for file in "$FOLDER"/*.rds; do
#   name=$(basename "$file" .rds)

#   echo "Running analysis for $name"

#   Rscript compDisMatrix.R -f "$file" -n 3

# done
#!/bin/bash

FOLDER="results/subsample_10_1000_manhattan_wardd2_3clus_shotgun"
MAX_JOBS=3

mkdir -p "$FOLDER/logs"

for file in "$FOLDER"/*.rds; do
  [ -e "$file" ] || continue

  name=$(basename "$file" .rds)

  echo "Running analysis for $name"

  Rscript compDisMatrix.R \
    -f "$file" \
    -n 3 \
    > "$FOLDER/logs/${name}.log" 2>&1 &

  while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
    sleep 1
  done
done

wait
echo "done"