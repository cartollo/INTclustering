#!/bin/bash

# FOLDER="results/subsample_10_1000_16S_fam"

# for file in "$FOLDER"/*.rds; do
#   name=$(basename "$file" .rds)

#   echo "Running analysis for $name"

#   Rscript compDisMatrix.R -f "$file" -n 3

# done
#!/bin/bash

FOLDER="results/subsample_10_1000_shotgun_5clus_correlation_average"
MAX_JOBS=10

mkdir -p "$FOLDER/logs"

find "$FOLDER" -name "*.rds" -print0 | \
xargs -0 -n 1 -P "$MAX_JOBS" -I {} bash -c '
file="$1"
name=$(basename "$file" .rds)
echo "Running analysis for $name"
Rscript compDisMatrix.R -f "$file" -n 5 > "'"$FOLDER"'/logs/${name}.log" 2>&1
' _ {}

echo "done"