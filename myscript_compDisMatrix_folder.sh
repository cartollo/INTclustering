#!/bin/bash

FOLDER="prova"

for file in "$FOLDER"/*.rds; do
  name=$(basename "$file" .rds)

  echo "Running analysis for $name"

  Rscript compDisMatrix.R -f "$file" -n 3

done