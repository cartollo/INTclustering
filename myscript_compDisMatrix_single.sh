#!/bin/bash

FOLDER="$1"

if [[ -z "$FOLDER" ]]; then
  echo "Uso: ./myscript.sh <cartella>"
  exit 1
fi

# Trova file
FILE_ISCORE=$(ls "$FOLDER"/*_iscore_*.rds 2>/dev/null)
FILE_RDS=$(ls "$FOLDER"/*.rds 2>/dev/null | grep -v "_iscore_")

# Controllo sicurezza
if [[ -z "$FILE_ISCORE" && -z "$FILE_RDS" ]]; then
  echo "Errore: file .rds non trovati correttamente nella cartella $FOLDER"
  exit 1
fi

echo "primo file : $FILE_RDS"
echo "seconod file: $FILE_ISCORE"
echo "===================================="

for k in {2..5}; do
  echo "Lancio k=$k"

  if [ -f "$FILE_RDS" ]; then
    Rscript compDisMatrix.R \
      -s "$FILE_ISCORE" \
      -f "$FILE_RDS" \
      -n "$k"
  else
    Rscript compDisMatrix.R \
      -s "$FILE_ISCORE" \
      -n "$k"
  fi
done

echo "Script completato."