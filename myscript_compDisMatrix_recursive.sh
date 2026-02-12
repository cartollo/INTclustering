#!/bin/bash

BASE_DIR=$(pwd)

echo "🚀 Avvio script"
echo "Base directory: $BASE_DIR"
echo "==================================="

for dir in results/*/; do

  # Nome pulito della sottocartella
  SUBDIR_NAME=$(basename "$dir")
  # - contengono NOWHITE
  # - NON contengono ispca
  if [[ "$SUBDIR_NAME" != *NOWHITE* || "$SUBDIR_NAME" == *ispca* ]]; then
    echo "⏭  Salto $SUBDIR_NAME (non valida)"
    echo "-----------------------------------"
    continue
  fi
  
  echo "➡️  Processando $SUBDIR_NAME"

  # File con path completo relativo alla directory iniziale
  FILE_ISCORE=$(ls "$dir"*iscore*.rds 2>/dev/null)
  FILE_RDS=$(ls "$dir"*.rds 2>/dev/null | grep -v iscore)

  if [[ -z "$FILE_ISCORE" || -z "$FILE_RDS" ]]; then
    echo "⚠️  File .rds non trovati in $SUBDIR_NAME, salto."
    echo "-----------------------------------"
    continue
  fi

  for n in {3..10}; do
    echo "   ▶️  Lancio n=$n"

    Rscript "$BASE_DIR/compDisMatrix.R" \
      -f "$FILE_RDS" \
      -c "$FILE_ISCORE" \
      -n "$n"
  done

  echo "✅ Completata $SUBDIR_NAME"
  echo "-----------------------------------"

done

echo "🎉 Script terminato"