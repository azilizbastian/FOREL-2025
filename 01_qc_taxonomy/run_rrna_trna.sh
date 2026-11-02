#!/usr/bin/env bash
# Barrnap + tRNAscan-SE sur tous les bins d'un dossier (MIMAG rRNA/tRNA).
# Usage : bash run_rrna_trna.sh <dossier_bins> <dossier_sortie> [nb_jobs]
# Variables optionnelles :
#   EXT=fa            extension des bins (fa, fasta, fna...) — fichiers NON compressés
#   GTDB_DIR=/chemin  dossier de sortie GTDB-Tk : les génomes listés dans
#                     gtdbtk.ar53.summary.tsv sont lancés en mode archée
set -euo pipefail

BINS_DIR=${1:?dossier des bins manquant}
OUT=${2:?dossier de sortie manquant}
JOBS=${3:-16}
EXT=${EXT:-fa}
GTDB_DIR=${GTDB_DIR:-}

mkdir -p "$OUT"/{barrnap,trnascan,logs,done}

# Liste des archées (le reste est traité comme bactérie)
ARC_LIST="$OUT/archaea.txt"
: > "$ARC_LIST"
if [[ -n "$GTDB_DIR" && -f "$GTDB_DIR/gtdbtk.ar53.summary.tsv" ]]; then
  tail -n +2 "$GTDB_DIR/gtdbtk.ar53.summary.tsv" | cut -f1 > "$ARC_LIST"
fi
echo "Archées détectées : $(wc -l < "$ARC_LIST")"

export OUT ARC_LIST EXT

run_one() {
  local f=$1
  local name
  name=$(basename "$f" ".$EXT")
  [[ -e "$OUT/done/$name" ]] && return 0   # déjà fait -> reprise possible

  local kb=bac tflag=-B
  if grep -qxF "$name" "$ARC_LIST"; then kb=arc; tflag=-A; fi

  barrnap --kingdom "$kb" --threads 1 --quiet "$f" \
    > "$OUT/barrnap/$name.gff" 2> "$OUT/logs/$name.barrnap.log"

  rm -f "$OUT/trnascan/$name.tsv"
  tRNAscan-SE "$tflag" -q --thread 1 --brief \
    -o "$OUT/trnascan/$name.tsv" "$f" > "$OUT/logs/$name.trnascan.log" 2>&1

  touch "$OUT/done/$name"
}
export -f run_one

ls "$BINS_DIR"/*."$EXT" | parallel -j "$JOBS" --bar --joblog "$OUT/logs/parallel.joblog" run_one {}

echo "Terminé : $(ls "$OUT/done" | wc -l) bins traités."
echo "Échecs éventuels (colonne Exitval != 0) : $OUT/logs/parallel.joblog"
