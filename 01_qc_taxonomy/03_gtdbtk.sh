#!/usr/bin/env bash
# Étape 3 : classification taxonomique des 450 MAGs (GTDB-Tk 2.8.0 + release R232)
# Lancer dans tmux :  bash 03_gtdbtk.sh 2>&1 | tee 03.log
set -euo pipefail

WD=/home/bastian/forel_func
REL=/mnt/gr-paoli/nas1/data/databases/GTDB/R232
ENV=gtdbtk_new
THREADS=${THREADS:-$(( $(nproc) / 2 ))}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"

# trouve le dossier de la release (celui qui contient taxonomy/)
if [ -d "$REL/taxonomy" ]; then DATA="$REL"
else DATA=$(dirname "$(find "$REL" -maxdepth 3 -type d -name taxonomy | head -1)")
fi
export GTDBTK_DATA_PATH="$DATA"
echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
gtdbtk --version
# gtdbtk check_install   # déjà validé : R232 OK avec GTDB-Tk 2.8.0

free -g | head -2      # pplacer demande beaucoup de RAM

gtdbtk classify_wf \
  --genome_dir "$WD/mags" -x fa \
  --out_dir "$WD/gtdbtk" \
  --cpus "$THREADS" --pplacer_cpus 1

# tableau unique bactéries + archées
cd "$WD/gtdbtk"
{ head -1 gtdbtk.bac120.summary.tsv
  tail -n +2 gtdbtk.bac120.summary.tsv
  [ -s gtdbtk.ar53.summary.tsv ] && tail -n +2 gtdbtk.ar53.summary.tsv
} > "$WD/meta/gtdbtk_summary.tsv"

echo "== Nouveauté taxonomique (rang le plus bas non assigné)"
awk -F'\t' 'NR>1{
  split($2,t,";"); lvl="species";
  if (t[7]=="s__") lvl="novel_species";
  if (t[6]=="g__") lvl="novel_genus";
  if (t[5]=="f__") lvl="novel_family";
  if (t[4]=="o__") lvl="novel_order";
  print $1"\t"lvl}' "$WD/meta/gtdbtk_summary.tsv" > "$WD/meta/mag_novelty.tsv"
cut -f2 "$WD/meta/mag_novelty.tsv" | sort | uniq -c
