#!/usr/bin/env bash
# Étape 1 : sélection des bins procaryotes + CheckM2 + filtre qualité (MIMAG ≥ medium)
# Lancer dans tmux :  bash 01_prep_checkm2.sh 2>&1 | tee 01.log
set -euo pipefail

P=/mnt/gr-paoli/nas1/data/seq/internal/FOREL/2025
WD=/home/bastian/forel_func
CHECKM2_DB=/mnt/gr-paoli/nas1/data/databases/CHECKM2/3.0/CheckM2_database/uniref100.KO.1.dmnd
THREADS=${THREADS:-$(( $(nproc) / 2 ))}
MIN_COMP=50
MAX_CONT=10

source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir -p "$WD"/{bins,meta}

echo "== 1. Tableau combiné + exclusion des bins majoritairement eucaryotes"
ALL="$WD/meta/bins_summary_all.tsv"
printf "sample\t%s\n" "$(head -1 $P/ANACONDA_INSITU_7M/ANACONDA_INSITU_7M-processed-bins/bins_summary.tsv)" > "$ALL"
for s in ANACONDA_INSITU_7M COBRA_ICEFIELD_1M; do
  awk -F'\t' -v s="$s" 'NR>1{print s"\t"$0}' "$P/$s/$s-processed-bins/bins_summary.tsv" >> "$ALL"
done
# colonnes : 1 sample, 3 new_bin_name, 8 n_prokaryote, 9 n_eukaryote
awk -F'\t' 'NR>1 && $9 > $8 {print $3}' "$ALL" > "$WD/meta/excluded_euk.txt"
awk -F'\t' 'NR>1 && $9 <= $8 {print $1"\t"$3}' "$ALL" > "$WD/meta/prok_bins.tsv"
echo "Total : $(($(wc -l < "$ALL")-1)) | exclus (euk) : $(wc -l < "$WD/meta/excluded_euk.txt") | gardés : $(wc -l < "$WD/meta/prok_bins.tsv")"

echo "== 2. Décompression des bins gardés"
while IFS=$'\t' read -r s b; do
  out="$WD/bins/$b.fa"
  [ -s "$out" ] || zcat "$P/$s/$s-processed-bins/renamed_bins/$b.fa.gz" > "$out"
done < "$WD/meta/prok_bins.tsv"
echo "$(ls "$WD/bins" | wc -l) fichiers dans $WD/bins"

echo "== 3. CheckM2 ($THREADS threads)"
conda activate checkm2
checkm2 --version
if [ ! -s "$WD/checkm2/quality_report.tsv" ]; then
  checkm2 predict -i "$WD/bins" -x fa -o "$WD/checkm2" -t "$THREADS" \
    --database_path "$CHECKM2_DB" --force
fi

echo "== 4. Filtre qualité (complétude ≥ $MIN_COMP %, contamination < $MAX_CONT %)"
Q="$WD/checkm2/quality_report.tsv"
awk -F'\t' -v c=$MIN_COMP -v k=$MAX_CONT 'NR>1 && $2>=c && $3<k {print $1}' "$Q" > "$WD/meta/mags_mq.txt"
awk -F'\t' 'NR>1 && $2>=90 && $3<5 {print $1}' "$Q" > "$WD/meta/mags_hq.txt"
mkdir -p "$WD/mags"
while read -r b; do ln -sf "$WD/bins/$b.fa" "$WD/mags/$b.fa"; done < "$WD/meta/mags_mq.txt"
echo "MAGs ≥ medium : $(wc -l < "$WD/meta/mags_mq.txt") | dont high quality (≥90/<5) : $(wc -l < "$WD/meta/mags_hq.txt")"
echo "== Terminé. MAGs retenus (liens) dans $WD/mags"
