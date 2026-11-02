#!/usr/bin/env bash
# Étape 2 : prédiction des gènes (Prodigal, paramètres OMDB v2) + catalogues non redondants (MMseqs2)
# Lancer dans tmux :  bash 02_prodigal_mmseqs.sh 2>&1 | tee 02.log
set -euo pipefail

WD=/home/bastian/forel_func
THREADS=${THREADS:-$(( $(nproc) / 2 ))}
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate funcannot

mkdir -p "$WD"/prodigal/{faa,fna,gff} "$WD/catalog"

echo "== 1. Prodigal (-c -q -m -p single) sur $(ls "$WD/mags" | wc -l) MAGs"
run_prodigal() {
  f=$1; WD=$2
  b=$(basename "$f" .fa)
  [ -s "$WD/prodigal/faa/$b.faa" ] && return 0
  prodigal -i "$f" -c -q -m -p single \
    -a "$WD/prodigal/faa/$b.faa.tmp" -d "$WD/prodigal/fna/$b.fna.tmp" \
    -o "$WD/prodigal/gff/$b.gff" -f gff
  # préfixe le nom du MAG à chaque gène pour des identifiants uniques : MAG|contig_N
  for ext in faa fna; do
    sed "s/^>\([^ ]*\).*/>$b|\1/" "$WD/prodigal/$ext/$b.$ext.tmp" > "$WD/prodigal/$ext/$b.$ext"
    rm "$WD/prodigal/$ext/$b.$ext.tmp"
  done
}
export -f run_prodigal
ls "$WD"/mags/*.fa | xargs -P "$THREADS" -I{} bash -c 'run_prodigal "$@"' _ {} "$WD"

echo "== 2. Catalogue complet"
cat "$WD"/prodigal/faa/*.faa | sed 's/\*$//' > "$WD/catalog/all_proteins.faa"
cat "$WD"/prodigal/fna/*.fna > "$WD/catalog/all_genes.fna"
grep '^>' "$WD/catalog/all_proteins.faa" | sed 's/^>//' | awk -F'|' '{print $0"\t"$1}' > "$WD/catalog/gene2mag.tsv"
seqkit stats -T "$WD/catalog/all_proteins.faa" "$WD/catalog/all_genes.fna"

echo "== 3. Clustering MMseqs2"
cd "$WD/catalog"
TMP="$WD/catalog/tmp"; mkdir -p "$TMP"
# AA : NR100, NR50, NR30 (couverture 80 %, cov-mode 0)
for id in 1.0 0.5 0.3; do
  n=$(awk -v i=$id 'BEGIN{printf "%d", i*100}')
  [ -s "AA_NR${n}_rep_seq.fasta" ] && continue
  mmseqs easy-cluster all_proteins.faa "AA_NR${n}" "$TMP" \
    --min-seq-id "$id" -c 0.8 --cov-mode 0 --threads "$THREADS" -v 1
done
# NT : NR100, NR95
for id in 1.0 0.95; do
  n=$(awk -v i=$id 'BEGIN{printf "%d", i*100}')
  [ -s "NT_NR${n}_rep_seq.fasta" ] && continue
  mmseqs easy-linclust all_genes.fna "NT_NR${n}" "$TMP" \
    --min-seq-id "$id" -c 0.8 --cov-mode 0 --threads "$THREADS" -v 1
done
rm -rf "$TMP"

echo "== 4. Résumé"
printf "catalogue\tn_sequences\n" > catalog_summary.tsv
printf "AA_all\t%s\n" "$(grep -c '^>' all_proteins.faa)" >> catalog_summary.tsv
printf "NT_all\t%s\n" "$(grep -c '^>' all_genes.fna)" >> catalog_summary.tsv
for f in AA_NR100 AA_NR50 AA_NR30 NT_NR100 NT_NR95; do
  printf "%s\t%s\n" "$f" "$(grep -c '^>' ${f}_rep_seq.fasta)" >> catalog_summary.tsv
done
column -t catalog_summary.tsv
echo "== Terminé. Catalogues dans $WD/catalog"
