#!/usr/bin/env bash
# Étape 5 : nouveauté des familles protéiques FOREL (AA NR30) par rapport au catalogue OMDB v2 (AA NR30)
#   tmux new -d -s novelty 'bash ~/05_novelty_omdb.sh 2>&1 | tee ~/05.log'
set -euo pipefail

WD=/home/bastian/forel_func
OMDB=/mnt/gr-paoli/nas1/data/databases/OMDB/2.0/catalog/clustering/OMDBv2.0_AA_G_NR30/OMDBv2.0_AA_G_NR30.faa.gz
DBDIR=/home/bastian/databases/omdb_nr30
OUT=$WD/novelty
THREADS=${THREADS:-$(( $(nproc) / 2 ))}
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate funcannot
mkdir -p "$DBDIR" "$OUT/tmp"

echo "== 1. Base MMseqs2 du catalogue OMDB v2 NR30 (une seule fois)"
[ -s "$DBDIR/omdb_nr30.dbtype" ] || mmseqs createdb "$OMDB" "$DBDIR/omdb_nr30" -v 1
echo "Séquences OMDB NR30 : $(wc -l < "$DBDIR/omdb_nr30.index")"

echo "== 2. Recherche des représentants FOREL NR30 contre OMDB NR30"
cd "$OUT"
[ -s forel_nr30.dbtype ] || mmseqs createdb "$WD/catalog/AA_NR30_rep_seq.fasta" forel_nr30 -v 1
if [ ! -s forel_vs_omdb.m8 ]; then
  mmseqs search forel_nr30 "$DBDIR/omdb_nr30" res tmp \
    -s 7.5 -e 1e-3 -c 0.8 --cov-mode 0 --max-seqs 10 --threads "$THREADS" -v 1
  mmseqs filterdb res res_best --extract-lines 1 -v 1   # meilleur hit par requête
  mmseqs convertalis forel_nr30 "$DBDIR/omdb_nr30" res_best forel_vs_omdb.m8 \
    --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits" -v 1
fi
rm -rf tmp

echo "== 3. Catégories par famille (représentant NR30)"
# novel = aucun hit (couverture ≥ 80 %, e ≤ 1e-3) ; sinon classe d'identité du meilleur hit
grep '^>' "$WD/catalog/AA_NR30_rep_seq.fasta" | sed 's/^>//; s/ .*//' > fam_ids.txt
awk -F'\t' 'NR==FNR{best[$1]=$3; next}
  { id=best[$1];
    if (id=="") c="novel_vs_OMDB";
    else if (id<0.5) c="distant_30-50";
    else if (id<0.9) c="related_50-90";
    else c="near_identical_90+";
    print $1"\t"(id==""?"NA":id)"\t"c }' forel_vs_omdb.m8 fam_ids.txt > family_novelty.tsv

echo "== 4. Propagation aux gènes (via les clusters NR30) et par MAG"
# AA_NR30_cluster.tsv : représentant <tab> membre
awk -F'\t' 'NR==FNR{cat[$1]=$3; next}{print $2"\t"$1"\t"cat[$1]}' \
  family_novelty.tsv "$WD/catalog/AA_NR30_cluster.tsv" > gene_novelty.tsv
awk -F'\t' '{split($1,a,"|"); n[a[1]]++; if($3=="novel_vs_OMDB") k[a[1]]++}
  END{print "mag\tn_genes\tn_genes_novel_fam\tfrac_novel";
      for(m in n) printf "%s\t%d\t%d\t%.4f\n", m, n[m], k[m]+0, (k[m]+0)/n[m]}' \
  gene_novelty.tsv > mag_novelty_fraction.tsv

echo "-- Familles (NR30) :"
cut -f3 family_novelty.tsv | sort | uniq -c
echo "-- Gènes :"
cut -f3 gene_novelty.tsv | sort | uniq -c
echo "== Terminé. Résultats dans $OUT"
