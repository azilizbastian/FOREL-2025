#!/usr/bin/env bash
# Vérification des chiffres du paragraphe "Metagenome assembly" (procaryotes + eucaryotes)
# Usage : bash ~/08_verif_counts.sh 2>&1 | tee ~/08_verif.txt
F=/mnt/gr-paoli/nas1/data/seq/internal/FOREL/2025
WD=~/forel_func
HQ=/mnt/gr-paoli/nas1/users/students/bastian/FOREL-2025-INPUT/mimag/HQ_MIMAG_103.tsv
SAMPLES="ANACONDA_INSITU_7M COBRA_ICEFIELD_1M"
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=== 1. Bins SemiBin2 (renamed_bins) et bins à majorité eucaryote (Whokaryote)"
for s in $SAMPLES; do
  B=$F/$s/$s-processed-bins
  n_fa=$(ls $B/renamed_bins | grep -c '\.fa')
  awk -F'\t' -v s=$s -v nfa=$n_fa 'NR>1{t++; if($8>$7) e++} END{printf "%-20s fichiers=%d  lignes_summary=%d  euk(n_euk>n_prok)=%d  prok=%d\n", s, nfa, t, e, t-e}' $B/bins_summary.tsv
done
echo "(colonnes : $(head -1 $F/ANACONDA_INSITU_7M/ANACONDA_INSITU_7M-processed-bins/bins_summary.tsv | tr '\t' ' '))"

echo; echo "=== 2. Procaryotes : qualité (CheckM2 + liste MIMAG HQ) par échantillon"
awk -F'\t' 'NR==FNR{if(FNR>1){g=$1; sub(/\.(fa|fna|fasta)(\.gz)?$/,"",g); hq[g]}; next}
 FNR==1{next}
 {s=$1; sub(/-BIN.*/,"",s)
  if($1 in hq) q="HQ"; else if($2>=50 && $3<10) q="MQ"; else if($2>=50) q="contam>=10"; else q="LQ<50"
  c[s"\t"q]++; tot[s]++}
 END{for(k in c) print k"\t"c[k]; for(s in tot) print s"\tTOTAL\t"tot[s]}' "$HQ" $WD/checkm2/quality_report.tsv | sort
echo "HQ dans la liste : $(tail -n +2 "$HQ" | wc -l)"
echo "MAG(s) exclus (pas de taxonomie au-delà du domaine) :"
awk -F'\t' 'NR>1 && $2 !~ /^d__[^;]*;p__[^;]+/ {print "  "$1"\t"$2}' $WD/meta/gtdbtk_summary.tsv

echo; echo "=== 3. Circularité des HQ (génomes complets)"
awk -F'\t' 'NR==FNR{if(FNR>1){g=$1; sub(/\.(fa|fna|fasta)(\.gz)?$/,"",g); hq[g]}; next}
 ($1 in hq){s=$1; sub(/-BIN.*/,"",s); c[s"\t"$2]++} END{for(k in c) print k"\t"c[k]}' "$HQ" $WD/meta/bin_circularity.tsv | sort

echo; echo "=== 4. Bins REMAG (eucaryotes)"
for s in $SAMPLES; do
  R=$F/$s/Medusa-filter-myloasm-v0.6.0-hifi-REMAG
  echo "$s : $(ls $R/bins 2>/dev/null | grep -cE '\.(fa|fasta|fna)(\.gz)?$') fichiers FASTA dans bins/ ; autres fichiers :"
  ls $R | sed 's/^/    /'
done

echo; echo "=== 5. Sorties EukCC trouvées"
find $F /mnt/gr-paoli/nas1/users/students/bastian ~ -maxdepth 6 -name "eukcc.csv" 2>/dev/null | while read f; do
  echo "--- $f ($(($(wc -l < "$f")-1)) bins)"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="completeness")C=i; if($i=="contamination")K=i}; next}
   {q=($C>=90&&$K<5)?"HQ":(($C>=50&&$K<10)?"MQ":"LQ"); n[q]++} END{printf "    HQ=%d MQ=%d LQ=%d\n", n["HQ"], n["MQ"], n["LQ"]}' "$f"
done

echo; echo "=== 6. Versions"
for e in "remag:remag --version" "semibin:SemiBin2 --version" "eukcc:eukcc --version" "eukcc:metaeuk version" \
         "checkm2:checkm2 --version" "gtdbtk_new:gtdbtk --version"; do
  env=${e%%:*}; cmd=${e#*:}
  printf "%-28s %s\n" "$cmd" "$(conda run -n $env $cmd 2>&1 | grep -m1 -iE '[0-9]+\.[0-9]+' || echo 'env/commande introuvable')"
done
for t in whokaryote.py tiara myloasm; do printf "%-28s %s\n" "$t" "$( (command -v $t >/dev/null && $t --version 2>&1 | head -1) || echo 'pas dans base - voir logs Medusa')"; done
grep -m1 -i version $F/ANACONDA_INSITU_7M/Medusa-filter-myloasm-v0.6.0-hifi/myloasm_*.log 2>/dev/null
