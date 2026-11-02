#!/usr/bin/env bash
# Étape 11 : déréplication des MAGs procaryotes (HQ+MQ) au niveau espèce (SGB, 95 % ANI)
#            et répartition des SGB entre les deux batches (exclusifs / partagés)
# Lancer :  tmux new -d -s drep 'bash ~/11_drep_sgb.sh 2>&1 | tee ~/11.log'
set -euo pipefail

WD=/home/bastian/forel_func
OUT=$WD/drep
THREADS=${THREADS:-$(( $(nproc) / 2 ))}
EXCL=COBRA_ICEFIELD_1M-BIN_00001680      # MAG sans taxonomie GTDB-Tk, exclu des analyses
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate drep_env
echo "dRep $(pip show drep 2>/dev/null | awk '/^Version/{print $2}')"

mkdir -p "$OUT/genomes"
# 1. liste des 449 MAGs + qualité CheckM2 (évite de relancer CheckM dans dRep)
echo "genome,completeness,contamination" > "$OUT/genomeInfo.csv"
awk -F'\t' -v ex=$EXCL 'NR==FNR{mq[$1];next} FNR>1 && ($1 in mq) && $1!=ex {print $1".fa,"$2","$3}' \
  $WD/meta/mags_mq.txt $WD/checkm2/quality_report.tsv >> "$OUT/genomeInfo.csv"
tail -n +2 "$OUT/genomeInfo.csv" | cut -d, -f1 | while read -r g; do
  ln -sf "$(readlink -f $WD/mags/$g)" "$OUT/genomes/$g"
done
echo "MAGs à dérépliquer : $(ls $OUT/genomes | wc -l)"

# 2. dRep : clusters primaires à 90 % (Mash), secondaires à 95 % ANI (niveau espèce), couverture ≥ 30 %
#    -comp 50 -con 10 pour garder tous les MAGs HQ+MQ (dRep filtre à 75 % par défaut)
run_drep() {
  dRep dereplicate "$OUT/dRep_out" -g "$OUT"/genomes/*.fa \
    --genomeInfo "$OUT/genomeInfo.csv" -comp 50 -con 10 \
    -pa 0.90 -sa 0.95 -nc 0.30 --S_algorithm "$1" -p "$THREADS"
}
if [ ! -s "$OUT/dRep_out/data_tables/Cdb.csv" ]; then
  run_drep skani || { echo "skani indisponible dans dRep -> fastANI"; rm -rf "$OUT/dRep_out"; run_drep fastANI; }
fi

# 3. SGB exclusifs / partagés entre batches
C=$OUT/dRep_out/data_tables/Cdb.csv
awk -F, 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
 {g=$c["genome"]; sub(/\.fa$/,"",g); s=g; sub(/-BIN.*/,"",s); k=$c["secondary_cluster"]
  n[k]++; has[k","s]=1; samp[s]=1; print g"\t"s"\t"k > "'$OUT'/mag2sgb.tsv"}
 END{
  for(k in n){ a=has[k",ANACONDA_INSITU_7M"]; b=has[k",COBRA_ICEFIELD_1M"]
    if(a&&b) sh++; else if(a) ea++; else eb++; tot++; if(n[k]>1) multi++ }
  printf "\nMAGs : %d\nSGB (95 %% ANI) : %d (dont %d avec ≥2 MAGs)\n", NR-1, tot, multi
  printf "Exclusifs ANACONDA : %d\nExclusifs COBRA    : %d\nPartagés           : %d (%.1f %%)\n", ea, eb, sh, 100*sh/tot }' "$C" \
  | tee "$OUT/sgb_summary.txt"

# 4. SGB partagés : représentant (Wdb) + taxonomie GTDB
W=$OUT/dRep_out/data_tables/Wdb.csv
awk -F, 'NR>1{g=$1; sub(/\.fa$/,"",g); print $2"\t"g}' "$W" | sort > "$OUT/sgb_representatives.tsv"
echo; echo "Fichiers : $OUT/mag2sgb.tsv, $OUT/sgb_representatives.tsv, $OUT/sgb_summary.txt"
