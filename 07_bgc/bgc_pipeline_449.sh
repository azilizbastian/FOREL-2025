#!/usr/bin/env bash
# =============================================================================
# BGC analysis on the 449 prokaryotic MAGs (FOREL 2025, POC manuscript)
# antiSMASH 8 (Prodigal) -> BiG-SCAPE 2 (cutoff 0.3) -> BiG-SLiCE query vs BiG-FAM (d <= 900)
# Run on grpaoli-haas011 inside tmux/screen:  bash bgc_pipeline_449.sh 2>&1 | tee bgc_449.log
# =============================================================================
set -euo pipefail

# ----------------------------- À ADAPTER -------------------------------------
BASE=/mnt/gr-paoli/nas1/users/students/bastian
MAG_DIR=/mnt/gr-paoli/nas1/users/students/bastian/BGC_449/mags_449
EXT=fa                                  # extension des MAGs
OUT=$BASE/BGC_449
PFAM=/home/bastian/ENTER/envs/antismash/lib/python3.11/site-packages/antismash/databases/pfam/35.0/Pfam-A.hmm
BIGFAM_DB=/mnt/gr-paoli/nas1/users/students/bastian/bigfam/bigfam_db/full_run_result
CPU_PER_JOB=4
N_JOBS=8                                # 8 x 4 = 32 coeurs
BIGSCAPE_CPU=32
ENV_ANTISMASH=antismash
ENV_BIGSCAPE=bigscape
ENV_BIGSLICE=bigslice_v1
# -----------------------------------------------------------------------------

mkdir -p "$OUT"/{renamed_mags,antismash,regions_gbk,bigscape,logs}
eval "$(conda shell.bash hook)"

echo "[$(date)] Versions :"
conda activate "$ENV_ANTISMASH"; antismash --version | tee "$OUT/logs/versions.txt"
conda activate "$ENV_BIGSCAPE";  bigscape --version 2>&1 | tee -a "$OUT/logs/versions.txt" || true
conda activate "$ENV_BIGSLICE";  bigslice --version 2>&1 | tee -a "$OUT/logs/versions.txt" || true

# --- 1. Préfixer les contigs par l'ID du MAG --------------------------------
# Évite les collisions de noms de contigs entre MAGs (BiG-SCAPE/BiG-SLiCE
# identifient les BGCs par nom de fichier .gbk = <contig>.regionXXX.gbk)
echo "[$(date)] Renommage des contigs"
for f in "$MAG_DIR"/*."$EXT"; do
  mag=$(basename "$f" ."$EXT")
  awk -v m="$mag" '/^>/{sub(/^>/,""); split($0,a," "); print ">" m "__" a[1]; next}{print}' \
      "$f" > "$OUT/renamed_mags/$mag.fna"
done
n=$(ls "$OUT"/renamed_mags/*.fna | wc -l)
echo "MAGs à traiter : $n"   # doit afficher 449

# --- 2. antiSMASH 8, gene finding Prodigal ----------------------------------
echo "[$(date)] antiSMASH"
conda activate "$ENV_ANTISMASH"
run_antismash() {
  f=$1; mag=$(basename "$f" .fna); o=$OUT/antismash/$mag
  [[ -f "$o/$mag.json" ]] && return 0          # reprise possible si interruption
  antismash "$f" \
    --taxon bacteria \
    --genefinding-tool prodigal \
    --cpus "$CPU_PER_JOB" \
    --output-dir "$o" \
    --output-basename "$mag" \
    --allow-long-headers \
    > "$OUT/logs/antismash_$mag.log" 2>&1 || echo "ECHEC $mag" >> "$OUT/logs/antismash_failed.txt"
}
export -f run_antismash; export OUT CPU_PER_JOB
ls "$OUT"/renamed_mags/*.fna | xargs -P "$N_JOBS" -I{} bash -c 'run_antismash {}'

# --- 3. Collecter les régions (1 .gbk = 1 BGC) -------------------------------
echo "[$(date)] Collecte des régions"
find "$OUT/antismash" -name "*.region*.gbk" -exec cp {} "$OUT/regions_gbk/" \;
echo "BGCs détectés : $(ls "$OUT"/regions_gbk | wc -l)"

# --- 4. BiG-SCAPE 2, cutoff 0.3 ----------------------------------------------
# Vérifier les options avec `bigscape cluster --help` (notamment --classify :
# 'legacy' = classes NRPS/PKSI/PKSother/RiPPs/Terpene/Others comme BiG-SCAPE 1)
echo "[$(date)] BiG-SCAPE"
conda activate "$ENV_BIGSCAPE"
bigscape cluster \
  -i "$OUT/antismash" \
  -o "$OUT/bigscape" \
  -p "$PFAM" \
  --gcf-cutoffs 0.3 \
  --classify legacy \
  --include-singletons \
  -c "$BIGSCAPE_CPU" \
  > "$OUT/logs/bigscape.log" 2>&1

# --- 5. BiG-SLiCE : query contre BiG-FAM -------------------------------------
# Le mode query écrit ses résultats dans $BIGFAM_DB/reports/<run_id>/data.db
echo "[$(date)] BiG-SLiCE query"
conda activate "$ENV_BIGSLICE"
bigslice --query "$OUT/antismash" --n_ranks 1 -t "$BIGSCAPE_CPU" "$BIGFAM_DB" \
  > "$OUT/logs/bigslice.log" 2>&1
ls -td "$BIGFAM_DB"/reports/*/ | head -1 | tee "$OUT/logs/bigslice_report_dir.txt"

echo "[$(date)] Terminé. Lancer ensuite : python summarise_bgc.py"
