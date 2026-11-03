#!/usr/bin/env bash
# BiG-SLiCE v1.1.0 query against BiG-FAM (replaces step 5 of bgc_pipeline_449.sh).
# BiG-SLiCE 1.1.0 only accepts antiSMASH 5/6 GenBank files: on copies of the antiSMASH 8
# region files, the version field is set to 6.1.1 (annotations are not modified).
set -euo pipefail
eval "$(conda shell.bash hook)"; conda activate bigslice_v1
# === À MODIFIER ===
B=/mnt/gr-paoli/nas1/users/students/bastian/BGC_449   # output of bgc_pipeline_449.sh
DB=~/bgc_449_local/full_run_result                     # BiG-FAM full_run_result (local copy)
# ==================
I=$B/bigslice_input; rm -rf "$I"; mkdir -p "$I"
for f in "$B"/regions_gbk/*.gbk; do
  sed 's/^\( *Version *:: *\)8\.[0-9.]*/\16.1.1/' "$f" > "$I/$(basename "$f")"
done
bigslice --query "$I" --n_ranks 1 -t 32 "$DB" > "$B/logs/bigslice.log" 2>&1
R=$(ls -td "$DB"/reports/*/ | head -1)
mkdir -p "$B/bigslice_report"; cp -r "$R" "$B/bigslice_report/"
echo "$B/bigslice_report/$(basename "$R")" > "$B/logs/bigslice_report_dir.txt"
