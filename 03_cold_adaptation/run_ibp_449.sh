#!/usr/bin/env bash
# Ice-binding proteins (Pfam PF11999) in the 449 MAGs vs OMDB v2 (polar / non-polar marine genomes).
# recA (KOfam K03553) is used as a detection control.
# Run from the working directory containing db/ (HMMs from 02_prepare_hmms.py) and results/.
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)                 # folder with the scripts
# === À MODIFIER ===
F=~/forel_func                                   # output of 01_qc_taxonomy / 02_functional_novelty
O=/mnt/gr-paoli/nas1/data/databases/OMDB/2.0     # OMDB v2.0
M=/mnt/gr-paoli/nas1/data/databases/OMDB/2.1/metadata
CPUS=32
# ==================
mkdir -p results449/omdb results449/figures results/omdb

# 0. Scan of the OMDB v2.0 protein catalogue (508,832,278 proteins; ~1 h, done once)
[ -s results/omdb/omdb_hits.tsv ] || pigz -dc $O/catalog/clustering/OMDBv2.0_AA_G_R/OMDBv2.0_AA_G_R.faa.gz | \
  python $S/10_hmm_scan_proteins.py --faa - --genome-from rep --pfam db/selected_pfam.hmm \
    --recA db/K03553.hmm --recA-threshold 188.03 --recA-score-type domain \
    --out results/omdb/omdb_hits.tsv --cpus $CPUS

# 1. Table of the 449 MAGs (CheckM2 quality + GTDB-Tk R232 taxonomy)
python - <<PY
import pandas as pd, os
F = os.path.expanduser("$F")
q = pd.read_csv(f"{F}/checkm2/quality_report.tsv", sep="\t")[["Name","Completeness","Contamination"]]
q.columns = ["genome_id","completeness","contamination"]
g = pd.read_csv(f"{F}/meta/gtdbtk_summary.tsv", sep="\t")[["user_genome","classification"]]
r = {}
for _, row in g.iterrows():
    d = dict(x.split("__", 1) for x in str(row.classification).split(";") if "__" in x)
    r[row.user_genome] = dict(phylum=d.get("p",""), family=d.get("f",""), species=d.get("s",""))
t = pd.DataFrame.from_dict(r, orient="index").reset_index().rename(columns={"index":"genome_id"})
m = q.merge(t, on="genome_id", how="inner")
m = m[m.genome_id != "COBRA_ICEFIELD_1M-BIN_00001680"]   # not classified beyond domain
m.to_csv("results449/genomes_mags.tsv", sep="\t", index=False)
print(len(m), "MAGs")
PY

# 2. Scan of the MAG proteins (Prodigal, headers MAG|protein)
python $S/10_hmm_scan_proteins.py --faa $F/catalog/all_proteins.faa --genome-from pipe \
  --pfam db/selected_pfam.hmm --recA db/K03553.hmm --recA-threshold 188.03 --recA-score-type domain \
  --out results449/omdb/mags_hits.tsv --cpus $CPUS

# 3. Enrichment (main: >=90% completeness; sensitivity: >=70%) and figures
for C in 90 70; do
  python $S/11_ibp_enrichment.py --omdb-hits results/omdb/omdb_hits.tsv --mag-hits results449/omdb/mags_hits.tsv \
    --omdb-genomes $M/OMDv2.0_genome-metadata-20250702.tsv.gz --omdb-samples $M/omdb21.samples.tsv.gz \
    --mags results449/genomes_mags.tsv --min-comp $C --max-cont 5 \
    --outdir results449/omdb/comp$C 2>&1 | tee results449/omdb/comp$C.log
  python $S/12_fig_ibp_enrichment.py --genomes results449/omdb/comp$C/ibp_enrichment_genomes.tsv \
    --level genus --out results449/figures/fig_ibp_enrichment_comp${C}_genus
done
