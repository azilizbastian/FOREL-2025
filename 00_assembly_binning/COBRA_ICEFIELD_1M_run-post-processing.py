#!/usr/bin/env python

import re
import gzip

from pathlib import Path
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
 
import polars as pl
 
# ---- config: edit and run --------------------------------------------
BINS_DIR = "Medusa-filter-myloasm-v0.6.0-hifi-Semibin2/output_bins"
BINS_INFO = "Medusa-filter-myloasm-v0.6.0-hifi-Semibin2/bins_info.tsv"
WHOKARYOTE = "Medusa-filter-myloasm-v0.6.0-hifi-whokaryote/whokaryote_predictions_T.tsv"
TIARA = "Medusa-filter-myloasm-v0.6.0-hifi-whokaryote/tiara_pred.txt"
SAMPLE_NAME = "COBRA_ICEFIELD_1M"
OUTDIR = "COBRA_ICEFIELD_1M-processed-bins"
MIN_BIN_SIZE = 500_000  # bp; bins smaller than this are dropped entirely
# ------------------------------------------------------------------------

HEADER_RE = re.compile(r'_len-(\d+)_circular-([^_]+)_depth-([\d-]+)_duplicated-(.+)$')
NUM_RE = re.compile(r'(\d+)$')
 
def clean_tiara(v):
    """Strip a leading 'mult=X.XX ' prefix from a Tiara stage-1 call, if present."""
    if not v:
        return v
    m = re.match(r'^mult=[\d.]+\s+(.*)$', v.strip())
    return m.group(1) if m else v.strip()
 
 
outdir = Path(OUTDIR)
bins_outdir = outdir / "renamed_bins"
bins_outdir.mkdir(parents=True, exist_ok=True)
 
# filter BINS_INFO to kept bins first; the glob below only processes members of this
bin_info_of = {
    Path(r["filename"]).name.removesuffix(".fa.gz"): r
    for r in pl.read_csv(BINS_INFO, separator="\t").filter(pl.col("nbps") >= MIN_BIN_SIZE).iter_rows(named=True)
}
 

# 1) one pass per bin: read, sort by length, rename, write, log scaffold rows
scaffold_rows = []
new_bin_of_orig = {}
for bf in sorted(Path(BINS_DIR).glob("*.fa.gz")):
    orig_bin = bf.name.removesuffix(".fa.gz")
    if orig_bin not in bin_info_of:
        continue
    new_bin = f"{SAMPLE_NAME}-BIN_{int(NUM_RE.search(orig_bin).group(1)):08d}"
    new_bin_of_orig[orig_bin] = new_bin
    print(f"Processing bin {orig_bin} into {new_bin}...")
 
    with gzip.open(bf, "rt") as fh:
        records = []
        for title, seq in SimpleFastaParser(fh):
            cid = title.split()[0]
            length, circular, depth, duplicated = HEADER_RE.search(cid).groups()
            records.append(dict(id=cid, seq=seq, length=len(seq),
                                 circular=circular, depth=depth, duplicated=duplicated))
    records.sort(key=lambda r: r["length"], reverse=True)
 
    with gzip.open(bins_outdir / f"{new_bin}.fa.gz", "wt", compresslevel=6) as out:
        for i, r in enumerate(records, 1):
            new_name = f"{new_bin}-scaffold_{i}"
            out.write(f">{new_name}\n{r['seq']}\n")
            scaffold_rows.append(dict(
                original_name=r["id"], new_name=new_name, size=r["length"],
                circularity=r["circular"], coverage=r["depth"], duplication=r["duplicated"],
                original_bin=orig_bin, new_bin=new_bin,
            ))
 
# 2) append Whokaryote + Tiara annotations to the scaffold rows
wk_df = pl.read_csv(WHOKARYOTE, separator="\t")
whokaryote_of = dict(zip(wk_df["contig"].to_list(), wk_df["predicted"].to_list()))
 
tiara_df = pl.read_csv(TIARA, separator="\t")
tiara_of = {r["sequence_id"]: clean_tiara(r["class_fst_stage"]) for r in tiara_df.iter_rows(named=True)}
 
for row in scaffold_rows:
    row["whokaryote"] = whokaryote_of.get(row["original_name"])
    row["tiara"] = tiara_of.get(row["original_name"])
 
pl.DataFrame(scaffold_rows).sort(["new_bin", "original_name"]).write_csv(
    outdir / "scaffolds_summary.tsv", separator="\t")
 
# 3) bin summary table: SemiBin2's own stats + prokaryote/eukaryote counts
wk_counts = {}
for row in scaffold_rows:
    wk_counts.setdefault(row["new_bin"], Counter())[row["whokaryote"]] += 1
 
bin_rows = [
    dict(original_bin_name=orig_bin, new_bin_name=new_bin,
         n_scaffolds=bin_info_of[orig_bin]["nr_contigs"], size=bin_info_of[orig_bin]["nbps"],
         N50=bin_info_of[orig_bin]["N50"], L50=bin_info_of[orig_bin]["L50"],
         n_prokaryote=wk_counts.get(new_bin, Counter()).get("prokaryote", 0),
         n_eukaryote=wk_counts.get(new_bin, Counter()).get("eukaryote", 0))
    for orig_bin, new_bin in new_bin_of_orig.items()
]
pl.DataFrame(bin_rows).sort("new_bin_name").write_csv(outdir / "bins_summary.tsv", separator="\t")
 
print(f"{len(bin_rows)} bins, {len(scaffold_rows)} binned scaffolds -> {outdir}/")
