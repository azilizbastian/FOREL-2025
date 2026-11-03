#!/usr/bin/env python3
"""Numbers reported in the BGC Results: BGCs per MAG and per phylum, completeness,
distribution of BiG-SCAPE classes across MAGs and phyla."""
import pandas as pd
# === À MODIFIER ===
BGC = "/mnt/gr-paoli/nas1/users/students/bastian/BGC_449/bgc_table.tsv"
TAXO = "/home/bastian/forel_func/phylo/mag_taxonomy.tsv"
# ==================
b = pd.read_csv(BGC, sep="\t")
t = pd.read_csv(TAXO, sep="\t")[["genome", "phylum"]]

n = b.groupby("mag").size().rename("n_bgc")
m = t.merge(n, left_on="genome", right_index=True, how="left").fillna({"n_bgc": 0})
pos = m[m.n_bgc > 0]
print(f"MAGs with >=1 BGC: {len(pos)}; median BGC/MAG: {pos.n_bgc.median()}; total: {int(pos.n_bgc.sum())}")
print(f"Complete BGCs (not on contig edge): {(b.contig_edge.astype(str) == 'False').sum()}")
print(pos.groupby("phylum").n_bgc.agg(["size", "median", "mean"]).sort_values("size", ascending=False).round(2).to_string())

b = b.merge(t, left_on="mag", right_on="genome", how="left")
print(pd.crosstab(b.products.str.contains("terpene-precursor"), b.bigscape_class).to_string())
g = b.groupby("bigscape_class").agg(BGCs=("bgc", "size"), MAGs=("mag", "nunique"), phyla=("phylum", "nunique"))
print(g.sort_values("MAGs", ascending=False).to_string())
