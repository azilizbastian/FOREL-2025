#!/usr/bin/env python
"""Complétude CheckM2 par phylum (449 MAGs HQ+MQ) : statistiques + boxplot
Usage : conda activate funcannot && python ~/10_completeness.py
Sorties : ~/forel_func/phylo/completeness_by_phylum.tsv + fig_completeness_by_phylum.png/.pdf
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal

WD = os.path.expanduser("~/forel_func")
OUT = f"{WD}/phylo"
TOP = ["Pseudomonadota", "Bacteroidota", "Actinomycetota", "Verrucomicrobiota"]
COL = dict(zip(TOP, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"

tax = pd.read_csv(f"{OUT}/mag_taxonomy.tsv", sep="\t")
q = pd.read_csv(f"{WD}/checkm2/quality_report.tsv", sep="\t")[["Name", "Completeness", "Contamination"]]
d = tax.merge(q, left_on="genome", right_on="Name")
print(f"MAGs : {len(d)}")

st = d.groupby("phylum").Completeness.agg(
    n="size", median="median", q1=lambda s: s.quantile(.25), q3=lambda s: s.quantile(.75),
    min="min", max="max").sort_values("n", ascending=False)
st["IQR"] = st.q3 - st.q1
st["pct_HQ"] = 100 * d.groupby("phylum").quality.apply(lambda s: (s == "HQ").mean())
st["pct_circular"] = 100 * d.groupby("phylum").structure.apply(lambda s: (s == "circular").mean())
st.to_csv(f"{OUT}/completeness_by_phylum.tsv", sep="\t", float_format="%.1f")
print(st.to_string(float_format="%.1f"))
print(f"\nToutes MAGs : médiane {d.Completeness.median():.1f} % "
      f"(IQR {d.Completeness.quantile(.25):.1f}-{d.Completeness.quantile(.75):.1f})")
groups = [d.loc[d.phylum == p, "Completeness"] for p in TOP]
H, p = kruskal(*groups)
print(f"Kruskal-Wallis entre les 4 phylums principaux : H = {H:.1f}, p = {p:.2e}")

# boxplot des 4 phylums principaux
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.edgecolor": TEXT2,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.labelcolor": TEXT,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
fig, ax = plt.subplots(figsize=(4.5, 3))
for i, ph in enumerate(TOP):
    v = d.loc[d.phylum == ph, "Completeness"]
    ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor=TEXT2), medianprops=dict(color=TEXT, linewidth=1.5),
               whiskerprops=dict(color=TEXT2), capprops=dict(color=TEXT2))
    jit = (pd.Series(range(len(v))).sample(frac=1, random_state=1).values / max(len(v), 1) - .5) * .35
    ax.scatter(i + jit, v, s=7, color=COL[ph], alpha=.65, linewidths=0, zorder=3)
ax.axhline(50, color=TEXT2, linestyle="--", linewidth=.8)
ax.axhline(90, color=TEXT2, linestyle=":", linewidth=.8)
ax.text(3.45, 50, "MQ", va="center", fontsize=6.5, color=TEXT2)
ax.text(3.45, 90, "HQ", va="center", fontsize=6.5, color=TEXT2)
ax.set_xticks(range(len(TOP)), [f"{p}\n(n={(d.phylum == p).sum()})" for p in TOP], fontsize=7)
ax.set_ylabel("CheckM2 completeness (%)")
ax.set_ylim(45, 101)
ax.yaxis.grid(True, color=GRID, linewidth=.5)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_completeness_by_phylum.png", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_completeness_by_phylum.pdf", bbox_inches="tight")
print(f"\nFigure : {OUT}/fig_completeness_by_phylum.png")
