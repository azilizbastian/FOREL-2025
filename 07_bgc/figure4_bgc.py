#!/usr/bin/env python3
"""Figure 4 - BGC abundance, distribution and novelty (449 MAGs, FOREL 2025).
(A) BGCs per MAG for the four major phyla (MAGs with >= 1 BGC)
(B) Average BGCs per MAG for each BGC class, by phylum
(C) Novelty vs BiG-FAM (BiG-SLiCE, novel if d > 900), BGCs and GCFs
Inputs : BGC_449/bgc_table.tsv (summarise_bgc.py + gcf_summary.py), GTDB-Tk taxonomy
Outputs: BGC_449/figure4/Figure4_BGC.{png,pdf,svg}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/mnt/gr-paoli/nas1/users/students/bastian/BGC_449"
TAXO = "/home/bastian/forel_func/phylo/mag_taxonomy.tsv"
FIG = f"{OUT}/figure4"
os.makedirs(FIG, exist_ok=True)

PHY = ["Pseudomonadota", "Bacteroidota", "Actinomycetota", "Verrucomicrobiota"]
SHORT = ["Pseudo.", "Bactero.", "Actino.", "Verruco."]
COL = ["#4A86C5", "#F2A541", "#5DB35E", "#A86BB8"]
CLASSES = ["Terpenes", "Terpene\nprecursors", "RiPPs", "NRPS", "PKS", "Others"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------- data ----------
b = pd.read_csv(f"{OUT}/bgc_table.tsv", sep="\t")
t = pd.read_csv(TAXO, sep="\t")[["genome", "phylum"]]
b = b.merge(t, left_on="mag", right_on="genome", how="left")
b["novel"] = b["novel"].astype(str) == "True"


def bgc_class(r):
    c = r.bigscape_class
    if c == "terpene":
        return "Terpenes"
    if c == "RiPP":
        return "RiPPs"
    if c == "NRPS":
        return "NRPS"
    if c in ("PKSI", "PKSother", "PKS-NRP_Hybrids"):
        return "PKS"
    if "terpene-precursor" in str(r.products):
        return "Terpene\nprecursors"
    return "Others"


b["cls"] = b.apply(bgc_class, axis=1)
per_mag = b.groupby(["mag", "phylum"]).size().rename("n").reset_index()

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8),
                         gridspec_kw={"width_ratios": [1, 1.6, 0.8], "wspace": 0.32})

# ---------- A : BGCs per MAG ----------
ax = axes[0]
rng = np.random.default_rng(1)
data = [per_mag.loc[per_mag.phylum == p, "n"].values for p in PHY]
bp = ax.boxplot(data, positions=range(4), widths=0.55, patch_artist=True,
                showfliers=False, medianprops=dict(color="black", linewidth=2))
for patch, c in zip(bp["boxes"], COL):
    patch.set_facecolor(c)
    patch.set_alpha(0.85)
for i, v in enumerate(data):
    ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v + rng.uniform(-0.15, 0.15, len(v)),
               s=6, color="0.25", alpha=0.6, zorder=3)
ax.set_xticks(range(4), [f"{s}\n(n = {len(v)})" for s, v in zip(SHORT, data)])
ax.set_ylabel("BGCs per MAG")
ax.set_title("BGC abundance per MAG")

# ---------- B : average BGCs per MAG per class ----------
ax = axes[1]
n_mag = per_mag.groupby("phylum").size()
counts = b.groupby(["phylum", "cls"]).size().unstack(fill_value=0)
avg = counts.reindex(index=PHY, columns=CLASSES, fill_value=0).div(n_mag.reindex(PHY), axis=0)
x = np.arange(len(CLASSES))
w = 0.2
for i, (p, c) in enumerate(zip(PHY, COL)):
    ax.bar(x + (i - 1.5) * w, avg.loc[p].values, w, color=c, label=p, edgecolor="none")
ax.set_xticks(x, CLASSES)
ax.set_xlabel("BGC class")
ax.set_ylabel("Average BGCs per MAG")
ax.set_title("Distribution of BGC classes")
leg = ax.legend(title="Phylum", frameon=False, loc="upper right")
for txt in leg.get_texts():
    txt.set_fontstyle("italic")

# ---------- C : novelty vs BiG-FAM ----------
ax = axes[2]
gnov = b.groupby("gcf")["novel"].all()
known = [int((~b.novel).sum()), int((~gnov).sum())]
novel = [int(b.novel.sum()), int(gnov.sum())]
tot = [k + n for k, n in zip(known, novel)]
xs = [0, 1]
ax.bar(xs, known, 0.6, color="0.85", edgecolor="black", linewidth=1.2, label="Known")
ax.bar(xs, novel, 0.6, bottom=known, color="0.45", edgecolor="black", linewidth=1.2, label="Novel")
for xi, k, n, tt in zip(xs, known, novel, tot):
    ax.text(xi, k / 2, f"{k}", ha="center", va="center", fontweight="bold")
    ax.text(xi, k + n / 2, f"{n}", ha="center", va="center", fontweight="bold", color="white")
    ax.text(xi, tt + max(tot) * 0.015, f"{tt}", ha="center", va="bottom", fontweight="bold", fontsize=10)
ax.set_xticks(xs, ["BGCs", "GCFs"])
ax.set_xlabel("vs. BiG-FAM")
ax.set_ylabel("Number of BGCs / GCFs")
ax.set_ylim(0, max(tot) * 1.1)
ax.set_title("Biosynthetic novelty")
ax.legend(frameon=True, loc="upper right", fontsize=8)

for ax, lab in zip(axes, "ABC"):
    ax.text(-0.14, 1.06, lab, transform=ax.transAxes, fontsize=15, fontweight="bold")

for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{FIG}/Figure4_BGC.{ext}", dpi=300, bbox_inches="tight")

# ---------- numbers for the legend ----------
print("A - MAGs with >=1 BGC per phylum and median:")
for p, v in zip(PHY, data):
    print(f"  {p}: n = {len(v)}, median = {np.median(v)}")
print("B - average BGCs per MAG:\n" + avg.round(2).to_string())
print(f"C - BGCs: {novel[0]}/{tot[0]} novel ; GCFs: {novel[1]}/{tot[1]} novel")
print(f"Figure written to {FIG}/Figure4_BGC.png/.pdf/.svg")
