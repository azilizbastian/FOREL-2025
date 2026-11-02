#!/usr/bin/env python3
import argparse

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BINS = [
    ("ANACONDA", "bin_5", 15.62, 0.62, "Dinobryaceae"),
    ("ANACONDA", "bin_13", 14.33, 10.67, "Eimeriorina"),
    ("ANACONDA", "bin_4", 26.8, 0.2, "Bathycoccaceae"),
    ("ANACONDA", "bin_1", 94.0, 7.4, "Micromonas"),
    ("ANACONDA", "bin_3", 64.09, 7.18, "MAST-4"),
    ("ANACONDA", "bin_2", 82.4, 2.2, "Micromonas"),
    ("ANACONDA", "bin_0", 75.2, 0.4, "Mamiellaceae"),
    ("ANACONDA", "bin_9", 9.2, 5.0, "Micromonas"),
    ("COBRA", "bin_6", 33.2, 16.6, "Mamiellaceae"),
    ("COBRA", "bin_17", 7.26, 3.23, "Pyramimonadaceae"),
    ("COBRA", "bin_5", 16.6, 5.0, "Mamiellaceae"),
    ("COBRA", "bin_40", 17.8, 5.6, "Micromonas"),
    ("COBRA", "bin_13", 2.8, 0.0, "Synurales"),
    ("COBRA", "bin_4", 38.2, 15.2, "Micromonas"),
    ("COBRA", "bin_39", 44.6, 5.8, "Micromonas"),
    ("COBRA", "bin_8", 21.69, 4.82, "Sar"),
    ("COBRA", "bin_1", 73.4, 7.0, "Micromonas"),
    ("COBRA", "bin_3", 33.4, 10.8, "Micromonas"),
    ("COBRA", "bin_2", 30.2, 6.6, "Mamiellaceae"),
    ("COBRA", "bin_0", 53.8, 0.4, "Mamiellaceae"),
    ("COBRA", "bin_26", 32.4, 0.2, "Micromonas"),
    ("COBRA", "bin_22", 1.6, 0.4, "Chrysochromulinaceae"),
    ("COBRA", "bin_9", 10.4, 9.6, "Micromonas"),
    ("COBRA", "bin_18", 5.6, 0.4, "Micromonas"),
]
GROUPS = {
    "Micromonas": ("Micromonas", "#1b7837"),
    "Mamiellaceae": ("Other Mamiellophyceae", "#a6dba0"),
    "Bathycoccaceae": ("Other Mamiellophyceae", "#a6dba0"),
}
OTHER = ("Other eukaryotes", "#9970ab")
MARKERS = {"ANACONDA": "o", "COBRA": "^"}


def n50(lengths):
    total, acc = sum(lengths), 0
    for rank, L in enumerate(sorted(lengths, reverse=True), 1):
        acc += L
        if acc >= total / 2:
            return L, rank


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="MAG shown in panel b")
    ap.add_argument("--label", default="MAG")
    ap.add_argument("--highlight", default="ANACONDA:bin_1")
    ap.add_argument("--out", default="overview")
    args = ap.parse_args()
    hl_ds, hl_bin = args.highlight.split(":")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8),
                                   gridspec_kw={"width_ratios": [1.15, 1]})

    ax1.axvspan(50, 100, ymin=0, ymax=10 / 18, color="#e8e8e8", zorder=0)
    ax1.axvspan(90, 100, ymin=0, ymax=5 / 18, color="#c8c8c8", zorder=0)
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 18)
    ax1.text(51, 9.6, "MIMAG medium\n(≥50%, <10%)",
             fontsize=7, va="top", color="#555555")
    ax1.text(91, 4.6, "high\n(>90%, <5%)", fontsize=7, va="top", color="#333333")
    n_medium = 0
    for ds, b, comp, cont, taxon in BINS:
        name, color = GROUPS.get(taxon, OTHER)
        ax1.scatter(comp, cont, s=55, marker=MARKERS[ds], c=color,
                    edgecolors="black", linewidths=0.5, zorder=3)
        if comp >= 50 and cont < 10:
            n_medium += 1
        if ds == hl_ds and b == hl_bin:
            ax1.scatter(comp, cont, s=260, marker="o", facecolors="none",
                        edgecolors="black", linewidths=1.4, zorder=4)
            ax1.annotate(args.label, (comp, cont), xytext=(-6, 52),
                         textcoords="offset points", ha="right", fontsize=8,
                         arrowprops=dict(arrowstyle="-", lw=0.7))
    ax1.set_xlabel("Completeness (%)")
    ax1.set_ylabel("Contamination (%)")
    handles = [Patch(fc=c, ec="black", lw=0.5, label=n) for n, c in
               [("Micromonas", "#1b7837"), ("Other Mamiellophyceae", "#a6dba0"),
                ("Other eukaryotes", "#9970ab")]]
    handles += [Line2D([], [], marker=m, ls="", mfc="white", mec="black", label=d)
                for d, m in MARKERS.items()]
    ax1.legend(handles=handles, loc="upper left", fontsize=7.5, frameon=False)
    ax1.set_title(f"a  Eukaryotic bins recovered (n = {len(BINS)})", loc="left",
                  fontsize=10, fontweight="bold")
    print(f"Bins with completeness >=50% and contamination <10%: {n_medium}/{len(BINS)}")

    lengths = [len(r.seq) for r in SeqIO.parse(args.fasta, "fasta")]
    lengths.sort(reverse=True)
    total = sum(lengths)
    cum = np.cumsum(lengths) / 1e6
    ranks = np.arange(1, len(lengths) + 1)
    N50, L50 = n50(lengths)
    ax2.plot(ranks, cum, color="#1b7837", lw=2)
    ax2.axhline(total / 2e6, color="grey", ls="--", lw=0.8)
    ax2.axvline(L50, color="grey", ls="--", lw=0.8)
    ax2.scatter([L50], [total / 2e6], color="black", zorder=3, s=25)
    ax2.set_xscale("log")
    ax2.set_xlabel("Contig rank (longest first, log scale)")
    ax2.set_ylabel("Cumulative length (Mb)")
    ax2.set_ylim(0, total / 1e6 * 1.05)
    txt = (f"{len(lengths)} contigs, {total / 1e6:.2f} Mb\n"
           f"N50 = {N50 / 1e3:,.0f} kb (L50 = {L50})\n"
           f"longest contig = {lengths[0] / 1e6:.2f} Mb")
    ax2.text(0.97, 0.05, txt, transform=ax2.transAxes, ha="right", va="bottom",
             fontsize=8.5, bbox=dict(fc="white", ec="#cccccc", boxstyle="round"))
    ax2.set_title(f"b  Contiguity of {args.label}", loc="left", fontsize=10,
                  fontweight="bold")
    print(f"{len(lengths)} contigs; total {total} bp; N50 {N50} bp; L50 {L50}; "
          f"longest {lengths[0]} bp")

    for a in (ax1, ax2):
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{args.out}.{ext}", dpi=300, bbox_inches="tight")
    print(f"Saved {args.out}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
