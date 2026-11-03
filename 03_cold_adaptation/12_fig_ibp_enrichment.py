#!/usr/bin/env python3
"""
Figure - IBP enrichment in Greenland MAGs vs the global ocean (OMDB v2.0).

  a  % genomes encoding >= 1 ice-binding protein (PF11999), log scale.
     Filled = all genomes, open = one point per species (a species counts as a
     carrier if >= 50 % of its genomes carry an IBP).
  b  Detection control: % genomes with recA (same groups, same filter).
  c  Same comparison inside the lineages that carry IBPs in Greenland
     (--level family or genus), so the enrichment is not just a taxonomic-composition effect.

Input: ibp_enrichment_genomes.tsv written by 11_ibp_enrichment.py
Also prints genus- and species-level tables for the carrier lineages.

Example:
  python scripts/12_fig_ibp_enrichment.py --genomes results/omdb/comp90/ibp_enrichment_genomes.tsv \
      --level family --out results/figures/fig_ibp_enrichment_comp90
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

GROUPS = ["Greenland (this study)", "OMDB polar", "OMDB non-polar"]
LABELS = {"Greenland (this study)": "Greenland MAGs", "OMDB polar": "Polar ocean\n(|lat| ≥ 60°)",
          "OMDB non-polar": "Rest of the ocean"}
COLORS = {"Greenland (this study)": "#eb6834", "OMDB polar": "#2a78d6", "OMDB non-polar": "#8a8984"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#d9d8d2"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(c - h, 0), min(c + h, 1)


def prev(x, col):
    k, n = int(x[col].sum()), len(x)
    lo, hi = wilson(k, n)
    return dict(k=k, n=n, pct=100 * k / n if n else np.nan, lo=100 * lo, hi=100 * hi)


def fmt_p(p):
    if p >= 0.01:
        return f"p = {p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"p = {m}×10$^{{{int(e)}}}$"


def table(d, by, col):
    rows = []
    for key, x in d.groupby(by):
        key = key if isinstance(key, tuple) else (key,)
        rows.append({**dict(zip(by, key)), **prev(x, col)})
    return pd.DataFrame(rows)


def species_derep(d):
    """One row per (group, species); a species is a carrier if >= 50 % of its genomes carry
    an IBP (majority rule: robust to one odd genome, unlike 'first' or 'any')."""
    agg = d.groupby(["group", "species"]).agg(has_ibp=("has_ibp", "mean"), has_recA=("has_recA", "mean"),
                                              family=("family", "first"), genus=("genus", "first"))
    agg["has_ibp"] = agg.has_ibp >= 0.5
    agg["has_recA"] = agg.has_recA >= 0.5
    return agg.reset_index()


def dot_panel(ax, rows, log, xlab, title, show_ylabels=True, counts=True):
    y = np.arange(len(GROUPS))[::-1]
    for yi, g in zip(y, GROUPS):
        for r, off, filled in rows[g]:
            if r is None or not r["n"]:
                continue
            col = COLORS[g]
            lo, hi, pct = r["lo"], r["hi"], r["pct"]
            if log:  # zero cannot be drawn on a log axis: draw at the axis floor
                lo, pct = max(lo, 0.05), max(pct, 0.05)
            ax.plot([lo, hi], [yi + off] * 2, color=col, lw=1.0, solid_capstyle="round", zorder=2)
            ax.plot(pct, yi + off, marker="o", ms=4.6, mew=0.9, mec=col, mfc=col if filled else "white",
                    ls="none", zorder=3)
            if filled and counts:
                ax.text(hi * (1.25 if log else 1) + (0 if log else 2), yi + off, f"{r['k']}/{r['n']}",
                        va="center", fontsize=5.8, color=INK)
    ax.set_yticks(y, [LABELS[g] for g in GROUPS] if show_ylabels else [""] * len(GROUPS))
    ax.set_ylim(-0.6, len(GROUPS) - 0.4)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlab)
    ax.set_title(title, fontsize=7, loc="left", color=INK)
    if log:
        ax.set_xscale("log")
        ax.set_xlim(0.05, 100)
        ticks = [0.1, 1, 10, 100]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True, help="ibp_enrichment_genomes.tsv from step 11")
    ap.add_argument("--level", choices=["family", "genus"], default="family")
    ap.add_argument("--out", required=True, help="output prefix (no extension)")
    ap.add_argument("--min-gl", type=int, default=2,
                    help="panel c: only lineages with >= this many Greenland genomes")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    d = pd.read_csv(a.genomes, sep="\t")
    d["species"] = d.species.fillna(d.genome).astype(str).str.replace(r"^s__", "", regex=True)
    d["family"] = d.family.fillna("").astype(str).str.replace(r"^f__", "", regex=True)
    d["genus"] = d.species.str.split().str[0]
    d.loc[~d.species.str.contains(" "), "genus"] = ""      # unnamed species (MOTU ids): no genus
    sp = species_derep(d)

    # ---- tables printed for the text -------------------------------------------------------
    carriers = d[(d.group == GROUPS[0]) & d.has_ibp]
    for lvl in ["family", "genus", "species"]:
        keys = [k for k in carriers[lvl].unique() if k]
        t = table(d[d[lvl].isin(keys)], [lvl, "group"], "has_ibp")
        print(f"\n== IBP prevalence in carrier {lvl} (all genomes) ==")
        print(t[[lvl, "group", "k", "n", "pct"]].round(1).to_string(index=False))

    # ---- figure ----------------------------------------------------------------------------
    n_gl = d[d.group == GROUPS[0]][a.level].value_counts()
    lin = [k for k in carriers[a.level].value_counts().index if k and n_gl.get(k, 0) >= a.min_gl]
    # sort by Greenland prevalence (highest on top)
    lin = sorted(lin, key=lambda L: -d[(d[a.level] == L) & (d.group == GROUPS[0])].has_ibp.mean())
    fig = plt.figure(figsize=(7.2, max(2.5, 0.55 * len(lin) + 0.9)))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 0.8, 1.45], wspace=0.12)
    axa, axb = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    ra = {g: [(prev(d[d.group == g], "has_ibp"), 0.13, True), (prev(sp[sp.group == g], "has_ibp"), -0.13, False)]
          for g in GROUPS}
    dot_panel(axa, ra, True, "Genomes with ≥ 1 IBP (%, log scale)", "a  Ice-binding proteins")
    rb = {g: [(prev(d[d.group == g], "has_recA"), 0.13, True), (prev(sp[sp.group == g], "has_recA"), -0.13, False)]
          for g in GROUPS}
    dot_panel(axb, rb, False, "Genomes with recA (%)", "b  Control: recA", show_ylabels=False, counts=False)
    axb.set_xlim(0, 115)
    axb.set_xticks([0, 50, 100])

    # c: one small row block per lineage
    axc = fig.add_subplot(gs[2])
    ypos, ylab = [], []
    y0 = 0
    for L in lin[::-1]:
        for j, g in enumerate(GROUPS[::-1]):
            r = prev(d[(d[a.level] == L) & (d.group == g)], "has_ibp")
            yy = y0 + j * 0.8
            if r["n"]:
                axc.plot([r["lo"], r["hi"]], [yy] * 2, color=COLORS[g], lw=1.0, solid_capstyle="round")
                axc.plot(r["pct"], yy, marker="o", ms=4.6, mfc=COLORS[g], mec=COLORS[g], ls="none")
                axc.text(min(r["hi"], 100) + 4, yy, f"{r['k']}/{r['n']}", va="center", fontsize=5.8, color=INK)
        gl = d[(d[a.level] == L) & (d.group == GROUPS[0])]
        po = d[(d[a.level] == L) & (d.group == GROUPS[1])]
        if len(po):
            _, pv = fisher_exact([[gl.has_ibp.sum(), (~gl.has_ibp).sum()], [po.has_ibp.sum(), (~po.has_ibp).sum()]])
            axc.annotate(f"{fmt_p(pv)} vs polar", xy=(1.02, y0 + 0.15), xycoords=("axes fraction", "data"),
                         fontsize=5.8, color=MUTED, va="center")
        ypos.append(y0 + 1.1)
        ylab.append(L)
        y0 += 3.2
    axc.set_yticks(ypos, ylab, fontstyle="italic" if a.level == "genus" else "normal")
    axc.yaxis.tick_right()
    axc.tick_params(axis="y", length=0)
    axc.set_xlim(0, 125)
    axc.set_xticks([0, 25, 50, 75, 100])
    axc.set_ylim(-0.7, y0 - 1.3)
    axc.xaxis.grid(True, color=GRID, lw=0.4)
    axc.set_axisbelow(True)
    axc.spines[["top", "left"]].set_visible(False)
    axc.set_xlabel("Genomes with ≥ 1 IBP (%)")
    axc.set_title(f"c  Within IBP-carrying {'families' if a.level == 'family' else 'genera'}", fontsize=7,
                  loc="left", color=INK)

    h = [plt.Line2D([], [], marker="o", ls="none", mfc=COLORS[g], mec=COLORS[g], ms=4.6, label=LABELS[g].replace("\n", " "))
         for g in GROUPS]
    h += [plt.Line2D([], [], marker="o", ls="none", mfc="white", mec=MUTED, ms=4.6, label="one point per species (a, b)")]
    fig.legend(handles=h, frameon=False, fontsize=6, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.02 / fig.get_figheight()),
               handletextpad=0.3, columnspacing=1.2)
    H = fig.get_figheight()
    fig.subplots_adjust(left=0.13, right=0.84, bottom=0.68 / H, top=1 - 0.25 / H)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300)

    src = pd.concat([table(d, ["group"], "has_ibp").assign(panel="a_all", gene="IBP"),
                     table(sp, ["group"], "has_ibp").assign(panel="a_species", gene="IBP"),
                     table(d, ["group"], "has_recA").assign(panel="b_all", gene="recA"),
                     table(sp, ["group"], "has_recA").assign(panel="b_species", gene="recA"),
                     table(d[d[a.level].isin(lin)], [a.level, "group"], "has_ibp").assign(panel="c", gene="IBP")])
    src.to_csv(f"{a.out}_source_data.tsv", sep="\t", index=False)
    print(f"\n[fig] {a.out}.pdf/.png + _source_data.tsv")


if __name__ == "__main__":
    main()
