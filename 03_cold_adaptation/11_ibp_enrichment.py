#!/usr/bin/env python3
"""
Step 11 - Are ice-binding proteins (IBPs) enriched in the Greenland MAGs compared
with the global ocean (OMDB v2.0)?

Groups (same quality filter for all: completeness >= --min-comp, contamination <= --max-cont):
  * Greenland  : our MAGs
  * OMDB polar : OMDB marine MAGs from samples with |latitude| >= --polar-lat
  * OMDB non-polar : other OMDB marine MAGs
Controls:
  * recA (near-universal gene) must be ~equally prevalent in all groups;
  * logistic regression  IBP ~ group + phylum + completeness  (taxonomy & completeness effects);
  * species-dereplicated version (one genome per species and group), because some
    OMDB species are sampled hundreds of times.

Outputs (in --outdir): ibp_enrichment_genomes.tsv, ibp_prevalence.tsv,
ibp_prevalence_by_phylum.tsv, ibp_logit.tsv, fig_ibp_enrichment.pdf/png

Example:
  python scripts/11_ibp_enrichment.py \
    --omdb-hits results/omdb/omdb_hits.tsv --mag-hits results/omdb/mags_hits.tsv \
    --omdb-genomes $M/OMDv2.0_genome-metadata-20250702.tsv.gz --omdb-samples $M/omdb21.samples.tsv.gz \
    --mags results/genomes_mags.tsv --exclude ctg45793 --outdir results/omdb
"""
import argparse
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

GROUPS = ["Greenland (this study)", "OMDB polar", "OMDB non-polar"]
COLORS = {"Greenland (this study)": "#eb6834", "OMDB polar": "#2a78d6", "OMDB non-polar": "#8a8984"}


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def load(a):
    hits_o = pd.read_csv(a.omdb_hits, sep="\t")
    hits_m = pd.read_csv(a.mag_hits, sep="\t")
    hits_m = hits_m[~hits_m.genome.isin(a.exclude) | (hits_m.profile != "Ice_binding")]

    g = pd.read_csv(a.omdb_genomes, sep="\t", low_memory=False)
    g = g.loc[:, ~g.columns.str.startswith("Unnamed")]
    s = pd.read_csv(a.omdb_samples, sep="\t", low_memory=False)
    s = s[["SUSHISAMPLE", "Ecosystem_type", "Ecosystem", "Latitude_degN"]]
    g = g.merge(s, on="SUSHISAMPLE", how="left")
    is_mag = g.IS_MAG.astype(str).str.lower().eq("true")
    marine = g.Ecosystem_type.fillna("").str.lower().str.contains("marine")
    g = g[is_mag & marine & g.Latitude_degN.notna()].copy()
    g["group"] = np.where(g.Latitude_degN.abs() >= a.polar_lat, "OMDB polar", "OMDB non-polar")
    o = pd.DataFrame(dict(genome=g.GENOME, group=g.group, completeness=g.COMPLETENESS,
                          contamination=g.CONTAMINATION, phylum=g.PHYLUM, family=g.FAMILY,
                          species=g.SPECIES.fillna(g.MOTU4.astype(str)), latitude=g.Latitude_degN))

    m = pd.read_csv(a.mags, sep="\t")
    m = pd.DataFrame(dict(genome=m.genome_id, group=GROUPS[0], completeness=m.completeness,
                          contamination=m.contamination, phylum=m.get("phylum"), family=m.get("family"),
                          species=m.get("species"), latitude=np.nan))
    m["species"] = m.species.fillna("").replace("", np.nan).fillna(m.genome)
    d = pd.concat([m, o], ignore_index=True)
    d["phylum"] = d.phylum.fillna("").replace("", "Unclassified")

    hits = pd.concat([hits_m, hits_o])
    for prof, col in [("Ice_binding", "has_ibp"), ("recA", "has_recA")]:
        d[col] = d.genome.isin(set(hits.loc[hits.profile == prof, "genome"]))
    n_ibp = hits[hits.profile == "Ice_binding"].groupby("genome").size()
    d["n_ibp"] = d.genome.map(n_ibp).fillna(0).astype(int)
    unmatched = set(hits_o.genome) - set(g.GENOME)
    print(f"[enrich] OMDB hits on genomes outside the marine-MAG set (ignored): {len(unmatched):,} genomes")
    return d


def prevalence(d, by):
    rows = []
    for key, x in d.groupby(by):
        key = key if isinstance(key, tuple) else (key,)
        r = dict(zip(by, key))
        r["n_genomes"] = len(x)
        for col in ["has_ibp", "has_recA"]:
            k = int(x[col].sum())
            lo, hi = wilson(k, len(x))
            r[col.replace("has_", "") + "_pct"] = 100 * k / len(x)
            r[col.replace("has_", "") + "_lo"], r[col.replace("has_", "") + "_hi"] = 100 * lo, 100 * hi
            r[col.replace("has_", "") + "_n"] = k
        rows.append(r)
    return pd.DataFrame(rows)


def logit(d, top_phyla):
    x = d.copy()
    x["phylum_m"] = np.where(x.phylum.isin(top_phyla), x.phylum, "Other")
    x["grp"] = pd.Categorical(x.group, categories=["OMDB non-polar", "OMDB polar", GROUPS[0]])
    x["y"] = x.has_ibp.astype(int)
    res = []
    for name, f in [("IBP ~ group", "y ~ grp"),
                    ("IBP ~ group + phylum + completeness", "y ~ grp + C(phylum_m) + completeness")]:
        try:
            fit = smf.glm(f, data=x, family=sm.families.Binomial()).fit()
            for k in [c for c in fit.params.index if c.startswith("grp[")]:
                lo, hi = np.exp(fit.conf_int().loc[k])
                res.append(dict(model=name, term=k.split("T.")[1].rstrip("]"), odds_ratio=np.exp(fit.params[k]),
                                ci_low=lo, ci_high=hi, p=fit.pvalues[k]))
        except Exception as e:  # e.g. perfect separation
            res.append(dict(model=name, term=f"failed: {e}"))
    return pd.DataFrame(res)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--omdb-hits", required=True)
    ap.add_argument("--mag-hits", required=True)
    ap.add_argument("--omdb-genomes", required=True)
    ap.add_argument("--omdb-samples", required=True)
    ap.add_argument("--mags", required=True)
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--min-comp", type=float, default=70)
    ap.add_argument("--max-cont", type=float, default=10)
    ap.add_argument("--polar-lat", type=float, default=60)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    d = load(a)
    d = d[(d.completeness >= a.min_comp) & (d.contamination <= a.max_cont)]
    d.to_csv(os.path.join(a.outdir, "ibp_enrichment_genomes.tsv"), sep="\t", index=False)

    # one row per species and group; carrier if >= 50 % of its genomes carry the gene (same rule as step 12)
    derep = (d.groupby(["group", "species"])
             .agg(has_ibp=("has_ibp", "mean"), has_recA=("has_recA", "mean"))
             .reset_index())
    derep["has_ibp"] = derep.has_ibp >= 0.5
    derep["has_recA"] = derep.has_recA >= 0.5
    p_all = prevalence(d, ["group"]).assign(level="all genomes")
    p_der = prevalence(derep, ["group"]).assign(level="one genome per species")
    prev = pd.concat([p_all, p_der])
    prev["group"] = pd.Categorical(prev.group, GROUPS)
    prev = prev.sort_values(["level", "group"])
    prev.to_csv(os.path.join(a.outdir, "ibp_prevalence.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(prev[["level", "group", "n_genomes", "ibp_n", "ibp_pct", "ibp_lo", "ibp_hi", "recA_pct"]].round(2).to_string(index=False))

    # Fisher tests Greenland vs each OMDB group
    for lvl, x in [("all", d), ("species-derep", derep)]:
        gl = x[x.group == GROUPS[0]]
        for other in GROUPS[1:]:
            ot = x[x.group == other]
            tab = [[gl.has_ibp.sum(), (~gl.has_ibp).sum()], [ot.has_ibp.sum(), (~ot.has_ibp).sum()]]
            orr, p = stats.fisher_exact(tab)
            print(f"[enrich] Fisher ({lvl}) Greenland vs {other}: OR = {orr:.2f}, p = {p:.2g}")

    top = d[d.group == GROUPS[0]].phylum.value_counts()
    top = [p for p in top.index if p != "Unclassified" and top[p] >= 3]
    bp = prevalence(d[d.phylum.isin(top)], ["phylum", "group"])
    bp.to_csv(os.path.join(a.outdir, "ibp_prevalence_by_phylum.tsv"), sep="\t", index=False)
    print(bp[["phylum", "group", "n_genomes", "ibp_n", "ibp_pct"]].round(2).to_string(index=False))

    lg = logit(d, top)
    lg.to_csv(os.path.join(a.outdir, "ibp_logit.tsv"), sep="\t", index=False)
    print(lg.round(4).to_string(index=False))
    plot(p_der, bp, os.path.join(a.outdir, "fig_ibp_enrichment"))


def plot(p, bp, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6})
    p = p.set_index("group").reindex(GROUPS)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(width_ratios=[1, 1, 1.6]))
    for ax, key, title in [(axes[0], "ibp", "Ice-binding protein (PF11999)"), (axes[1], "recA", "Control: recA")]:
        y = np.arange(len(GROUPS))[::-1]
        for yi, gname in zip(y, GROUPS):
            r = p.loc[gname]
            ax.barh(yi, r[f"{key}_pct"], color=COLORS[gname], height=0.62)
            ax.plot([r[f"{key}_lo"], r[f"{key}_hi"]], [yi, yi], color="#1f1f1e", lw=0.8)
            ax.text(r[f"{key}_hi"] + 1, yi, f"{int(r[f'{key}_n'])}/{int(r.n_genomes)}", va="center", fontsize=6)
        ax.set_yticks(y, GROUPS if ax is axes[0] else [""] * len(GROUPS))
        ax.set_xlabel("% of genomes (one per species)")
        ax.set_title(title, fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[1].set_xlim(0, 110)
    ax = axes[2]
    phyla = list(dict.fromkeys(bp.phylum))
    w = 0.26
    for i, gname in enumerate(GROUPS):
        s = bp[bp.group == gname].set_index("phylum").reindex(phyla)
        ax.bar(np.arange(len(phyla)) + (i - 1) * w, s.ibp_pct, width=w, color=COLORS[gname], label=gname)
    ax.set_xticks(np.arange(len(phyla)), phyla, rotation=30, ha="right")
    ax.set_ylabel("% genomes with IBP")
    ax.set_title("By phylum (all genomes)", fontsize=7)
    ax.legend(frameon=False, fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{path}.{ext}", dpi=300)


if __name__ == "__main__":
    main()
