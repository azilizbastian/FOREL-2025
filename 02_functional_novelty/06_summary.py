#!/usr/bin/env python
"""Étape 6 : synthèse « Novelty in the functional potential: overall protein-level annotations »

Croise, gène par gène puis par MAG :
  - annotations Pfam (pyHMMER, GA), eggNOG-mapper, KofamScan
  - nouveauté des familles NR30 par rapport à OMDB v2
  - qualité CheckM2 et nouveauté taxonomique GTDB-Tk

Usage (env funcannot + pandas/matplotlib/scipy) :
  conda install -n funcannot -c conda-forge -y pandas matplotlib scipy
  conda activate funcannot && python ~/06_summary.py
Sorties : ~/forel_func/summary/
"""
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WD = os.path.expanduser("~/forel_func")
A = f"{WD}/annotation"
OUT = f"{WD}/summary"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- palette (thème clair, figures papier)
TEXT, TEXT2, GRID, SURF = "#0b0b0b", "#52514e", "#e1e0d9", "#ffffff"
BLUE = {"150": "#b7d3f6", "300": "#6da7ec", "450": "#2a78d6", "650": "#104281"}
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.edgecolor": TEXT2,
    "axes.labelcolor": TEXT, "xtick.color": TEXT2, "ytick.color": TEXT2,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300,
})

def log(msg):
    print(msg, flush=True)

# ---------------------------------------------------------------- 1. gènes
genes = pd.read_csv(f"{WD}/catalog/gene2mag.tsv", sep="\t", names=["gene", "mag"])

# Exclusion (même critère que les tables S1/S2) : pas de classification GTDB-Tk au-delà du domaine
if os.path.exists(f"{WD}/meta/gtdbtk_summary.tsv"):
    _gt = pd.read_csv(f"{WD}/meta/gtdbtk_summary.tsv", sep="\t")
    _ok = _gt.classification.astype(str).str.match(r"^d__[^;]*;p__[^;]+")
    EXCLUDED = set(_gt.loc[~_ok, "user_genome"])
    if EXCLUDED:
        log(f"MAGs exclus (taxonomie non définie ou domaine seul) : {', '.join(sorted(EXCLUDED))}")
        genes = genes[~genes.mag.isin(EXCLUDED)]
log(f"Gènes : {len(genes):,} | MAGs : {genes.mag.nunique()}")

# Pfam
pf = pd.read_csv(f"{A}/pfam_hits.tsv", sep="\t")
pf["is_duf"] = pf.pfam_name.str.startswith("DUF")
pf_g = pf.groupby("gene").agg(pfam_any=("pfam_acc", "size"),
                              pfam_known=("is_duf", lambda s: (~s).any()))
genes = genes.merge(pf_g, on="gene", how="left")
genes["pfam_any"] = genes.pfam_any.notna()
genes["pfam_known"] = genes.pfam_known.fillna(False).astype(bool)

# eggNOG-mapper
with open(f"{A}/forel.emapper.annotations") as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("##")]
eg = pd.DataFrame(rows[1:], columns=[c.lstrip("#") for c in rows[0]])
eg["eggnog_known"] = (eg.Description != "-") & (~eg.COG_category.isin(["S", "-"]))
eg["eggnog_ko"] = eg.KEGG_ko != "-"
eg = eg.rename(columns={"query": "gene"})[["gene", "COG_category", "Description",
                                           "eggnog_known", "eggnog_ko"]]
genes = genes.merge(eg, on="gene", how="left")
genes["eggnog_any"] = genes.Description.notna()
for c in ["eggnog_known", "eggnog_ko"]:
    genes[c] = genes[c].fillna(False).astype(bool)

# KofamScan (optionnel tant que la base n'est pas prête)
kofam_ok = os.path.exists(f"{A}/kofam_best.tsv")
if kofam_ok:
    ko = pd.read_csv(f"{A}/kofam_best.tsv", sep="\t",
                     names=["gene", "KO", "ko_score", "ko_evalue", "ko_def"])
    genes = genes.merge(ko[["gene", "KO"]], on="gene", how="left")
    genes["kofam"] = genes.KO.notna()
else:
    log("ATTENTION : kofam_best.tsv absent - KofamScan exclu de cette version")
    genes["kofam"] = False

# Statut fonctionnel : connu / homologue sans fonction / aucun hit
known = genes.kofam | genes.eggnog_known | genes.pfam_known
anyhit = genes.pfam_any | genes.eggnog_any | genes.kofam
genes["annot_status"] = "No hit"
genes.loc[anyhit, "annot_status"] = "Homolog, unknown function"
genes.loc[known, "annot_status"] = "Function assigned"

# Nouveauté vs OMDB v2
nov = pd.read_csv(f"{WD}/novelty/gene_novelty.tsv", sep="\t",
                  names=["gene", "nr30_rep", "omdb_class"])
genes = genes.merge(nov, on="gene", how="left")
genes.to_csv(f"{OUT}/gene_table.tsv.gz", sep="\t", index=False)

# ---------------------------------------------------------------- 2. résumé global
n = len(genes)
pct = lambda s: 100 * s.sum() / n
glob = pd.DataFrame({
    "metric": ["Pfam (any)", "Pfam (non-DUF)", "eggNOG (any hit)", "eggNOG (informative)",
               "KO (KofamScan)", "KO (eggNOG)", "Any annotation", "Function assigned",
               "No hit"],
    "pct_genes": [pct(genes.pfam_any), pct(genes.pfam_known), pct(genes.eggnog_any),
                  pct(genes.eggnog_known), pct(genes.kofam), pct(genes.eggnog_ko),
                  pct(anyhit), pct(known), pct(~anyhit)],
})
glob.to_csv(f"{OUT}/annotation_overall.tsv", sep="\t", index=False, float_format="%.2f")
log("\n" + glob.to_string(index=False, float_format="%.1f"))

fam = pd.read_csv(f"{WD}/novelty/family_novelty.tsv", sep="\t",
                  names=["nr30_rep", "best_id", "omdb_class"])
order = ["novel_vs_OMDB", "distant_30-50", "related_50-90", "near_identical_90+"]
fam_pct = fam.omdb_class.value_counts(normalize=True).reindex(order).fillna(0) * 100
gen_pct = genes.omdb_class.value_counts(normalize=True).reindex(order).fillna(0) * 100
nov_tab = pd.DataFrame({"families_NR30_pct": fam_pct, "genes_pct": gen_pct,
                        "families_n": fam.omdb_class.value_counts().reindex(order).fillna(0).astype(int)})
nov_tab.to_csv(f"{OUT}/novelty_vs_omdb.tsv", sep="\t", float_format="%.2f")
log("\n" + nov_tab.to_string(float_format="%.1f"))

# Croisement : les familles nouvelles sont-elles annotables ?
cross = pd.crosstab(genes.omdb_class, genes.annot_status, normalize="index") * 100
cross = cross.reindex(order)
cross.to_csv(f"{OUT}/novelty_x_annotation.tsv", sep="\t", float_format="%.2f")
log("\n" + cross.to_string(float_format="%.1f"))

# ---------------------------------------------------------------- 3. par MAG
mag = genes.groupby("mag").agg(
    n_genes=("gene", "size"),
    frac_function=("annot_status", lambda s: (s == "Function assigned").mean()),
    frac_unknown=("annot_status", lambda s: (s == "Homolog, unknown function").mean()),
    frac_nohit=("annot_status", lambda s: (s == "No hit").mean()),
    frac_novel_omdb=("omdb_class", lambda s: (s == "novel_vs_OMDB").mean()),
).reset_index()
mag["sample"] = mag.mag.str.split("-BIN").str[0]
q = pd.read_csv(f"{WD}/checkm2/quality_report.tsv", sep="\t")
mag = mag.merge(q[["Name", "Completeness", "Contamination", "Genome_Size"]],
                left_on="mag", right_on="Name", how="left").drop(columns="Name")
gt_ok = os.path.exists(f"{WD}/meta/mag_novelty.tsv")
if gt_ok:
    tx = pd.read_csv(f"{WD}/meta/mag_novelty.tsv", sep="\t", names=["mag", "tax_novelty"])
    gs = pd.read_csv(f"{WD}/meta/gtdbtk_summary.tsv", sep="\t")[["user_genome", "classification"]]
    mag = mag.merge(tx, on="mag", how="left").merge(gs, left_on="mag",
                                                    right_on="user_genome", how="left").drop(columns="user_genome")
else:
    log("ATTENTION : sorties GTDB-Tk absentes - pas de ventilation taxonomique")
mag.to_csv(f"{OUT}/mag_summary.tsv", sep="\t", index=False, float_format="%.4f")

# ---------------------------------------------------------------- 4. figures
def stacked_row(ax, y, values, colors, labels):
    left = 0
    for v, c, lab in zip(values, colors, labels):
        ax.barh(y, v, left=left, color=c, edgecolor=SURF, linewidth=1.5, height=0.6, label=lab)
        if v >= 6:
            ink = SURF if c in (BLUE["450"], BLUE["650"]) else TEXT
            ax.text(left + v / 2, y, f"{v:.0f}%", ha="center", va="center", color=ink, fontsize=7)
        left += v

fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.5, 3.6), gridspec_kw={"hspace": 1.4})
st = ["Function assigned", "Homolog, unknown function", "No hit"]
v = genes.annot_status.value_counts(normalize=True).reindex(st).fillna(0) * 100
stacked_row(a1, 0, v.values, [BLUE["650"], BLUE["300"], BLUE["150"]], st)
a1.set_yticks([0], [f"Genes (n={n:,})"])
a1.set_title("a  Functional annotation (Pfam, eggNOG, KOfam)" + ("" if kofam_ok else " - KOfam pending"),
             loc="left", fontsize=8, color=TEXT, fontweight="bold")
a1.legend(ncol=3, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.55), fontsize=7)

labs = ["No OMDB v2 homolog", "30-50% identity", "50-90% identity", "≥90% identity"]
cols = [BLUE["650"], BLUE["450"], BLUE["300"], BLUE["150"]]
stacked_row(a2, 1, fam_pct.values, cols, labs)
stacked_row(a2, 0, gen_pct.values, cols, [None] * 4)
a2.set_yticks([1, 0], [f"Families NR30 (n={len(fam):,})", "Genes"])
a2.set_title("b  Similarity to the OMDB v2 protein catalog (NR30)", loc="left",
             fontsize=8, color=TEXT, fontweight="bold")
handles, hl = a2.get_legend_handles_labels()
a2.legend(handles[:4], hl[:4], ncol=4, frameon=False, loc="upper left",
          bbox_to_anchor=(0, -0.55), fontsize=7)
for ax in (a1, a2):
    ax.set_xlim(0, 100)
    ax.set_xlabel("%", color=TEXT2)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
fig.savefig(f"{OUT}/fig_annotation_novelty.png", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_annotation_novelty.pdf", bbox_inches="tight")

if gt_ok:
    lv = ["species", "novel_species", "novel_genus", "novel_family", "novel_order"]
    lv = [l for l in lv if l in set(mag.tax_novelty)]
    names = {"species": "Known\nspecies", "novel_species": "Novel\nspecies",
             "novel_genus": "Novel\ngenus", "novel_family": "Novel\nfamily", "novel_order": "Novel\norder"}
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.8), sharex=True)
    for ax, col, ttl in [(axs[0], "frac_novel_omdb", "c  Genes in families absent from OMDB v2"),
                         (axs[1], "frac_nohit", "d  Genes without any annotation")]:
        for i, l in enumerate(lv):
            vals = mag.loc[mag.tax_novelty == l, col] * 100
            jitter = (pd.Series(range(len(vals))).sample(frac=1, random_state=1).values / max(len(vals), 1) - 0.5) * 0.5
            ax.scatter(i + jitter, vals, s=10, color=BLUE["450"], alpha=0.6, linewidths=0)
            ax.hlines(vals.median(), i - 0.3, i + 0.3, color=TEXT, linewidth=2)
        ns = [(mag.tax_novelty == l).sum() for l in lv]
        ax.set_xticks(range(len(lv)), [f"{names[l]}\n(n={k})" for l, k in zip(lv, ns)], fontsize=7)
        ax.set_ylabel("% of genes per MAG")
        ax.set_title(ttl, loc="left", fontsize=8, color=TEXT, fontweight="bold")
        ax.yaxis.grid(True, color=GRID, linewidth=0.5)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_novelty_by_taxrank.png", bbox_inches="tight")
    fig.savefig(f"{OUT}/fig_novelty_by_taxrank.pdf", bbox_inches="tight")
    try:
        from scipy.stats import kruskal, spearmanr
        groups = [mag.loc[mag.tax_novelty == l, "frac_novel_omdb"] for l in lv]
        rank = mag.tax_novelty.map({l: i for i, l in enumerate(lv)})
        log(f"\nKruskal-Wallis (frac novel vs rang) : H={kruskal(*groups).statistic:.1f}, p={kruskal(*groups).pvalue:.2e}")
        r = spearmanr(rank, mag.frac_novel_omdb, nan_policy="omit")
        log(f"Spearman (rang de nouveauté vs frac novel) : rho={r.statistic:.2f}, p={r.pvalue:.2e}")
        log("Médianes par rang :\n" + mag.groupby("tax_novelty")[["frac_novel_omdb", "frac_nohit"]].median().reindex(lv).to_string(float_format="%.3f"))
    except ImportError:
        log("scipy absent : tests non calculés")

log(f"\nTerminé. Fichiers dans {OUT}")
