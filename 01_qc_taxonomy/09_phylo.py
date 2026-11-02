#!/usr/bin/env python
"""Étape 9 : distribution phylogénétique des MAGs procaryotes (HQ + MQ)

1. Tableau taxonomique par MAG (GTDB-Tk R232) + qualité MIMAG, échantillon, circularité, nouveauté
2. Comptes par phylum x échantillon x qualité + figure en barres
3. Arbres phylogénomiques (FastTree, WAG+Gamma) sur les alignements de marqueurs GTDB-Tk
   (bac120 pour les bactéries, ar53 pour les archées) + fichiers d'annotation iTOL

Usage :
  conda activate funcannot && python ~/09_phylo.py
Sorties : ~/forel_func/phylo/
Arbres : à visualiser sur https://itol.embl.de (glisser le .nwk puis les fichiers itol_*.txt)
"""
import os
import re
import gzip
import glob
import shutil
import subprocess
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WD = os.path.expanduser("~/forel_func")
OUT = f"{WD}/phylo"
HQ_FILE = os.environ.get("HQ_FILE") or "/mnt/gr-paoli/nas1/users/students/bastian/FOREL-2025-INPUT/mimag/HQ_MIMAG_103.tsv"
GTDBTK_ENV = "gtdbtk_new"   # contient FastTree
os.makedirs(OUT, exist_ok=True)

# palette (thème clair) : 8 teintes catégorielles fixes + gris pour "Other"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#b3b2ab"
HQ_C, MQ_C = "#1c5cab", "#9ec5f4"
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"

# ------------------------------------------------------------------ 1. tableau taxonomique
gt = pd.read_csv(f"{WD}/meta/gtdbtk_summary.tsv", sep="\t")
gt = gt.rename(columns={"user_genome": "genome"})
mq = {l.strip() for l in open(f"{WD}/meta/mags_mq.txt") if l.strip()}
gt = gt[gt.genome.isin(mq)].copy()
ok = gt.classification.astype(str).str.match(r"^d__[^;]*;p__[^;]+")
print(f"MAGs HQ+MQ : {len(gt)} | exclus (pas de taxonomie au-delà du domaine) : "
      f"{', '.join(gt.loc[~ok, 'genome']) or 'aucun'}")
gt = gt[ok].copy()

ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
def split_tax(c):
    parts = [p.split("__", 1)[1] for p in c.split(";")]
    return [p if p else None for p in (parts + [None] * 7)[:7]]
gt[ranks] = pd.DataFrame(gt.classification.map(split_tax).tolist(), index=gt.index)

hq = set()
for l in open(HQ_FILE):
    g = re.sub(r"\.(fa|fna|fasta)(\.gz)?$", "", l.split("\t")[0].strip())
    hq.add(g)
gt["quality"] = gt.genome.map(lambda g: "HQ" if g in hq else "MQ")
gt["sample"] = gt.genome.str.split("-BIN").str[0].str.split("_").str[0]
circ = pd.read_csv(f"{WD}/meta/bin_circularity.tsv", sep="\t", names=["genome", "structure"])
gt = gt.merge(circ, on="genome", how="left")
gt["novel_species"] = gt.species.isna()

cols = ["genome", "sample", "quality", "structure", "novel_species"] + ranks + ["classification"]
gt[cols].sort_values(ranks[:6] + ["genome"]).to_csv(f"{OUT}/mag_taxonomy.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 2. comptes par phylum
cnt = gt.pivot_table(index="phylum", columns=["sample", "quality"], values="genome",
                     aggfunc="count", fill_value=0)
cnt.columns = [f"{s}_{q}" for s, q in cnt.columns]
for _s in sorted(gt["sample"].unique()):
    for _q in ("HQ", "MQ"):
        if f"{_s}_{_q}" not in cnt.columns:
            cnt[f"{_s}_{_q}"] = 0
cnt = cnt[sorted(cnt.columns)]
cnt["total"] = cnt.sum(axis=1)
cnt["pct"] = 100 * cnt.total / cnt.total.sum()
cnt["novel_species"] = gt.groupby("phylum").novel_species.sum()
cnt["circular"] = gt.assign(c=gt.structure.eq("circular")).groupby("phylum").c.sum()
cnt = cnt.sort_values("total", ascending=False)
cnt.to_csv(f"{OUT}/phylum_counts.tsv", sep="\t", float_format="%.1f")
print(f"\n{len(cnt)} phylums, {gt['class'].nunique()} classes, {gt.order.nunique()} ordres, "
      f"{gt.family.nunique()} familles, {gt.genus.nunique()} genres")
print(cnt.to_string(float_format="%.1f"))

top_cls = (gt.groupby(["phylum", "class"]).size().rename("n").reset_index()
           .sort_values("n", ascending=False))
top_cls.to_csv(f"{OUT}/class_counts.tsv", sep="\t", index=False)
print("\nClasses les plus représentées :\n" + top_cls.head(10).to_string(index=False))

# ------------------------------------------------------------------ figure : barres par phylum
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.edgecolor": TEXT2,
                     "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.labelcolor": TEXT,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
samples = sorted(gt["sample"].unique())
order = cnt.index[::-1]
fig, axes = plt.subplots(1, len(samples), figsize=(6.5, 0.22 * len(order) + 1.2), sharey=True)
axes = [axes] if len(samples) == 1 else axes
xmax = max((cnt[f"{s}_HQ"] + cnt[f"{s}_MQ"]).max() for s in samples) * 1.15
for ax, s in zip(axes, samples):
    h = cnt.reindex(order)[f"{s}_HQ"].values
    m = cnt.reindex(order)[f"{s}_MQ"].values
    y = range(len(order))
    ax.barh(y, h, color=HQ_C, height=0.7, edgecolor="white", linewidth=1, label="High quality")
    ax.barh(y, m, left=h, color=MQ_C, height=0.7, edgecolor="white", linewidth=1, label="Medium quality")
    for yi, t in zip(y, h + m):
        if t:
            ax.text(t + xmax * 0.01, yi, str(int(t)), va="center", fontsize=6.5, color=TEXT2)
    ax.set_xlim(0, xmax)
    ax.set_title(f"{s} (n={int((h + m).sum())})", loc="left", fontsize=8, color=TEXT, fontweight="bold")
    ax.set_xlabel("MAGs")
    ax.xaxis.grid(True, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
axes[0].set_yticks(range(len(order)), order)
axes[0].legend(frameon=False, loc="lower right", fontsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_phylum_distribution.png", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_phylum_distribution.pdf", bbox_inches="tight")

# ------------------------------------------------------------------ 3. arbres + iTOL
phyla_rank = list(cnt.index)
pcol = {p: (CAT[i] if i < len(CAT) else OTHER) for i, p in enumerate(phyla_rank)}

def itol_colorstrip(path, label, mapping, colors):
    with open(path, "w") as f:
        present = set(mapping.values())
        names = {}
        for v, c in colors.items():          # ordre du dictionnaire = ordre de la légende
            if v in present:
                names.setdefault(c, v if c != OTHER else "Other")
        legend = [c for c in names if c != OTHER] + ([OTHER] if OTHER in names else [])
        f.write(f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t#000000\n"
                f"STRIP_WIDTH\t25\nLEGEND_TITLE\t{label}\n"
                f"LEGEND_SHAPES\t" + "\t".join("1" for _ in legend) + "\n"
                f"LEGEND_COLORS\t" + "\t".join(legend) + "\n"
                f"LEGEND_LABELS\t" + "\t".join(names[c] for c in legend) + "\nDATA\n")
        for g, v in mapping.items():
            f.write(f"{g}\t{colors[v]}\t{v}\n")

def itol_binary(path, label, members, color):
    with open(path, "w") as f:
        f.write(f"DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t{color}\n"
                f"FIELD_SHAPES\t2\nFIELD_LABELS\t{label}\nFIELD_COLORS\t{color}\nDATA\n")
        for g, v in members.items():
            f.write(f"{g}\t{1 if v else -1}\n")

g_idx = gt.set_index("genome")
itol_colorstrip(f"{OUT}/itol_phylum.txt", "Phylum", g_idx.phylum.to_dict(), pcol)
itol_colorstrip(f"{OUT}/itol_quality.txt", "MIMAG quality", g_idx.quality.to_dict(), {"HQ": HQ_C, "MQ": MQ_C})
itol_colorstrip(f"{OUT}/itol_sample.txt", "Sample", g_idx["sample"].to_dict(),
                dict(zip(samples, ["#52514e", "#c3c2b7", "#8a8984"])))
itol_binary(f"{OUT}/itol_circular.txt", "Circular", g_idx.structure.eq("circular").to_dict(), "#0b0b0b")
itol_binary(f"{OUT}/itol_novel_species.txt", "Novel species", g_idx.novel_species.to_dict(), "#e34948")
with open(f"{OUT}/itol_labels.txt", "w") as f:
    f.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for g, r in g_idx.iterrows():
        v = lambda k: r[k] if isinstance(r[k], str) and r[k] else None
        lab = v("species") or (f"{v('genus')} sp." if v("genus") else
                               f"{v('family') or v('order') or v('class') or v('phylum')} bacterium")
        f.write(f"{g}\t{lab} ({g.split('-BIN_')[0].split('_')[0][:3]}{int(g.split('_')[-1])})\n")

keep = set(gt.genome)
for dom, pat in [("bacteria", "*bac120*user_msa.fasta*"), ("archaea", "*ar53*user_msa.fasta*")]:
    hits = glob.glob(f"{WD}/gtdbtk/**/{pat}", recursive=True)
    if not hits:
        print(f"\n[{dom}] alignement GTDB-Tk introuvable ({pat}) - arbre non construit")
        continue
    msa = hits[0]
    opener = gzip.open if msa.endswith(".gz") else open
    aln, n, name, seq = f"{OUT}/{dom}_msa.faa", 0, None, []
    with opener(msa, "rt") as fi, open(aln, "w") as fo:
        def flush():
            global n
            if name in keep:
                fo.write(f">{name}\n{''.join(seq)}\n")
                n += 1
        for line in fi:
            if line.startswith(">"):
                if name: flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        if name: flush()
    print(f"\n[{dom}] {n} MAGs dans l'alignement ({os.path.basename(msa)})")
    if n < 3:
        continue
    tree = f"{OUT}/{dom}_tree.nwk"
    if shutil.which("FastTree"):
        cmd = ["FastTree", "-wag", "-gamma", "-quiet", "-out", tree, aln]
    else:
        cmd = ["conda", "run", "-n", GTDBTK_ENV, "FastTree", "-wag", "-gamma", "-quiet", "-out", tree, aln]
    try:
        subprocess.run(cmd, check=True)
        print(f"[{dom}] arbre : {tree}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[{dom}] FastTree a échoué : {e}")

print(f"\nTerminé. Fichiers dans {OUT}")
