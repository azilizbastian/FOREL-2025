#!/usr/bin/env python
"""Étape 12 : nouveauté des SGB FOREL vs catalogue arctique (Royo-Llonch et al. 2021) et OMDB v2.1

Références :
  - ARCTIC : tous les génomes OMDB des 952 échantillons arctiques (arctic-subset.list, incl. Tara Polar Circle)
  - OMDB   : génomes représentatifs des 31 058 mOTUs (OMDB 2.1 / mOTUs 4.1, IS_REP == True)
Comparaison : skani dist, représentant de chaque SGB vs références.
Un SGB est "connu" s'il a une référence à ANI >= 95 % et fraction alignée >= 30 % (query ou ref).

Usage :  conda activate funcannot && python ~/12_sgb_novelty.py
Sorties : ~/forel_func/novelty_sgb/
"""
import os
import re
import sys
import glob
import subprocess
import pandas as pd

WD = os.path.expanduser("~/forel_func")
OUT = f"{WD}/novelty_sgb"
os.makedirs(OUT, exist_ok=True)
THREADS = str(max(1, os.cpu_count() // 2))
ANI_MIN, AF_MIN = 95.0, 30.0
ARCTIC_PAT = os.environ.get("ARCTIC_PAT", r"s41564-021-00979")   # DOI Royo-Llonch 2021

MOTUS = next((m for m in ["/mnt/gr-paoli/nfs/data/databases/MOTUS/4.1",
                          "/mnt/gr-paoli/nas1/data/databases/MOTUS/4.1"] if os.path.isdir(m)), None)
if MOTUS is None:
    sys.exit("Dossier MOTUS/4.1 introuvable")
META = f"{MOTUS}/metadata/OMDB2.1/omdb21.genomes.tsv.gz"
print(f"MOTUS : {MOTUS}")

# ------------------------------------------------------------------ 1. références
g = pd.read_csv(META, sep="\t", low_memory=False, dtype=str)
g["IS_REP"] = g.IS_REP.str.lower().isin(["true", "1"])
reps = g[g.IS_REP]
# Référence arctique : tous les génomes OMDB issus des échantillons du sous-ensemble arctique
# (952 échantillons, dont Tara Oceans Polar Circle via l'étude TOPC)
ARCTIC_LIST = "/mnt/gr-paoli/nas1/data/databases/OMDB/2.0/arctic-subset/arctic-subset.list"
arc_samples = pd.read_csv(ARCTIC_LIST, sep="\t", header=None, names=["sample", "study"])
arc = g[g.SUSHISAMPLE.isin(set(arc_samples["sample"]))]
print(f"Génomes OMDB 2.1 : {len(g):,} | représentants mOTU : {len(reps):,}")
print(f"Échantillons arctiques : {len(arc_samples)} | génomes arctiques : {len(arc):,} "
      f"| dont TOPC (Tara Polar Circle) : {(arc.SUSHISTUDY == 'TOPC').sum():,}")
print(arc.SUSHISTUDY.value_counts().head(10).to_string())
if len(arc) == 0:
    sys.exit("Aucun génome arctique trouvé : vérifier la correspondance SUSHISAMPLE / arctic-subset.list")

# ------------------------------------------------------------------ 2. chemins locaux des génomes
def candidates(r):
    gf = r.get("GENOME_FILE")
    gf = gf.strip() if isinstance(gf, str) else ""
    rel = re.sub(r"^.*?/genomes/genomes/", "", gf) if "/genomes/genomes/" in gf else gf.lstrip("/")
    st, sa, ge = r.get("SUSHISTUDY"), r.get("SUSHISAMPLE"), r.get("GENOME")
    return [f"{MOTUS}/genomes/{rel}", f"{MOTUS}/{rel}",
            f"{MOTUS}/genomes/{st}/{sa}/{ge}/{ge}.fa.gz", f"{MOTUS}/genomes/{st}/{sa}/{ge}.fa.gz",
            f"{MOTUS}/genomes/{st}/{sa}/{ge}.fa", f"{MOTUS}/genomes/{st}/{sa}/{ge}/{ge}.fa"]

# repère le motif de chemin valide sur un échantillon, puis l'applique à tous
probe = pd.concat([reps.head(15), arc.head(15)])
hits = [next((i for i, p in enumerate(candidates(r)) if os.path.exists(p)), None)
        for _, r in probe.iterrows()]
found = [h for h in hits if h is not None]
if not found:
    r = probe.iloc[0]
    print("Aucun chemin trouvé. Exemple GENOME_FILE :", r.get("GENOME_FILE"))
    print("Contenu :", glob.glob(f"{MOTUS}/genomes/{r.SUSHISTUDY}/*")[:5])
    sys.exit("Indique-moi la structure du dossier genomes/")
pat = max(set(found), key=found.count)
print(f"Motif de chemin retenu : #{pat} ({found.count(pat)}/{len(probe)} génomes de test trouvés)")

def write_list(df, name):
    paths = [candidates(r)[pat] for _, r in df.iterrows()]
    ok = [p for p in paths if os.path.exists(p)]
    with open(f"{OUT}/{name}_list.txt", "w") as f:
        f.write("\n".join(ok) + "\n")
    print(f"{name} : {len(ok):,}/{len(paths):,} fichiers trouvés")
    return f"{OUT}/{name}_list.txt"

ref_lists = {"ARCTIC": write_list(arc, "ref_arctic"), "OMDB": write_list(reps, "ref_omdb")}

# ------------------------------------------------------------------ 3. requêtes : représentants des SGB
drep = f"{WD}/drep/dRep_out/dereplicated_genomes"
q = sorted(glob.glob(f"{drep}/*.fa"))
with open(f"{OUT}/query_list.txt", "w") as f:
    f.write("\n".join(os.path.realpath(p) for p in q) + "\n")
print(f"Représentants de SGB : {len(q)}")

# ------------------------------------------------------------------ 4. skani
best = {}
for ref, lst in ref_lists.items():
    out = f"{OUT}/skani_{ref}.tsv"
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        cmd = ["skani", "dist", "--ql", f"{OUT}/query_list.txt", "--rl", lst, "-t", THREADS, "-o", out]
        if not subprocess.run(["which", "skani"], capture_output=True).returncode == 0:
            cmd = ["conda", "run", "-n", "skani_env"] + cmd
        print(f"[{ref}] skani dist ...", flush=True)
        subprocess.run(cmd, check=True)
    s = pd.read_csv(out, sep="\t")
    s["q"] = s.Query_file.map(lambda p: os.path.basename(p).replace(".fa", ""))
    s["AF"] = s[["Align_fraction_query", "Align_fraction_ref"]].max(axis=1)
    s = s[s.AF >= AF_MIN].sort_values("ANI", ascending=False).drop_duplicates("q")
    best[ref] = s.set_index("q")[["ANI", "AF", "Ref_name"]]

# ------------------------------------------------------------------ 5. synthèse
m2s = pd.read_csv(f"{WD}/drep/mag2sgb.tsv", sep="\t", names=["genome", "sample", "sgb"])
rep = pd.read_csv(f"{WD}/drep/sgb_representatives.tsv", sep="\t", names=["sgb", "rep"])
tax = pd.read_csv(f"{WD}/phylo/mag_taxonomy.tsv", sep="\t")
cat = m2s.groupby("sgb")["sample"].agg(lambda s: "shared" if s.nunique() > 1 else
                                       ("excl_" + s.iloc[0].split("_")[0])).rename("category")
t = rep.merge(cat, on="sgb").merge(tax[["genome", "phylum", "genus", "species", "novel_species"]],
                                   left_on="rep", right_on="genome", how="left")
for ref in best:
    b = best[ref]
    t[f"ANI_{ref}"] = t.rep.map(b.ANI)
    t[f"novel_vs_{ref}"] = ~(t[f"ANI_{ref}"] >= ANI_MIN)
t["novel_vs_both"] = t.novel_vs_ARCTIC & t.novel_vs_OMDB
t.to_csv(f"{OUT}/sgb_novelty.tsv", sep="\t", index=False, float_format="%.2f")

n = len(t)
print(f"\nSGB : {n}  (seuils : ANI >= {ANI_MIN} %, AF >= {AF_MIN} %)")
for c in ["novel_vs_ARCTIC", "novel_vs_OMDB", "novel_vs_both"]:
    print(f"  {c:16s}: {t[c].sum():4d} ({100 * t[c].mean():.1f} %)")
print("\nPar catégorie :")
print(t.groupby("category")[["novel_vs_ARCTIC", "novel_vs_OMDB", "novel_vs_both"]]
      .agg(["sum", "mean"]).round(3).to_string())
print("\nCroisement avec GTDB (SGB sans espèce GTDB) :")
print(pd.crosstab(t.novel_species, t.novel_vs_OMDB, margins=True).to_string())
print(f"\nTableau : {OUT}/sgb_novelty.tsv")
