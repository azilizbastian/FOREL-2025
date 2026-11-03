#!/usr/bin/env python3
"""
Summarise BGC results for the 449 MAGs:
- number of BGCs, per antiSMASH product / BiG-SCAPE class, per MAG
- number of GCFs (BiG-SCAPE, cutoff 0.3)
- novelty vs BiG-FAM (BiG-SLiCE, novel if best distance > 900)
Outputs: bgc_table.tsv (one row per BGC) + summary printed to screen.
"""
import glob, os, sqlite3, sys
import pandas as pd
from Bio import SeqIO

OUT = "/mnt/gr-paoli/nas1/users/students/bastian/BGC_449"   # = $OUT du script bash
THRESHOLD = 900

# ---------- 1. BGCs antiSMASH ----------
rows = []
for gbk in glob.glob(f"{OUT}/regions_gbk/*.region*.gbk"):
    rec = next(SeqIO.parse(gbk, "genbank"))
    region = [f for f in rec.features if f.type == "region"][0]
    q = region.qualifiers
    rows.append({
        "bgc": os.path.basename(gbk).replace(".gbk", ""),
        "mag": rec.id.split("__")[0],
        "products": ";".join(q.get("product", [])),
        "contig_edge": q.get("contig_edge", ["NA"])[0],
        "length_bp": len(rec.seq),
    })
df = pd.DataFrame(rows)

# ---------- 2. GCFs BiG-SCAPE 2 ----------
# Le format de sortie peut varier selon la version : on cherche les fichiers
# *_clustering_c0.30.tsv (colonnes 'Record'/'BGC' et 'GCF number'/'Family')
cl = glob.glob(f"{OUT}/bigscape/output_files/*/mix/*clustering_c0.30*.tsv") \
   + glob.glob(f"{OUT}/bigscape/**/mix*clustering*0.3*.tsv", recursive=True)
if cl:
    c = pd.read_csv(cl[0], sep="\t")
    bgc_col = [x for x in c.columns if x.lower() in ("record", "bgc", "bgc name", "#bgc name")][0]
    gcf_col = [x for x in c.columns if "family" in x.lower() or "gcf" in x.lower()][0]
    c["bgc"] = c[bgc_col].astype(str).str.replace(".gbk", "", regex=False)
    df = df.merge(c[["bgc", gcf_col]].rename(columns={gcf_col: "gcf"}), on="bgc", how="left")
else:
    print("! Fichier de clustering BiG-SCAPE introuvable, vérifier output_files/", file=sys.stderr)

# ---------- 3. Nouveauté BiG-SLiCE ----------
report = open(f"{OUT}/logs/bigslice_report_dir.txt").read().strip()
con = sqlite3.connect(os.path.join(report, "data.db"))
q = """SELECT b.orig_filename AS file, m.membership_value AS dist
       FROM gcf_membership m JOIN bgc b ON b.id = m.bgc_id
       WHERE m.rank = 0"""
s = pd.read_sql(q, con)
s["bgc"] = s["file"].apply(lambda x: os.path.basename(x).replace(".gbk", ""))
df = df.merge(s[["bgc", "dist"]], on="bgc", how="left")
df["novel"] = df["dist"] > THRESHOLD

df.to_csv(f"{OUT}/bgc_table.tsv", sep="\t", index=False)

# ---------- Résumé ----------
n = len(df)
print(f"BGCs : {n} dans {df['mag'].nunique()} MAGs (sur 449)")
print(f"BGCs en bord de contig : {(df['contig_edge']=='True').sum()} ({(df['contig_edge']=='True').mean():.1%})")
if "gcf" in df:
    print(f"GCFs : {df['gcf'].nunique()} (dont singletons : {(df['gcf'].value_counts()==1).sum()})")
print(f"Nouveaux vs BiG-FAM (d > {THRESHOLD}) : {df['novel'].sum()} ({df['novel'].mean():.1%})")
print(f"Sans distance BiG-SLiCE : {df['dist'].isna().sum()}")
print("\nProduits antiSMASH :")
print(df["products"].str.split(";").explode().value_counts().to_string())
