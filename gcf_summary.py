#!/usr/bin/env python3
"""Ajoute classe BiG-SCAPE et GCF (cutoff 0.3) à bgc_table.tsv et résume."""
import glob, os
import pandas as pd

OUT = "/mnt/gr-paoli/nas1/users/students/bastian/BGC_449"
files = glob.glob(f"{OUT}/bigscape/output_files/*_c0.3/*/*_clustering_c0.3.tsv")
parts = []
for f in files:
    cls = os.path.basename(os.path.dirname(f))
    c = pd.read_csv(f, sep="\t")
    bcol = "GBK" if "GBK" in c.columns else c.columns[0]
    fcol = [x for x in c.columns if "family" in x.lower() or "gcf" in x.lower()][0]
    parts.append(pd.DataFrame({
        "bgc": c[bcol].astype(str).str.replace(".gbk", "", regex=False),
        "bigscape_class": cls,
        "gcf": cls + "_" + c[fcol].astype(str)}))
    print(f"{cls}: colonnes {list(c.columns)}") if f == files[0] else None
cl = pd.concat(parts).drop_duplicates("bgc")

df = pd.read_csv(f"{OUT}/bgc_table.tsv", sep="\t")
df = df.drop(columns=[x for x in ("gcf", "bigscape_class") if x in df.columns]).merge(cl, on="bgc", how="left")
df.to_csv(f"{OUT}/bgc_table.tsv", sep="\t", index=False)

print(f"\nBGCs appariés à un GCF : {df['gcf'].notna().sum()} / {len(df)}")
gs = df["gcf"].value_counts()
print(f"GCFs : {df['gcf'].nunique()} (singletons : {(gs == 1).sum()})")
gnov = df.groupby("gcf")["novel"].all()
print(f"GCFs entièrement nouveaux vs BiG-FAM : {gnov.sum()} ({gnov.mean():.1%})")
print("\nPar classe BiG-SCAPE :")
t = df.groupby("bigscape_class").agg(BGCs=("bgc", "size"), GCFs=("gcf", "nunique"),
                                     novel_pct=("novel", "mean"))
t["novel_pct"] = (100 * t["novel_pct"]).round(1)
print(t.sort_values("BGCs", ascending=False).to_string())
