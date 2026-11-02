#!/usr/bin/env python
"""Étape 7 : tables supplémentaires des MAGs FOREL (une ligne par génome)

  Table S1 - identité, qualité (CheckM2), assemblage, taxonomie (GTDB-Tk R232)
  Table S2 - identité, fonction (Pfam, eggNOG, KofamScan, nouveauté vs OMDB v2), données

Usage : conda activate funcannot && python ~/07_tables.py
Sortie : ~/forel_func/summary/FOREL_MAG_tables.xlsx (+ un .tsv par table)
Relancer après KofamScan : la colonne KO KofamScan se remplit toute seule.
"""
import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

WD = os.path.expanduser("~/forel_func")
OUT = f"{WD}/summary"
XLSX = f"{OUT}/FOREL_MAG_tables.xlsx"

# ------------------------------------------------------------------ sources
mags = [l.strip() for l in open(f"{WD}/meta/mags_mq.txt") if l.strip()]
hq = {l.strip() for l in open(f"{WD}/meta/mags_hq.txt") if l.strip()}

bins = pd.read_csv(f"{WD}/meta/bins_summary_all.tsv", sep="\t")
bins = bins.rename(columns={"new_bin_name": "Genome", "original_bin_name": "Original bin",
                            "sample": "Sample"})[["Genome", "Sample", "Original bin"]]

q = pd.read_csv(f"{WD}/checkm2/quality_report.tsv", sep="\t").rename(columns={"Name": "Genome"})

gt = pd.read_csv(f"{WD}/meta/gtdbtk_summary.tsv", sep="\t").rename(columns={"user_genome": "Genome"})
nov = pd.read_csv(f"{WD}/meta/mag_novelty.tsv", sep="\t", names=["Genome", "Taxonomic novelty"])

genes = pd.read_csv(f"{OUT}/gene_table.tsv.gz", sep="\t")

# ------------------------------------------------------------------ Table S1
s1 = pd.DataFrame({"Genome": mags}).merge(bins, on="Genome", how="left").merge(q, on="Genome", how="left")
s1["MIMAG quality"] = s1.Genome.map(lambda g: "High" if g in hq else "Medium")

s1 = s1.merge(gt, on="Genome", how="left").merge(nov, on="Genome", how="left")
ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

def split_tax(c):
    if not isinstance(c, str) or not c.startswith("d__"):
        return [None] * 7
    parts = [p.split("__", 1)[1] if "__" in p else p for p in c.split(";")]
    return [p if p else None for p in (parts + [None] * 7)[:7]]

s1[ranks] = pd.DataFrame(s1["classification"].map(split_tax).tolist(), index=s1.index)
s1["Taxonomic novelty"] = s1["Taxonomic novelty"].str.replace("_", " ").str.capitalize().replace({"Species": "Known species"})

def col(df, *names):
    """premier nom de colonne présent (les noms GTDB-Tk changent selon les versions)"""
    for n in names:
        if n in df.columns:
            return df[n]
    return pd.Series([None] * len(df), index=df.index)

s1["Closest reference"] = col(s1, "closest_genome_reference", "fastani_reference")
s1["ANI to closest reference (%)"] = pd.to_numeric(col(s1, "closest_genome_ani", "fastani_ani"), errors="coerce")
s1["AF to closest reference"] = pd.to_numeric(col(s1, "closest_genome_af", "fastani_af"), errors="coerce")
s1["Classification method"] = col(s1, "classification_method")

# Exclusion : génomes sans classification GTDB-Tk, ou classés au domaine seulement
# (probables eucaryotes / bins non procaryotes passés au travers du filtre whokaryote)
excl_mask = s1["Phylum"].isna()
EXCLUDED = set(s1.loc[excl_mask, "Genome"])
excl_tab = s1.loc[excl_mask, ["Genome", "Sample", "Completeness", "Contamination", "classification"]]
excl_tab.to_csv(f"{OUT}/excluded_unclassified_MAGs.tsv", sep="\t", index=False)
print(f"Exclus de la Table S1 (taxonomie non définie ou domaine seul) : {len(EXCLUDED)}")
if len(EXCLUDED):
    print(excl_tab.to_string(index=False))
s1 = s1.loc[~excl_mask].reset_index(drop=True)

S1_COLS = {
    "Genome": "Genome identifier",
    "Sample": "Sample of origin",
    "Original bin": "Bin name in the SemiBin2 output",
    "Completeness": "CheckM2 completeness (%)",
    "Contamination": "CheckM2 contamination (%)",
    "MIMAG quality": "High: ≥90% completeness and <5% contamination; Medium: ≥50% and <10%",
    "Completeness_Model_Used": "CheckM2 completeness model",
    "Genome_Size": "Genome size (bp)",
    "Total_Contigs": "Number of contigs",
    "Contig_N50": "Contig N50 (bp)",
    "Max_Contig_Length": "Longest contig (bp)",
    "GC_Content": "GC content (fraction)",
    "Coding_Density": "Coding density (fraction, CheckM2)",
    "classification": "Full GTDB-Tk classification (R232)",
    **{r: f"GTDB {r.lower()}" for r in ranks},
    "Taxonomic novelty": "Lowest GTDB rank without an assigned taxon",
    "Closest reference": "Closest GTDB reference genome",
    "ANI to closest reference (%)": "Average nucleotide identity to the closest reference (skani)",
    "AF to closest reference": "Alignment fraction to the closest reference",
    "Classification method": "GTDB-Tk classification method",
}
s1 = s1[list(S1_COLS)].rename(columns={
    "Completeness": "Completeness (%)", "Contamination": "Contamination (%)",
    "Completeness_Model_Used": "CheckM2 model", "Genome_Size": "Genome size (bp)",
    "Total_Contigs": "Contigs", "Contig_N50": "N50 (bp)", "Max_Contig_Length": "Longest contig (bp)",
    "GC_Content": "GC", "Coding_Density": "Coding density", "classification": "GTDB classification"})
S1_DESC = dict(zip(s1.columns, S1_COLS.values()))

# ------------------------------------------------------------------ Table S2
g = genes.copy()
kofam_done = os.path.exists(f"{WD}/annotation/kofam_best.tsv") and "kofam" in g.columns and g["kofam"].any()
fr = lambda s: s.mean()
agg = g.groupby("mag").agg(
    n_genes=("gene", "size"),
    n_fam=("nr30_rep", "nunique"),
    f_function=("annot_status", lambda s: (s == "Function assigned").mean()),
    f_unknown=("annot_status", lambda s: (s == "Homolog, unknown function").mean()),
    f_nohit=("annot_status", lambda s: (s == "No hit").mean()),
    f_pfam=("pfam_any", fr),
    f_eggnog=("eggnog_any", fr),
    f_ko_kofam=("kofam", fr),
    f_ko_eggnog=("eggnog_ko", fr),
    f_novel=("omdb_class", lambda s: (s == "novel_vs_OMDB").mean()),
    f_30_50=("omdb_class", lambda s: (s == "distant_30-50").mean()),
    f_50_90=("omdb_class", lambda s: (s == "related_50-90").mean()),
    f_90=("omdb_class", lambda s: (s == "near_identical_90+").mean()),
).reset_index().rename(columns={"mag": "Genome"})
novfam = (g[g.omdb_class == "novel_vs_OMDB"].groupby("mag").nr30_rep.nunique()
          .rename("n_novel_fam").reset_index().rename(columns={"mag": "Genome"}))
agg = agg.merge(novfam, on="Genome", how="left")
agg["n_novel_fam"] = agg.n_novel_fam.fillna(0).astype(int)
if not kofam_done:
    agg["f_ko_kofam"] = None

s2 = pd.DataFrame({"Genome": [m for m in mags if m not in EXCLUDED]}) \
    .merge(bins[["Genome", "Sample"]], on="Genome", how="left").merge(agg, on="Genome", how="left")
s2["Genome file"] = s2.Genome + ".fa"
s2["Protein file"] = s2.Genome + ".faa"
s2["Accession"] = None

S2_COLS = {
    "Genome": ("Genome", "Genome identifier"),
    "Sample": ("Sample", "Sample of origin"),
    "n_genes": ("Protein-coding genes", "Genes predicted by Prodigal v2.6.3 (-c -q -m -p single)"),
    "n_fam": ("Protein families (NR30)", "Distinct AA NR30 clusters the genome's genes belong to"),
    "f_function": ("Function assigned", "Genes with a KO, a non-DUF Pfam domain or an informative eggNOG annotation (COG ≠ S)"),
    "f_unknown": ("Homolog, unknown function", "Genes with only uninformative hits (DUF, COG S, no description)"),
    "f_nohit": ("No hit", "Genes without a hit in Pfam, eggNOG or KOfam"),
    "f_pfam": ("Pfam", "Genes with ≥1 Pfam v38.2 domain (pyHMMER, gathering thresholds)"),
    "f_eggnog": ("eggNOG", "Genes with an eggNOG v5.0.2 ortholog (eggNOG-mapper v2.1.13)"),
    "f_ko_kofam": ("KO (KofamScan)", "Genes with a KO above the KofamScan v1.3.0 adaptive threshold" + ("" if kofam_done else " - pending")),
    "f_ko_eggnog": ("KO (eggNOG)", "Genes with a KO transferred by eggNOG-mapper"),
    "f_novel": ("Absent from OMDB v2", "Genes whose NR30 family has no OMDB v2 NR30 hit (MMseqs2, e ≤ 1e-3, coverage ≥ 80%)"),
    "f_30_50": ("OMDB v2 30-50% id.", "Genes whose family's best OMDB v2 hit has 30-50% identity"),
    "f_50_90": ("OMDB v2 50-90% id.", "Genes whose family's best OMDB v2 hit has 50-90% identity"),
    "f_90": ("OMDB v2 ≥90% id.", "Genes whose family's best OMDB v2 hit has ≥90% identity"),
    "n_novel_fam": ("Families absent from OMDB v2", "Number of NR30 families absent from OMDB v2"),
    "Genome file": ("Genome file", "Nucleotide FASTA file name"),
    "Protein file": ("Protein file", "Predicted protein FASTA file name"),
    "Accession": ("Accession", "Public repository accession - to fill in at submission"),
}
s2 = s2[list(S2_COLS)].rename(columns={k: v[0] for k, v in S2_COLS.items()})
S2_DESC = {v[0]: v[1] for v in S2_COLS.values()}
PCT_S2 = [S2_COLS[k][0] for k in S2_COLS if k.startswith("f_")]

# ------------------------------------------------------------------ export
s1.to_csv(f"{OUT}/Table_S1_MAG_quality_taxonomy.tsv", sep="\t", index=False)
s2.to_csv(f"{OUT}/Table_S2_MAG_function.tsv", sep="\t", index=False)

readme = pd.DataFrame(
    [("Table S1", c, d) for c, d in S1_DESC.items()] +
    [("Table S2", c, d + (" (fraction of the genome's genes)" if c in PCT_S2 else "")) for c, d in S2_DESC.items()],
    columns=["Sheet", "Column", "Description"])

with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
    readme.to_excel(xw, sheet_name="README", index=False, startrow=3)
    s1.to_excel(xw, sheet_name="Table S1", index=False)
    s2.to_excel(xw, sheet_name="Table S2", index=False)

wb = load_workbook(XLSX)
ARIAL, BOLD = Font(name="Arial", size=10), Font(name="Arial", size=10, bold=True)
HEAD = PatternFill("solid", fgColor="D9E2F3")
FILL = PatternFill("solid", fgColor="FFF2CC")
fmt1 = {"Completeness (%)": "0.00", "Contamination (%)": "0.00", "GC": "0.000",
        "Coding density": "0.000", "ANI to closest reference (%)": "0.00",
        "AF to closest reference": "0.000", "Genome size (bp)": "#,##0", "N50 (bp)": "#,##0",
        "Longest contig (bp)": "#,##0", "Protein-coding genes": "#,##0",
        "Protein families (NR30)": "#,##0", "Families absent from OMDB v2": "#,##0",
        **{c: "0.0%" for c in PCT_S2}}

for ws in wb.worksheets:
    hrow = 4 if ws.title == "README" else 1
    for row in ws.iter_rows():
        for c in row:
            c.font = ARIAL
    for c in ws[hrow]:
        c.font, c.fill = BOLD, HEAD
        c.alignment = Alignment(wrap_text=True, vertical="center")
    for i, c in enumerate(ws[hrow], 1):
        f = fmt1.get(c.value)
        width = 14
        if c.value in ("GTDB classification",):
            width = 60
        elif c.value in ("Genome", "Description"):
            width = 34 if c.value == "Genome" else 95
        elif c.value in ("Closest reference", "Species", "Original bin", "Sample"):
            width = 22
        ws.column_dimensions[get_column_letter(i)].width = width
        if f:
            for r in range(hrow + 1, ws.max_row + 1):
                ws.cell(r, i).number_format = f
        if c.value == "Accession":
            for r in range(hrow + 1, ws.max_row + 1):
                ws.cell(r, i).fill = FILL
    ws.freeze_panes = ws.cell(hrow + 1, 2 if ws.title != "README" else 1)
    if ws.title != "README":
        ws.auto_filter.ref = ws.dimensions

rd = wb["README"]
rd["A1"] = "FOREL 2025 - prokaryotic MAGs (CheckM2 ≥50% completeness, <10% contamination)"
rd["A1"].font = Font(name="Arial", size=12, bold=True)
rd["A2"] = ("Table S1: genome quality, assembly and taxonomy. Table S2: functional annotation and novelty. "
            "Genomes without a GTDB-Tk classification beyond domain level were excluded. "
            "Percent columns in Table S2 are fractions of each genome's protein-coding genes. "
            "Yellow cells (Accession) are to be filled in at data submission.")
rd["A2"].font = ARIAL

wb.save(XLSX)

print(f"Table S1 : {s1.shape[0]} génomes x {s1.shape[1]} colonnes")
print(f"Table S2 : {s2.shape[0]} génomes x {s2.shape[1]} colonnes")
print("KofamScan :", "inclus" if kofam_done else "en attente (colonne vide)")
print(f"Écrit : {XLSX}")
