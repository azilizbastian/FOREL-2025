#!/usr/bin/env python3
"""Fusionne CheckM2 + Barrnap + tRNAscan-SE et attribue la classe MIMAG.

Usage :
  python mimag_merge.py <quality_report.tsv> <dossier_sortie_run> <sortie.tsv> [--full-rrna-only]

MIMAG (Bowers et al. 2017) :
  High   : compl. > 90 %, contam. < 5 %, 23S + 16S + 5S présents, >= 18 tRNA (acides aminés distincts)
  Medium : compl. >= 50 %, contam. < 10 %
  Low    : compl. < 50 %, contam. < 10 %
Les bins >90/<5 sans les rRNA/tRNA requis tombent en Medium.
--full-rrna-only : ne compter que les rRNA non marqués "(partial)" par Barrnap.
"""
import sys
from pathlib import Path
import pandas as pd

STANDARD_AA = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
               "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
ALIASES = {"iMet": "Met", "fMet": "Met", "Ile2": "Ile"}


def parse_barrnap(gff, full_only):
    counts = {"5S": 0, "16S": 0, "23S": 0}
    partial = {"5S": 0, "16S": 0, "23S": 0}
    if not gff.exists():
        return counts, partial
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        attrs = line.split("\t")[8]
        for r in counts:
            if f"Name={r}_rRNA" in attrs:
                if "partial" in attrs:
                    partial[r] += 1
                    if full_only:
                        continue
                counts[r] += 1
    return counts, partial


def parse_trnascan(tsv):
    aa, n_total = set(), 0
    if not tsv.exists():
        return aa, n_total
    for line in tsv.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 9 or cols[0].startswith(("Sequence", "Name", "-----")):
            continue
        note = cols[9].lower() if len(cols) > 9 else ""
        if "pseudo" in note:
            continue
        t = ALIASES.get(cols[4].strip(), cols[4].strip())
        if t in STANDARD_AA:
            aa.add(t)
            n_total += 1
    return aa, n_total


def mimag_class(row):
    c, k = row["Completeness"], row["Contamination"]
    has_rrna = row["n_16S"] > 0 and row["n_23S"] > 0 and row["n_5S"] > 0
    if c > 90 and k < 5 and has_rrna and row["n_tRNA_aa"] >= 18:
        return "High"
    if c >= 50 and k < 10:
        return "Medium"
    if k < 10:
        return "Low"
    return "Fail (contam >= 10)"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    full_only = "--full-rrna-only" in sys.argv
    checkm2, run_dir, out = Path(args[0]), Path(args[1]), Path(args[2])

    q = pd.read_csv(checkm2, sep="\t")[["Name", "Completeness", "Contamination"]]
    rows = []
    for name in q["Name"]:
        counts, partial = parse_barrnap(run_dir / "barrnap" / f"{name}.gff", full_only)
        aa, n_trna = parse_trnascan(run_dir / "trnascan" / f"{name}.tsv")
        done = (run_dir / "done" / name).exists()
        rows.append({"Name": name, "processed": done,
                     "n_16S": counts["16S"], "n_23S": counts["23S"], "n_5S": counts["5S"],
                     "partial_16S": partial["16S"], "partial_23S": partial["23S"],
                     "partial_5S": partial["5S"],
                     "n_tRNA_total": n_trna, "n_tRNA_aa": len(aa),
                     "missing_aa": ",".join(sorted(STANDARD_AA - aa))})
    df = q.merge(pd.DataFrame(rows), on="Name")
    df["CheckM2_HQ"] = (df["Completeness"] >= 90) & (df["Contamination"] < 5)
    df["MIMAG"] = df.apply(mimag_class, axis=1)
    df.to_csv(out, sep="\t", index=False)

    missing = (~df["processed"]).sum()
    if missing:
        print(f"ATTENTION : {missing} bins du rapport CheckM2 sans sortie Barrnap/tRNAscan.")
    print(df["MIMAG"].value_counts().to_string())
    hq = df[df["CheckM2_HQ"]]
    print(f"\nHQ CheckM2 seul : {len(hq)} | HQ MIMAG complet : {(df['MIMAG'] == 'High').sum()}")
    lost = hq[hq["MIMAG"] != "High"]
    if len(lost):
        print("Rétrogradés (HQ CheckM2 -> Medium MIMAG) :")
        for _, r in lost.iterrows():
            why = []
            for x in ("16S", "23S", "5S"):
                if r[f"n_{x}"] == 0:
                    why.append(f"no {x}")
            if r["n_tRNA_aa"] < 18:
                why.append(f"{r['n_tRNA_aa']} aa tRNA")
            print(f"  {r['Name']}: {', '.join(why)}")


if __name__ == "__main__":
    main()
