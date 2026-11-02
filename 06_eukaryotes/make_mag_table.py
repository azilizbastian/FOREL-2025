#!/usr/bin/env python3
import argparse
import os

# NCBI taxids found in the EukCC lineages of the Lucas bins (names checked in the
# EukCC terminal_name table; 211026 is still to be verified).
TAXID_NAMES = {
    "38832": "Micromonas",
    "41873": "Mamiellaceae",
    "1525212": "Bathycoccaceae",
    "1735725": "Stramenopiles MAST-4",
    "704141": "Dinobryaceae",
    "423054": "Eimeriorina",
    "211026": "Ciliophora? (to verify)",
}


def contig_lengths(path):
    lengths, cur = [], 0
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    return lengths


def n50_l50(lengths):
    total, acc = sum(lengths), 0
    for rank, L in enumerate(sorted(lengths, reverse=True), 1):
        acc += L
        if acc >= total / 2:
            return L, rank


def read_eukcc(path):
    with open(path) as fh:
        head = next(fh).rstrip("\n").split("\t")
        rows = []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                d = dict(zip(head, c))
                rows.append(d)
    return rows


def mimag(comp, cont):
    if comp > 90 and cont < 5:
        return "high*"      # needs rRNA/tRNA evidence for a formal MIMAG 'high'
    if comp >= 50 and cont < 10:
        return "medium"
    return "low"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs=3, action="append", required=True,
                    metavar=("NAME", "EUKCC_TSV", "BINS_DIR"))
    ap.add_argument("--taxa")
    ap.add_argument("--out", default="mag_table")
    args = ap.parse_args()

    taxa = {}
    if args.taxa:
        with open(args.taxa) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 3 and c[0] != "dataset":
                    taxa[(c[0], c[1])] = c[2]

    cols = ["dataset", "bin", "taxon", "completeness_%", "contamination_%", "MIMAG",
            "size_Mb", "contigs", "N50_kb", "L50", "longest_contig_Mb"]
    rows = []
    for name, tsv, bdir in args.set:
        for d in read_eukcc(tsv):
            b = d["bin"]
            fa = os.path.join(bdir, b)
            if not os.path.exists(fa):
                print(f"WARNING: {fa} not found, bin skipped")
                continue
            L = contig_lengths(fa)
            n50, l50 = n50_l50(L)
            comp, cont = float(d["completeness"]), float(d["contamination"])
            lin = d.get("ncbi_lng", "")
            tid = lin.split("-")[-1] if lin else ""
            taxon = taxa.get((name, b)) or TAXID_NAMES.get(tid) or (f"taxid {tid}" if tid else "NA")
            rows.append([name, b.replace(".fa", ""), taxon, round(comp, 1), round(cont, 1),
                         mimag(comp, cont), round(sum(L) / 1e6, 2), len(L),
                         round(n50 / 1e3), l50, round(max(L) / 1e6, 2)])
    if not rows:
        raise SystemExit("No bin found: check the EukCC table and the bins directory.")
    rows.sort(key=lambda r: (r[0], -r[3]))

    with open(f"{args.out}.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    print(f"Wrote {args.out}.tsv ({len(rows)} bins)")

    try:
        from openpyxl import Workbook
        wb = Workbook()
        ws = wb.active
        ws.title = "MAGs"
        ws.append(cols)
        for r in rows:
            ws.append(r)
        wb.save(f"{args.out}.xlsx")
        print(f"Wrote {args.out}.xlsx")
    except ImportError:
        print("openpyxl not installed: only the .tsv was written (open it in Excel).")

    n_med = sum(1 for r in rows if r[5] in ("medium", "high*"))
    print(f"\n{len(rows)} bins; {n_med} with completeness >= 50% and contamination < 10%")
    print("* 'high' = completeness > 90% and contamination < 5% only; rRNA/tRNA not assessed.")


if __name__ == "__main__":
    main()
