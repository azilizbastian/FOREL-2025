#!/usr/bin/env python3
import argparse
import os
import statistics as st

import matplotlib.pyplot as plt

# Duncan et al. 2022 (Microbiome 10:67), 21 eukaryotic MAGs (bins < 50% complete discarded).
# The article calls the contiguity value "L50" but reports it in bp (N50-type length).
REF = dict(completeness=67.82, contamination=2.82, size_mb=24.31, n50_kb=5.229)
MARKERS = ["o", "^", "s", "D"]


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
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            c = line.split()
            if len(c) >= 3:
                rows.append((c[0], float(c[1]), float(c[2])))
    return rows


def summarize(label, rows):
    print(f"\n{label}")
    if not rows:
        print("  no bin")
        return
    comp = [r["completeness"] for r in rows]
    cont = [r["contamination"] for r in rows]
    size = [r["size_mb"] for r in rows]
    n50 = [r["n50_kb"] for r in rows]
    print(f"  n = {len(rows)}")
    print(f"  completeness  mean {st.mean(comp):.1f}% (range {min(comp):.1f}-{max(comp):.1f})"
          f"   | Duncan et al.: mean {REF['completeness']}%")
    print(f"  contamination mean {st.mean(cont):.2f}% (range {min(cont):.2f}-{max(cont):.2f})"
          f"   | Duncan et al.: mean {REF['contamination']}%")
    print(f"  size          median {st.median(size):.2f} Mb (range {min(size):.2f}-{max(size):.2f})"
          f"   | Duncan et al.: median {REF['size_mb']} Mb")
    print(f"  N50           median {st.median(n50):,.0f} kb (range {min(n50):,.0f}-{max(n50):,.0f})"
          f"   | Duncan et al.: median {REF['n50_kb']} kb")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs=3, action="append", required=True,
                    metavar=("NAME", "EUKCC_TSV", "BINS_DIR"))
    ap.add_argument("--min-comp", type=float, default=50.0)
    ap.add_argument("--max-cont", type=float, default=10.0)
    ap.add_argument("--out", default="comparison")
    args = ap.parse_args()

    bins = []
    for name, tsv, bdir in args.set:
        for b, comp, cont in read_eukcc(tsv):
            fa = os.path.join(bdir, b)
            if not os.path.exists(fa):
                print(f"WARNING: {fa} not found, bin skipped")
                continue
            L = contig_lengths(fa)
            n50, l50 = n50_l50(L)
            bins.append(dict(dataset=name, bin=b, completeness=comp, contamination=cont,
                             size_mb=sum(L) / 1e6, n50_kb=n50 / 1e3, l50=l50,
                             contigs=len(L), longest_mb=max(L) / 1e6))
    if not bins:
        raise SystemExit("No bin found: check the EukCC table and the bins directory.")

    with open(f"{args.out}_bins.tsv", "w") as fo:
        keys = list(bins[0].keys())
        fo.write("\t".join(keys) + "\n")
        for r in bins:
            fo.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k])
                               for k in keys) + "\n")

    sel = [r for r in bins if r["completeness"] >= args.min_comp
           and r["contamination"] < args.max_cont]
    summarize(f"Bins with completeness >= {args.min_comp:g}% and contamination < "
              f"{args.max_cont:g}% (comparable to Duncan et al., who kept bins >= 50% complete)", sel)
    summarize("All bins with an EukCC result (NOT comparable: includes low-quality bins)", bins)
    print(f"\nWrote {args.out}_bins.tsv")

    datasets = [s[0] for s in args.set]
    mk = {d: MARKERS[i % len(MARKERS)] for i, d in enumerate(datasets)}
    panels = [("completeness", "Completeness (%)", REF["completeness"], "mean", False),
              ("contamination", "Contamination (%)", REF["contamination"], "mean", False),
              ("size_mb", "Bin size (Mb)", REF["size_mb"], "median", False),
              ("n50_kb", "N50 (kb, log scale)", REF["n50_kb"], "median", True)]
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.8))
    for ax, (key, ylabel, ref, kind, logy) in zip(axes, panels):
        for i, r in enumerate(sel):
            jitter = ((i * 37) % 11 - 5) * 0.02
            ax.scatter(jitter, r[key], marker=mk[r["dataset"]], s=50, c="#1b7837",
                       edgecolors="black", linewidths=0.5, zorder=3)
        ax.axhline(ref, color="#b2182b", ls="--", lw=1)
        ax.text(0.5, ref, f" Duncan et al. 2022\n ({kind} {ref:g})", color="#b2182b",
                fontsize=7, va="bottom", transform=ax.get_yaxis_transform())
        if sel:
            vals = [r[key] for r in sel]
            own = st.mean(vals) if kind == "mean" else st.median(vals)
            ax.scatter(0.0, own, marker="_", s=900, c="black", zorder=4)
        ax.set_xlim(-0.5, 1.1)
        ax.set_xticks([])
        ax.set_ylabel(ylabel)
        if logy:
            ax.set_yscale("log")
        ax.spines[["top", "right", "bottom"]].set_visible(False)
    handles = [plt.Line2D([], [], marker=mk[d], ls="", mfc="#1b7837", mec="black", label=d)
               for d in datasets]
    handles.append(plt.Line2D([], [], marker="_", ls="", mec="black", ms=14,
                              label="this study (mean / median)"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               fontsize=8)
    fig.suptitle(f"Bins with ≥{args.min_comp:g}% completeness and <{args.max_cont:g}% "
                 f"contamination (n = {len(sel)}) vs eukaryotic MAGs of Duncan et al. 2022 (n = 21)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{args.out}.{ext}", dpi=300, bbox_inches="tight")
    print(f"Wrote {args.out}.png/.pdf/.svg")


if __name__ == "__main__":
    main()
