#!/usr/bin/env python3
"""
Step 10 - Scan a (very large) protein FASTA for ice-binding proteins (Pfam PF11999,
gathering threshold) and for the control gene recA (KOfam K03553, KOfam threshold),
streaming the file in blocks so that nothing is decompressed to disk.

Works for the OMDB v2.0 protein catalogue (header: '>OMDBv2.0_AA_G_R_... rep=<GENOME>-<scaffold>_<n>;size=1')
and for our MAGs (header: '>ctgXXXX|<protein_id>').

Output: one line per hit (genome, protein, profile, score).

Examples:
  # OMDB (use pigz if available: much faster decompression)
  pigz -dc $OMDB/catalog/clustering/OMDBv2.0_AA_G_R/OMDBv2.0_AA_G_R.faa.gz | \
    python scripts/10_hmm_scan_proteins.py --faa - --genome-from rep \
      --pfam db/selected_pfam.hmm --recA db/K03553.hmm --recA-threshold <from ko_list> \
      --out results/omdb/omdb_hits.tsv --cpus 32
  # our MAGs
  python scripts/10_hmm_scan_proteins.py --faa results/omdb/mags_all.faa --genome-from pipe ...
"""
import argparse
import re
import sys
import time

import pyhmmer

REP = re.compile(rb"rep=([^;\s]+)")


def genome_of(name, desc, mode):
    if mode == "pipe":
        return name.split(b"|", 1)[0].decode()
    m = REP.search(desc or b"")
    gene = m.group(1) if m else name
    return gene.rsplit(b"-", 1)[0].decode()   # '<GENOME>-scaffold_1_1' -> '<GENOME>'


def b(x):
    return x if isinstance(x, (bytes, bytearray)) else (x or "").encode()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", required=True, help="protein FASTA (plain); '-' = stdin")
    ap.add_argument("--genome-from", choices=["rep", "pipe"], required=True)
    ap.add_argument("--pfam", required=True, help="HMM file containing Ice_binding (PF11999)")
    ap.add_argument("--recA", required=True, help="KOfam K03553 HMM")
    ap.add_argument("--recA-threshold", type=float, required=True)
    ap.add_argument("--recA-score-type", choices=["full", "domain"], default="full")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=16)
    ap.add_argument("--block", type=int, default=500_000, help="sequences per block")
    a = ap.parse_args()

    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(a.pfam) as hf:
        ibp = [h for h in hf if b(h.name) == b"Ice_binding"]
    with pyhmmer.plan7.HMMFile(a.recA) as hf:
        reca = list(hf)
    if len(ibp) != 1 or len(reca) != 1:
        sys.exit(f"expected 1 Ice_binding and 1 recA HMM, got {len(ibp)} and {len(reca)}")

    src = sys.stdin.buffer if a.faa == "-" else a.faa
    n_seq, t0 = 0, time.time()
    seqfile = pyhmmer.easel.SequenceFile(src, digital=True, alphabet=alpha, format="fasta")
    with open(a.out, "w") as out:
        out.write("genome\tprotein\tprofile\tscore\tbest_domain_score\n")
        while True:
            block = seqfile.read_block(sequences=a.block)
            if not block:
                break
            n_seq += len(block)
            for hits in pyhmmer.hmmsearch(ibp, block, cpus=a.cpus, bit_cutoffs="gathering"):
                for h in hits:
                    if h.included:
                        out.write(f"{genome_of(b(h.name), b(h.description), a.genome_from)}\t{b(h.name).decode()}\t"
                                  f"Ice_binding\t{h.score:.1f}\t{h.best_domain.score:.1f}\n")
            for hits in pyhmmer.hmmsearch(reca, block, cpus=a.cpus, E=1e-5):
                for h in hits:
                    sc = h.score if a.recA_score_type == "full" else h.best_domain.score
                    if sc >= a.recA_threshold:
                        out.write(f"{genome_of(b(h.name), b(h.description), a.genome_from)}\t{b(h.name).decode()}\t"
                                  f"recA\t{h.score:.1f}\t{h.best_domain.score:.1f}\n")
            out.flush()
            print(f"[scan] {n_seq:,} proteins, {(time.time() - t0) / 60:.1f} min", file=sys.stderr, flush=True)
    seqfile.close()
    print(f"[scan] done: {n_seq:,} proteins -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
