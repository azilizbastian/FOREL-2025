#!/usr/bin/env python3
"""
Step 2 - Download Pfam-A and KOfam, extract only the profiles listed in
config/categories.tsv, and check that every ID matches the expected gene/domain
name. Read the printed check table before going further.

Outputs (in --dbdir):
  selected_pfam.hmm, selected_kofam.hmm, ko_thresholds.tsv, config_check.tsv

If the lab already has Pfam-A.hmm / KOfam locally, pass --pfam-hmm, --kofam-tar
and --ko-list to skip the downloads.

Example:
  python scripts/02_prepare_hmms.py --categories config/categories.tsv --dbdir db
"""
import argparse
import gzip
import os
import shutil
import tarfile
import urllib.request

import pandas as pd
import pyhmmer

PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
KOFAM_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
KOLIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"


def fetch(url, dest):
    if os.path.exists(dest):
        return dest
    print(f"[hmms] downloading {url}")
    with urllib.request.urlopen(url) as r, open(dest + ".part", "wb") as out:
        shutil.copyfileobj(r, out, length=1 << 20)
    os.rename(dest + ".part", dest)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--categories", required=True)
    ap.add_argument("--dbdir", required=True)
    ap.add_argument("--pfam-hmm", help="existing Pfam-A.hmm(.gz)")
    ap.add_argument("--kofam-tar", help="existing KOfam profiles.tar.gz")
    ap.add_argument("--kofam-profiles", help="existing extracted KOfam profiles/ directory (K*.hmm)")
    ap.add_argument("--ko-list", help="existing ko_list(.gz)")
    a = ap.parse_args()
    os.makedirs(a.dbdir, exist_ok=True)
    cat = pd.read_csv(a.categories, sep="\t", dtype=str)
    checks = []

    # ---------------- Pfam ----------------
    pf = cat[cat.db == "pfam"]
    if len(pf):
        src = a.pfam_hmm or fetch(PFAM_URL, os.path.join(a.dbdir, "Pfam-A.hmm.gz"))
        wanted = set(pf.id)
        found = {}
        with pyhmmer.plan7.HMMFile(src) as hf, open(os.path.join(a.dbdir, "selected_pfam.hmm"), "wb") as out:
            for hmm in hf:
                acc = hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession
                base = acc.split(".")[0]
                if base in wanted:
                    hmm.write(out)
                    name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                    found[base] = name
        for _, r in pf.iterrows():
            name = found.get(r.id, "")
            checks.append(dict(category=r.category, db="pfam", id=r.id, expected=r.expected,
                               found=name or "NOT FOUND",
                               ok=bool(name) and r.expected.lower() == name.lower()))

    # ---------------- KOfam ----------------
    ko = cat[cat.db == "ko"]
    if len(ko):
        kl = a.ko_list or fetch(KOLIST_URL, os.path.join(a.dbdir, "ko_list.gz"))
        kolist = pd.read_csv(kl, sep="\t", dtype=str)
        kolist = kolist[kolist.knum.isin(ko.id)]
        kolist[["knum", "threshold", "score_type", "definition"]].to_csv(
            os.path.join(a.dbdir, "ko_thresholds.tsv"), sep="\t", index=False)
        defs = dict(zip(kolist.knum, kolist.definition))
        wanted = {f"{k}.hmm" for k in ko.id}
        n = 0
        with open(os.path.join(a.dbdir, "selected_kofam.hmm"), "wb") as out:
            if a.kofam_profiles:
                for w in sorted(wanted):
                    p = os.path.join(a.kofam_profiles, w)
                    if os.path.exists(p):
                        with open(p, "rb") as fh:
                            out.write(fh.read())
                        n += 1
            else:
                tar = a.kofam_tar or fetch(KOFAM_URL, os.path.join(a.dbdir, "profiles.tar.gz"))
                with tarfile.open(tar, "r:gz") as t:
                    for m in t:
                        if os.path.basename(m.name) in wanted:
                            out.write(t.extractfile(m).read())
                            n += 1
        print(f"[hmms] {n}/{len(wanted)} KOfam profiles extracted")
        for _, r in ko.iterrows():
            d = defs.get(r.id, "")
            symbols = d.split(";")[0].lower().replace(" ", "").split(",") if d else []
            checks.append(dict(category=r.category, db="ko", id=r.id, expected=r.expected,
                               found=d or "NOT FOUND", ok=r.expected.lower() in symbols))

    chk = pd.DataFrame(checks)
    chk.to_csv(os.path.join(a.dbdir, "config_check.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200, "display.max_colwidth", 80, "display.max_rows", 200)
    print(chk.to_string(index=False))
    bad = chk[~chk.ok]
    if len(bad):
        print(f"\n[hmms] {len(bad)} IDs do not match the expected name: check them and fix "
              "config/categories.tsv before step 3 (a mismatch may only be a naming difference).")
    else:
        print("\n[hmms] all IDs match their expected names.")


if __name__ == "__main__":
    main()
