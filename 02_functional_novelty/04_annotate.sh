#!/usr/bin/env bash
# Étape 4 : annotation fonctionnelle des protéines (approche OMDB v2, KofamScan à la place de KEGG)
# Un outil par session, en parallèle :
#   tmux new -d -s pfam   'bash ~/04_annotate.sh pfam   2>&1 | tee ~/04_pfam.log'
#   tmux new -d -s eggnog 'bash ~/04_annotate.sh eggnog 2>&1 | tee ~/04_eggnog.log'
#   tmux new -d -s kofam  'bash ~/04_annotate.sh kofam  2>&1 | tee ~/04_kofam.log'
set -euo pipefail

TOOL=${1:?usage: 04_annotate.sh pfam|eggnog|kofam}
WD=/home/bastian/forel_func
DBDIR=/home/bastian/databases
FAA=$WD/catalog/all_proteins.faa
OUT=$WD/annotation
THREADS=${THREADS:-$(( $(nproc) / 3 ))}
mkdir -p "$OUT"
source "$(conda info --base)/etc/profile.d/conda.sh"
echo "== $TOOL | $(grep -c '^>' "$FAA") protéines | $THREADS threads | $(date)"

case "$TOOL" in

pfam)
  conda activate funcannot
  cat "$DBDIR/pfam/Pfam.version"
  python - "$FAA" "$DBDIR/pfam/Pfam-A.hmm" "$OUT/pfam_hits.tsv" "$THREADS" <<'EOF'
import sys, pyhmmer
faa, hmm, out, cpus = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
txt = lambda x: x.decode() if isinstance(x, bytes) else (x or "")
alph = pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alph) as sf:
    seqs = sf.read_block()
n = 0
with pyhmmer.plan7.HMMFile(hmm) as hf, open(out, "w") as fo:
    fo.write("gene\tpfam_acc\tpfam_name\tevalue\tbitscore\tali_from\tali_to\n")
    # seuils de gathering (GA) de chaque famille, comme OMDB v2
    for hits in pyhmmer.hmmsearch(hf, seqs, cpus=cpus, bit_cutoffs="gathering"):
        acc, name = txt(hits.query.accession), txt(hits.query.name)
        for hit in hits.included:
            for d in hit.domains.included:
                fo.write(f"{txt(hit.name)}\t{acc}\t{name}\t{d.i_evalue:.2e}\t{d.score:.1f}\t{d.alignment.target_from}\t{d.alignment.target_to}\n")
                n += 1
print(f"{n} domaines Pfam retenus")
EOF
  echo "Protéines avec ≥1 Pfam : $(tail -n +2 "$OUT/pfam_hits.tsv" | cut -f1 | sort -u | wc -l)"
  ;;

eggnog)
  conda activate funcannot
  emapper.py -i "$FAA" --itype proteins \
    --data_dir "$DBDIR/eggnog" \
    -m diamond --sensmode more-sensitive \
    --dbmem --cpu "$THREADS" \
    --output forel --output_dir "$OUT" --temp_dir "$OUT" --override
  echo "Protéines annotées eggNOG : $(grep -vc '^#' "$OUT/forel.emapper.annotations")"
  ;;

kofam)
  conda activate kofam
  KO=$DBDIR/kofam
  echo "Profils KOfam : $(ls "$KO/profiles" | wc -l)"
  mkdir -p "$OUT/kofam_tmp"
  exec_annotation -o "$OUT/kofam_detail.tsv" -f detail-tsv \
    -p "$KO/profiles/prokaryote.hal" -k "$KO/ko_list" \
    --cpu "$THREADS" --tmp-dir "$OUT/kofam_tmp" "$FAA"
  rm -rf "$OUT/kofam_tmp"
  # hits significatifs (seuil adaptatif KOfam, colonne '*') ; meilleur KO par gène
  awk -F'\t' '$1=="*"{print $2"\t"$3"\t"$5"\t"$6"\t"$7}' "$OUT/kofam_detail.tsv" \
    | sort -t$'\t' -k1,1 -k3,3gr | awk -F'\t' '!seen[$1]++' > "$OUT/kofam_best.tsv"
  echo "Protéines avec un KO significatif : $(wc -l < "$OUT/kofam_best.tsv")"
  ;;

*) echo "outil inconnu : $TOOL"; exit 1 ;;
esac
echo "== Terminé $TOOL | $(date)"
