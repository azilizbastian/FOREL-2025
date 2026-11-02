# 00 – Read preparation, assembly and binning

Scripts are run once per sample (`ANACONDA_INSITU_7M`, `COBRA_ICEFIELD_1M`), from the sample directory.
Output folders carry the historical prefix `Medusa-filter-myloasm-v0.6.0-hifi`.

| Order | Script | Step | Software |
|---|---|---|---|
| 0 | `00_prepare_reads.sh` | Deduplicated HiFi BAM → FASTQ | pbtk [X] |
| 1 | `<SAMPLE>_run-myloasm.sh` | Metagenome assembly (`--hifi`) | myloasm v0.6.0 |
| 2 | `<SAMPLE>_run-semibin.sh` | Read mapping, then prokaryotic binning (self-supervised, long-read) | minimap2 [X], samtools [X], SemiBin2 [X] |
| 3 | `<SAMPLE>_run-whokaryote.sh` | Contig classification eukaryote/prokaryote | Whokaryote [X], Tiara [X] |
| 4 | `<SAMPLE>_run-post-processing.py` | Bin renaming, removal of bins < 500 kb, eukaryotic annotation → `<SAMPLE>-processed-bins/` | Python, polars, Biopython |
| 5 | `<SAMPLE>_run-remag.sh` | Eukaryotic binning (metagenomics mode) | REMAG [X] |

The BAM from step 2 is reused by REMAG (step 5).
Downstream analyses (`01_qc_taxonomy/`) start from `<SAMPLE>-processed-bins/`.
