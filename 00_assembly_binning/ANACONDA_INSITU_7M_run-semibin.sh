#!/usr/bin/bash

set -euo pipefail

minimap2 -ax map-hifi -t 110 --secondary=no Medusa-filter-myloasm-v0.6.0-hifi/assembly_primary.fa deduplicated.ccs.fastq.gz \
  | samtools sort -@ 16 -m 4G -o Medusa-filter-myloasm-v0.6.0-hifi.sorted.bam -

SemiBin2 single_easy_bin \
    --input-fasta Medusa-filter-myloasm-v0.6.0-hifi/assembly_primary.fa \
    --input-bam Medusa-filter-myloasm-v0.6.0-hifi.sorted.bam \
    --self-supervised \
    --sequencing-type=long_read \
    --threads 126 \
    --output Medusa-filter-myloasm-v0.6.0-hifi-Semibin2 

