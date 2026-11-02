#!/usr/bin/bash
# HiFi reads deduplicated by the sequencing facility (University of Bern, AmpliFi protocol)
# Input: deduplicated.ccs.bam (PCR duplicates removed by the facility)
set -euo pipefail
pbindex deduplicated.ccs.bam
bam2fastq -o deduplicated.ccs deduplicated.ccs.bam   # -> deduplicated.ccs.fastq.gz
