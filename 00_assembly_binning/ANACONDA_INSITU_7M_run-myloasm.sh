#!/usr/bin/bash

set -euo pipefail

myloasm deduplicated.ccs.fastq.gz -o Medusa-filter-myloasm-v0.6.0-hifi --hifi --threads 64
