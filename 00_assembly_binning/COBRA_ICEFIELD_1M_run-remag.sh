#!/usr/bin/bash

set -euo pipefail

remag Medusa-filter-myloasm-v0.6.0-hifi/assembly_primary.fa -c Medusa-filter-myloasm-v0.6.0-hifi.sorted.bam -o Medusa-filter-myloasm-v0.6.0-hifi-REMAG -m metagenomics -t 110 

