#!/usr/bin/bash

set -euo pipefail

whokaryote.py --contigs Medusa-filter-myloasm-v0.6.0-hifi/assembly_primary.fa --outdir Medusa-filter-myloasm-v0.6.0-hifi-whokaryote --minsize 100000 --threads 64 
