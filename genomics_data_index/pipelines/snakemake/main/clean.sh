#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
rm -rf ${SCRIPT_DIR}/{assemblies,reads,gdi-input.fofn,logs,mlst,mlst.tsv,reference,.snakemake,snpeff_db}
