#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
rm -rf ${SCRIPT_DIR}/{align,consensus,gdi-input.fofn,logs,mlst,mlst.tsv,reference,sketch,variant,.snakemake,snpeff_db}
