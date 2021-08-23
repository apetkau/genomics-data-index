# Snakemake analysis pipeline

This Snakemake workflow takes assembled genomes as input and generates VCF variant calls and a BED file of unknown/missing regions 
to be used as input for indexing. This works by:

1. Taking as input a config.yaml file defining the inputs (see `config/config.yaml` for an example)
    1. `reference`: The reference input file.
    2. `samples`: The input sample nams/files (as assemblies or reads genomes).
2. For each assembly
    1. Aligns samples to the reference genome using `minimap2`.
    2. Uses `bcftools` to identify variants (output VCF files).
    3. Used `bedtools` to identify unknown/missing regions (output BED files).
    4. (Optional) Annotates VCF using snpeff.
3. For each set of reads
    1. Identifies variants using [snippy](https://github.com/tseemann/snippy)
    2. Uses `bedtools` and the `.aligned.fa` and `snps.vcf.gz` files from snippy to generate a BED file of unknown/missing regions.
4. (Optional) If requested, uses [sourmash](https://github.com/sourmash-bio/sourmash) to generate k-mer sketches.
5. (Optional) If requested, uses [mlst](https://github.com/tseemann/mlst) to generate MLST results (only for assemblies).
5. Outputs a file of file names to be used as input to the larger genomics index application. 

# Test

## Basic assemblies

To test this pipeline independent of the Python software you can run the following from within this directory:

```bash
snakemake -j 1
```

## snpeff assemblies

To test the pipeline with snpeff annotations you can run:

```bash
snakemake --configfile config/config-snpeff.yaml -j 1
```

## Basic reads

```bash
# Need to --use-conda since snippy dependencies (e.g., samtools) requires specific versions
snakemake --use-conda --configfile config/config-reads.yaml -j 1
```

# Cleanup

To cleanup, you can run:

```bash
./clean.sh
```
