# Assembly input

This Snakemake workflow takes assembled genomes as input and generates VCF variant calls and a masked consenus sequences to be used as input for indexing. This works by:

1. Taking as input a config.yaml file defining the inputs (see `config/config.yaml` for an example)
    1. `reference`: The reference input file.
    2. `samples`: The input sample nams/files (as assembled genomes).
2. Aligns samples to the reference genome using `minimap2`.
3. Uses `htsbox` to identify variants (VCF files) and a masked consenus sequence for each sample.
4. Outputs a file of file names to be used as input to the larger genomics index application. 
