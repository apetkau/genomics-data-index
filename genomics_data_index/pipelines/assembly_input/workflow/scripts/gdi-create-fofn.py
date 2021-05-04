import os

# Generates a file of file names (fofn) from the given samples and VCF/consensus files
# This file of file names is then used for input to the genomic indexer (gdi)

# TODO: Instead of parallel lists I should probably have a dictionary keyed by the sample name
sample_vcfs = snakemake.input.sample_vcfs
sample_consensus = snakemake.input.sample_consensus

if len(sample_vcfs) != len(sample_consensus):
    raise Exception(f'sample_vcfs length [{len(sample_vcfs)}] '
                    f'does not equal sample_consensus length [{len(sample_consensus)}]')

with open(snakemake.output[0], 'w') as fh:
    fh.write('Sample\tVCF\tMask File\n')
    for idx in range(0,len(sample_consensus)):
        vcf = os.path.abspath(sample_vcfs[idx])
        consensus = os.path.abspath(sample_consensus[idx])
        fh.write(f'\t{vcf}\t{consensus}\n')
