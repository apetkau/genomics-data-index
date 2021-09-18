import os


# Generates a file of file names (fofn) from the given samples and VCF/consensus files
# This file of file names is then used for input to the genomic indexer (gdi)


# From https://stackoverflow.com/a/51779415
def remove_suffix(s, suffix):
    return s[:-len(suffix)] if s.endswith(suffix) else s


# TODO: Instead of parallel lists I should probably have a dictionary keyed by the sample name
sample_vcfs = snakemake.input.sample_vcfs
sample_masks = snakemake.input.sample_masks
sample_sketches = snakemake.input.sample_sketches
include_kmer = snakemake.config["include_kmer"]
output_file = snakemake.output[0]

if len(sample_vcfs) != len(sample_masks):
    raise Exception(f'sample_vcfs length [{len(sample_vcfs)}] '
                    f'does not equal sample_masks length [{len(sample_masks)}]')

if include_kmer and len(sample_vcfs) != len(sample_sketches):
    raise Exception(f'sample_vcfs length [{len(sample_vcfs)}] '
                    f'does not equal sample_sketches length [{len(sample_sketches)}]')

with open(output_file, 'w') as fh:
    if include_kmer:
        fh.write('Sample\tVCF\tMask File\tSketch File\n')
    else:
        fh.write('Sample\tVCF\tMask File\n')

    for idx in range(0, len(sample_masks)):
        vcf = os.path.abspath(sample_vcfs[idx])
        mask = os.path.abspath(sample_masks[idx])

        # I check both sample names to guarantee that files don't get mixed up
        # I don't know enough about snakemake to solve this differently. Ideally I'd like to just pass a dictionary
        # to this script which maps the sample name to the vcf and consensus file
        sample_name_vcf = remove_suffix(os.path.basename(vcf), '.vcf.gz')
        sample_name_mask = remove_suffix(os.path.basename(mask), '.bed.gz')

        if sample_name_mask != sample_name_vcf:
            raise Exception(f'VCF=[{vcf}] has a different sample name from mask=[{mask}]')

        if include_kmer:
            sketch = os.path.abspath(sample_sketches[idx])
            sample_name_sketch = remove_suffix(os.path.basename(sketch), '.sig.gz')

            if sample_name_vcf != sample_name_sketch:
                raise Exception(f'VCF=[{vcf}] has a different sample name from sketch=[{sketch}]')

            fh.write(f'{sample_name_vcf}\t{vcf}\t{mask}\t{sketch}\n')
        else:
            fh.write(f'{sample_name_vcf}\t{vcf}\t{mask}\n')
