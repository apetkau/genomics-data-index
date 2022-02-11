from Bio import SeqIO

# Converts sequences in the reference fasta file to a genome BED file
# <https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format>
# Used to get the complement of all regions in the VCF file with reference/variants (i.e.,
# used to construct a BED file of all unknown/missing positions)
# Outputs a format like:
## SequenceID    length
# seq1  100
# seq2 500
# ...

input_fasta = snakemake.input.reference
output_file = snakemake.output.genome_bedtools

with open(output_file, 'w') as oh:
    with open(input_fasta, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            oh.write(f'{record.id}\t{len(record)}\n')
