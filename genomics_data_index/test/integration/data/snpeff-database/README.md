To generate example VCF file to apply annotations to:

1. Extracted sequence from `NC_011083-5000.gbk.gz` to `ref.fasta.gz`
2. Copied `ref.fasta.gz` to `input/SampleA.fasta.gz`
3. Manually modified `input/SampleA.fasta.gz` to create mutations at specific codons (lower-case letters in sequence)
4. Ran `gdi analysis assembly` on these files to create a VCF file containing the mutations I inserted.
5. The VCF file from the index was copied to `SmapleA.vcf.gz` which acts as the input file for running SnpEff and validating the correct effects get inserted into the VCF.
