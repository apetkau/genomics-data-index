# snpeff test data

This data is for testing parsing the `ANN` field of [snpeff](http://pcingola.github.io/SnpEff/). To prepare the data I did the following:

1. Downloaded sequence data for `SH10-014`, `SH14-001`, and `SH14-014` from ENA.
2. Downloaded `NC_011083.gbk` from Genbank.
3. Converted `NC_011083.gbk` to fasta format `NC_011083.fasta`.
4. Ran snippy with `NC_011083.fasta` as the reference genome.
5. Created snpeff database with `NC_011083.gbk`: `snpEff build -genbank -v NC_011083`. This also involved creating a custom snpeff config file.
6. Ran snpeff on snippy output VCF files: `snpEff ann -c ../../snpEff/snpEff.config NC_011083 snps.vcf.gz | bgzip > SH14-014.vcf.gz`
7. Copied data over.

I ran snpeff myself instead of through snippy so that the default snpeff annotations show up in `ANN` (so I can make sure I handle these cases).

## missing/unknown BED files

I downloaded the `salmonella-project.zip` archive from <https://github.com/apetkau/genomics-data-index-examples>, which contains processed files (by snippy). I copied over the corresponding `*.bed.gz` files containing unknown/missing positions from this archive.

