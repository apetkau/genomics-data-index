The `Sample*/` folders contain data generated by snippy (when run on genomes with simulated mutations).

Generated `Sample*/snps.fill-tags.vcf.gz` from snippy output file `snps.vcf.gz` by running:

```
bcftools plugin fill-tags snps.vcf.gz -O z -o snps.fill-tags.vcf.gz -- -t TYPE
bcftools index snps.fill-tags.vcf.gz
```

Generated `Sample*/mutations-dataframe.snps.tsv` with:

```
echo -e "Mutation\tSequence\tPosition\tDeletion\tInsertion" > mutations-dataframe.snps.tsv
zgrep -v '^#' snps.fill-tags.vcf.gz | grep SNP | cut -f 1,2,4,5 | tr '\t' ':' | perl -ne 'chomp;@a=split(/:/, $_);print($_, "\t", $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3], "\n");' >> mutations-dataframe.snps.tsv
```

Generated `Sample*/mutations-dataframe.all.tsv` with:

```
echo -e "Mutation\tSequence\tPosition\tDeletion\tInsertion" > mutations-dataframe.all.tsv
zgrep -v '^#' snps.fill-tags.vcf.gz | cut -f 1,2,4,5 | tr '\t' ':' | perl -ne 'chomp;@a=split(/:/, $_);print($_, "\t", $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3], "\n");' >> mutations-dataframe.all.tsv
```