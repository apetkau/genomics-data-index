from pathlib import Path

from genomics_data_index.storage.io.mutation.SnpEffDatabase import SnpEffDatabase

# Annotate VCF with snpeff

snpeff_database_conf = Path(snakemake.input.snpeff_db_conf)
snpeff_database = SnpEffDatabase.create_from_config(snpeff_database_conf)
input_vcf = Path(snakemake.input.vcf).absolute()
output_vcf = Path(snakemake.output.snpeff_vcf).absolute()

snpeff_database.annotate(input_vcf_file=input_vcf, output_vcf_file=output_vcf)
