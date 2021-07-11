from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

# Prepares a SnpEff database from a reference genome

input_reference = Path(snakemake.input.reference)
output_db_conf = Path(snakemake.output.snpeff_db_conf)
snpeff_db_dir = output_db_conf.parent.absolute()

sequence_file = SequenceFile(input_reference)
snpeff_database = sequence_file.create_snpeff_database(database_dir=snpeff_db_dir)
snpeff_database.write_database_config(output_db_conf)
