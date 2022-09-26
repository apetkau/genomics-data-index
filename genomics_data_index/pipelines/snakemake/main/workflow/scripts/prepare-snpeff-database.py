from pathlib import Path
from packaging import version

import genomics_data_index
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

gdi_version = version.parse(genomics_data_index.__version__)

# Prepares a SnpEff database from a reference genome

input_reference = Path(snakemake.input.reference)
output_db_conf = Path(snakemake.output.snpeff_db_conf)
snpeff_db_dir = output_db_conf.parent.absolute()

sequence_file = SequenceFile(input_reference)

# Because of the use of conda environments and different versions of GDI
# I have to have two different calls to create the SnpEff database here
# To account for different versions of the code
if gdi_version <= version.parse("0.7.0"):
    snpeff_database = sequence_file.create_snpeff_database(database_dir=snpeff_db_dir)
else:
    snpeff_database = sequence_file.create_snpeff_database(database_dir=snpeff_db_dir,
            no_check_protein=snakemake.config["snpeff"]["no_check_protein"],
            no_check_cds=snakemake.config["snpeff"]["no_check_cds"])

snpeff_database.write_database_config(output_db_conf)
