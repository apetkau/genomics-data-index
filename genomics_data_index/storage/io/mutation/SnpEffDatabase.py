from __future__ import annotations

import configparser
import logging
import os
import tempfile
from pathlib import Path

from genomics_data_index.storage.util import execute_command

logger = logging.getLogger(__name__)


class SnpEffDatabase:

    def __init__(self, snpeff_config: Path, database_dir: Path, genome_name: str):
        self._snpeff_config = snpeff_config
        self._database_dir = database_dir
        self._genome_name = genome_name

    def write_database_config(self, file: Path) -> None:
        database_config = configparser.ConfigParser()
        database_config['snpeff'] = {
            'snpeff_config': str(self._snpeff_config),
            'database_dir': str(self._database_dir),
            'genome_name': str(self._genome_name),
        }

        logger.debug(f'Writing snpeff database config to [{file}]')
        with open(file, 'w') as fh:
            database_config.write(fh)

    @property
    def config(self) -> Path:
        return self._snpeff_config

    @property
    def genome_name(self) -> str:
        return self._genome_name

    @property
    def database_dir(self) -> Path:
        return self._database_dir

    def annotate(self, input_vcf_file: Path, output_vcf_file: Path) -> Path:
        """
        Annotates the given VCF file with the stored snpeff database.
        :param input_vcf_file: The VCF file to annotate.
        :param output_vcf_file: The output VCF file to write to.
        :return: The annotated VCF file.
        """
        # Create temporary directory to avoid cluttering up location of main VCF file with additional outputs (if any)
        # And to provide a place to output a VCF file and then bgzip it
        with tempfile.TemporaryDirectory() as out_dir:
            tmp_out_dir = Path(out_dir)
            tmp_vcf_path = tmp_out_dir / input_vcf_file.name
            tmp_vcf_out = tmp_out_dir / 'output.vcf'
            os.symlink(input_vcf_file.absolute(), tmp_vcf_path)

            snpeff_command = ['snpEff', 'ann', '-v', '-c', str(self._snpeff_config),
                              '-hgvs', '-hgvs1LetterAa', '-no-upstream', '-no-downstream', '-i', 'vcf', '-o', 'vcf',
                              '-noStats', '-sequenceOntology', '-nodownload', '-noLog',
                              self._genome_name, str(tmp_vcf_path)]

            # Execute snpeff and write data to output_vcf_h
            with open(tmp_vcf_out, 'w') as out_vcf_h:
                execute_command(snpeff_command, stdout=out_vcf_h)

            # compress file with bgzip
            bgzip_command = ['bgzip', '--stdout', str(tmp_vcf_out)]
            with open(output_vcf_file, 'w') as bgzip_out_h:
                execute_command(bgzip_command, stdout=bgzip_out_h)

        return output_vcf_file

    @classmethod
    def create_from_config(cls, config_file: Path) -> SnpEffDatabase:
        database_config = configparser.ConfigParser()
        database_config.read(config_file)

        snpeff_config = database_config['snpeff']['snpeff_config']
        database_dir = database_config['snpeff']['database_dir']
        genome_name = database_config['snpeff']['genome_name']

        if snpeff_config is None or database_dir is None or genome_name is None:
            raise Exception(f'Invalid snpeff database config in [{config_file}]')
        else:
            snpeff_config = Path(snpeff_config)
            database_dir = Path(database_dir)

        return SnpEffDatabase(snpeff_config=snpeff_config,
                              database_dir=database_dir,
                              genome_name=genome_name)
