from pathlib import Path
import tempfile
import os
import shutil

from genomics_data_index.storage.util import execute_command


class SnpEffDatabase:

    def __init__(self, snpeff_config: Path, database_dir: Path, genome_name: str):
        self._snpeff_config = snpeff_config
        self._database_dir = database_dir
        self._genome_name = genome_name

    @property
    def config(self) -> Path:
        return self._snpeff_config

    @property
    def genome_name(self) -> str:
        return self._genome_name

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
            os.symlink(input_vcf_file, tmp_vcf_path)

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
