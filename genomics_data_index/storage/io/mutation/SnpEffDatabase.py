from pathlib import Path


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

    def annotate(self, input_vcf_file: Path, output_vcf_file: Path = None) -> Path:
        """
        Annotates the given VCF file with the stored snpeff database.
        :param input_vcf_file: The VCF file to annotate.
        :param output_vcf_file: The output VCF file to write to. If left as None will write to {input_vcf_file}.ann
        :return: The annotated VCF file.
        """
        if output_vcf_file is None:
            output_vcf_file = input_vcf_file.with_name(str(input_vcf_file.name) + '.ann')

        return output_vcf_file
