import pytest
from typing import Dict
import tempfile
from pathlib import Path
import re

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file_5000_snpeff, snpeff_vcf_file


def test_annotate_vcf_file():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        output_vcf_file = database_dir / 'output.vcf.gz'

        snpeff_database = SequenceFile(reference_file_5000_snpeff).create_snpeff_database(database_dir)

        returned_output = snpeff_database.annotate(input_vcf_file=snpeff_vcf_file, output_vcf_file=output_vcf_file)

        assert output_vcf_file == returned_output

