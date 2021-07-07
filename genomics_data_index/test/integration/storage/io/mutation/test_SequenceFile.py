import tempfile
from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file_5000_snpeff


def test_create_snpeff_database():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        sequence_file = SequenceFile(reference_file_5000_snpeff)
        snpeff_config = sequence_file.create_snpeff_database(database_dir)

        assert snpeff_config.name == 'snpEff.config'
        assert snpeff_config.exists()
        assert snpeff_config.parent == database_dir
