from typing import Dict
import tempfile
from pathlib import Path
import re

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file_5000_snpeff


def parse_snpeff_config(config_file: Path) -> Dict[str, str]:
    key_values = {}
    with open(config_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith('#'):
                tokens = re.split(r'[:=]', line)

                if len(tokens) == 2:
                    key_values[tokens[0].strip()] = tokens[1].strip()

    return key_values


def test_create_snpeff_database():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        snpeff_database_dir = database_dir / 'db'
        sequence_file = SequenceFile(reference_file_5000_snpeff)
        snpeff_config = sequence_file.create_snpeff_database(database_dir)

        assert snpeff_config.name == 'snpEff.config'
        assert snpeff_config.exists()
        assert snpeff_config.parent == database_dir

        entries = parse_snpeff_config(snpeff_config)

        print(entries.keys())
        assert entries['data.dir'] == f'{snpeff_database_dir}'
