from typing import Dict
import tempfile
from pathlib import Path
import re

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file_5000_snpeff, reference_file_5000_snpeff_2, reference_file


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


def test_read_fasta_file():
    sequence_file = SequenceFile(reference_file)
    name, records = sequence_file.parse_sequence_file()

    assert 'genome' == name
    assert 1 == len(records)
    record = records[0]
    assert 'reference' == record.id
    assert 5180 == len(record)


def test_read_genbank_file():
    sequence_file = SequenceFile(reference_file_5000_snpeff)
    name, records = sequence_file.parse_sequence_file()

    assert 'NC_011083-5000' == name
    assert 1 == len(records)
    record = records[0]
    assert 'NC_011083.1' == record.id
    assert 5000 == len(record)


def test_read_genbank_file_multiple_sequences():
    sequence_file = SequenceFile(reference_file_5000_snpeff_2)
    name, records = sequence_file.parse_sequence_file()

    assert 'NC_011083_CP001602-5000' == name
    assert 2 == len(records)

    record = records[0]
    assert 'NC_011083.1' == record.id
    assert 5000 == len(record)

    record = records[1]
    assert 'CP001602.2' == record.id
    assert 5000 == len(record)


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

        assert entries['data.dir'] == f'{snpeff_database_dir}'
        assert 'NC_011083-5000.genome' in entries
        assert entries['NC_011083-5000.chromosomes'] == 'NC_011083.1'
        assert entries['NC_011083-5000.NC_011083.1.codonTable'] == 'Bacterial_and_Plant_Plastid'
