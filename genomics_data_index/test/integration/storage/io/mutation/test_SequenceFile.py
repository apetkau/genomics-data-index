import re
import tempfile
from pathlib import Path
from typing import Dict

from Bio import SeqIO

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file_5000_snpeff, reference_file_5000_snpeff_2, \
    reference_file


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


def test_write():
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir, 'output.fasta')
        sequence_file = SequenceFile(reference_file_5000_snpeff)
        sequence_file.write(out_file, 'fasta')

        assert out_file.exists()

        with open(out_file, 'r') as fh:
            sequences_out = list(SeqIO.parse(fh, 'fasta'))

        ref_name, sequences_in = sequence_file.parse_sequence_file()

        assert 1 == len(sequences_out)
        assert 1 == len(sequences_in)

        assert sequences_in[0].id == sequences_out[0].id
        assert sequences_in[0].seq == sequences_out[0].seq


def test_write_multiple_sequences():
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir, 'output.fasta')
        sequence_file = SequenceFile(reference_file_5000_snpeff_2)
        sequence_file.write(out_file, 'fasta')

        assert out_file.exists()

        with open(out_file, 'r') as fh:
            sequences_out = list(SeqIO.parse(fh, 'fasta'))

        ref_name, sequences_in = sequence_file.parse_sequence_file()

        assert 2 == len(sequences_out)
        assert 2 == len(sequences_in)

        assert sequences_in[0].id == sequences_out[0].id
        assert sequences_in[0].seq == sequences_out[0].seq
        assert sequences_in[1].id == sequences_out[1].id
        assert sequences_in[1].seq == sequences_out[1].seq


def test_create_snpeff_database():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        snpeff_database_dir = database_dir / 'db'
        sequence_file = SequenceFile(reference_file_5000_snpeff)
        snpeff_database = sequence_file.create_snpeff_database(database_dir)

        expected_genbank_path = snpeff_database_dir / 'NC_011083-5000' / 'genes.gbk.gz'
        expected_snpeff_bin_file = snpeff_database_dir / 'NC_011083-5000' / 'snpEffectPredictor.bin'

        assert snpeff_database.config.name == 'snpEff.config'
        assert snpeff_database.config.exists()
        assert snpeff_database.config.parent == database_dir

        entries = parse_snpeff_config(snpeff_database.config)

        assert entries['data.dir'] == f'{snpeff_database_dir}'
        assert 'NC_011083-5000.genome' in entries
        assert entries['NC_011083-5000.chromosomes'] == 'NC_011083.1'
        assert entries['NC_011083-5000.NC_011083.1.codonTable'] == 'Standard'

        assert expected_genbank_path.exists()
        assert expected_genbank_path.is_symlink()
        assert expected_snpeff_bin_file.exists()


def test_can_use_snpeff():
    assert SequenceFile(reference_file_5000_snpeff).can_use_snpeff()
    assert SequenceFile(reference_file_5000_snpeff_2).can_use_snpeff()
    assert not SequenceFile(reference_file).can_use_snpeff()


def test_create_snpeff_database_multiple_contigs():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        snpeff_database_dir = database_dir / 'db'
        sequence_file = SequenceFile(reference_file_5000_snpeff_2)

        snpeff_database = sequence_file.create_snpeff_database(database_dir, codon_type='Bacterial_and_Plant_Plastid')

        expected_genbank_path = snpeff_database_dir / 'NC_011083_CP001602-5000' / 'genes.gbk.gz'
        expected_snpeff_bin_file = snpeff_database_dir / 'NC_011083_CP001602-5000' / 'snpEffectPredictor.bin'

        assert snpeff_database.config.name == 'snpEff.config'
        assert snpeff_database.config.exists()
        assert snpeff_database.config.parent == database_dir

        entries = parse_snpeff_config(snpeff_database.config)

        assert entries['data.dir'] == f'{snpeff_database_dir}'
        assert 'NC_011083_CP001602-5000.genome' in entries
        assert entries['NC_011083_CP001602-5000.chromosomes'] == 'NC_011083.1,CP001602.2'
        assert entries['NC_011083_CP001602-5000.NC_011083.1.codonTable'] == 'Bacterial_and_Plant_Plastid'
        assert entries['NC_011083_CP001602-5000.CP001602.2.codonTable'] == 'Bacterial_and_Plant_Plastid'

        assert expected_genbank_path.exists()
        assert expected_genbank_path.is_symlink()
        assert expected_snpeff_bin_file.exists()
