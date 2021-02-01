from typing import Dict, Any
import pytest
import tempfile
import gzip
from Bio import SeqIO
from pathlib import Path

from storage.variant.model import Reference
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.test.integration.variant import reference_file


@pytest.fixture
def setup() -> Dict[str, Any]:
    seq_repo_root = Path(tempfile.mkdtemp(prefix='index-test'))
    database = DatabaseConnection('sqlite:///:memory:')

    reference_service = ReferenceService(database, seq_repo_root)

    return dict(database=database, reference_service=reference_service)


def test_insert_reference_genome(setup):
    assert 0 == setup['database'].get_session().query(Reference).count(), 'Database should be empty initially'
    setup['reference_service'].create_reference_genome(reference_file)
    assert 1 == setup['database'].get_session().query(Reference).count(), 'Database should have one entry'
    assert 'genome' == setup['database'].get_session().query(Reference).all()[0].name, 'Name should match'


def test_find_reference_genome(setup):
    assert 0 == setup['database'].get_session().query(Reference).count(), 'Database should be empty initially'
    setup['reference_service'].create_reference_genome(reference_file)

    reference = setup['reference_service'].find_reference_genome('genome')
    assert 'genome' == reference.name, 'Reference name should match'
    assert 5180 == reference.length, 'Reference length should match'


def test_add_get_reference_genome(setup):
    reference_service = setup['reference_service']

    reference_service.add_reference_genome(reference_file)
    seq_record = reference_service.get_sequence('reference')

    assert seq_record is not None, 'No matching seq record'
    assert seq_record.id == 'reference', 'Incorrect record id'

    with gzip.open(reference_file, mode='rt', encoding='ascii') as f:
        records = list(SeqIO.parse(f, 'fasta'))
        assert records[0].seq == seq_record.seq, 'Incorrect sequence'


def test_get_nonexistent_sequence(setup):
    reference_service = setup['reference_service']

    with pytest.raises(KeyError) as execinfo:
        reference_service.get_sequence('does_not_exist')
    assert 'Alias does_not_exist' in str(execinfo.value)