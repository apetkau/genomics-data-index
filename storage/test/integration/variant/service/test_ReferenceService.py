import gzip

import pytest
from Bio import SeqIO

from storage.test.integration.variant import reference_file
from storage.variant.model import Reference


def test_insert_reference_genome(database, reference_service):
    assert 0 == database.get_session().query(Reference).count(), 'Database should be empty initially'
    reference_service.create_reference_genome(reference_file)
    assert 1 == database.get_session().query(Reference).count(), 'Database should have one entry'
    assert 'genome' == database.get_session().query(Reference).all()[0].name, 'Name should match'


def test_find_reference_genome(database, reference_service):
    assert 0 == database.get_session().query(Reference).count(), 'Database should be empty initially'
    reference_service.create_reference_genome(reference_file)

    reference = reference_service.find_reference_genome('genome')
    assert 'genome' == reference.name, 'Reference name should match'
    assert 5180 == reference.length, 'Reference length should match'


def test_add_get_reference_genome(reference_service):
    reference_service.add_reference_genome(reference_file)
    seq_record = reference_service.get_sequence('reference')

    assert seq_record is not None, 'No matching seq record'
    assert seq_record.id == 'reference', 'Incorrect record id'

    with gzip.open(reference_file, mode='rt', encoding='ascii') as f:
        records = list(SeqIO.parse(f, 'fasta'))
        assert records[0].seq == seq_record.seq, 'Incorrect sequence'


def test_get_nonexistent_sequence(reference_service):
    with pytest.raises(KeyError) as execinfo:
        reference_service.get_sequence('does_not_exist')
    assert 'Alias does_not_exist' in str(execinfo.value)
