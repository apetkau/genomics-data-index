import gzip

import pytest
from Bio import SeqIO
from ete3 import Tree

from storage.test.integration.variant import reference_file
from storage.variant.model import Reference
from storage.variant.service import EntityExistsError
from storage.test.integration.variant import tree_file


@pytest.fixture
def example_tree() -> Tree:
    return Tree(str(tree_file))


def test_insert_reference_genome(database, reference_service):
    assert 0 == database.get_session().query(Reference).count(), 'Database should be empty initially'
    reference_service.add_reference_genome(reference_file)
    assert 1 == database.get_session().query(Reference).count(), 'Database should have one entry'
    assert 'genome' == database.get_session().query(Reference).all()[0].name, 'Name should match'

    seq_record = reference_service.get_sequence('reference')

    assert seq_record is not None, 'No matching seq record'
    assert seq_record.id == 'reference', 'Incorrect record id'

    with gzip.open(reference_file, mode='rt', encoding='ascii') as f:
        records = list(SeqIO.parse(f, 'fasta'))
        assert records[0].seq == seq_record.seq, 'Incorrect sequence'


def test_double_insert_reference_genome(database, reference_service):
    assert 0 == database.get_session().query(Reference).count(), 'Database should be empty initially'
    reference_service.add_reference_genome(reference_file)
    assert 1 == database.get_session().query(Reference).count(), 'Database should have one entry'
    assert 'genome' == database.get_session().query(Reference).all()[0].name, 'Name should match'

    with pytest.raises(EntityExistsError) as execinfo:
        reference_service.add_reference_genome(reference_file)
    assert 'Reference genome [genome] already exists in database' in str(execinfo.value)
    assert 1 == database.get_session().query(Reference).count(), 'Database should have one entry'


def test_find_reference_genome(database, reference_service):
    assert 0 == database.get_session().query(Reference).count(), 'Database should be empty initially'
    reference_service.add_reference_genome(reference_file)

    reference = reference_service.find_reference_genome('genome')
    assert 'genome' == reference.name, 'Reference name should match'
    assert 5180 == reference.length, 'Reference length should match'


def test_get_nonexistent_sequence(reference_service):
    with pytest.raises(KeyError) as execinfo:
        reference_service.get_sequence('does_not_exist')
    assert 'Alias does_not_exist' in str(execinfo.value)


def test_get_reference_genomes(reference_service):
    reference_service.add_reference_genome(reference_file)
    assert {'genome'} == {genome.name for genome in reference_service.get_reference_genomes()}

def test_get_reference_genomes_empty(reference_service):
    assert [] == reference_service.get_reference_genomes()


def test_update_tree(database, reference_service, example_tree):
    reference_service.add_reference_genome(reference_file)

    reference_genome = database.get_session().query(Reference).filter(Reference.name == 'genome').one()
    assert reference_genome is not None
    with pytest.raises(Exception) as execinfo:
        reference_genome.tree()
    assert 'Cannot convert an empty tree' in str(execinfo.value)

    reference_service.update_tree(reference_name='genome', tree=example_tree)
    reference_genome = database.get_session().query(Reference).filter(Reference.name == 'genome').one()
    assert reference_genome.tree.write() == example_tree.write()


def test_find_references_for_sample(reference_service_with_data, variation_service):
    found_references = reference_service_with_data.find_references_for_sample('SampleA')
    assert len(found_references) == 1
    assert 'genome' == found_references[0].name

    found_references = reference_service_with_data.find_references_for_sample('SampleB')
    assert len(found_references) == 1
    assert 'genome' == found_references[0].name

    found_references = reference_service_with_data.find_references_for_sample('SampleC')
    assert len(found_references) == 1
    assert 'genome' == found_references[0].name


def test_find_references_for_sample_not_exist(reference_service_with_data, variation_service):
    found_references = reference_service_with_data.find_references_for_sample('not_exist')
    assert len(found_references) == 0
