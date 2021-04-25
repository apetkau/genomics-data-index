import gzip

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from ete3 import Tree
from sqlalchemy.orm.exc import NoResultFound

from genomics_data_index.storage.model.db import Reference
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.test.integration import reference_file
from genomics_data_index.test.integration import tree_file


@pytest.fixture
def example_tree() -> Tree:
    return Tree(str(tree_file))


def test_count_references(reference_service):
    assert 0 == reference_service.count_reference_genomes()


def count_references_with_data(reference_service_with_data):
    assert 1 == reference_service_with_data.count_reference_genomes()


def test_translate_spdi(reference_service_with_data):
    spdi_ids = {
        'reference:3762:1:G',
        'reference:374:1:ATTCTAGGGTAGACGCT',
        'reference:3576:4:T',
        'reference:3513:1:G',
        'reference:1984:7:TTGA',
        'reference:1483:21:A',
    }

    translated_ids = reference_service_with_data.translate_spdi(spdi_ids, to='spdi_ref')
    assert translated_ids['reference:3762:1:G'] == 'reference:3762:A:G'
    assert translated_ids['reference:374:1:ATTCTAGGGTAGACGCT'] == 'reference:374:A:ATTCTAGGGTAGACGCT'
    assert translated_ids['reference:3576:4:T'] == 'reference:3576:TTTC:T'
    assert translated_ids['reference:3513:1:G'] == 'reference:3513:T:G'
    assert translated_ids['reference:1984:7:TTGA'] == 'reference:1984:GTGATTG:TTGA'
    assert translated_ids['reference:1483:21:A'] == 'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A'

    spdi_ids_ref = {
        'reference:3762:A:G',
        'reference:374:A:ATTCTAGGGTAGACGCT',
        'reference:3576:TTTC:T',
        'reference:3513:T:G',
        'reference:1984:GTGATTG:TTGA',
        'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A'
    }

    translated_ids = reference_service_with_data.translate_spdi(spdi_ids_ref, to='spdi_ref')
    assert translated_ids['reference:3762:A:G'] == 'reference:3762:A:G'
    assert translated_ids['reference:374:A:ATTCTAGGGTAGACGCT'] == 'reference:374:A:ATTCTAGGGTAGACGCT'
    assert translated_ids['reference:3576:TTTC:T'] == 'reference:3576:TTTC:T'
    assert translated_ids['reference:3513:T:G'] == 'reference:3513:T:G'
    assert translated_ids['reference:1984:GTGATTG:TTGA'] == 'reference:1984:GTGATTG:TTGA'
    assert translated_ids['reference:1483:AAAGAGGGGCTGCTGGAGCCG:A'] == 'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A'


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
        reference_genome.tree
    assert 'Cannot convert an empty tree' in str(execinfo.value)

    reference_service.update_tree(reference_name='genome', tree=example_tree, alignment_length=1000)
    reference_genome = database.get_session().query(Reference).filter(Reference.name == 'genome').one()
    assert reference_genome.tree.write() == example_tree.write()
    assert 1000 == reference_genome.tree_alignment_length


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


def test_get_reference_sequences(reference_service_with_data, variation_service):
    reference_sequences = reference_service_with_data.get_reference_sequences('genome')

    assert {'reference'} == set(reference_sequences.keys())


def test_find_reference_for_sequence(reference_service_with_data, variation_service):
    reference = reference_service_with_data.find_reference_for_sequence('reference')
    assert 'genome' == reference.name


def test_find_reference_for_sequence_not_exist(reference_service_with_data, variation_service):
    with pytest.raises(NoResultFound) as execinfo:
        reference_service_with_data.find_reference_for_sequence('not_exist')
    assert 'No row was found' in str(execinfo.value)


def test_get_sequence(reference_service_with_data):
    seq_record = reference_service_with_data.get_sequence('reference')

    assert 'reference' == seq_record.id
    assert isinstance(seq_record.seq, Seq)


def test_get_reference_genome_records(reference_service_with_data):
    records = reference_service_with_data.get_reference_genome_records('genome')
    assert 1 == len(records)

    assert 'reference' == records[0].id
    assert 5180 == len(records[0])
