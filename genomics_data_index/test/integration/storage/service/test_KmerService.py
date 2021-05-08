import pytest

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.KmerService import KmerService


def test_find_matches_all(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=1.0)

    assert {sampleA.id, sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_all_in_sub_universe(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    samples_universe = SampleSet([sampleB.id, sampleC.id])

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=1.0,
                                                                 samples_universe=samples_universe)

    assert {sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_all_in_empty_universe(database: DatabaseConnection, kmer_service_with_data: KmerService):
    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=1.0,
                                                                 samples_universe=SampleSet.create_empty())
    assert set() == set(all_matches_set)


def test_find_matches_all_lower_threshold(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=0.522)

    assert {sampleA.id, sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_two(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=0.521)

    assert {sampleA.id, sampleC.id} == set(all_matches_set)


def test_find_matches_one(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=31,
                                                                 distance_threshold=0.49)

    assert {sampleA.id} == set(all_matches_set)


def test_find_matches_all_lower_different_genome(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleC'], kmer_size=31,
                                                                 distance_threshold=0.5)

    assert {sampleA.id, sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_one_different_genome(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleC'], kmer_size=31,
                                                                 distance_threshold=0.49)

    assert {sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_two_different_genomes(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA', 'SampleC'], kmer_size=31,
                                                                 distance_threshold=0.49)

    assert {sampleA.id, sampleB.id, sampleC.id} == set(all_matches_set)


def test_find_matches_different_kmer(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=21,
                                                                 distance_threshold=0.4)

    assert {sampleA.id, sampleC.id} == set(all_matches_set)


def test_find_matches_nonindexed_kmer_size(database: DatabaseConnection, kmer_service_with_data: KmerService):
    with pytest.raises(Exception) as execinfo:
        kmer_service_with_data.find_matches_within(['SampleA'], kmer_size=11,
                                                   distance_threshold=1.0)
    assert 'Could not run' in str(execinfo.value)
