import math

from genomics_data_index.storage.model.db import Sample
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.KmerService import KmerService


def test_find_matches_all(database: DatabaseConnection, kmer_service_with_data: KmerService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    all_matches_set = kmer_service_with_data.find_matches_within(['SampleA'], distance_threshold=1.0)

    assert {sampleA.id, sampleB.id, sampleC.id} == set(all_matches_set)
