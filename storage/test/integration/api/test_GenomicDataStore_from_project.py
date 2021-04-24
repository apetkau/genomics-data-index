from storage.api.GenomicDataStore import GenomicDataStore
from storage.variant.model.db import Sample
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation


def test_query_single_mutation(loaded_data_store_from_project_dir: GenomicDataStore):
    ds = loaded_data_store_from_project_dir
    db = ds.connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = ds.samples_query().has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
