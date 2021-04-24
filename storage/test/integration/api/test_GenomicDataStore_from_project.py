from storage.api.GenomicDataStore import GenomicDataStore
from storage.variant.model.db import Sample
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST


def test_query_single_mutation(loaded_data_store_from_project_dir: GenomicDataStore):
    ds = loaded_data_store_from_project_dir
    db = ds.connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = ds.samples_query().has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_alleles(loaded_data_store_from_project_dir: GenomicDataStore):
    ds = loaded_data_store_from_project_dir
    db = ds.connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = ds.samples_query().has(
        QueryFeatureMLST('lmonocytogenes:abcZ:1')).has(
        QueryFeatureMLST('lmonocytogenes:lhkA:4'))
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
