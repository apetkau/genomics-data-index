import pytest

from storage.variant.service.MLSTQueryService import MLSTQueryService, QueryFeatureMLST
from storage.variant.model import Sample


@pytest.fixture
def mlst_query_service(mlst_service_loaded, sample_service) -> MLSTQueryService:
    return MLSTQueryService(mlst_service=mlst_service_loaded,
                                sample_service=sample_service)


def test_find_by_features(database, mlst_query_service: MLSTQueryService):
    sample = database.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()

    matches_df = mlst_query_service.find_by_features([QueryFeatureMLST('campylobacter:aspA:2')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'2014D-0067'} == set(matches_df['Sample Name'].tolist())
    assert {sample.id} == set(matches_df['Sample ID'].tolist())
    assert {'campylobacter:aspA:2'} == set(matches_df['Feature'].tolist())
    assert {'mlst'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 1
