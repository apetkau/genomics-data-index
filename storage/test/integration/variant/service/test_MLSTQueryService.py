import math

import pandas as pd
import pytest

from storage.variant.model import Sample
from storage.variant.service.MLSTQueryService import MLSTQueryService, QueryFeatureMLST


@pytest.fixture
def mlst_query_service(mlst_service_loaded, sample_service) -> MLSTQueryService:
    return MLSTQueryService(mlst_service=mlst_service_loaded,
                            sample_service=sample_service)

@pytest.fixture
def mlst_query_service_unknown(mlst_service_loaded_unknown, sample_service) -> MLSTQueryService:
    return MLSTQueryService(mlst_service=mlst_service_loaded_unknown,
                            sample_service=sample_service)


def test_find_by_features(database, mlst_query_service: MLSTQueryService):
    sample1 = database.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == '2014D-0068').one()

    matches_df = mlst_query_service.find_by_features([QueryFeatureMLST('campylobacter:aspA:2')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'2014D-0067', '2014D-0068'} == set(matches_df['Sample Name'].tolist())
    assert {sample1.id, sample2.id} == set(matches_df['Sample ID'].tolist())
    assert {'campylobacter:aspA:2'} == set(matches_df['Feature'].tolist())
    assert {'mlst'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 2


def test_find_by_features_unknown(database, mlst_query_service_unknown: MLSTQueryService):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    matches_df = mlst_query_service_unknown.find_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1')],
                                                     include_unknown=True)

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'CFSAN002349', 'CFSAN023463'} == set(matches_df['Sample Name'].tolist())
    assert {sample1.id, sample2.id} == set(matches_df['Sample ID'].tolist())
    assert {'lmonocytogenes:abcZ:1'} == set(matches_df['Feature'].tolist())
    assert {'mlst'} == set(matches_df['Type'].tolist())
    assert {'Present', 'Unknown'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 2


def test_find_by_features_empty(mlst_query_service: MLSTQueryService):
    matches_df = mlst_query_service.find_by_features([QueryFeatureMLST('campylobacter:abcZ:20')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())
    assert len(matches_df) == 0


def test_find_by_features_unknown_two_features(database, mlst_query_service_unknown: MLSTQueryService):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    matches_df = mlst_query_service_unknown.find_by_features([
            QueryFeatureMLST('lmonocytogenes:abcZ:1'),
            QueryFeatureMLST('lmonocytogenes:bglA:51'),
        ],
        include_unknown=True)

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert ['CFSAN002349', 'CFSAN023463', 'CFSAN002349', 'CFSAN023463'] == list(matches_df['Sample Name'].tolist())
    assert [sample1.id, sample2.id, sample1.id, sample2.id] == list(matches_df['Sample ID'].tolist())
    assert ['lmonocytogenes:abcZ:1', 'lmonocytogenes:abcZ:1',
            'lmonocytogenes:bglA:51', 'lmonocytogenes:bglA:51'] == list(matches_df['Feature'].tolist())
    assert {'mlst'} == set(matches_df['Type'].tolist())
    assert ['Present', 'Unknown', 'Unknown', 'Unknown'] == list(matches_df['Status'].tolist())
    assert len(matches_df) == 4


def test_find_by_features_two_features(database, mlst_query_service: MLSTQueryService):
    sample1 = database.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == '2014D-0068').one()
    sample3 = database.get_session().query(Sample).filter(Sample.name == '2014C-3598').one()
    sample4 = database.get_session().query(Sample).filter(Sample.name == '2014C-3599').one()

    matches_df = mlst_query_service.find_by_features([QueryFeatureMLST('campylobacter:aspA:2'),
                                                      QueryFeatureMLST('ecoli:adk:100')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'2014D-0067', '2014D-0068', '2014C-3598', '2014C-3599'} == set(matches_df['Sample Name'].tolist())
    assert {sample1.id, sample2.id, sample3.id, sample4.id} == set(matches_df['Sample ID'].tolist())
    assert {'campylobacter:aspA:2', 'ecoli:adk:100'} == set(matches_df['Feature'].tolist())
    assert {'mlst'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 4


def test_count_by_features(mlst_query_service: MLSTQueryService):
    matches_df = mlst_query_service.count_by_features([QueryFeatureMLST('campylobacter:aspA:2')],
                                                      include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['campylobacter:aspA:2'] == list(matches_df['Feature'].tolist())
    assert [2] == list(matches_df['Present'].tolist())
    assert [0] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert [2] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 2 / 2, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])


def test_count_by_features_empty(mlst_query_service: MLSTQueryService):
    matches_df = mlst_query_service.count_by_features([QueryFeatureMLST('campylobacter:aspA:20')],
                                                      include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['campylobacter:aspA:20'] == list(matches_df['Feature'].tolist())
    assert [0] == list(matches_df['Present'].tolist())
    assert [2] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert [2] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 0 / 2, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 2 / 2, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])


def test_count_by_features_unknown(mlst_query_service_unknown: MLSTQueryService):
    matches_df = mlst_query_service_unknown.count_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1')],
                                                      include_unknown=True)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['lmonocytogenes:abcZ:1'] == list(matches_df['Feature'].tolist())
    assert [1] == list(matches_df['Present'].tolist())
    assert [0] == list(matches_df['Absent'].tolist())
    assert [1] == list(matches_df['Unknown'].tolist())
    assert [2] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 1 / 2, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 1 / 2, matches_df['% Unknown'].tolist()[0])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[0])


def test_count_by_features_two_results(mlst_query_service: MLSTQueryService):
    matches_df = mlst_query_service.count_by_features([QueryFeatureMLST('campylobacter:aspA:2'),
                                                       QueryFeatureMLST('ecoli:adk:100')],
                                                      include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['campylobacter:aspA:2', 'ecoli:adk:100'] == list(matches_df['Feature'].tolist())
    assert [2, 2] == list(matches_df['Present'].tolist())
    assert [0, 0] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert pd.isna(matches_df['Unknown'].tolist()[1])
    assert [2, 2] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 2 / 2, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])

    assert math.isclose(100 * 2 / 2, matches_df['% Present'].tolist()[1])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[1])
    assert pd.isna(matches_df['% Unknown'].tolist()[1])


def test_count_by_features_two_results_same_scheme(mlst_query_service: MLSTQueryService):
    matches_df = mlst_query_service.count_by_features([QueryFeatureMLST('campylobacter:aspA:2'),
                                                       QueryFeatureMLST('campylobacter:glyA:3')],
                                                      include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['campylobacter:aspA:2', 'campylobacter:glyA:3'] == list(matches_df['Feature'].tolist())
    assert [2, 2] == list(matches_df['Present'].tolist())
    assert [0, 0] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert pd.isna(matches_df['Unknown'].tolist()[1])
    assert [2, 2] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 2 / 2, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])

    assert math.isclose(100 * 2 / 2, matches_df['% Present'].tolist()[1])
    assert math.isclose(100 * 0 / 2, matches_df['% Absent'].tolist()[1])
    assert pd.isna(matches_df['% Unknown'].tolist()[1])
