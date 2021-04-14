import math

import pandas as pd
import pytest

from storage.variant.model import Sample
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.service.MutationQueryService import MutationQueryService


@pytest.fixture
def mutation_query_service(tree_service_with_tree_stored, reference_service_with_data, sample_service,
                           variation_service) -> MutationQueryService:
    return MutationQueryService(reference_service=reference_service_with_data,
                                sample_service=sample_service,
                                tree_service=tree_service_with_tree_stored)


def test_find_matchesC(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.find_matches(['SampleC'])

    assert ['Type', 'Reference Genome', 'Query', 'Match', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleC'} == set(matches_df['Query'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['genome', 'SampleB', 'SampleA'] == matches_df['Match'].tolist()


def test_find_matchesB(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.find_matches(['SampleB'])

    assert ['Type', 'Reference Genome', 'Query', 'Match', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleB'} == set(matches_df['Query'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['SampleC', 'genome', 'SampleA'] == matches_df['Match'].tolist()


def test_find_matchesA(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.find_matches(['SampleA'])

    assert ['Type', 'Reference Genome', 'Query', 'Match', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleA'} == set(matches_df['Query'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['genome', 'SampleC', 'SampleB'] == matches_df['Match'].tolist()


def test_find_matchesAB(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.find_matches(['SampleA', 'SampleB'])

    assert ['Type', 'Reference Genome', 'Query', 'Match', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert ['SampleA', 'SampleA', 'SampleA',
            'SampleB', 'SampleB', 'SampleB'] == matches_df['Query'].tolist()
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['genome', 'SampleC', 'SampleB',
            'SampleC', 'genome', 'SampleA'] == matches_df['Match'].tolist()


def test_find_by_features(database, mutation_query_service: MutationQueryService):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:5061:G:A')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'SampleB'} == set(matches_df['Sample Name'].tolist())
    assert {sampleB.id} == set(matches_df['Sample ID'].tolist())
    assert {'reference:5061:G:A'} == set(matches_df['Feature'].tolist())
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 1


def test_count_by_features(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.count_by_features([QueryFeatureMutation('reference:5061:G:A')],
                                                          include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['reference:5061:G:A'] == list(matches_df['Feature'].tolist())
    assert [1] == list(matches_df['Present'].tolist())
    assert [2] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert [3] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 1 / 3, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 2 / 3, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])


def test_find_by_features_2_results(database, mutation_query_service: MutationQueryService):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:3063:A:ATGCAGC')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'SampleB', 'SampleC'} == set(matches_df['Sample Name'].tolist())
    assert {sampleB.id, sampleC.id} == set(matches_df['Sample ID'].tolist())
    assert {'reference:3063:A:ATGCAGC'} == set(matches_df['Feature'].tolist())
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 2


def test_count_by_features_2_results(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.count_by_features([QueryFeatureMutation('reference:3063:A:ATGCAGC')],
                                                          include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['reference:3063:A:ATGCAGC'] == list(matches_df['Feature'].tolist())
    assert [2] == list(matches_df['Present'].tolist())
    assert [1] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert [3] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 2 / 3, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 1 / 3, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])


def test_find_by_features_2_features(database, mutation_query_service: MutationQueryService):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:5061:G:A'),
                                                          QueryFeatureMutation('reference:3063:A:ATGCAGC')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert ['SampleB', 'SampleC', 'SampleB'] == matches_df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id, sampleB.id] == matches_df['Sample ID'].tolist()
    assert ['reference:3063:A:ATGCAGC', 'reference:3063:A:ATGCAGC', 'reference:5061:G:A'] == matches_df[
        'Feature'].tolist()
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert ['Present', 'Present', 'Present'] == matches_df['Status'].tolist()
    assert len(matches_df) == 3


def test_count_by_features_2_features(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.count_by_features([QueryFeatureMutation('reference:5061:G:A'),
                                                           QueryFeatureMutation('reference:3063:A:ATGCAGC')],
                                                          include_unknown=False)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['reference:3063:A:ATGCAGC', 'reference:5061:G:A'] == list(matches_df['Feature'].tolist())
    assert [2, 1] == list(matches_df['Present'].tolist())
    assert [1, 2] == list(matches_df['Absent'].tolist())
    assert pd.isna(matches_df['Unknown'].tolist()[0])
    assert pd.isna(matches_df['Unknown'].tolist()[1])
    assert [3, 3] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 2 / 3, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 1 / 3, matches_df['% Absent'].tolist()[0])
    assert pd.isna(matches_df['% Unknown'].tolist()[0])

    assert math.isclose(100 * 1 / 3, matches_df['% Present'].tolist()[1])
    assert math.isclose(100 * 2 / 3, matches_df['% Absent'].tolist()[1])
    assert pd.isna(matches_df['% Unknown'].tolist()[1])


def test_find_by_features_unknown(database, mutation_query_service: MutationQueryService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:190:A:G')],
                                                         include_unknown=True)

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert ['SampleA', 'SampleB', 'SampleC'] == matches_df['Sample Name'].tolist()
    assert [sampleA.id, sampleB.id, sampleC.id] == matches_df['Sample ID'].tolist()
    assert {'reference:190:A:G'} == set(matches_df['Feature'].tolist())
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert ['Unknown', 'Present', 'Unknown'] == matches_df['Status'].tolist()
    assert len(matches_df) == 3

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:190:A:G')],
                                                         include_unknown=False)
    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert {'SampleB'} == set(matches_df['Sample Name'].tolist())
    assert {sampleB.id} == set(matches_df['Sample ID'].tolist())
    assert {'reference:190:A:G'} == set(matches_df['Feature'].tolist())
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert {'Present'} == set(matches_df['Status'].tolist())
    assert len(matches_df) == 1


def test_count_by_features_unknown(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.count_by_features([QueryFeatureMutation('reference:190:A:G')],
                                                          include_unknown=True)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['reference:190:A:G'] == list(matches_df['Feature'].tolist())
    assert [1] == list(matches_df['Present'].tolist())
    assert [0] == list(matches_df['Absent'].tolist())
    assert [2] == list(matches_df['Unknown'].tolist())
    assert [3] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 1 / 3, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 0 / 3, matches_df['% Absent'].tolist()[0])
    assert math.isclose(100 * 2 / 3, matches_df['% Unknown'].tolist()[0])


def test_find_by_features_found_unknown(database, mutation_query_service: MutationQueryService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    matches_df = mutation_query_service.find_by_features([QueryFeatureMutation('reference:5061:G:A'),
                                                          QueryFeatureMutation('reference:190:A:G')],
                                                         include_unknown=True)

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID', 'Status'] == list(matches_df.columns.tolist())

    assert ['SampleA', 'SampleB', 'SampleC',
            'SampleA', 'SampleB'] == matches_df['Sample Name'].tolist()
    assert [sampleA.id, sampleB.id, sampleC.id,
            sampleA.id, sampleB.id] == matches_df['Sample ID'].tolist()
    assert ['reference:190:A:G', 'reference:190:A:G', 'reference:190:A:G',
            'reference:5061:G:A', 'reference:5061:G:A'] == matches_df['Feature'].tolist()
    assert {'mutation'} == set(matches_df['Type'].tolist())
    assert ['Unknown', 'Present', 'Unknown',
            'Unknown', 'Present'] == matches_df['Status'].tolist()
    assert len(matches_df) == 5


def test_count_by_features_found_unknown(mutation_query_service: MutationQueryService):
    matches_df = mutation_query_service.count_by_features([QueryFeatureMutation('reference:5061:G:A'),
                                                           QueryFeatureMutation('reference:190:A:G')],
                                                          include_unknown=True)

    assert ['Type', 'Feature', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == list(matches_df.columns.tolist())

    assert ['reference:190:A:G', 'reference:5061:G:A'] == list(matches_df['Feature'].tolist())
    assert [1, 1] == list(matches_df['Present'].tolist())
    assert [0, 1] == list(matches_df['Absent'].tolist())
    assert [2, 1] == list(matches_df['Unknown'].tolist())
    assert [3, 3] == list(matches_df['Total'].tolist())

    assert math.isclose(100 * 1 / 3, matches_df['% Present'].tolist()[0])
    assert math.isclose(100 * 0 / 3, matches_df['% Absent'].tolist()[0])
    assert math.isclose(100 * 2 / 3, matches_df['% Unknown'].tolist()[0])

    assert math.isclose(100 * 1 / 3, matches_df['% Present'].tolist()[1])
    assert math.isclose(100 * 1 / 3, matches_df['% Absent'].tolist()[1])
    assert math.isclose(100 * 1 / 3, matches_df['% Unknown'].tolist()[1])
