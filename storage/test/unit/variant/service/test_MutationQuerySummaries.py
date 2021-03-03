import pytest
import pandas as pd
import math

from storage.variant.service.MutationQueryService import MutationQuerySummaries


@pytest.fixture
def mutation_query_summaries() -> MutationQuerySummaries:
    return MutationQuerySummaries()


def test_find_by_features_summary(mutation_query_summaries: MutationQuerySummaries):
    df = pd.DataFrame({
        'Type': ['mutation'],
        'Feature': ['seq:1:A:T'],
        'Sample Name': ['A'],
        'Sample ID': [1],
        'Status': ['Present'],
    })
    sample_counts = {
        'seq:1:A:T': 1,
    }

    summary_df = mutation_query_summaries.find_by_features_summary(sample_counts=sample_counts,
                                                                   find_by_features_results=df,
                                                                   )

    assert ['Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == summary_df.columns.tolist()
    assert ['seq:1:A:T'] == summary_df.index.tolist()
    assert [1] == summary_df['Present'].tolist()
    assert [0] == summary_df['Absent'].tolist()
    assert [0] == summary_df['Unknown'].tolist()
    assert [1] == summary_df['Total'].tolist()
    assert math.isclose(100, summary_df.loc['seq:1:A:T', '% Present'])
    assert math.isclose(0, summary_df.loc['seq:1:A:T', '% Absent'])
    assert math.isclose(0, summary_df.loc['seq:1:A:T', '% Unknown'])


def test_find_by_features_summary_2(mutation_query_summaries: MutationQuerySummaries):
    df = pd.DataFrame({
        'Type': ['mutation', 'mutation'],
        'Feature': ['seq:1:A:T', 'seq:1:A:T'],
        'Sample Name': ['A', 'B'],
        'Sample ID': [1,2],
        'Status': ['Present', 'Present'],
    })
    sample_counts = {
        'seq:1:A:T': 3,
    }

    summary_df = mutation_query_summaries.find_by_features_summary(sample_counts=sample_counts,
                                                                   find_by_features_results=df,
                                                                   )

    assert ['Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == summary_df.columns.tolist()
    assert ['seq:1:A:T'] == summary_df.index.tolist()
    assert [2] == summary_df['Present'].tolist()
    assert [1] == summary_df['Absent'].tolist()
    assert [0] == summary_df['Unknown'].tolist()
    assert [3] == summary_df['Total'].tolist()
    assert math.isclose(100*2/3, summary_df.loc['seq:1:A:T', '% Present'])
    assert math.isclose(100*1/3, summary_df.loc['seq:1:A:T', '% Absent'])
    assert math.isclose(0, summary_df.loc['seq:1:A:T', '% Unknown'])


def test_find_by_features_summary_3(mutation_query_summaries: MutationQuerySummaries):
    df = pd.DataFrame({
        'Type': ['mutation', 'mutation', 'mutation'],
        'Feature': ['seq:1:A:T', 'seq:1:A:T', 'seq:1:A:T'],
        'Sample Name': ['A', 'B', 'C'],
        'Sample ID': [1,2,3],
        'Status': ['Unknown', 'Present', 'Unknown'],
    })
    sample_counts = {
        'seq:1:A:T': 5,
    }

    summary_df = mutation_query_summaries.find_by_features_summary(sample_counts=sample_counts,
                                                                   find_by_features_results=df,
                                                                   )

    assert ['Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == summary_df.columns.tolist()
    assert ['seq:1:A:T'] == summary_df.index.tolist()
    assert [1] == summary_df['Present'].tolist()
    assert [2] == summary_df['Absent'].tolist()
    assert [2] == summary_df['Unknown'].tolist()
    assert [5] == summary_df['Total'].tolist()
    assert math.isclose(100*1/5, summary_df.loc['seq:1:A:T', '% Present'])
    assert math.isclose(100*2/5, summary_df.loc['seq:1:A:T', '% Absent'])
    assert math.isclose(100*2/5, summary_df.loc['seq:1:A:T', '% Unknown'])


def test_find_by_features_summary_two_mutations(mutation_query_summaries: MutationQuerySummaries):
    df = pd.DataFrame({
        'Type': ['mutation', 'mutation', 'mutation', 'mutation'],
        'Feature': ['seq:1:A:T', 'seq:1:A:T', 'seq:2:C:G', 'seq:2:C:G'],
        'Sample Name': ['A', 'B', 'A', 'C'],
        'Sample ID': [1, 2, 1, 3],
        'Status': ['Present', 'Present', 'Unknown', 'Present'],
    })
    sample_counts = {
        'seq:1:A:T': 5,
        'seq:2:C:G': 5,
    }

    summary_df = mutation_query_summaries.find_by_features_summary(sample_counts=sample_counts,
                                                                   find_by_features_results=df,
                                                                   )

    assert ['Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == summary_df.columns.tolist()
    assert ['seq:1:A:T', 'seq:2:C:G'] == summary_df.index.tolist()
    assert [2, 1] == summary_df['Present'].tolist()
    assert [3, 3] == summary_df['Absent'].tolist()
    assert [0, 1] == summary_df['Unknown'].tolist()
    assert [5, 5] == summary_df['Total'].tolist()

    assert math.isclose(100*2/5, summary_df.loc['seq:1:A:T', '% Present'])
    assert math.isclose(100*3/5, summary_df.loc['seq:1:A:T', '% Absent'])
    assert math.isclose(100*0/5, summary_df.loc['seq:1:A:T', '% Unknown'])

    assert math.isclose(100*1/5, summary_df.loc['seq:2:C:G', '% Present'])
    assert math.isclose(100*3/5, summary_df.loc['seq:2:C:G', '% Absent'])
    assert math.isclose(100*1/5, summary_df.loc['seq:2:C:G', '% Unknown'])
