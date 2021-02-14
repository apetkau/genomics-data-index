import pytest

from storage.variant.service.SampleQueryService import SampleQueryService


@pytest.fixture
def sample_query_service(tree_service_with_tree_stored,
                         reference_service_with_data, variation_service) -> SampleQueryService:
    return SampleQueryService(tree_service=tree_service_with_tree_stored,
                              reference_service=reference_service_with_data)


def test_find_matchesC(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches('SampleC')

    assert ['Reference Genome', 'Sample A', 'Sample B', 'Distance (subs/site)',
            'Distance (subs)', 'Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['Alignment Length'].tolist())
    assert {'SampleC'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['reference', 'SampleB', 'SampleA'] == matches_df['Sample B'].tolist()


def test_find_matchesB(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches('SampleB')

    assert ['Reference Genome', 'Sample A', 'Sample B', 'Distance (subs/site)',
            'Distance (subs)', 'Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['Alignment Length'].tolist())
    assert {'SampleB'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['SampleC', 'reference', 'SampleA'] == matches_df['Sample B'].tolist()


def test_find_matchesA(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches('SampleA')

    assert ['Reference Genome', 'Sample A', 'Sample B', 'Distance (subs/site)',
            'Distance (subs)', 'Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['Alignment Length'].tolist())
    assert {'SampleA'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['reference', 'SampleC', 'SampleB'] == matches_df['Sample B'].tolist()