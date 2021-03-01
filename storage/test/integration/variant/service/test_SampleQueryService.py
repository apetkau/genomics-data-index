import pytest

from storage.variant.service.SampleQueryService import SampleQueryService
from storage.variant.service.SampleQueryService import QueryFeatureMutation


@pytest.fixture
def sample_query_service(tree_service_with_tree_stored,
                         reference_service_with_data, sample_service, variation_service) -> SampleQueryService:
    return SampleQueryService(tree_service=tree_service_with_tree_stored,
                              reference_service=reference_service_with_data,
                              sample_service=sample_service)


def test_find_matchesC(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches(['SampleC'])

    assert ['Type', 'Reference Genome', 'Sample A', 'Sample B', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleC'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['reference', 'SampleB', 'SampleA'] == matches_df['Sample B'].tolist()


def test_find_matchesB(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches(['SampleB'])

    assert ['Type', 'Reference Genome', 'Sample A', 'Sample B', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleB'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['SampleC', 'reference', 'SampleA'] == matches_df['Sample B'].tolist()


def test_find_matchesA(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches(['SampleA'])

    assert ['Type', 'Reference Genome', 'Sample A', 'Sample B', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert {'SampleA'} == set(matches_df['Sample A'])
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['reference', 'SampleC', 'SampleB'] == matches_df['Sample B'].tolist()


def test_find_matchesAB(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_matches(['SampleA', 'SampleB'])

    assert ['Type', 'Reference Genome', 'Sample A', 'Sample B', 'Distance', 'Distance (subs/site)',
            'SNV Alignment Length'] == list(matches_df.columns)
    assert {58} == set(matches_df['SNV Alignment Length'].tolist())
    assert ['SampleA', 'SampleA', 'SampleA',
            'SampleB', 'SampleB', 'SampleB'] == matches_df['Sample A'].tolist()
    assert {'genome'} == set(matches_df['Reference Genome'])
    assert ['reference', 'SampleC', 'SampleB',
            'SampleC', 'reference', 'SampleA'] == matches_df['Sample B'].tolist()


def test_find_by_features(sample_query_service: SampleQueryService):
    matches_df = sample_query_service.find_by_features([QueryFeatureMutation('reference:5061:G:A')])

    assert ['Type', 'Feature', 'Sample Name', 'Sample ID'] == list(matches_df.columns.tolist())

    assert {'SampleB'} == set(matches_df['Sample Name'].tolist())
    assert {2} == set(matches_df['Sample ID'].tolist())
    assert {'reference:5061:G:A'} == set(matches_df['Feature'].tolist())
    assert {'SNV'} == set(matches_df['Type'].tolist())
    assert len(matches_df) == 1
