import pytest

from storage.variant.service.TreeService import TreeService


@pytest.fixture
def tree_service(database) -> TreeService:
    return TreeService(database)


def test_build_tree(tree_service, core_alignment_service):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree = tree_service.build_tree(alignment)

    assert {'SampleA', 'SampleB', 'SampleC', 'reference'} == set(tree.get_leaf_names())


def test_build_tree_two_samples(tree_service, core_alignment_service):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB'])

    tree = tree_service.build_tree(alignment)

    assert {'SampleA', 'SampleB', 'reference'} == set(tree.get_leaf_names())
