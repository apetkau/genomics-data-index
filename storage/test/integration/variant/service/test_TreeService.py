import pytest
from ete3 import Tree

from storage.test.integration.variant import tree_file
from storage.variant.service.TreeService import TreeService


@pytest.fixture
def tree_service(database) -> TreeService:
    return TreeService(database)


@pytest.fixture
def expected_tree() -> Tree:
    return Tree(str(tree_file))


def test_build_tree_core_fasttree(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='fasttree')

    assert {'SampleA', 'SampleB', 'SampleC', 'reference'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0


def test_build_tree_core_iqtree(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree',
                                        extra_params='--seed 42 -m GTR+ASC')

    assert {'SampleA', 'SampleB', 'SampleC', 'reference'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * tree.get_distance('SampleA', 'reference')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 5
    print(actual_distance)


def test_build_tree_core_iqtree_2cores(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree', num_cores=2)

    assert {'SampleA', 'SampleB', 'SampleC', 'reference'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0


def test_build_tree_full(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'], align_type='full')

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree',
                                        extra_params='--seed 42 -m GTR')

    assert {'SampleA', 'SampleB', 'SampleC', 'reference'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * tree.get_distance('SampleA', 'reference')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 5


def test_build_tree_two_samples(tree_service, core_alignment_service):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB'])

    tree, out = tree_service.build_tree(alignment)

    assert {'SampleA', 'SampleB', 'reference'} == set(tree.get_leaf_names())
