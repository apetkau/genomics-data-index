import pytest
from ete3 import Tree

from genomics_data_index.test.integration import tree_file


@pytest.fixture
def expected_tree() -> Tree:
    return Tree(str(tree_file))


def test_build_tree_core_fasttree(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='fasttree')

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0


def test_build_tree_core_iqtree(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree',
                                        extra_params='--seed 42 -m GTR+ASC')

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 6


def test_build_tree_core_iqtree_2cores(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'])

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree', num_cores=2)

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0


def test_build_tree_full(tree_service, core_alignment_service, expected_tree):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB', 'SampleC'], align_type='full')

    tree, out = tree_service.build_tree(alignment, tree_build_type='iqtree',
                                        extra_params='--seed 42 -m GTR')

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 6


def test_build_tree_two_samples(tree_service, core_alignment_service):
    alignment = core_alignment_service.construct_alignment(
        reference_name='genome', samples=['SampleA', 'SampleB'], align_type='full')

    tree, out = tree_service.build_tree(alignment)

    assert {'SampleA', 'SampleB', 'genome'} == set(tree.get_leaf_names())


def test_rebuild_tree(tree_service, reference_service_with_data, expected_tree):
    tree_service.rebuild_tree(reference_name='genome')
    reference_genome = reference_service_with_data.find_reference_genome('genome')
    tree = reference_genome.tree

    assert 58 == reference_genome.tree_alignment_length
    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(tree.get_leaf_names())

    tree_comparison = expected_tree.compare(tree, unrooted=True)
    assert tree_comparison['rf'] == 0
