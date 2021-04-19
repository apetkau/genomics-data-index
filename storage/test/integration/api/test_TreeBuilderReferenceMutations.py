from ete3 import Tree

from storage.test.integration import tree_file
from storage.api.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations

expected_tree = Tree(str(tree_file))


def test_build_tree(loaded_database_connection):
    tree_builder = TreeBuilderReferenceMutations(database_connection=loaded_database_connection,
                                                 reference_name='genome')

    actual_tree = tree_builder.build(samples_set=['SampleA', 'SampleB', 'SampleC'],
                                     method='iqtree',
                                     align_type='full',
                                     extra_params='--seed 42 -m GTR'
                                     )

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(actual_tree.get_leaf_names())

    tree_comparison = expected_tree.compare(actual_tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * actual_tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 5


def test_build_tree_core(loaded_database_connection):
    tree_builder = TreeBuilderReferenceMutations(database_connection=loaded_database_connection,
                                                 reference_name='genome')

    actual_tree = tree_builder.build(samples_set=['SampleA', 'SampleB', 'SampleC'],
                                     method='iqtree',
                                     align_type='core',
                                     extra_params='--seed 42 -m GTR+ASC'
                                     )

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(actual_tree.get_leaf_names())

    tree_comparison = expected_tree.compare(actual_tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * actual_tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 5
