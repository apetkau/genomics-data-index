from ete3 import Tree

from genomics_data_index.api.query.impl.TreeBuilderReferenceMutations import TreeBuilderReferenceMutations
from genomics_data_index.test.integration import tree_file

expected_tree = Tree(str(tree_file))


def test_build_tree(loaded_database_connection):
    tree_builder = TreeBuilderReferenceMutations(database_connection=loaded_database_connection,
                                                 reference_name='genome')

    actual_tree, alignment_length, tree_set = tree_builder.build(samples_set=['SampleA', 'SampleB', 'SampleC'],
                                                                 method='iqtree',
                                                                 align_type='full',
                                                                 include_reference=True,
                                                                 extra_params='--seed 42 -m GTR'
                                                                 )

    assert 5180 == alignment_length

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(actual_tree.get_leaf_names())

    tree_comparison = expected_tree.compare(actual_tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * actual_tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 6


def test_build_tree_core(loaded_database_connection):
    tree_builder = TreeBuilderReferenceMutations(database_connection=loaded_database_connection,
                                                 reference_name='genome')

    actual_tree, alignment_length, tree_set = tree_builder.build(samples_set=['SampleA', 'SampleB', 'SampleC'],
                                                                 method='iqtree',
                                                                 align_type='core',
                                                                 include_reference=True,
                                                                 extra_params='--seed 42 -m GTR+ASC'
                                                                 )

    # See test_CoreAlignmentService.py for where my test alignment is coming from and how I get the length
    assert 58 == alignment_length

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(actual_tree.get_leaf_names())

    tree_comparison = expected_tree.compare(actual_tree, unrooted=True)
    assert tree_comparison['rf'] == 0

    actual_distance = 5180 * actual_tree.get_distance('SampleA', 'genome')
    expected_distance = 26
    assert abs(expected_distance - actual_distance) < 6


def test_build_tree_exclude_reference(loaded_database_connection):
    tree_builder = TreeBuilderReferenceMutations(database_connection=loaded_database_connection,
                                                 reference_name='genome')

    actual_tree, alignment_length, tree_set = tree_builder.build(samples_set=['SampleA', 'SampleB', 'SampleC'],
                                                                 method='iqtree',
                                                                 align_type='full',
                                                                 include_reference=False,
                                                                 extra_params='--seed 42 -m GTR'
                                                                 )

    assert {'SampleA', 'SampleB', 'SampleC'} == set(actual_tree.get_leaf_names())
