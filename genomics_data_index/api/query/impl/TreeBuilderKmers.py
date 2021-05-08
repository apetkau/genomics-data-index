from typing import Tuple

import ete3
import scipy.cluster.hierarchy as sch
import skbio.tree
from ete3 import Tree
from scipy.spatial.distance import squareform

from genomics_data_index.api.query.TreeBuilder import TreeBuilder
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class TreeBuilderKmers(TreeBuilder):
    BUILD_METHODS = ['single-linkage']

    def __init__(self, database_connection: DataIndexConnection):
        super().__init__()
        self._database_connection = database_connection

    def _build_kmer_tree(self, samples_set: SampleSet, linkage_method: str, kmer_size: int,
                         ncores: int = 1) -> Tuple[
        Tree, int, SampleSet]:
        distance_matrix, labels = self._database_connection.kmer_service.get_distance_matrix(
            sample_ids=samples_set,
            kmer_size=kmer_size,
            ncores=ncores
        )

        compressed_matrix = squareform(distance_matrix)
        linkage_matrix = sch.linkage(compressed_matrix, method=linkage_method)

        skb_tree = skbio.tree.TreeNode.from_linkage_matrix(linkage_matrix, id_list=labels)
        skb_newick = str(skb_tree).strip().replace("'", "")
        ete_tree = ete3.ClusterNode(newick=skb_newick)

        return ete_tree, -1, samples_set

    def build(self, samples_set: SampleSet, method: str = 'single-linkage', **kwargs) -> Tuple[Tree, int, SampleSet]:
        if method == 'single-linkage':
            return self._build_kmer_tree(samples_set=samples_set, linkage_method='single',
                                         **kwargs)
        else:
            raise Exception(f'Invalid method=[{method}]. Must be one of {self.BUILD_METHODS}.')
