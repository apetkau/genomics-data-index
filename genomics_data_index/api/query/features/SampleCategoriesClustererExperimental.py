import ete3
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import skbio
import skbio.tree
from ete3 import Tree
from scipy.spatial.distance import jensenshannon
from scipy.spatial.distance import squareform


class SampleCategoriesClustererExperimental:
    """
    This class contains experimental methods for clustering sample categories (e.g., lineage or other
    metadata) by the features found in these categories. This is called experimental since I have
    not had time to properly test these out yet. Eventually they may make it into the rest of the API.
    """

    def __init__(self):
        pass


def cluster_features_comparison(comparison_df: pd.DataFrame) -> Tree:
    """
    Clusters vectors of feature proportions in particular categories by their Jensen-Shannon distance.
    That is, takes as input a dataframe generated from "SamplesQuery.features_comparison(..., unit='proportion')"
    (which compares the proportion of genomes with a particular feature in some category such as a linage). Produces a
    pair-wise distance matrix from the distances between these feature vectors and clusters this distance matrix
    (using single-linkage clustering).
    :param comparison_df: The dataframe produced as output from "SamplesQuery.features_comparison(..., unit='proportion')"
    :return: An ete3.ClusterTree which stores the hierarchical cluster of each of the sample categories.
    """
    categories = [c[:-len('_proportion')] for c in comparison_df]

    # Create 2d-array of distances
    distances = []
    for i1, c1 in enumerate(categories):
        row = []
        for i2, c2 in enumerate(categories):
            vec1 = comparison_df[f'{c1}_proportion'].to_numpy()
            vec2 = comparison_df[f'{c2}_proportion'].to_numpy()
            row.append(jensenshannon(vec1, vec2))
        distances.append(row)

    distance_matrix = np.array(distances)
    compressed_matrix = squareform(distance_matrix)
    linkage_matrix = sch.linkage(compressed_matrix, method='single')
    skb_tree = skbio.tree.TreeNode.from_linkage_matrix(linkage_matrix, id_list=categories)
    skb_newick = str(skb_tree).strip().replace("'", "")

    return ete3.ClusterNode(newick=skb_newick)
