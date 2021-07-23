import logging
from typing import Union, Iterable, Dict

from ete3 import Tree, TreeStyle, NodeStyle

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual

logger = logging.getLogger(__name__)


class HighlightTreeSamplesVisual(TreeSamplesVisual):
    """
    Class used to highlight the set of samples with a particular color/node style.
    """

    def __init__(self, samples: Union[SamplesQuery, Iterable[str]],
                 present_node_style: NodeStyle,
                 unknown_node_style: NodeStyle,
                 include_unknown: bool,
                 legend_color: str, legend_label: str,
                 legend_nodesize: int,
                 legend_fontsize: int,
                 legend_columns: Dict[str, int]):
        """
        Builds a new HighlightTreeSamplesVisual with the given information.
        :param samples: The set of samples to highlight.
        :param present_node_style: The style of the nodes in the tree for the present set of samples.
        :param unknown_node_style: The style of the ndoes in the tree for the unknown set of samples.
        :param include_unknown: Whether or not unknown samples should be included.
        :param legend_color: The color of the legend item.
        :param legend_label: The legend label.
        :param legend_nodesize: The size of the legend node.
        :param legend_fontsize: The size of the legend font.
        :param legend_columns: A dict mapping legend item types to column numbers.
        :return: A new HighlightTreeSamplesVisual.
        """
        super().__init__(samples=samples, legend_nodesize=legend_nodesize,
                         legend_fontsize=legend_fontsize, legend_columns=legend_columns)
        self._present_node_style = present_node_style
        self._unknown_node_style = unknown_node_style
        self._include_unknown = include_unknown
        self._legend_color = legend_color
        self._legend_label = legend_label
        self._unknown_legend_color = self._unknown_node_style['bgcolor']

    def _get_and_set_style(self, sample_name: str, node_style: NodeStyle, tree: Tree) -> None:
        nodes = tree.get_leaves_by_name(sample_name)
        if len(nodes) == 0:
            logger.warning(f'Could not find sample=[{sample_name}] in tree. Not highlighting.')
        elif len(nodes) > 1:
            raise Exception(f'More than one node in the tree matched sample=[{sample_name}]')
        else:
            node = nodes[0]
            node.set_style(node_style)

    def apply_visual(self, tree: Tree, tree_style: TreeStyle) -> None:
        has_unknowns = self._include_unknown and len(self.unknown_sample_names) > 0
        self._add_legend_entry(self._legend_label, legend_color=self._legend_color,
                               unknown_color=self._unknown_legend_color,
                               include_unknown=self._include_unknown,
                               has_unknowns=has_unknowns,
                               kind='rect', tree_style=tree_style)

        for name in self.present_sample_names:
            self._get_and_set_style(name, node_style=self._present_node_style, tree=tree)

        if self._include_unknown:
            for name in self.unknown_sample_names:
                self._get_and_set_style(name, node_style=self._unknown_node_style, tree=tree)
