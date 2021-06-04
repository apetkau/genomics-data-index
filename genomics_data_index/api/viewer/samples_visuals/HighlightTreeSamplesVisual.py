import logging
from typing import Union, Iterable

from ete3 import Tree, TreeStyle, NodeStyle

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual

logger = logging.getLogger(__name__)


class HighlightTreeSamplesVisual(TreeSamplesVisual):
    """
    Class used to highlight the set of samples with a particular color/node style.
    """

    def __init__(self, samples: Union[SamplesQuery, Iterable[str]],
                 node_style: NodeStyle,
                 legend_color: str, legend_label: str,
                 legend_nodesize: int,
                 legend_fontsize: int):
        """
        Builds a new HighlightTreeSamplesVisual with the given information.
        :param samples: The set of samples to highlight.
        :param node_style: The style of the nodes in the tree for the set of samples.
        :param legend_color: The color of the legend item.
        :param legend_label: The legend label.
        :param legend_nodesize: The size of the legend node.
        :param legend_fontsize: The size of the legend font.
        :return: A new HighlightTreeSamplesVisual.
        """
        super().__init__(samples=samples, legend_nodesize=legend_nodesize,
                         legend_fontsize=legend_fontsize)
        self._node_style = node_style
        self._legend_color = legend_color
        self._legend_label = legend_label

    def apply_visual(self, tree: Tree, tree_style: TreeStyle) -> None:
        self._add_legend_entry(self._legend_label, legend_color=self._legend_color,
                               kind='rect', tree_style=tree_style)

        for name in self.sample_names:
            nodes = tree.get_leaves_by_name(name)
            if len(nodes) == 0:
                logger.warning(f'Could not find sample=[{name}] in tree. Not highlighting.')
            elif len(nodes) > 1:
                raise Exception(f'More than one node in the tree matched sample=[{name}]')
            else:
                node = nodes[0]
                node.set_style(self._node_style)
