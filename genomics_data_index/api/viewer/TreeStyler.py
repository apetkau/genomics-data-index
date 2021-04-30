from __future__ import annotations

import copy
import logging
from typing import List, Dict, Any, Union, Iterable

from ete3 import Tree, NodeStyle, TreeStyle, CircleFace, TextFace

from genomics_data_index.api.query.SamplesQuery import SamplesQuery

logger = logging.getLogger(__name__)

fg_color1 = '#41ae76'
fg_color2 = '#ef6548'
fg_color3 = '#8c6bb1'

nstyle1 = NodeStyle()
nstyle1['fgcolor'] = fg_color1
nstyle1['bgcolor'] = '#e5f5f9'
nstyle1['size'] = 10

style1 = {
    'nstyle': nstyle1,
    'legend_color': fg_color1
}

nstyle2 = NodeStyle()
nstyle2['fgcolor'] = fg_color2
nstyle2['bgcolor'] = '#fee8c8'
nstyle2['size'] = 10

style2 = {
    'nstyle': nstyle2,
    'legend_color': fg_color2
}

nstyle3 = NodeStyle()
nstyle3['fgcolor'] = fg_color3
nstyle3['bgcolor'] = '#e0ecf4'
nstyle3['size'] = 10

style3 = {
    'nstyle': nstyle3,
    'legend_color': fg_color3
}

DEFAULT_HIGHLIGHT_STYLES = [style1, style2, style3]


class TreeStyler:

    def __init__(self, tree: Tree, default_highlight_styles: List[Dict[str, Any]], tree_style: TreeStyle,
                 legend_nsize: int = 10, legend_fsize: int = 11):
        self._tree = tree
        self._default_highlight_styles = default_highlight_styles
        self._tree_style = tree_style
        self._legend_nsize = legend_nsize
        self._legend_fsize = legend_fsize

    def highlight(self, samples: Union[SamplesQuery, Iterable[str]],
                  nstyle: NodeStyle = None, legend_color: str = None,
                  legend_label: str = None) -> TreeStyler:
        if nstyle is None and legend_color is None:
            nstyle = self._default_highlight_styles[0]['nstyle']
            legend_color = self._default_highlight_styles[0]['legend_color']

            # Shift default styles by 1
            new_default_styles = copy.copy(self._default_highlight_styles)
            new_default_styles.append(new_default_styles.pop(0))
        else:
            new_default_styles = self._default_highlight_styles

        if isinstance(samples, SamplesQuery):
            sample_names = samples.tolist(names=True)
            query_expression = samples.query_expression()
        else:
            sample_names = set(samples)
            query_expression = f'set({len(sample_names)} samples)'

        # Add legend item
        ts = copy.deepcopy(self._tree_style)
        if legend_label is not None:
            ts.legend.add_face(CircleFace(radius=self._legend_nsize / 2, color=legend_color), column=0)
            ts.legend.add_face(TextFace(legend_label, fsize=self._legend_fsize), column=1)

        # Highlight nodes
        tree = copy.deepcopy(self._tree)
        for name in sample_names:
            nodes = tree.get_leaves_by_name(name)
            if len(nodes) == 0:
                logger.warning(f'Could not find sample=[{name}] in tree. Not highlighting.')
            elif len(nodes) > 1:
                raise Exception(f'More than one node in the tree matched sample=[{name}]')
            else:
                node = nodes[0]
                node.set_style(nstyle)

        return TreeStyler(tree, default_highlight_styles=new_default_styles,
                          tree_style=ts, legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize)

    def render(self, file_name: str = '%%inline', w: int = 400, h: int = 300,
               tree_style: TreeStyle = None):
        if tree_style is None:
            tree_style = self._tree_style
        return self._tree.render(file_name=file_name, w=w, tree_style=tree_style)

    @property
    def tree(self) -> Tree:
        return copy.deepcopy(self._tree)

    @property
    def tree_style(self) -> TreeStyle:
        return copy.deepcopy(self._tree_style)
