from __future__ import annotations

import copy
import logging

from ete3 import Tree, NodeStyle

from storage.api.SamplesQuery import SamplesQuery

logger = logging.getLogger(__name__)


class TreeStyler:

    def __init__(self, tree: Tree):
        self._tree = tree

    def highlight(self, samples_query: SamplesQuery, nstyle: NodeStyle) -> TreeStyler:
        sample_names = samples_query.tolist(names=True)

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

        return TreeStyler(tree)

    def inlineplot(self, w: int = 400):
        return self._tree.render("%%inline", w=w)

    @property
    def tree(self) -> Tree:
        return copy.deepcopy(self._tree)
