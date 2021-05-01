from __future__ import annotations

import copy
import logging
from typing import List, Dict, Any, Union, Iterable

from ete3 import Tree, NodeStyle, TreeStyle, CircleFace, TextFace, RectFace, Face

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

    ANNOTATE_KINDS = ['circle', 'rect', 'rectangle']

    def __init__(self, tree: Tree, default_highlight_styles: List[Dict[str, Any]],  annotate_column: int,
                 tree_style: TreeStyle,
                 legend_nsize: int = 10, legend_fsize: int = 11,
                 annotate_color_present: str = '#41ae76',
                 annotate_color_absent: str = 'white',
                 annotate_outline_color: str = 'black',
                 annotate_kind: str = 'rect'):
        self._tree = tree
        self._default_highlight_styles = default_highlight_styles
        self._tree_style = tree_style
        self._legend_nsize = legend_nsize
        self._legend_fsize = legend_fsize
        self._annotate_outline_color = annotate_outline_color
        self._annotate_color_present = annotate_color_present
        self._annotate_color_absent = annotate_color_absent
        self._annotate_column = annotate_column

        if annotate_kind not in self.ANNOTATE_KINDS:
            raise Exception(f'Invalid value for annotate_kind={annotate_kind}.'
                            f' Must be one of {self.ANNOTATE_KINDS}')
        self._annotate_kind = annotate_kind

    def _build_annotate_face(self, width: int, height: int, fgcolor: str, bgcolor: str) -> Face:
        if self._annotate_kind == 'rect' or self._annotate_kind == 'rectangle':
            return RectFace(width=width, height=height, fgcolor=fgcolor, bgcolor=bgcolor)
        elif self._annotate_kind == 'circle':
            # Make circle radius such that it fits in bounding box defined by width and height
            min_dimension = min(width, height)
            radius = min_dimension/2

            return CircleFace(radius=radius, color=bgcolor)
        else:
            raise Exception(f'Invalid value for annotate_kind={self._annotate_kind}.'
                            f' Must be one of {self.ANNOTATE_KINDS}')

    def annotate(self, samples: Union[SamplesQuery, Iterable[str]],
                 annotate_color_present: str = None, annotate_color_absent: str = None) -> TreeStyler:
        """
        Adds an annotation column beside the tree showing the which samples are in the passed set.
        :param samples: The samples to show as being present.
        :param annotate_color_present: The color to use when a sample is present in this set (defaults class-defined color).
        :param annotate_color_absent: The color to use when a sample is absent (defaults to class-defined color).
        :return: A new TreeStyler object which contains the completed annotation column.
        """
        # Add legend item
        # ts = copy.deepcopy(self._tree_style)
        # if legend_label is not None:
        #     ts.legend.add_face(CircleFace(radius=self._legend_nsize / 2, color=legend_color), column=0)
        #     ts.legend.add_face(TextFace(legend_label, fsize=self._legend_fsize), column=1)

        if annotate_color_absent is None:
            annotate_color_absent = self._annotate_color_absent
        if annotate_color_present is None:
            annotate_color_present = self._annotate_color_present

        if isinstance(samples, SamplesQuery):
            sample_names = set(samples.tolist(names=True))
        else:
            sample_names = set(samples)

        face_width = 40
        face_height = 30

        # Annotate nodes
        tree = self._tree.copy(method='deepcopy')
        for leaf in tree.iter_leaves():
            if leaf.name in sample_names:
                annotate_face = self._build_annotate_face(width=face_width, height=face_height,
                                                          fgcolor=self._annotate_outline_color,
                                                          bgcolor=annotate_color_present)
            else:
                annotate_face = self._build_annotate_face(width=face_width, height=face_height,
                                                          fgcolor=self._annotate_outline_color,
                                                          bgcolor=annotate_color_absent)

            leaf.add_face(annotate_face, column=self._annotate_column, position='aligned')

        return TreeStyler(tree, default_highlight_styles=self._default_highlight_styles,
                          tree_style=self._tree_style, legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize,
                          annotate_column=self._annotate_column + 1,
                          annotate_color_present=self._annotate_color_present,
                          annotate_color_absent=self._annotate_color_absent,
                          annotate_outline_color=self._annotate_outline_color,
                          annotate_kind=self._annotate_kind)


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
                          tree_style=ts, legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize,
                          annotate_column=self._annotate_column,
                          annotate_color_present=self._annotate_color_present,
                          annotate_color_absent=self._annotate_color_absent,
                          annotate_outline_color=self._annotate_outline_color,
                          annotate_kind=self._annotate_kind)

    def render(self, file_name: str = '%%inline', w: int = None, h: int = None,
               tree_style: TreeStyle = None, units: str = 'px', dpi: int = 90):
        if tree_style is None:
            tree_style = self._tree_style

        # Set default width if no width or height specified
        # Do this here instead of as a default method value since at least one of
        # width or height must be None to preserve the correct aspect ratio (so I only
        # want a default value for width if neither width nor height are set.
        if w is None and h is None:
            w = 400

        return self._tree.render(file_name=file_name, w=w, h=h, tree_style=tree_style,
                                 units=units, dpi=dpi)

    @property
    def tree(self) -> Tree:
        return copy.deepcopy(self._tree)

    @property
    def tree_style(self) -> TreeStyle:
        return copy.deepcopy(self._tree_style)
