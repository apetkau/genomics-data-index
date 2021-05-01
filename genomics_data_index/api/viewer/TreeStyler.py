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
                 annotate_color_present: str = '#66c2a4',
                 annotate_color_absent: str = 'white',
                 annotate_border_color: str = 'black',
                 annotate_kind: str = 'rect',
                 annotate_box_width: int = 30,
                 annotate_box_height: int = 30,
                 annotate_border_width: int = 1,
                 annotate_margin: int = 0):
        self._tree = tree
        self._default_highlight_styles = default_highlight_styles
        self._tree_style = tree_style
        self._legend_nsize = legend_nsize
        self._legend_fsize = legend_fsize
        self._annotate_border_color = annotate_border_color
        self._annotate_color_present = annotate_color_present
        self._annotate_color_absent = annotate_color_absent
        self._annotate_column = annotate_column
        self._annotate_box_width = annotate_box_width
        self._annotate_box_height = annotate_box_height
        self._annotate_border_width = annotate_border_width
        self._annotate_margin = annotate_margin

        if annotate_kind not in self.ANNOTATE_KINDS:
            raise Exception(f'Invalid value for annotate_kind={annotate_kind}.'
                            f' Must be one of {self.ANNOTATE_KINDS}')
        self._annotate_kind = annotate_kind

    def _build_annotate_face(self, width: int, height: int, border_color: str, bgcolor: str,
                             label: Union[str, Dict[str, Any]] = None) -> Face:
        if self._annotate_kind == 'rect' or self._annotate_kind == 'rectangle':
            rf = RectFace(width=width, height=height, fgcolor=None, bgcolor=bgcolor, label=label)
            rf.border.width = self._annotate_border_width
            rf.margin_top = self._annotate_margin
            rf.margin_bottom = self._annotate_margin
            rf.margin_left = self._annotate_margin
            rf.margin_right = self._annotate_margin
            rf.border.color = border_color
            return rf
        elif self._annotate_kind == 'circle':
            # Make circle radius such that it fits in bounding box defined by width and height
            # Shrink a bit since I noticed it was being clipped slightly
            min_dimension = min(width, height)
            radius = min_dimension/2

            cf = CircleFace(radius=radius, color=bgcolor, label=label)
            cf.border.width = self._annotate_border_width
            cf.margin_top = self._annotate_margin
            cf.margin_bottom = self._annotate_margin
            cf.margin_left = self._annotate_margin
            cf.margin_right = self._annotate_margin
            cf.border.color = border_color
            return cf
        else:
            raise Exception(f'Invalid value for annotate_kind={self._annotate_kind}.'
                            f' Must be one of {self.ANNOTATE_KINDS}')

    def annotate(self, samples: Union[SamplesQuery, Iterable[str]],
                 label_present: Union[str, Dict[str, Any]] = None, legend_label: str = None,
                 box_width: int = None, box_height: int = None,
                 color_present: str = None, color_absent: str = None) -> TreeStyler:
        """
        Adds an annotation column beside the tree showing the which samples are in the passed set.
        :param samples: The samples to show as being present.
        :param label_present: An optional label to display for any present items. Can be text or dict
                              with  attributes text, font, color, and fontsize (this is passed to the underlying ete3 Face)
        :param legend_label: A label to use for a legend item. The color should match the color for this annotated set.
        :param box_width: The width of the bounding box (defaults to class variable annotate_box_width).
        :param box_height: The height of the bounding box (defaults to class variable annotate_box_height).
        :param color_present: The color to use when a sample is present in this set (defaults class-defined color).
        :param color_absent: The color to use when a sample is absent (defaults to class-defined color).
        :return: A new TreeStyler object which contains the completed annotation column.
        """
        if color_absent is None:
            color_absent = self._annotate_color_absent
        if color_present is None:
            color_present = self._annotate_color_present
        if isinstance(label_present, str):
            # Pick default color since ete3 by default colors the same as what I'm using for the fill color
            label_present = {'text': label_present, 'color': 'black'}

        if isinstance(samples, SamplesQuery):
            sample_names = set(samples.tolist(names=True))
        else:
            sample_names = set(samples)

        if box_width is None:
            face_width = self._annotate_box_width
        else:
            face_width = box_width

        if box_height is None:
            face_height = self._annotate_box_height
        else:
            face_height = box_height

        # Annotate nodes
        tree = self._tree.copy(method='deepcopy')
        for leaf in tree.iter_leaves():
            if leaf.name in sample_names:
                annotate_face = self._build_annotate_face(width=face_width, height=face_height,
                                                          border_color=self._annotate_border_color,
                                                          bgcolor=color_present, label=label_present)
            else:
                annotate_face = self._build_annotate_face(width=face_width, height=face_height,
                                                          border_color=self._annotate_border_color,
                                                          bgcolor=color_absent, label=None)

            leaf.add_face(annotate_face, column=self._annotate_column, position='aligned')

        # Add legend item
        if legend_label is not None:
            ts = copy.deepcopy(self._tree_style)
            ts.legend.add_face(CircleFace(radius=self._legend_nsize / 2, color=color_present), column=0)
            tf = TextFace(legend_label, fsize=self._legend_fsize)
            tf.margin_left = 10
            tf.margin_right = 10
            ts.legend.add_face(tf, column=1)
        else:
            ts = self._tree_style

        return TreeStyler(tree, default_highlight_styles=self._default_highlight_styles,
                          tree_style=ts, legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize,
                          annotate_column=self._annotate_column + 1,
                          annotate_color_present=self._annotate_color_present,
                          annotate_color_absent=self._annotate_color_absent,
                          annotate_border_color=self._annotate_border_color,
                          annotate_kind=self._annotate_kind,
                          annotate_box_width=self._annotate_box_width,
                          annotate_box_height=self._annotate_box_height,
                          annotate_border_width=self._annotate_border_width,
                          annotate_margin=self._annotate_margin)


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
            tf = TextFace(legend_label, fsize=self._legend_fsize)
            tf.margin_left = 10
            tf.margin_left = 10
            ts.legend.add_face(tf, column=1)

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
                          annotate_border_color=self._annotate_border_color,
                          annotate_kind=self._annotate_kind,
                          annotate_box_height=self._annotate_box_height,
                          annotate_box_width=self._annotate_box_width,
                          annotate_border_width=self._annotate_border_width,
                          annotate_margin=self._annotate_margin)

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
