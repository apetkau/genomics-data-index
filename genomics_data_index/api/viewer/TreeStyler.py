from __future__ import annotations

import copy
import logging
from typing import List, Dict, Any, Union, Iterable, Tuple

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


class TreeStyler:
    MODES = ['r', 'c']
    DEFAULT_HIGHLIGHT_STYLES = [style1, style2, style3]
    ANNOTATE_KINDS = ['circle', 'rect', 'rectangle']

    def __init__(self, tree: Tree, default_highlight_styles: List[Dict[str, Any]], annotate_column: int,
                 tree_style: TreeStyle,
                 legend_nsize: int = 10, legend_fsize: int = 11,
                 annotate_color_present: str = '#66c2a4',
                 annotate_color_absent: str = 'white',
                 annotate_border_color: str = 'black',
                 annotate_kind: str = 'rect',
                 annotate_box_width: int = 30,
                 annotate_box_height: int = 30,
                 annotate_border_width: int = 1,
                 annotate_margin: int = 0,
                 annotate_show_box_label: bool = False,
                 annotate_box_label_color: str = 'white'):
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
        self._annotate_show_box_label = annotate_show_box_label
        self._annotate_box_label_color = annotate_box_label_color

        if annotate_kind not in self.ANNOTATE_KINDS:
            raise Exception(f'Invalid value for annotate_kind={annotate_kind}.'
                            f' Must be one of {self.ANNOTATE_KINDS}')
        self._annotate_kind = annotate_kind

    def _build_legend_item(self, color: str, legend_label: str) -> Tuple[Face, Face]:
        cf = CircleFace(radius=self._legend_nsize / 2, color=color)
        cf.hz_align=2
        tf = TextFace(legend_label, fsize=self._legend_fsize)
        tf.margin_left = 10
        tf.margin_right = 10
        return cf, tf

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
            rf.hz_align = 1
            rf.vt_align = 1
            return rf
        elif self._annotate_kind == 'circle':
            # Make circle radius such that it fits in bounding box defined by width and height
            # Shrink a bit since I noticed it was being clipped slightly
            min_dimension = min(width, height)
            radius = min_dimension / 2

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
                 label: Union[str, Dict[str, Any]] = None,
                 annotate_show_box_label: bool = None, annotate_box_label_color: str = None,
                 legend_label: str = None,
                 box_width: int = None, box_height: int = None,
                 color_present: str = None, color_absent: str = None) -> TreeStyler:
        """
        Adds an annotation column beside the tree showing the which samples are in the passed set.
        :param samples: The samples to show as being present.
        :param label: A label for the column. Can be text or dict with attributes text, font, color, and fontsize
                      (this is passed to the underlying ete3 Face).
        :param annotate_show_box_label: Whether or not the label should be added to every present item in the figure.
        :param annotate_box_label_color: If 'annotate_show_box_label' is set, this defines the color (overrides class setting).
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
        if isinstance(label, str):
            # Pick default color since ete3 by default colors the same as what I'm using for the fill color
            label = {'text': label, 'color': 'black', 'font': 'Verdana', 'fontsize': 14}

        if (label is not None) and (self._annotate_show_box_label or annotate_show_box_label):
            label_present = copy.deepcopy(label)

            if annotate_box_label_color is not None:
                label_present['color'] = annotate_box_label_color
            else:
                label_present['color'] = self._annotate_box_label_color
        else:
            label_present = None

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

        ts = copy.deepcopy(self._tree_style)
        if label is not None:
            text = label.get('text', None)
            fsize = label.get('fontsize', 12)
            ftype = label.get('font', 'Verdana')
            color = label.get('color', 'black')
            tf = TextFace(text, fsize=fsize, ftype=ftype, fgcolor=color)
            tf.margin_bottom = 10
            tf.hz_align = 1
            ts.aligned_header.add_face(tf, self._annotate_column)

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
            color_face, text_face = self._build_legend_item(color=color_present, legend_label=legend_label)
            ts.legend.add_face(color_face, column=0)
            ts.legend.add_face(text_face, column=1)

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
                          annotate_margin=self._annotate_margin,
                          annotate_show_box_label=self._annotate_show_box_label,
                          annotate_box_label_color=self._annotate_box_label_color)

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
        if legend_label is not None:
            ts = copy.deepcopy(self._tree_style)
            color_face, text_face = self._build_legend_item(color=legend_color, legend_label=legend_label)
            ts.legend.add_face(color_face, column=0)
            ts.legend.add_face(text_face, column=1)
        else:
            ts = self._tree_style

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
                          annotate_margin=self._annotate_margin,
                          annotate_show_box_label=self._annotate_show_box_label,
                          annotate_box_label_color=self._annotate_box_label_color)

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

    @classmethod
    def create(cls, tree: Tree,
               initial_style: TreeStyle = None,
               mode='r',
               highlight_styles=None,
               legend_nsize: int = 10, legend_fsize: int = 11,
               annotate_color_present: str = 'black',
               annotate_color_absent: str = 'white',
               annotate_border_color: str = 'black',
               annotate_kind: str = 'rect',
               annotate_box_width: int = None,
               annotate_box_height: int = None,
               annotate_border_width: int = 1,
               annotate_margin: int = 0,
               annotate_guiding_lines: bool = True,
               annotate_guiding_lines_color: str = 'gray',
               figure_margin: int = None,
               show_border: bool = True,
               title: str = None,
               title_fsize: int = 16,
               legend_title: str = None,
               annotate_show_box_label: bool = False,
               annotate_box_label_color: str = 'white') -> TreeStyler:
        if initial_style is not None:
            tree_style_elements = {'mode': mode, 'annotate_guiding_lines': annotate_guiding_lines,
                                   'figure_margin': figure_margin, 'show_border': show_border,
                                   'title': title, 'legend_title': legend_title, 'title_fsize': title_fsize}
            tree_style_elements_none = [tree_style_elements[k] is None for k in tree_style_elements]

            if any(tree_style_elements_none):
                logger.warning(f'Both initial_style=[{initial_style}] and one of parameters {tree_style_elements.keys()}'
                               f' are set. Will ignore these listed parameters.')
            ts = initial_style
        else:
            ts = TreeStyle()

            if mode is not None and mode not in cls.MODES:
                raise Exception(f'Invalid value mode=[{mode}]. Must be one of {cls.MODES}')
            elif mode is not None:
                ts.mode = mode

            ts.draw_guiding_lines = annotate_guiding_lines
            ts.guiding_lines_color = annotate_guiding_lines_color
            ts.show_border = show_border

            if figure_margin is None:
                figure_margin = 10

            ts.margin_top = figure_margin
            ts.margin_bottom = figure_margin
            ts.margin_left = figure_margin
            ts.margin_right = figure_margin

        # Setup default box width/height here
        if annotate_box_height is None and annotate_box_width is None:
            annotate_box_height = 30
            annotate_box_width = 30
        elif annotate_box_height is None:
            annotate_box_height = annotate_box_width
        elif annotate_box_width is None:
            annotate_box_width = annotate_box_height

        if highlight_styles is None:
            highlight_styles = cls.DEFAULT_HIGHLIGHT_STYLES

        if legend_title is not None:
            margin_bottom = 10
            cf = RectFace(width=10, height=10, fgcolor=None, bgcolor=None)
            cf.margin_bottom = margin_bottom
            ts.legend.add_face(cf, column=0)
            tf = TextFace(legend_title, fsize=legend_fsize + 3)
            tf.margin_bottom = margin_bottom
            ts.legend.add_face(tf, column=1)

        if title is not None:
            tf = TextFace(title, fsize=title_fsize)
            ts.title.add_face(tf, column=0)

        return TreeStyler(tree=tree,
                          default_highlight_styles=highlight_styles,
                          tree_style=ts,
                          legend_nsize=legend_nsize,
                          legend_fsize=legend_fsize,
                          annotate_column=1,
                          annotate_color_present=annotate_color_present,
                          annotate_color_absent=annotate_color_absent,
                          annotate_border_color=annotate_border_color,
                          annotate_kind=annotate_kind,
                          annotate_box_width=annotate_box_width,
                          annotate_box_height=annotate_box_height,
                          annotate_border_width=annotate_border_width,
                          annotate_margin=annotate_margin,
                          annotate_show_box_label=annotate_show_box_label,
                          annotate_box_label_color=annotate_box_label_color)
