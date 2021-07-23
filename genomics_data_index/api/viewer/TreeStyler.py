from __future__ import annotations

import copy
import logging
from typing import List, Dict, Any, Union, Iterable

from ete3 import Tree, NodeStyle, TreeStyle, TextFace, RectFace

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual
from genomics_data_index.api.viewer.samples_visuals.AnnotateTreeSamplesVisual import AnnotateTreeSamplesVisual
from genomics_data_index.api.viewer.samples_visuals.HighlightTreeSamplesVisual import HighlightTreeSamplesVisual

logger = logging.getLogger(__name__)


class TreeStyler:
    """
    A class used to style and render a tree and sets of samples on this tree.
    """

    MODES = ['r', 'c']

    def __init__(self, tree: Tree, default_highlight_styles: HighlightStyle, annotate_column: int,
                 tree_style: TreeStyle,
                 legend_columns: Dict[str, int],
                 samples_styles_list: List[TreeSamplesVisual],
                 legend_nsize: int = 20, legend_fsize: int = 11,
                 annotate_color_present: str = '#66c2a4',
                 annotate_color_absent: str = 'white',
                 annotate_color_unknown: str = 'lightgray',
                 annotate_opacity_present: float = 1.0,
                 annotate_opacity_absent: float = 0.0,
                 annotate_opacity_unknown: float = 0.0,
                 annotate_border_color: str = 'black',
                 annotate_kind: str = 'rect',
                 annotate_box_width: int = 30,
                 annotate_box_height: int = 30,
                 annotate_border_width: int = 1,
                 annotate_margin: int = 0,
                 include_unknown: bool = True,
                 annotate_show_box_label: bool = False,
                 annotate_box_label_color: str = 'white',
                 annotate_label_fontsize: int = 12):
        """
        Creates a new TreeStyler with the given default settings.
        :param tree: The tree to style and render.
        :param highlight_style: A style used to define how the highlight() method should behave.
                                Can either be one of the named highlight styles ['light', 'light_hn', 'pastel', 'medium', dark']
                                or an instance of a :py:class:`genomics_data_index.api.viewer.TreeStyler.HighlightStyle`.
        :param annotate_column: Which column in the annotate table this annotate visual belongs to (starting from column 0).
        :param tree_style: The ete3.TreeStyle object defining the style of the tree.
        :param legend_columns: A dictionary mapping legend column names to integers defining the column numbers.
        :param samples_styles_list: A list of :py:class:`genomics_data_index.api.viewer.TreeSamplesVisual` objects
                                   that define the different styles to apply to sets of samples (e.g., highlight or annotate).
        :param legend_nsize: The legend node size.
        :param legend_fsize: The legend font size.
        :param annotate_color_present: The default color of samples which are present in the set for the annotate() method.
        :param annotate_color_absent: The default color of samples which are absent in the set for the annotate() method.
        :param annotate_color_unknown: The default color of samples which are unknown if present/absent in the set for the annotate() method.
        :param annotate_opacity_present: The default opacity of samples which are present in the set for the annotate() method.
        :param annotate_opacity_absent: The default opacity of samples which are absent in the set for the annotate() method.
        :param annotate_opacity_unknown: The default opacity of samples which are unknown if present/absent in the set for the annotate() method.
        :param annotate_border_color: The default border color of the drawn annotations.
        :param annotate_kind: The default kind color of the drawn annotations (either 'circle' or 'rectangle').
        :param annotate_box_width: The width of the boxes for the drawn annotations.
        :param annotate_box_height: The height of the boxes for the drawn annotations.
        :param annotate_border_width: The width of the border for the boxes for the drawn annotations.
        :param annotate_margin: The margin width of the boxes for the drawn annotations.
        :param include_unknown: Whether or not to include unknowns in highlight/annotations.
        :param annotate_show_box_label: True if labels should be shown in the annotation boxes, False otherwise.
        :param annotate_box_label_color: The color of the labels in the annotation boxes.
        :param annotate_label_fontsize: The font size of the annotation labels.
        :return: A new TreeStyler object.
        """
        self._tree = tree
        self._default_highlight_styles = default_highlight_styles
        self._tree_style = tree_style
        self._legend_columns = legend_columns
        self._legend_nsize = legend_nsize
        self._legend_fsize = legend_fsize
        self._annotate_border_color = annotate_border_color
        self._annotate_color_present = annotate_color_present
        self._annotate_color_absent = annotate_color_absent
        self._annotate_color_unknown = annotate_color_unknown
        self._annotate_opacity_present = annotate_opacity_present
        self._annotate_opacity_absent = annotate_opacity_absent
        self._annotate_opacity_unknown = annotate_opacity_unknown
        self._annotate_column = annotate_column
        self._annotate_box_width = annotate_box_width
        self._annotate_box_height = annotate_box_height
        self._annotate_border_width = annotate_border_width
        self._annotate_margin = annotate_margin
        self._include_unknown = include_unknown
        self._annotate_show_box_label = annotate_show_box_label
        self._annotate_box_label_color = annotate_box_label_color
        self._annotate_box_label_unknown_color = 'black'
        self._annotate_label_fontsize = annotate_label_fontsize

        if annotate_kind not in AnnotateTreeSamplesVisual.ANNOTATE_KINDS:
            raise Exception(f'Invalid value for annotate_kind={annotate_kind}.'
                            f' Must be one of {AnnotateTreeSamplesVisual.ANNOTATE_KINDS}')
        self._annotate_kind = annotate_kind

        self._samples_styles_list = samples_styles_list

    def annotate(self, samples: Union[SamplesQuery, Iterable[str]],
                 box_label: Union[str, Dict[str, Any]] = None,
                 box_label_unknown: Union[str, Dict[str, Any]] = None,
                 annotate_show_box_label: bool = None, annotate_box_label_color: str = None,
                 legend_label: str = None,
                 box_width: int = None, box_height: int = None,
                 color_present: str = None, color_absent: str = None, color_unknown: str = None) -> TreeStyler:
        """
        Adds an annotation column beside the tree showing the which samples are in the passed set.
        :param samples: The samples to show as being present.
        :param box_label: A label for the column. Can be text or dict with attributes text, font, color, and fontsize
                      (this is passed to the underlying ete3 Face).
        :param box_label_unknown: A label for any samples with unknown status. Can be text or dict with attributes text,
                                  font, color, and fontsize (this is passed to the underlying ete3 Face).
        :param annotate_show_box_label: Whether or not the label should be added to every present item in the figure.
        :param annotate_box_label_color: If 'annotate_show_box_label' is set, this defines the color (overrides class setting).
        :param legend_label: A label to use for a legend item. The color should match the color for this annotated set.
        :param box_width: The width of the bounding box (defaults to class variable annotate_box_width).
        :param box_height: The height of the bounding box (defaults to class variable annotate_box_height).
        :param color_present: The color to use when a sample is present in this set (defaults class-defined color).
        :param color_absent: The color to use when a sample is absent (defaults to class-defined color).
        :param color_unknown: The color to use when a sample is unknown if present/absent (defaults to class-defined color).
        :return: A new TreeStyler object which contains the completed annotation column.
        """
        if box_width is None:
            box_width = self._annotate_box_width
        if box_height is None:
            box_height = self._annotate_box_height
        if annotate_show_box_label is None:
            annotate_show_box_label = self._annotate_show_box_label
        if annotate_box_label_color is None:
            annotate_box_label_color = self._annotate_box_label_color
        if color_present is None:
            color_present = self._annotate_color_present
        if color_absent is None:
            color_absent = self._annotate_color_absent
        if color_unknown is None:
            color_unknown = self._annotate_color_unknown

        if isinstance(box_label, str):
            # Pick default color since ete3 by default colors the same as what I'm using for the fill color
            box_label = {'text': box_label, 'color': 'black', 'font': 'Verdana',
                         'fontsize': self._annotate_label_fontsize}

        if box_label_unknown is None and box_label is not None:
            box_label_unknown = 'U'

        if isinstance(box_label_unknown, str):
            # Pick default color since ete3 by default colors the same as what I'm using for the fill color
            box_label_unknown = {'text': box_label_unknown, 'color': 'black', 'font': 'Verdana',
                                 'fontsize': self._annotate_label_fontsize}

        if box_label is not None:
            label_present = copy.deepcopy(box_label)
            if 'fontsize' not in label_present:
                label_present['fontsize'] = self._annotate_label_fontsize
                label_present['color'] = self._annotate_box_label_color
        else:
            label_present = None

        box_label_color_unknown = self._annotate_box_label_unknown_color
        if box_label_unknown is not None:
            label_unknown = copy.deepcopy(box_label_unknown)
            if 'fontsize' not in label_unknown:
                label_unknown['fontsize'] = self._annotate_label_fontsize
            if 'color' not in label_unknown:
                label_unknown['color'] = self._annotate_box_label_unknown_color
            else:
                box_label_color_unknown = label_unknown['color']
        else:
            label_unknown = None

        samples_visual = AnnotateTreeSamplesVisual(samples=samples,
                                                   label=label_present,
                                                   label_unknown=label_unknown,
                                                   annotate_show_box_label=annotate_show_box_label,
                                                   annotate_box_label_color=annotate_box_label_color,
                                                   annotate_box_label_unknown_color=box_label_color_unknown,
                                                   annotate_label_fontsize=self._annotate_label_fontsize,
                                                   legend_label=legend_label,
                                                   box_width=box_width,
                                                   box_height=box_height,
                                                   color_present=color_present,
                                                   color_absent=color_absent,
                                                   color_unknown=color_unknown,
                                                   legend_nodesize=self._legend_nsize,
                                                   legend_fontsize=self._legend_fsize,
                                                   annotate_column=self._annotate_column,
                                                   annotate_kind=self._annotate_kind,
                                                   annotate_border_color=self._annotate_border_color,
                                                   annotate_opacity_present=self._annotate_opacity_present,
                                                   annotate_opacity_absent=self._annotate_opacity_absent,
                                                   annotate_opacity_unknown=self._annotate_opacity_unknown,
                                                   border_width=self._annotate_border_width,
                                                   margin=self._annotate_margin,
                                                   include_unknown=self._include_unknown,
                                                   legend_columns=self._legend_columns)
        samples_styles_list_new = copy.copy(self._samples_styles_list)
        samples_styles_list_new.append(samples_visual)

        return TreeStyler(self._tree, default_highlight_styles=self._default_highlight_styles,
                          tree_style=self._tree_style,
                          samples_styles_list=samples_styles_list_new,
                          legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize,
                          annotate_column=self._annotate_column + 1,
                          annotate_color_present=self._annotate_color_present,
                          annotate_color_absent=self._annotate_color_absent,
                          annotate_color_unknown=self._annotate_color_unknown,
                          annotate_opacity_present=self._annotate_opacity_present,
                          annotate_opacity_absent=self._annotate_opacity_absent,
                          annotate_opacity_unknown=self._annotate_opacity_unknown,
                          annotate_border_color=self._annotate_border_color,
                          annotate_kind=self._annotate_kind,
                          annotate_box_width=self._annotate_box_width,
                          annotate_box_height=self._annotate_box_height,
                          annotate_border_width=self._annotate_border_width,
                          annotate_margin=self._annotate_margin,
                          include_unknown=self._include_unknown,
                          annotate_show_box_label=self._annotate_show_box_label,
                          annotate_box_label_color=self._annotate_box_label_color,
                          annotate_label_fontsize=self._annotate_label_fontsize,
                          legend_columns=self._legend_columns)

    def highlight(self, samples: Union[SamplesQuery, Iterable[str]],
                  present_node_style: NodeStyle = None, unknown_node_style: NodeStyle = None,
                  legend_color: str = None, legend_label: str = None) -> TreeStyler:
        """
        Highlights the selected samples in the tree with a particular style defined by the passed default_highlight_styles
        on creation of this TreeStyler.
        :param samples: The set of samples to highlight.
        :param present_node_style: Overrides the node style for present samples
                                  corresponding to leaf nodes in the tree defined by the highlight styles.
        :param unknown_node_style: Overrides the node style for unknown samples
                                  corresponding to leaf nodes in the tree defined by the highlight styles.
        :param legend_color: Overrides the color of the legend item.
        :param legend_label: Defines the label of the legend item.
        :return: A new TreeStyler with the given samples highlighted.
        """

        if present_node_style is None and legend_color is None:
            present_node_style = self._default_highlight_styles.present_node_style
            legend_color = self._default_highlight_styles.legend_color

            if unknown_node_style is None:
                unknown_node_style = self._default_highlight_styles.unknown_node_style

            new_default_styles = self._default_highlight_styles.next()
        else:
            new_default_styles = self._default_highlight_styles
            if unknown_node_style is None:
                unknown_node_style = self._default_highlight_styles.unknown_node_style

        samples_visual = HighlightTreeSamplesVisual(samples=samples,
                                                    present_node_style=present_node_style,
                                                    unknown_node_style=unknown_node_style,
                                                    include_unknown=self._include_unknown,
                                                    legend_color=legend_color,
                                                    legend_label=legend_label,
                                                    legend_nodesize=self._legend_nsize,
                                                    legend_fontsize=self._legend_fsize,
                                                    legend_columns=self._legend_columns)
        samples_styles_list_new = copy.copy(self._samples_styles_list)
        samples_styles_list_new.append(samples_visual)

        return TreeStyler(self._tree, default_highlight_styles=new_default_styles,
                          tree_style=self._tree_style,
                          samples_styles_list=samples_styles_list_new,
                          legend_fsize=self._legend_fsize, legend_nsize=self._legend_nsize,
                          annotate_column=self._annotate_column,
                          annotate_color_present=self._annotate_color_present,
                          annotate_color_absent=self._annotate_color_absent,
                          annotate_color_unknown=self._annotate_color_unknown,
                          annotate_opacity_present=self._annotate_opacity_present,
                          annotate_opacity_absent=self._annotate_opacity_absent,
                          annotate_opacity_unknown=self._annotate_opacity_unknown,
                          annotate_border_color=self._annotate_border_color,
                          annotate_kind=self._annotate_kind,
                          annotate_box_height=self._annotate_box_height,
                          annotate_box_width=self._annotate_box_width,
                          annotate_border_width=self._annotate_border_width,
                          annotate_margin=self._annotate_margin,
                          include_unknown=self._include_unknown,
                          annotate_show_box_label=self._annotate_show_box_label,
                          annotate_box_label_color=self._annotate_box_label_color,
                          annotate_label_fontsize=self._annotate_label_fontsize,
                          legend_columns=self._legend_columns)

    def _apply_samples_styles(self, tree: Tree, tree_style: TreeStyle) -> None:
        for samples_style in self._samples_styles_list:
            samples_style.apply_visual(tree=tree, tree_style=tree_style)

    def render(self, file_name: str = '%%inline', w: int = None, h: int = None,
               tree_style: TreeStyle = None, units: str = 'px', dpi: int = 90):
        """
        Renders the tree with the styles defined by this TreeStyler object to an image.
        :param file_name: The image file name to save, use '%%inline' to draw inline in Jupyter notebooks (default).
        :param w: The width of the image.
        :param h: The height of the image.
        :param tree_style: The ete3.TreeStyle object (overrides the object defined by this TreeStyler).
        :param units: The units of the width/height (default in pixels).
        :param dpi: The dots per inch.
        :return: An image/image data of the drawn tree styled according to this TreeStyler.
        """
        if tree_style is None:
            tree_style = self._tree_style

        # Set default width if no width or height specified
        # Do this here instead of as a default method value since at least one of
        # width or height must be None to preserve the correct aspect ratio (so I only
        # want a default value for width if neither width nor height are set.
        if w is None and h is None:
            w = 400

        tree = self._tree.copy('newick')
        tree_style = copy.deepcopy(tree_style)
        self._apply_samples_styles(tree=tree, tree_style=tree_style)

        return tree.render(file_name=file_name, w=w, h=h, tree_style=tree_style,
                           units=units, dpi=dpi)

    @property
    def tree(self) -> Tree:
        """
        Gets a copy of the tree this TreeStyler is applied to.
        Note that any highlight() or annotate() methods only get applied when running render().
        So the returned tree will *not* be styled according to the defined highlight() or annotate() objects.
        This is in part because for large trees copying all the tree style/node style objects from ete3 was causing
        recursion errors.

        :return: The tree this TreeStyler is applied to.
        """
        return self._tree.copy('newick')

    @property
    def tree_style(self) -> TreeStyle:
        """
        Gets a copy of the ete3.TreeStyle object being used to draw this tree.
        :return: The ete3.TreeStyle object.
        """
        return copy.deepcopy(self._tree_style)

    @classmethod
    def create(cls, tree: Tree,
               initial_style: TreeStyle = None,
               mode='r',
               highlight_style: Union[str, HighlightStyle] = 'light',
               legend_nsize: int = 20, legend_fsize: int = 11,
               annotate_color_present: str = 'black',
               annotate_color_absent: str = 'white',
               annotate_color_unknown: str = 'lightgray',
               annotate_opacity_present: float = 1.0,
               annotate_opacity_absent: float = 0.0,
               annotate_opacity_unknown: float = 1.0,
               annotate_border_color: str = 'black',
               annotate_kind: str = 'rect',
               annotate_box_width: int = 30,
               annotate_box_height: int = 30,
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
               annotate_box_label_color: str = 'white',
               annotate_arc_span: int = 350,
               annotate_label_fontsize: int = 12,
               show_leaf_names: bool = True,
               include_unknown: bool = True,
               show_legend_type_labels: bool = True,
               legend_type_label_present: str = 'Pr.',
               legend_type_label_unknown: str = 'Un.',
               tree_scale: float = None) -> TreeStyler:
        """
        Constructs a new TreeStyler object used to style and visualize trees.
        All parameters listed below are optional except for tree.
        :param tree: The tree to style.
        :param initial_style: The initial ete3.TreeStyle to start with.
        :param mode: Either 'r' (rectangular) or 'c' (circular).
        :param highlight_style: A style used to define how the highlight() method should behave.
                                Can either be one of the named highlight styles ['light', 'light_hn', 'pastel', 'medium', dark']
                                or an instance of a :py:class:`genomics_data_index.api.viewer.TreeStyler.HighlightStyle`.
        :param legend_nsize: The legend node size.
        :param legend_fsize: The legend font size.
        :param annotate_color_present: The default color of samples which are present in the set for the annotate() method.
        :param annotate_color_absent: The default color of samples which are absent in the set for the annotate() method.
        :param annotate_color_unknown: The default color of samples which are unknown if present/absent in the set for the annotate() method.
        :param annotate_opacity_present: The default opacity of samples which are present in the set for the annotate() method.
        :param annotate_opacity_absent: The default opacity of samples which are absent in the set for the annotate() method.
        :param annotate_opacity_unknown: The default opacity of samples which are unknown if present/absent in the set for the annotate() method.
        :param annotate_border_color: The default border color of the drawn annotations.
        :param annotate_kind: The default kind color of the drawn annotations (either 'circle' or 'rectangle').
        :param annotate_box_width: The width of the boxes for the drawn annotations.
        :param annotate_box_height: The height of the boxes for the drawn annotations.
        :param annotate_border_width: The width of the border for the boxes for the drawn annotations.
        :param annotate_margin: The margin width of the boxes for the drawn annotations.
        :param annotate_guiding_lines: True if guiding lines should be drawn that matches sample names to the annotation boxes,
                                       False otherwise.
        :param annotate_guiding_lines_color:  The color of the annotate guiding lines.
        :param figure_margin: The margin spacing (used for all of top, bottom, left, and right) for the overall figure.
        :param show_border: True if a border should be shown around the overall figure, False otherwise.
        :param title: A title for the figure.
        :param title_fsize: The font size of the figure title.
        :param legend_title: The title of the legend.
        :param annotate_show_box_label: True if labels should be shown in the annotation boxes, False otherwise.
        :param annotate_box_label_color: The color of the labels in the annotation boxes.
        :param annotate_arc_span: For mode='c' (circular) the degrees the circular tree should span.
        :param annotate_label_fontsize: The font size of the annotation labels.
        :param include_unknown: Whether or not to include unknown samples in highlight/annotation boxes.
        :param show_leaf_names: True if leaf names should be shown on the tree, False otherwise.
        :param show_legend_type_labels: Whether or not to show labels for legend types/categories (present or unknown).
        :param legend_type_label_present: Text to show above legend color for present items.
        :param legend_type_label_unknown: Text to show above legend color for unknown items.
        :param tree_scale: A scale factor for the tree.
        :return: A new TreeStyler object used to style and visualize trees.
        """
        if initial_style is not None:
            tree_style_elements = {'mode': mode, 'annotate_guiding_lines': annotate_guiding_lines,
                                   'figure_margin': figure_margin, 'show_border': show_border,
                                   'title': title, 'legend_title': legend_title, 'title_fsize': title_fsize,
                                   'annotate_arc_span': annotate_arc_span, 'show_leaf_names': show_leaf_names,
                                   'tree_scale': tree_scale}
            tree_style_elements_none = [tree_style_elements[k] is None for k in tree_style_elements]

            if any(tree_style_elements_none):
                logger.warning(
                    f'Both initial_style=[{initial_style}] and one of parameters {tree_style_elements.keys()}'
                    f' are set. Will ignore these listed parameters.')
            ts = initial_style
        else:
            ts = TreeStyle()
            ts.arc_span = annotate_arc_span
            ts.show_leaf_name = show_leaf_names

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
            ts.scale = tree_scale

        if include_unknown:
            legend_columns = {
                'present': 0,
                'unknown': 1,
                'text': 2
            }
        else:
            legend_columns = {
                'present': 0,
                'text': 1
            }

        if highlight_style is None:
            highlight_style = HighlightStyle.create('pastel')
        elif isinstance(highlight_style, str):
            highlight_style = HighlightStyle.create(highlight_style)

        if legend_title is not None:
            margin_bottom = 10
            cf = RectFace(width=10, height=10, fgcolor=None, bgcolor=None)
            cf.margin_bottom = margin_bottom
            tf = TextFace(legend_title, fsize=legend_fsize + 3)
            tf.margin_bottom = margin_bottom

            if include_unknown:
                ts.legend.add_face(cf, column=legend_columns['present'])
                ts.legend.add_face(cf, column=legend_columns['unknown'])
                ts.legend.add_face(tf, column=legend_columns['text'])
            else:
                ts.legend.add_face(cf, column=legend_columns['present'])
                ts.legend.add_face(tf, column=legend_columns['text'])

        if show_legend_type_labels:
            margin = 10
            item_ptf = TextFace(legend_type_label_present, fsize=legend_fsize)
            item_ptf.margin_bottom = margin
            item_ptf.margin_left = margin
            item_utf = TextFace(legend_type_label_unknown, fsize=legend_fsize)
            item_utf.margin_bottom = margin
            item_utf.margin_left = margin
            item_cf = RectFace(width=10, height=10, fgcolor=None, bgcolor=None)
            item_cf.margin_bottom = margin
            item_cf.margin_left = margin

            if include_unknown:
                ts.legend.add_face(item_ptf, column=legend_columns['present'])
                ts.legend.add_face(item_utf, column=legend_columns['unknown'])
                ts.legend.add_face(item_cf, column=legend_columns['text'])
            else:
                ts.legend.add_face(item_ptf, column=legend_columns['present'])
                ts.legend.add_face(item_cf, column=legend_columns['text'])

        if title is not None:
            tf = TextFace(title, fsize=title_fsize)
            ts.title.add_face(tf, column=0)

        return TreeStyler(tree=tree,
                          default_highlight_styles=highlight_style,
                          tree_style=ts,
                          legend_columns=legend_columns,
                          samples_styles_list=[],
                          legend_nsize=legend_nsize,
                          legend_fsize=legend_fsize,
                          annotate_column=1,
                          annotate_color_present=annotate_color_present,
                          annotate_color_absent=annotate_color_absent,
                          annotate_color_unknown=annotate_color_unknown,
                          annotate_opacity_present=annotate_opacity_present,
                          annotate_opacity_absent=annotate_opacity_absent,
                          annotate_opacity_unknown=annotate_opacity_unknown,
                          annotate_border_color=annotate_border_color,
                          annotate_kind=annotate_kind,
                          annotate_box_width=annotate_box_width,
                          annotate_box_height=annotate_box_height,
                          annotate_border_width=annotate_border_width,
                          annotate_margin=annotate_margin,
                          include_unknown=include_unknown,
                          annotate_show_box_label=annotate_show_box_label,
                          annotate_box_label_color=annotate_box_label_color,
                          annotate_label_fontsize=annotate_label_fontsize)


class HighlightStyle:
    """
    Defines a list of styles used to highlight different sets of samples in a tree.
    """

    THEMES = ['light', 'light_hn', 'pastel', 'medium', 'dark']

    def __init__(self, node_styles: List[Dict[str, Union[str, NodeStyle]]],
                 index: int):
        """
        Builds a new HighlightStyle with the given information.
        :param node_styles: A list of NodeStyle and legend color objects
                           (a list of dictionaries with two keys: 'nstyle' and 'legend_color').
        :param index: The index for the current element in the list of styles (the one that will be used next).
        :return: A new HighlightStyle.
        """
        self._node_styles = node_styles
        self._index = index

    @property
    def present_node_style(self) -> NodeStyle:
        """
        Gets the current node style from the HighlightStyle for present samples.
        :return: The node style for present samples
        """
        return self._node_styles[self._index]['present_nstyle']

    @property
    def unknown_node_style(self) -> NodeStyle:
        """
        Gets the current node style from the HighlightStyle for unknown samples.
        :return: The node style for unknown samples
        """
        return self._node_styles[self._index]['unknown_nstyle']

    @property
    def legend_color(self) -> NodeStyle:
        """
        Gets the current legend color from the HighlightStyle.
        :return: The legend color.
        """
        return self._node_styles[self._index]['legend_color']

    def next(self) -> HighlightStyle:
        """
        Advances the current highlight style by 1. Will wrap around to the beginning of the list of styles
        if the end is reached. This is used so that every execution of TreeStyler.highlight() will use a different highlight style.
        :return: A new HighlightStyle object with the next style selected in the list.
        """
        # Advance by 1, wrapping around if we reach the end of the highlight styles
        new_index = (self._index + 1) % len(self._node_styles)
        return HighlightStyle(self._node_styles, index=new_index)

    @classmethod
    def create(cls, kind: str) -> HighlightStyle:
        """
        Creates a new pre-defined HighlightStyle.
        :param kind: The kind (name) of the pre-defined HighlightStyle to create.
        :return: A new HighlightStyle.
        """
        unknown_fg_color = 'lightgray'
        unknown_bg_color = 'lightgray'

        if kind == 'light':
            fg_colors = ['#e5f5f9', '#fee8c8', '#e0ecf4', '#deebf7']
            bg_colors = fg_colors
            return cls._create_highlights(fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'light_hn':
            fg_colors = ['#41ae76', '#ef6548', '#8c6bb1', '#4292c6']
            bg_colors = ['#e5f5f9', '#fee8c8', '#e0ecf4', '#deebf7']
            return cls._create_highlights(fg_colors=fg_colors, bg_colors=bg_colors, nsize=20,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'pastel':
            fg_colors = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6']
            bg_colors = fg_colors
            return cls._create_highlights(fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'medium':
            fg_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
            bg_colors = fg_colors
            return cls._create_highlights(fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'dark':
            fg_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']
            bg_colors = fg_colors
            return cls._create_highlights(fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        else:
            raise Exception(f'kind=[{kind}] must be one of {cls.THEMES}')

    @classmethod
    def _create_highlights(cls, fg_colors: List[str], bg_colors: List[str],
                           unknown_fg_color: str, unknown_bg_color: str,
                           nsize: int = None) -> HighlightStyle:
        if len(fg_colors) != len(bg_colors):
            raise Exception(f'fg_colors={fg_colors} and bg_colors={bg_colors} do not have the same length')

        styles = []
        for i in range(0, len(fg_colors)):
            nstyle = NodeStyle()
            nstyle['fgcolor'] = fg_colors[i]
            nstyle['bgcolor'] = bg_colors[i]
            if nsize is not None:
                nstyle['size'] = nsize

            unknown_nstyle = NodeStyle()
            unknown_nstyle['fgcolor'] = unknown_fg_color
            unknown_nstyle['bgcolor'] = unknown_bg_color
            if nsize is not None:
                unknown_nstyle['size'] = nsize

            styles.append({
                'present_nstyle': nstyle,
                'unknown_nstyle': unknown_nstyle,
                'legend_color': fg_colors[i]
            })

        return HighlightStyle(styles, index=0)
