from __future__ import annotations

import copy
import logging
from typing import List, Dict, Any, Union, Iterable, Tuple, Callable

import pandas as pd
from ete3 import Tree, TreeNode

from genomics_data_index.ete3treeview import NodeStyle, TreeStyle, TextFace, RectFace
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual
from genomics_data_index.api.viewer.samples_visuals.AnnotateTreeSamplesVisual import AnnotateTreeSamplesVisual
from genomics_data_index.api.viewer.samples_visuals.AnnotationSpacingTreeSamplesVisual import \
    AnnotationSpacingTreeSamplesVisual
from genomics_data_index.api.viewer.samples_visuals.HighlightTreeSamplesVisual import HighlightTreeSamplesVisual

logger = logging.getLogger(__name__)

try:
    DEFAULT_NODE_STYLE = NodeStyle()
    DEFAULT_NODE_STYLE['size'] = 0
except Exception as e:
    logger.warning(e)
    DEFAULT_NODE_STYLE = None


class TreeStyler:
    """
    A class used to style and render a tree and sets of samples on this tree.
    """

    MODES = ['r', 'c']

    def __init__(self, tree: Tree, default_highlight_styles: HighlightStyle, annotate_column: int,
                 tree_style: TreeStyle,
                 node_style: NodeStyle,
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
                 annotate_label_fontsize: int = 12,
                 show_leaf_names: bool = True,
                 leaf_name_fontsize: int = 12,
                 leaf_name_func: Callable[[str, pd.Series], str] = None,
                 sample_metadata: pd.DataFrame = None):
        """
        Creates a new TreeStyler with the given default settings.
        :param tree: The tree to style and render.
        :param highlight_style: A style used to define how the highlight() method should behave.
                                Can either be one of the named highlight styles ['light', 'light_hn', 'pastel', 'medium', dark']
                                or an instance of a :py:class:`genomics_data_index.api.viewer.TreeStyler.HighlightStyle`.
        :param annotate_column: Which column in the annotate table this annotate visual belongs to (starting from column 0).
        :param tree_style: The ete3.TreeStyle object defining the style of the tree.
        :param node_style: The default ete3.NodeStyle object for styling the nodes of the tree.
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
        :param show_leaf_names: Whether or not to show leaf names.
        :param leaf_name_fontsize: The font size of leaf names.
        :param leaf_name_func: A function which lets you create custom leaf names to display.
                               The function should look like: func(name: str, metadata: pd.Series) -> str.
                               That is it takes as input the leaf name and a pandas.Series of metadata for the
                               particular sample (derived from the query.toframe() set of data).
                               For example: leaf_name_func=lambda name, metadata: f'{name}_{metadata["Location"]}'
                               This would display a label like "SampleX_Canada" for each sample (assumes that "Location"
                               is a column name in the table produced by query.toframe()).
        :param sample_metadata: A dataframe of sample metadata indexed by Sample Name.
                                This needs to be set if leaf_name_func is set.
        :return: A new TreeStyler object.
        """
        self._tree = tree
        self._default_highlight_styles = default_highlight_styles
        self._tree_style = tree_style
        self._node_style = node_style
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

        self._show_leaf_names = show_leaf_names
        self._leaf_name_fontsize = leaf_name_fontsize
        self._leaf_name_func = leaf_name_func

        if self._leaf_name_func is not None and sample_metadata is None:
            raise Exception(f'leaf_name_func is not None but sample_metadata is None')

        self._sample_metadata = sample_metadata

    def add_spacing(self, width: int = None, height: int = None, color: str = None,
                    border_width: int = None, border_color: str = None) -> TreeStyler:
        """
        Adds an empty column in the annotation table as additional space.
        :param width: The width of the spacing (default: annotation column width).
        :param height: The height of each annotation block used for spacing (default: annotation column block height).
        :param color: The color of the spacing (default: no color).
        :param border_width: Specify the width of the border (default: None for no border).
        :param border_color: Specify the border color (default: default border color of other annotation elements).
        :return: A new TreeStyler object with the spacing added in.
        """
        if width is None:
            width = self._annotate_box_width
        if height is None:
            height = self._annotate_box_height
        if border_color is None:
            border_color = self._annotate_border_color

        samples_visual = AnnotationSpacingTreeSamplesVisual(width=width,
                                                            height=height,
                                                            annotate_column=self._annotate_column,
                                                            color=color,
                                                            border_width=border_width,
                                                            border_color=border_color)

        samples_styles_list_new = copy.copy(self._samples_styles_list)
        samples_styles_list_new.append(samples_visual)

        return TreeStyler(self._tree, default_highlight_styles=self._default_highlight_styles,
                          tree_style=self._tree_style,
                          node_style=self._node_style,
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
                          legend_columns=self._legend_columns,
                          show_leaf_names=self._show_leaf_names,
                          leaf_name_fontsize=self._leaf_name_fontsize,
                          leaf_name_func=self._leaf_name_func,
                          sample_metadata=self._sample_metadata)

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
                          node_style=self._node_style,
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
                          legend_columns=self._legend_columns,
                          show_leaf_names=self._show_leaf_names,
                          leaf_name_fontsize=self._leaf_name_fontsize,
                          leaf_name_func=self._leaf_name_func,
                          sample_metadata=self._sample_metadata)

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
                          node_style=self._node_style,
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
                          legend_columns=self._legend_columns,
                          show_leaf_names=self._show_leaf_names,
                          leaf_name_fontsize=self._leaf_name_fontsize,
                          leaf_name_func=self._leaf_name_func,
                          sample_metadata=self._sample_metadata)

    def _apply_samples_styles(self, tree: Tree, tree_style: TreeStyle) -> None:
        for samples_style in self._samples_styles_list:
            samples_style.apply_visual(tree=tree, tree_style=tree_style)

    def render(self, file_name: str = '%%inline', w: int = None, h: int = None,
               tree_style: TreeStyle = None, units: str = 'px', dpi: int = 90,
               ladderize: bool = False):
        """
        Renders the tree with the styles defined by this TreeStyler object to an image.
        :param file_name: The image file name to save, use '%%inline' to draw inline in Jupyter notebooks (default).
        :param w: The width of the image.
        :param h: The height of the image.
        :param tree_style: The ete3.TreeStyle object (overrides the object defined by this TreeStyler).
        :param units: The units of the width/height (default in pixels).
        :param dpi: The dots per inch.
        :param ladderize: Runs ete3.Tree.ladderize() to sort branches of a given tree according to size.
        :return: An image/image data of the drawn tree styled according to this TreeStyler.
        """

        # Set default width if no width or height specified
        # Do this here instead of as a default method value since at least one of
        # width or height must be None to preserve the correct aspect ratio (so I only
        # want a default value for width if neither width nor height are set.
        if w is None and h is None:
            w = 400

        tree, tree_style = self.prerender(tree_style=tree_style, ladderize=ladderize)

        return self.draw(tree=tree, tree_style=tree_style,
                         file_name=file_name, w=w, h=h,
                         units=units, dpi=dpi)

    def _add_node_name(self, node: TreeNode) -> None:
        if self._leaf_name_func is not None and node.name in self._sample_metadata.index:
            sample_metadata_s = self._sample_metadata.loc[node.name]
            node_name = self._leaf_name_func(node.name, sample_metadata_s)
        else:
            node_name = node.name
        tf = TextFace(node_name, fsize=self._leaf_name_fontsize)
        node.add_face(tf, 0, position='branch-right')

    def prerender(self, tree_style: TreeStyle = None, ladderize: bool = False) -> Tuple[Tree, TreeStyle]:
        """
        Pre-renders the tree, returning an ete3.Tree and ete3.TreeStyle object which can be used to draw the tree.
        :param tree_style: The ete3.TreeStyle object (overrides the object defined by this TreeStyler).
        :param ladderize: Runs ete3.Tree.ladderize() to sort branches of a given tree according to size.
        :return: A tuple of <ete3.Tree, ete3.TreeStyle> with highlights/annotations added.
        """
        if tree_style is None:
            tree_style = self._tree_style

        tree = self._tree.copy('newick')
        tree_style = copy.deepcopy(tree_style)

        # Sets default node style for all nodes in the tree and adds the leaf label
        for n in tree.traverse():
            n.set_style(self._node_style)

            if self._show_leaf_names and n.is_leaf():
                self._add_node_name(n)

        self._apply_samples_styles(tree=tree, tree_style=tree_style)

        if ladderize:
            tree.ladderize()

        return tree, tree_style

    @classmethod
    def draw(cls, tree: Tree, tree_style: TreeStyle,
             file_name: str = '%%inline', w: int = None, h: int = None,
             units: str = 'px', dpi: int = 90):
        """
        Draws the passed tree with the passed styles to an image.
        :param tree: The tree to draw.
        :param tree_style: The ete3.TreeStyle object (overrides the object defined by this TreeStyler).
        :param file_name: The image file name to save, use '%%inline' to draw inline in Jupyter notebooks (default).
        :param w: The width of the image.
        :param h: The height of the image.
        :param units: The units of the width/height (default in pixels).
        :param dpi: The dots per inch.
        :return: An image/image data of the drawn tree styled according to this TreeStyler.
        """
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
               mode: str = 'r',
               highlight_style: Union[str, HighlightStyle] = 'light',
               node_style: NodeStyle = None,
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
               leaf_name_fontsize: int = 12,
               leaf_name_func: Callable[[str, pd.Series], str] = None,
               sample_metadata: pd.DataFrame = None,
               include_unknown: bool = True,
               show_legend_type_labels: bool = True,
               legend_type_label_present: str = 'Pr.',
               legend_type_label_unknown: str = 'Un.',
               rotation: float = 0,
               allow_face_overlap: bool = False,
               show_branch_length: bool = False,
               show_branch_support: bool = False,
               tree_line_width: int = None,
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
        :param node_style: The default ete3.NodeStyle object for styling the nodes of the tree.
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
        :param leaf_name_fontsize: The font size of leaf names.
        :param leaf_name_func: A function which lets you create custom leaf names to display.
                               The function should look like: func(name: str, metadata: pd.Series) -> str.
                               That is it takes as input the leaf name and a pandas.Series of metadata for the
                               particular sample (derived from the query.toframe() set of data).
                               For example: leaf_name_func=lambda name, metadata: f'{name}_{metadata["Location"]}'
                               This would display a label like "SampleX_Canada" for each sample (assumes that "Location"
                               is a column name in the table produced by query.toframe()).
        :param sample_metadata: A dataframe of sample metadata indexed by Sample Name.
                                This needs to be set if leaf_name_func is set.
        :param show_legend_type_labels: Whether or not to show labels for legend types/categories (present or unknown).
        :param legend_type_label_present: Text to show above legend color for present items.
        :param legend_type_label_unknown: Text to show above legend color for unknown items.
        :param rotation: The rotation of the tree in degrees (for circular mode).
        :param allow_face_overlap: Allow overlap in node faces for circular images.
        :param show_branch_length: Show branch lengths.
        :param show_branch_support: Show branch supports.
        :param tree_line_width: The line width for the tree. Overrides 'hz_line_width' and 'vt_line_width'
                                in node_style. Default: None (no overriding of line width in node_style).
        :param tree_scale: A scale factor for the tree.
        :return: A new TreeStyler object used to style and visualize trees.
        """
        if initial_style is not None:
            tree_style_elements = {'mode': mode, 'annotate_guiding_lines': annotate_guiding_lines,
                                   'figure_margin': figure_margin, 'show_border': show_border,
                                   'title': title, 'legend_title': legend_title, 'title_fsize': title_fsize,
                                   'annotate_arc_span': annotate_arc_span,
                                   'tree_scale': tree_scale, 'rotation': rotation,
                                   'allow_face_overlap': allow_face_overlap, 'show_branch_length': show_branch_length,
                                   'show_branch_support': show_branch_support}
            tree_style_elements_none = [tree_style_elements[k] is None for k in tree_style_elements]

            if any(tree_style_elements_none):
                logger.warning(
                    f'Both initial_style=[{initial_style}] and one of parameters {tree_style_elements.keys()}'
                    f' are set. Will ignore these listed parameters.')
            ts = copy.deepcopy(initial_style)
            ts.show_leaf_name = False
        else:
            ts = TreeStyle()
            ts.arc_span = annotate_arc_span
            ts.rotation = rotation
            ts.allow_face_overlap = allow_face_overlap
            ts.show_branch_length = show_branch_length
            ts.show_branch_support = show_branch_support

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

            # ete3 includes a built-in method to show leaf names, but I need more control over the text to display
            # and font sizes. So instead of using the default method I implement my own (and so always set
            # show_leaf_names to False).
            ts.show_leaf_name = False

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

        if node_style is None:
            node_style = DEFAULT_NODE_STYLE

        if tree_line_width is not None:
            node_style = copy.deepcopy(node_style)
            node_style['hz_line_width'] = tree_line_width
            node_style['vt_line_width'] = tree_line_width

        if highlight_style is None:
            highlight_style = HighlightStyle.create('pastel', base_node_style=node_style)
        elif isinstance(highlight_style, str):
            highlight_style = HighlightStyle.create(highlight_style, base_node_style=node_style)
        elif tree_line_width is not None:
            highlight_style = highlight_style.copy_and_update(tree_line_width=tree_line_width)

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

        if leaf_name_func is not None and sample_metadata is None:
            raise Exception(f'leaf_name_func is not None but sample_metadata is None')

        return TreeStyler(tree=tree,
                          default_highlight_styles=highlight_style,
                          tree_style=ts,
                          node_style=node_style,
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
                          annotate_label_fontsize=annotate_label_fontsize,
                          show_leaf_names=show_leaf_names,
                          leaf_name_fontsize=leaf_name_fontsize,
                          leaf_name_func=leaf_name_func,
                          sample_metadata=sample_metadata)


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

    def copy_and_update(self, tree_line_width: int = None) -> HighlightStyle:
        """
        Copies the entire list of highlight styles and updates with the given parameters. Used to update
        certain style elements after the highlight styles are created.
        :param tree_line_width: The width of the lines used to draw the tree. Default: None (do not update).
        :return: A new HighlightStyle with the passed style elements updated.
        """
        node_styles_copy = copy.deepcopy(self._node_styles)

        if tree_line_width is not None:
            for style_dict in node_styles_copy:
                style_dict['present_nstyle']['hz_line_width'] = tree_line_width
                style_dict['present_nstyle']['vt_line_width'] = tree_line_width
                style_dict['unknown_nstyle']['hz_line_width'] = tree_line_width
                style_dict['unknown_nstyle']['vt_line_width'] = tree_line_width

        return HighlightStyle(node_styles=node_styles_copy, index=self._index)

    @classmethod
    def create(cls, kind: str = None, colors: List[str] = None, unknown_bg_color: str = 'lightgray',
               base_node_style: NodeStyle = None) -> HighlightStyle:
        """
        Creates a new pre-defined HighlightStyle.
        :param kind: The kind (name) of the pre-defined HighlightStyle to create.
        :param colors: A list of colors to use for highlights (unused if kind is set).
        :param unknown_bg_color: The color for unknown nodes.
        :param base_node_style: The base node style to use.
        :return: A new HighlightStyle.
        """
        if base_node_style is None:
            base_node_style = DEFAULT_NODE_STYLE

        default_fg_color = base_node_style['fgcolor']
        unknown_fg_color = default_fg_color

        if kind is None and colors is not None:
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=[default_fg_color] * len(colors),
                                          bg_colors=colors, unknown_bg_color=unknown_bg_color,
                                          unknown_fg_color=unknown_fg_color)
        elif kind is None:
            kind = 'light'

        if kind == 'light':
            fg_colors = [default_fg_color] * 4
            bg_colors = ['#e5f5f9', '#fee8c8', '#e0ecf4', '#deebf7']
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'light_hn':
            fg_colors = ['#41ae76', '#ef6548', '#8c6bb1', '#4292c6']
            bg_colors = ['#e5f5f9', '#fee8c8', '#e0ecf4', '#deebf7']
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=fg_colors, bg_colors=bg_colors,
                                          nsize=20,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'pastel':
            fg_colors = [default_fg_color] * 5
            bg_colors = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6']
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'medium':
            fg_colors = [default_fg_color] * 5
            bg_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        elif kind == 'dark':
            fg_colors = [default_fg_color] * 5
            bg_colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e']
            return cls._create_highlights(base_node_style=base_node_style, fg_colors=fg_colors, bg_colors=bg_colors,
                                          unknown_bg_color=unknown_bg_color, unknown_fg_color=unknown_fg_color)
        else:
            raise Exception(f'kind=[{kind}] must be one of {cls.THEMES}')

    @classmethod
    def _create_highlights(cls, base_node_style: NodeStyle, fg_colors: List[str], bg_colors: List[str],
                           unknown_fg_color: str, unknown_bg_color: str,
                           nsize: int = None) -> HighlightStyle:
        if len(fg_colors) != len(bg_colors):
            raise Exception(f'fg_colors={fg_colors} and bg_colors={bg_colors} do not have the same length')

        styles = []
        for i in range(0, len(bg_colors)):
            nstyle = copy.deepcopy(base_node_style)
            nstyle['fgcolor'] = fg_colors[i]
            nstyle['bgcolor'] = bg_colors[i]
            if nsize is not None:
                nstyle['size'] = nsize

            unknown_nstyle = copy.deepcopy(base_node_style)
            unknown_nstyle['fgcolor'] = unknown_fg_color
            unknown_nstyle['bgcolor'] = unknown_bg_color
            if nsize is not None:
                unknown_nstyle['size'] = nsize

            styles.append({
                'present_nstyle': nstyle,
                'unknown_nstyle': unknown_nstyle,
                'legend_color': bg_colors[i]
            })

        return HighlightStyle(styles, index=0)
