import copy
import logging
from typing import Union, Iterable, Dict, Any

from ete3 import Tree, TreeStyle, TextFace, RectFace, Face, CircleFace

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual

logger = logging.getLogger(__name__)


class AnnotateTreeSamplesVisual(TreeSamplesVisual):
    """
    Creates a visual style which adds boxes which are colored based on the presence/absence of samples in the set.
    """

    ANNOTATE_KINDS = TreeSamplesVisual.LEGEND_FACE_KINDS

    def __init__(self,
                 samples: Union[SamplesQuery, Iterable[str]],
                 label: Union[str, Dict[str, Any]],
                 label_unknown: Union[str, Dict[str, Any]],
                 annotate_show_box_label: bool,
                 annotate_box_label_color: str,
                 annotate_box_label_unknown_color: str,
                 annotate_label_fontsize: int,
                 legend_label: str,
                 box_width: int, box_height: int,
                 color_present: str, color_absent: str, color_unknown: str,
                 legend_nodesize: int,
                 legend_fontsize: int,
                 annotate_column: int,
                 annotate_kind: str,
                 annotate_border_color: str,
                 annotate_opacity_present: float,
                 annotate_opacity_absent: float,
                 annotate_opacity_unknown: float,
                 border_width: int,
                 margin: int,
                 include_unknown: bool,
                 legend_columns: Dict[str, int]):
        """
        Creates a new AnnotateTreeSamplesVisual with the given information. This is not intended to be created directly
        but instead is created via methods in the :py:class:`genomics_data_index.api.viewer.TreeStyler` class.
        :param samples: The set of samples this visual style applies to.
        :param label: A label for the this set of samples.
        :param label_unknown: A label for the unknown set of samples.
        :param annotate_show_box_label: Whether or not the label should be printed in each annotate box.
        :param annotate_box_label_color: The color of the annotate box label if it is to be displayed.
        :param annotate_box_label_unknown_color: The color of the annotate box label (for unknowns) if it is to be displayed.
        :param annotate_label_fontsize: The font size of the annotate box label if it is to be displayed.
        :param legend_label: The legend label.
        :param box_width: The width of the annotate box.
        :param box_height: The height of the annotate box.
        :param color_present: The color of the samples present in this set.
        :param color_absent: The color of the samples absent in this set.
        :param color_unknown: The color of the samples which are unknown if present/absent in this set.
        :param legend_nodesize: The node size for the legend.
        :param legend_fontsize: The legend font size.
        :param annotate_column: Which column in the annotate table this annotate visual belongs to (starting from column 0).
        :param annotate_kind: The type of annotation to display for present samples ['circle', 'rectangle'].
        :param annotate_border_color: The border color of teh annotate box.
        :param annotate_opacity_present: The opacity of the present samples.
        :param annotate_opacity_absent: The opacity of the absent samples.
        :param annotate_opacity_unknown: The opacity of the unknown samples.
        :param border_width: The width of the border.
        :param margin: The width of the margin in each annotate box.
        :param include_unknown: Whether or not unknowns should be included/represented.
        :param legend_columns: A dict mapping legend item types to column numbers.
        :return: A new AnnotateTreeSamplesVisual.
        """
        super().__init__(samples, legend_nodesize=legend_nodesize,
                         legend_fontsize=legend_fontsize, legend_columns=legend_columns)
        self._annotate_show_box_label = annotate_show_box_label
        self._annotate_label_fontsize = annotate_label_fontsize
        self._label = label
        self._legend_label = legend_label
        self._box_width = box_width
        self._box_height = box_height
        self._color_present = color_present
        self._color_absent = color_absent
        self._color_unknown = color_unknown
        self._annotate_column = annotate_column

        if annotate_kind not in self.ANNOTATE_KINDS:
            raise Exception(f'annotate_kind=[{annotate_kind}] is invalid. Must be one of {self.ANNOTATE_KINDS}')
        else:
            self._annotate_kind = annotate_kind

        self._annotate_border_color = annotate_border_color
        self._annotate_opacity_present = annotate_opacity_present
        self._annotate_opacity_absent = annotate_opacity_absent
        self._annotate_opacity_unknown = annotate_opacity_unknown
        self._annotate_border_width = border_width
        self._annotate_margin = margin
        self._include_unknown = include_unknown

        if self._annotate_show_box_label:
            self._box_label_present = copy.deepcopy(label)
            self._box_label_unknown = copy.deepcopy(label_unknown)
            if self._box_label_present is not None and annotate_box_label_color is not None:
                self._box_label_present['color'] = annotate_box_label_color
            if self._box_label_unknown is not None and annotate_box_label_unknown_color is not None:
                self._box_label_unknown['color'] = annotate_box_label_unknown_color
        else:
            self._box_label_present = None
            self._box_label_unknown = None

    def apply_visual(self, tree: Tree, tree_style: TreeStyle) -> None:
        if self._label is not None:
            text = self._label.get('text', None)
            fsize = self._label.get('fontsize', self._annotate_label_fontsize)
            ftype = self._label.get('font', 'Verdana')
            color = self._label.get('color', 'black')
            tf = TextFace(text, fsize=fsize, ftype=ftype, fgcolor=color)
            tf.margin_bottom = 10
            tf.margin_left = 10
            tf.margin_right = 10
            tf.margin_top = 10
            tf.hz_align = 1
            tree_style.aligned_header.add_face(tf, self._annotate_column)
            tree_style.aligned_foot.add_face(tf, self._annotate_column)

        # Annotate nodes
        unknown_samples_count = 0
        for leaf in tree.iter_leaves():
            if leaf.name in self.present_sample_names:
                annotate_face = self._build_annotate_face(width=self._box_width, height=self._box_height,
                                                          border_color=self._annotate_border_color,
                                                          bgcolor=self._color_present,
                                                          opacity=self._annotate_opacity_present,
                                                          label=self._box_label_present)
            elif self._include_unknown and (leaf.name in self.unknown_sample_names):
                annotate_face = self._build_annotate_face(width=self._box_width, height=self._box_height,
                                                          border_color=self._annotate_border_color,
                                                          bgcolor=self._color_unknown,
                                                          opacity=self._annotate_opacity_unknown,
                                                          label=self._box_label_unknown)
                unknown_samples_count = unknown_samples_count + 1
            else:
                annotate_face = self._build_annotate_face(width=self._box_width, height=self._box_height,
                                                          border_color=self._annotate_border_color,
                                                          bgcolor=self._color_absent,
                                                          opacity=self._annotate_opacity_absent,
                                                          label=None)

            leaf.add_face(annotate_face, column=self._annotate_column, position='aligned')

        has_unknowns = unknown_samples_count > 0
        self._add_legend_entry(self._legend_label, legend_color=self._color_present,
                               unknown_color=self._color_unknown, include_unknown=self._include_unknown,
                               has_unknowns=has_unknowns,
                               kind=self._annotate_kind, tree_style=tree_style)

    def _build_annotate_face(self, width: int, height: int, border_color: str, bgcolor: str,
                             opacity: float, label: Union[str, Dict[str, Any]] = None) -> Face:
        if self._annotate_kind == 'rect' or self._annotate_kind == 'rectangle' or self._annotate_kind == 'r':
            rf = RectFace(width=width, height=height, fgcolor=None, bgcolor=bgcolor, label=label)
            rf.border.width = self._annotate_border_width
            rf.margin_top = self._annotate_margin
            rf.margin_bottom = self._annotate_margin
            rf.margin_left = self._annotate_margin
            rf.margin_right = self._annotate_margin
            rf.border.color = border_color
            rf.background.color = bgcolor
            rf.opacity = opacity
            rf.hz_align = 1
            rf.vt_align = 1
            return rf
        elif self._annotate_kind == 'circle' or self._annotate_kind == 'circ' or self._annotate_kind == 'c':
            # Make circle radius such that it fits in bounding box defined by width and height
            min_dimension = min(width, height)
            radius = min_dimension / 2

            cf = CircleFace(radius=radius, color=bgcolor, label=label)
            cf.border.width = self._annotate_border_width
            cf.margin_top = self._annotate_margin
            cf.margin_bottom = self._annotate_margin
            cf.margin_left = self._annotate_margin
            cf.margin_right = self._annotate_margin
            cf.border.color = border_color
            cf.opacity = opacity
            cf.hz_align = 1
            cf.vt_align = 1
            return cf
        else:
            raise Exception(f'Invalid value for annotate_kind={self._annotate_kind}.'
                            f' Must be one of {self.LEGEND_FACE_KINDS}')
