import logging

from ete3 import Tree

from genomics_data_index.ete3treeview import TreeStyle
from genomics_data_index.api.viewer.samples_visuals.AbstractAnnotationTreeSamplesVisual import \
    AbstractAnnotationTreeSamplesVisual

logger = logging.getLogger(__name__)


class AnnotationSpacingTreeSamplesVisual(AbstractAnnotationTreeSamplesVisual):
    """
    Creates a visual style which adds empty spaces in annotation table.
    """

    def __init__(self, width: int, height: int, annotate_column: int, color: str = None,
                 border_width: int = None, border_color: str = None):
        super().__init__(samples=set(),
                         legend_nodesize=0,
                         legend_fontsize=0,
                         legend_columns={},
                         annotate_kind='r',
                         annotate_border_width=border_width,
                         annotate_margin=0)
        self._box_width = width
        self._box_height = height
        self._annotate_column = annotate_column
        self._color = color
        self._annotation_opacity = 0.0
        self._annotate_border_color = border_color

    def apply_visual(self, tree: Tree, tree_style: TreeStyle) -> None:
        for leaf in tree.iter_leaves():
            annotate_face = self._build_annotate_face(width=self._box_width, height=self._box_height,
                                                      border_color=self._annotate_border_color,
                                                      bgcolor=self._color,
                                                      opacity=self._annotation_opacity,
                                                      label=None)

            leaf.add_face(annotate_face, column=self._annotate_column, position='aligned')
