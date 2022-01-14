import abc
import logging
from typing import Union, Iterable, Dict, Any, Optional

from ete3 import Tree, TreeStyle, TextFace, RectFace, Face, CircleFace, TreeNode

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.viewer.TreeSamplesVisual import TreeSamplesVisual

logger = logging.getLogger(__name__)


class AbstractAnnotationTreeSamplesVisual(TreeSamplesVisual, abc.ABC):
    """
    A tree style for adding annotations (columns next to leaves of a tree).
    """

    def __init__(self, samples: Union[SamplesQuery, Iterable[str]],
                 legend_nodesize: int,
                 legend_fontsize: int,
                 legend_columns: Dict[str, int],
                 annotate_kind: str,
                 annotate_border_width: Optional[int],
                 annotate_margin: Optional[int],
                 ):
        super().__init__(samples=samples,
                         legend_nodesize=legend_nodesize,
                         legend_fontsize=legend_fontsize,
                         legend_columns=legend_columns)

        self._annotate_kind = annotate_kind
        self._annotate_border_width = annotate_border_width,
        self._annotate_margin = annotate_margin

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
