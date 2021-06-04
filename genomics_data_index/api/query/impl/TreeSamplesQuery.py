from __future__ import annotations

import abc
import logging
from typing import Union, List, cast

from ete3 import Tree, TreeStyle

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.api.viewer.TreeStyler import TreeStyler, HighlightStyle
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class TreeSamplesQuery(WrappedSamplesQuery, abc.ABC):
    """
    A class to handle queries that are joined to a tree. Subclasses are used for different implementations
    depending on what sort of data was used to build the tree (e.g., mutations or kmers).
    """

    BUILD_TREE_KINDS = ['mutation']
    ISIN_TREE_TYPES = ['distance', 'mrca']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree):
        """
        Builds a new TreeSamplesQuery from the given information. In most normal operations SamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery (e.g., build_tree() or join_tree()).

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around and join to the passed tree.
        :param tree: The tree to join to this query.
        :return: A new TreeSamplesQuery object.
        """
        super().__init__(connection=connection, wrapped_query=wrapped_query)
        self._tree = tree

    def _within_mrca(self, data: Union[str, List[str], SamplesQuery, SampleSet]) -> SamplesQuery:
        if isinstance(data, str):
            data = [data]

        sample_names, query_infix = self._get_sample_names_query_infix_from_data(data)
        if len(sample_names) == 0:
            found_samples = SampleSet.create_empty()
        else:
            sample_name_ids_self = self._get_sample_name_ids()
            sample_leaves_list = []
            for name in sample_names:
                sample_leaves = self._tree.get_leaves_by_name(name)
                if len(sample_leaves) != 1:
                    raise Exception(
                        f'Invalid number of matching leaves for sample [{name}], leaves {sample_leaves}')
                else:
                    sample_leaves_list.append(sample_leaves[0])

            if len(sample_leaves_list) == 0:
                raise Exception(f'Should at least have some leaves in the tree matching data={data}')
            elif len(sample_leaves_list) == 1:
                found_sample_names = [sample_leaves_list[0].name]
            else:
                # According to https://github.com/etetoolkit/ete/issues/503 'get_common_ancestor' does *not*
                # Include the first node calling the function so I have to leave it in the list.
                first_sample_leaf = sample_leaves_list[0]
                ancestor_node = first_sample_leaf.get_common_ancestor(sample_leaves_list)
                found_sample_names = ancestor_node.get_leaf_names()

            found_samples_list = []
            for name in found_sample_names:
                if name in sample_name_ids_self:
                    found_samples_list.append(sample_name_ids_self[name])
            found_samples = SampleSet(found_samples_list)
        return self.intersect(found_samples, f'within(mrca of {query_infix})')

    def _isin_internal(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str, **kwargs) -> SamplesQuery:
        if kind == 'distance':
            return self._within_distance(data=data, **kwargs)
        elif kind == 'mrca':
            return self._within_mrca(data=data)
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self._isin_kinds()}')

    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        samples_query = self._wrapped_query.build_tree(kind=kind, **kwargs)

        if isinstance(samples_query, TreeSamplesQuery):
            samples_query = cast(TreeSamplesQuery, samples_query)
            return samples_query._wrap_create(self, self.universe_set)
        else:
            raise Exception(f'Build tree is not of type {TreeSamplesQuery.__class__}')

    def join_tree(self, tree: Tree, kind='mutation', **kwargs) -> SamplesQuery:
        samples_query = self._wrapped_query.join_tree(tree=tree, kind=kind, **kwargs)

        if isinstance(samples_query, TreeSamplesQuery):
            samples_query = cast(TreeSamplesQuery, samples_query)
            return samples_query._wrap_create(self, self.universe_set)
        else:
            raise Exception(f'Build tree is not of type {TreeSamplesQuery.__class__}')

    def _within_distance(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float,
                         units: str, **kwargs) -> SamplesQuery:
        if self._can_handle_distance_units(units):
            return self._within_distance_internal(data=data,
                                                  distance=distance, units=units)
        else:
            return self._wrap_create(self._wrapped_query._within_distance(data=data,
                                                                          distance=distance,
                                                                          units=units,
                                                                          **kwargs))

    @abc.abstractmethod
    def _within_distance_internal(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float,
                                  units: str) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _can_handle_distance_units(self, units: str) -> bool:
        pass

    def _isin_kinds(self) -> List[str]:
        return list(set(super()._isin_kinds() + self.ISIN_TREE_TYPES))

    def _can_handle_isin_kind(self, kind: str) -> bool:
        return kind in self.ISIN_TREE_TYPES

    def tree_styler(self,
                    initial_style: TreeStyle = None,
                    mode='r',
                    highlight_style: Union[str, HighlightStyle] = 'light',
                    legend_nsize: int = 20, legend_fsize: int = 11,
                    annotate_color_present: str = 'black',
                    annotate_color_absent: str = 'white',
                    annotate_opacity_present: float = 1.0,
                    annotate_opacity_absent: float = 0.0,
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
                    tree_scale: float = None) -> TreeStyler:
        """
        Constructs a new :py:class:`genomics_data_index.api.viewer.TreeStyler` object used to style and visualize trees.
        All parameters listed below are optional.
        :param initial_style: The initial ete3.TreeStyle to start with.
        :param mode: Either 'r' (rectangular) or 'c' (circular).
        :param highlight_style: A style used to define how the highlight() method should behave.
                                Can either be one of the named highlight styles ['light', 'light_hn', 'pastel', 'medium', dark']
                                or an instance of a :py:class:`genomics_data_index.api.viewer.TreeStyler.HighlightStyle`.
        :param legend_nsize: The legend node size.
        :param legend_fsize: The legend font size.
        :param annotate_color_present: The default color of samples which are present in the set for the annotate() method.
        :param annotate_color_absent: The default color of samples which are absent in the set for the annotate() method.
        :param annotate_opacity_present: The default opacity of samples which are present in the set for the annotate() method.
        :param annotate_opacity_absent: The default opacity of samples which are absent in the set for the annotate() method.
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
        :param show_leaf_names: True if leaf names should be shown on the tree, False otherwise.
        :param tree_scale: A scale factor for the tree.
        :return: A new :py:class:`genomics_data_index.api.viewer.TreeStyler` object used to style and visualize trees.
        """
        return TreeStyler.create(tree=self._tree,
                                 initial_style=initial_style,
                                 mode=mode,
                                 highlight_style=highlight_style,
                                 legend_nsize=legend_nsize,
                                 legend_fsize=legend_fsize,
                                 annotate_color_present=annotate_color_present,
                                 annotate_color_absent=annotate_color_absent,
                                 annotate_opacity_present=annotate_opacity_present,
                                 annotate_opacity_absent=annotate_opacity_absent,
                                 annotate_border_color=annotate_border_color,
                                 annotate_kind=annotate_kind,
                                 annotate_box_width=annotate_box_width,
                                 annotate_box_height=annotate_box_height,
                                 annotate_border_width=annotate_border_width,
                                 annotate_margin=annotate_margin,
                                 annotate_guiding_lines=annotate_guiding_lines,
                                 annotate_guiding_lines_color=annotate_guiding_lines_color,
                                 figure_margin=figure_margin,
                                 show_border=show_border,
                                 title=title,
                                 title_fsize=title_fsize,
                                 legend_title=legend_title,
                                 annotate_show_box_label=annotate_show_box_label,
                                 annotate_box_label_color=annotate_box_label_color,
                                 annotate_arc_span=annotate_arc_span,
                                 annotate_label_fontsize=annotate_label_fontsize,
                                 show_leaf_names=show_leaf_names,
                                 tree_scale=tree_scale)

    @property
    def tree(self):
        return self._tree

    def has_tree(self) -> bool:
        return True
