from __future__ import annotations

import abc
import logging
from typing import Union, List, cast

import pandas as pd
from ete3 import Tree, TreeStyle

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.api.viewer.TreeStyler import TreeStyler, HighlightStyle
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class TreeSamplesQuery(WrappedSamplesQuery, abc.ABC):
    BUILD_TREE_KINDS = ['mutation']
    ISIN_TREE_TYPES = ['distance', 'mrca']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery, tree: Tree):
        super().__init__(connection=connection, wrapped_query=wrapped_query)
        self._tree = tree

    def _within_mrca(self, sample_names: Union[str, List[str]]) -> SamplesQuery:
        if isinstance(sample_names, str):
            sample_names = [sample_names]

        sample_name_ids = self._get_sample_name_ids()
        sample_leaves_list = []
        for name in sample_names:
            sample_leaves = self._tree.get_leaves_by_name(name)
            if len(sample_leaves) != 1:
                raise Exception(
                    f'Invalid number of matching leaves for sample [{name}], leaves {sample_leaves}')
            else:
                sample_leaves_list.append(sample_leaves[0])

        if len(sample_leaves_list) == 0:
            raise Exception(f'Should at least have some leaves in the tree matching sample_names={sample_names}')
        elif len(sample_leaves_list) == 1:
            found_sample_names = sample_names
        else:
            first_sample_leaf = sample_leaves_list.pop()
            ancestor_node = first_sample_leaf.get_common_ancestor(sample_leaves_list)
            found_sample_names = ancestor_node.get_leaf_names()

        found_samples_list = []
        for name in found_sample_names:
            if name in sample_name_ids:
                found_samples_list.append(sample_name_ids[name])
        found_samples = SampleSet(found_samples_list)
        return self.intersect(found_samples, f'within(mrca of {sample_names})')

    def _isin_internal(self, data: Union[str, List[str], pd.Series], kind: str, **kwargs) -> SamplesQuery:
        if kind == 'distance':
            return self._within_distance(sample_names=data, **kwargs)
        elif kind == 'mrca':
            return self._within_mrca(sample_names=data)
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self._isin_kinds()}')

    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        samples_query = self._wrapped_query.build_tree(kind=kind, **kwargs)

        if isinstance(samples_query, TreeSamplesQuery):
            samples_query = cast(TreeSamplesQuery, samples_query)
            return samples_query._wrap_create(self, self.universe_set)
        else:
            raise Exception(f'Build tree is not of type {TreeSamplesQuery.__class__}')

    def _within_distance(self, sample_names: Union[str, List[str]], distance: float,
                         units: str, **kwargs) -> SamplesQuery:
        if self._can_handle_distance_units(units):
            return self._within_distance_internal(sample_names=sample_names,
                                                  distance=distance, units=units)
        else:
            return self._wrap_create(self._wrapped_query._within_distance(sample_names=sample_names,
                                                                          distance=distance,
                                                                          units=units,
                                                                          **kwargs))

    @abc.abstractmethod
    def _within_distance_internal(self, sample_names: Union[str, List[str]], distance: float,
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

        return TreeStyler.create(tree=self._tree.copy(method='cpickle'),
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
