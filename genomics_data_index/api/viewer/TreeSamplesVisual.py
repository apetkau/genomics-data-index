import abc
from typing import Union, Iterable, Tuple, Set, Dict

from ete3 import Tree, TreeStyle, Face, RectFace, CircleFace, TextFace

from genomics_data_index.api.query.SamplesQuery import SamplesQuery


class TreeSamplesVisual(abc.ABC):
    """
    Defines a visual element on a tree which applies to a set of samples (e.g., highlight or annotate).
    """

    LEGEND_FACE_KINDS = ['circle', 'circ', 'c',
                         'rectangle', 'rect', 'r']

    def __init__(self, samples: Union[SamplesQuery, Iterable[str]],
                 legend_nodesize: int,
                 legend_fontsize: int,
                 legend_columns: Dict[str, int]):
        """
        Creates a new TreeSamplesVisual with the given information.
        :param samples: The set of samples to apply this visual style to.
        :param legend_nodesize: The legend node size for this set of samples.
        :param legend_fontsize: The legend font size for this set of samples.
        :param legend_columns: A dict mapping legend item types to column numbers.
        :return: A new TreeSamplesVisual.
        """
        self._samples = samples
        self._legend_nodesize = legend_nodesize
        self._legend_fontsize = legend_fontsize
        self._present_sample_names = None
        self._unknown_sample_names = None
        self._legend_columns = legend_columns

    @property
    def present_sample_names(self) -> Set[str]:
        """
        Gets the names of the samples present in the selected set this visual style applies to.
        :return: A set of the names of present samples this visual style applies to.
        """
        if self._present_sample_names is None:
            if isinstance(self._samples, SamplesQuery):
                self._present_sample_names = set(self._samples.tolist(names=True))
            else:
                self._present_sample_names = set(self._samples)
        return self._present_sample_names

    @property
    def unknown_sample_names(self) -> Set[str]:
        """
        Gets the names of the samples unknown in the selected set this visual style applies to.
        :return: A set of the names of unknown samples this visual style applies to.
        """
        if self._unknown_sample_names is None:
            if isinstance(self._samples, SamplesQuery):
                self._unknown_sample_names = set(self._samples.select_unknown().tolist(names=True))
            else:
                # If we don't have a SamplesQuery, there are no unknown samples
                self._unknown_sample_names = set()
        return self._unknown_sample_names

    @abc.abstractmethod
    def apply_visual(self, tree: Tree, tree_style: TreeStyle) -> None:
        """
        Applies a visual style defined by the passed samples to the Tree and TreeStyle.
        Will modify tree and tree_style in-place (this avoids having to copy large trees
        which was causing me issues).

        :param tree: The Tree to apply visual information to.
        :param tree_style: The TreeStyle to apply visual information to.
        :returns: None. Modifies tree and tree_style in-place.
        """
        pass

    def _add_legend_entry(self, legend_label: str, legend_color: str, unknown_color: str,
                          include_unknown: bool, has_unknowns: bool, kind: str, tree_style: TreeStyle) -> None:
        if legend_label is not None:
            color_face, unknown_face, text_face = self._build_legend_item(color=legend_color,
                                                                          unknown_color=unknown_color,
                                                                          include_unknown=include_unknown,
                                                                          has_unknowns=has_unknowns,
                                                                          legend_label=legend_label,
                                                                          kind=kind)

            tree_style.legend.add_face(color_face, column=self._legend_columns['present'])
            if include_unknown:
                tree_style.legend.add_face(unknown_face, column=self._legend_columns['unknown'])
            tree_style.legend.add_face(text_face, column=self._legend_columns['text'])

    def _build_legend_item(self, color: str, unknown_color: str, include_unknown: bool, has_unknowns: bool,
                           legend_label: str, kind: str) -> Tuple[Face, Face, Face]:
        if kind == 'rect' or kind == 'rectangle' or kind == 'r':
            cf = RectFace(width=self._legend_nodesize, height=self._legend_nodesize, bgcolor=color,
                          fgcolor='black')
            if include_unknown:
                if has_unknowns:
                    ucf = RectFace(width=self._legend_nodesize, height=self._legend_nodesize, bgcolor=unknown_color,
                                   fgcolor='black')
                else:
                    ucf = RectFace(width=self._legend_nodesize, height=self._legend_nodesize, bgcolor=None,
                                   fgcolor=None)
            else:
                ucf = None
        elif kind == 'circle' or kind == 'circ' or kind == 'c':
            cf = CircleFace(radius=self._legend_nodesize / 2, color=color)
            if include_unknown:
                if has_unknowns:
                    ucf = CircleFace(radius=self._legend_nodesize / 2, color=unknown_color)
                else:
                    ucf = CircleFace(radius=self._legend_nodesize / 2, color=None)
            else:
                ucf = None
        else:
            raise Exception(f'kind=[{kind}] must be one of {self.LEGEND_FACE_KINDS}')
        cf.hz_align = 2
        cf.margin_left = 3
        cf.margin_right = 3
        cf.margin_top = 3
        cf.margin_bottom = 3
        if ucf is not None:
            ucf.hz_align = 2
            ucf.margin_left = 3
            ucf.margin_right = 3
            ucf.margin_top = 3
            ucf.margin_bottom = 3
        tf = TextFace(legend_label, fsize=self._legend_fontsize)
        tf.margin_left = 10
        tf.margin_right = 10
        return cf, ucf, tf
