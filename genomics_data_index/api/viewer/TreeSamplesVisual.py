import abc
from typing import Union, Iterable, Tuple, Set

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
                 legend_fontsize: int):
        """
        Creates a new TreeSamplesVisual with the given information.
        :param samples: The set of samples to apply this visual style to.
        :param legend_nodesize: The legend node size for this set of samples.
        :param legend_fontsize: The legend font size for this set of samples.
        :return: A new TreeSamplesVisual.
        """
        self._samples = samples
        self._legend_nodesize = legend_nodesize
        self._legend_fontsize = legend_fontsize
        self._present_sample_names = None
        self._unknown_sample_names = None

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

    def _add_legend_entry(self, legend_label: str, legend_color: str, kind: str, tree_style: TreeStyle) -> None:
        if legend_label is not None:
            color_face, text_face = self._build_legend_item(color=legend_color,
                                                            legend_label=legend_label,
                                                            kind=kind)
            tree_style.legend.add_face(color_face, column=0)
            tree_style.legend.add_face(text_face, column=1)

    def _build_legend_item(self, color: str, legend_label: str, kind: str) -> Tuple[Face, Face]:
        if kind == 'rect' or kind == 'rectangle' or kind == 'r':
            cf = RectFace(width=self._legend_nodesize, height=self._legend_nodesize, bgcolor=color,
                          fgcolor='black')
        elif kind == 'circle' or kind == 'circ' or kind == 'c':
            cf = CircleFace(radius=self._legend_nodesize / 2, color=color)
        else:
            raise Exception(f'kind=[{kind}] must be one of {self.LEGEND_FACE_KINDS}')
        cf.hz_align = 2
        cf.margin_left = 3
        cf.margin_right = 3
        cf.margin_top = 3
        cf.margin_bottom = 3
        tf = TextFace(legend_label, fsize=self._legend_fontsize)
        tf.margin_left = 10
        tf.margin_right = 10
        return cf, tf
