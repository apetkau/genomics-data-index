import abc
from typing import Tuple

from ete3 import Tree

from genomics_data_index.storage.SampleSet import SampleSet


class TreeBuilder(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def build(self, samples_set: SampleSet, method: str, **kwargs) -> Tuple[Tree, int, SampleSet]:
        pass
