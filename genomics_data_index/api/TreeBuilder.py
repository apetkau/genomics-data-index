from typing import Tuple

from ete3 import Tree

from genomics_data_index.variant.SampleSet import SampleSet


class TreeBuilder:

    def __init__(self):
        pass

    def build(self, samples_set: SampleSet, method: str, **kwargs) -> Tuple[Tree, int, SampleSet]:
        pass
