from typing import List

import abc


class DataRepresentation(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def get_alleles(self, samples: List[str]):
        pass
