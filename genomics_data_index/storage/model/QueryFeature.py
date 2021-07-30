from __future__ import annotations

import abc
from typing import List, Optional


class QueryFeature(abc.ABC):
    WILD = '*'
    SPLIT_CHAR = ':'

    def __init__(self):
        pass

    @property
    def id(self) -> Optional[str]:
        return f'{self.prefix}{self.id_no_prefix}'

    @property
    @abc.abstractmethod
    def id_no_prefix(self) -> Optional[str]:
        pass

    def is_wild(self) -> bool:
        return self.WILD in self.id

    @abc.abstractmethod
    def is_unknown(self) -> bool:
        pass

    @property
    @abc.abstractmethod
    def scope(self) -> str:
        pass

    @abc.abstractmethod
    def to_unknown(self) -> QueryFeature:
        """
        Converts this given QueryFeature to one representing an unknown (e.g., unknown mutation or MLST allele).
        :return: The equivalent of this feature but representing an unknown.
        """
        pass

    @abc.abstractmethod
    def to_unknown_explode(self) -> List[QueryFeature]:
        """
        Converts this given QueryFeature to one representing an unknown (e.g., unknown mutation or MLST allele).
        This will explode into separate unknown features (e.g., in cases of an indel/complex mutation).
        :return: A list of unknown features, where each feature corresponds to a single position.
        """
        pass

    @property
    def prefix(self) -> str:
        return ''

    def __str__(self) -> str:
        return self.id
