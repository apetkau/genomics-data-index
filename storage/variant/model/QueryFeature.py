from __future__ import annotations
import abc


class QueryFeature(abc.ABC):
    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def id(self):
        pass

    @property
    @abc.abstractmethod
    def scope(self):
        pass

    @abc.abstractmethod
    def to_unknown(self) -> QueryFeature:
        """
        Converts this given QueryFeature to one representing an unknown (e.g., unknown mutation or MLST allele).
        :return: The equivalent of this feature but representing an unknown.
        """
        pass