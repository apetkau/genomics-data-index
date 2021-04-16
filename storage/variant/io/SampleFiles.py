from __future__ import annotations
import abc


class SampleFiles(abc.ABC):

    def __init__(self, sample_name: str):
        self._name = sample_name

    @property
    def sample_name(self) -> str:
        return self._name

    @abc.abstractmethod
    def is_preprocessed(self) -> bool:
        pass

    def preprocess(self) -> SampleFiles:
        if not self.is_preprocessed():
            return self._do_preprocess()
        else:
            return self

    @abc.abstractmethod
    def _do_preprocess(self) -> SampleFiles:
        pass
