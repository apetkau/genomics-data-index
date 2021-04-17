from __future__ import annotations

import abc
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class SampleFiles(abc.ABC):

    def __init__(self, sample_name: str):
        self._name = sample_name

    @property
    def sample_name(self) -> str:
        return self._name

    @abc.abstractmethod
    def is_preprocessed(self) -> bool:
        pass

    def preprocess(self, output_dir: Path) -> SampleFiles:
        return self.persist(output_dir)

    def persist(self, output_dir: Path) -> SampleFiles:
        if not self.is_preprocessed():
            logger.debug(f'Processing sample [{self.sample_name}] and saving to [{output_dir}]')
            return self._do_preprocess_and_persist(output_dir)
        else:
            return self._do_persist(output_dir)

    @abc.abstractmethod
    def _do_preprocess_and_persist(self, output_dir: Path) -> SampleFiles:
        pass

    @abc.abstractmethod
    def _do_persist(self, output_dir: Path) -> SampleFiles:
        pass
