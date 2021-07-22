from __future__ import annotations

import abc
import logging
import uuid
from pathlib import Path

from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)

persistence_namespace = uuid.UUID('260bb0a7-6573-4282-bcdb-f4008aa42af8')


def generate_sample_persistence_name(sample_name: str) -> str:
    return uuid.uuid3(namespace=persistence_namespace, name=sample_name).hex


class SampleData(abc.ABC):

    def __init__(self, sample_name: str):
        self._name = sample_name
        self._sample_name_persistence = generate_sample_persistence_name(sample_name)

    @property
    def sample_name(self) -> str:
        return self._name

    @property
    def sample_name_persistence(self) -> str:
        return self._sample_name_persistence

    @abc.abstractmethod
    def is_preprocessed(self) -> bool:
        pass

    def preprocess(self, output_dir: Path) -> SampleData:
        return self.persist(output_dir)

    def persist(self, output_dir: Path) -> SampleData:
        if not self.is_preprocessed():
            logger.log(TRACE_LEVEL, f'Processing sample [{self.sample_name}] and saving to [{output_dir}]')
            return self._do_preprocess_and_persist(output_dir)
        else:
            return self._do_persist(output_dir)

    @abc.abstractmethod
    def _do_preprocess_and_persist(self, output_dir: Path) -> SampleData:
        pass

    @abc.abstractmethod
    def _do_persist(self, output_dir: Path) -> SampleData:
        pass
