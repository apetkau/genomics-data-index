import abc
import logging
from pathlib import Path
from typing import List, Set, Any, Dict, Optional, Generator

import pandas as pd

from genomics_data_index.storage.io.FeaturesReader import FeaturesReader
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.SampleService import SampleService

logger = logging.getLogger(__name__)

AUTO_SCOPE = '__AUTO_SCOPE__'


class FeatureService(abc.ABC):

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService,
                 max_insert_batch_size: int = 500):
        self._connection = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir
        self._max_insert_batch_size = max_insert_batch_size
        self._min_insert_batch_size = 5

    @abc.abstractmethod
    def get_correct_data_package(self) -> Any:
        pass

    @abc.abstractmethod
    def get_correct_sample_data(self) -> Any:
        pass

    @abc.abstractmethod
    def _create_feature_identifier(self, features_df: pd.DataFrame) -> str:
        pass

    @abc.abstractmethod
    def _get_sample_id_series(self, features_df: pd.DataFrame, sample_name_ids: Dict[str, int]) -> pd.Series:
        pass

    @abc.abstractmethod
    def _create_feature_object(self, features_df: pd.DataFrame):
        pass

    @abc.abstractmethod
    def aggregate_feature_column(self) -> Dict[str, Any]:
        pass

    def _modify_df_types(self, features_df: pd.DataFrame) -> pd.DataFrame:
        return features_df

    def _create_feature_objects(self, features_df: pd.DataFrame, sample_names: Set[str]) -> List[Any]:
        sample_name_ids = self._sample_service.find_sample_name_ids(sample_names)

        features_df['_FEATURE_ID'] = features_df.apply(self._create_feature_identifier, axis='columns')
        features_df['_SAMPLE_ID'] = self._get_sample_id_series(features_df, sample_name_ids)

        index_df = features_df.groupby('_FEATURE_ID').agg(self.aggregate_feature_column())
        index_df = self._modify_df_types(index_df).reset_index()

        return index_df.apply(self._create_feature_object, axis='columns').tolist()

    @abc.abstractmethod
    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        """
        Checks if any of the passed sample names already have features.
        :param sample_names: The set of sample names.
        :param feature_scope_name: A name defining the feature scope (i.e., reference genome name).
        :return: True if any of the samples have features, false otherwise.
        """
        pass

    def _verify_correct_data_package(self, data_package: SampleDataPackage) -> None:
        if not isinstance(data_package, self.get_correct_data_package()):
            raise Exception(f'data_package=[{data_package}] is not of type'
                            f' {self.get_correct_data_package().__name__}')

    def _verify_correct_sample_data(self, sample_data: SampleData) -> None:
        if not isinstance(sample_data, self.get_correct_sample_data()):
            raise Exception(f'sample_data=[{sample_data}] is not of type'
                            f' {self.get_correct_sample_data().__name__}')

    def _get_or_create_sample(self, sample_name: str) -> Sample:
        if self._sample_service.exists(sample_name):
            return self._sample_service.get_sample(sample_name)
        else:
            return Sample(name=sample_name)

    def _get_or_create_samples(self, sample_names: List[str]) -> Dict[str, Sample]:
        existing_samples = self._sample_service.get_existing_samples_by_names(sample_names)
        samples_dict = {sample.name: sample for sample in existing_samples}

        for name in sample_names:
            if name not in samples_dict:
                samples_dict[name] = Sample(name=name)

        return samples_dict

    def _verify_correct_feature_scope(self, feature_scope_name: str) -> None:
        if feature_scope_name is None:
            raise Exception('feature_scope_name cannot be None')

    def log_progress(self, number: int, total: int) -> None:
        logger.info(f'Proccessed {number / total * 100:0.0f}% ({number}/{total}) samples')

    def _set_batch_size(self, num_samples: int) -> int:
        batch_size = max(self._min_insert_batch_size, int(num_samples / 50))
        batch_size = min(self._max_insert_batch_size, batch_size)
        return batch_size

    def _sample_data_iter_batch(self, data_package: SampleDataPackage,
                                batch_size: int) -> Generator[Dict[str, SampleData], None, None]:
        sample_data_batch = {}
        for sample_data in data_package.iter_sample_data():
            sample_data_batch[sample_data.sample_name] = sample_data

            if len(sample_data_batch) >= batch_size:
                yield sample_data_batch
                sample_data_batch = {}

        if len(sample_data_batch) > 0:
            yield sample_data_batch

    def insert(self, data_package: SampleDataPackage, feature_scope_name: str = AUTO_SCOPE) -> None:
        self._verify_correct_data_package(data_package=data_package)
        self._verify_correct_feature_scope(feature_scope_name)

        sample_names = data_package.sample_names()
        num_samples = len(sample_names)

        if self.check_samples_have_features(sample_names, feature_scope_name):
            raise EntityExistsError(f'Passed samples already have features for feature scope [{feature_scope_name}], '
                                    f'will not insert any new features')

        batch_size = self._set_batch_size(num_samples)
        logger.debug(f'Batch size {batch_size}')
        processed_samples = 0
        persisted_sample_data_dict = {}
        self.log_progress(processed_samples, total=num_samples)
        for sample_data_batch in self._sample_data_iter_batch(data_package, batch_size=batch_size):
            persisted_sample_data = self._handle_batch(sample_data_batch, feature_scope_name)
            persisted_sample_data_dict.update(persisted_sample_data)
            processed_samples += len(sample_data_batch)
            self.log_progress(processed_samples, total=num_samples)

        self._connection.get_session().commit()
        logger.info(f'Finished processing {num_samples} samples')

        persisted_features_reader = self._create_persisted_features_reader(sample_data_dict=persisted_sample_data_dict,
                                                                           data_package=data_package)
        self.index_features(features_reader=persisted_features_reader, feature_scope_name=feature_scope_name)

    def _handle_batch(self, sample_data_batch: Dict[str, SampleData], feature_scope_name: str) -> Dict[str, SampleData]:
        samples_dict = self._get_or_create_samples(list(sample_data_batch.keys()))
        persisted_sample_data_dict = {}
        for sample_name in samples_dict:
            logger.debug(f'Persisting sample {sample_name}')

            sample = samples_dict[sample_name]
            persisted_sample_data = self._persist_sample_data(sample_data_batch[sample.name])
            sample_feature_object = self.build_sample_feature_object(sample=sample, sample_data=persisted_sample_data,
                                                                     feature_scope_name=feature_scope_name)

            persisted_sample_data_dict[persisted_sample_data.sample_name] = persisted_sample_data
            self._connection.get_session().add(sample_feature_object)

        return persisted_sample_data_dict

    def _update_scope(self, features_df: pd.DataFrame, feature_scope_name: str) -> pd.DataFrame:
        return features_df

    def index_features(self, features_reader: FeaturesReader, feature_scope_name: str) -> None:
        logger.info('Indexing features from all samples')
        features_df = features_reader.get_features_table()
        features_df = self._update_scope(features_df, feature_scope_name)
        sample_names = features_reader.samples_set()
        self._connection.get_session().bulk_save_objects(self._create_feature_objects(features_df, sample_names))
        self._connection.get_session().commit()
        logger.info('Finished indexing features from all samples')

    @abc.abstractmethod
    def build_sample_feature_object(self, sample: Sample, sample_data: SampleData, feature_scope_name: str) -> Any:
        pass

    def _persist_sample_data(self, sample_data: SampleData) -> Optional[SampleData]:
        if sample_data is not None:
            return sample_data.persist(self._features_dir)
        else:
            return None

    @abc.abstractmethod
    def _create_persisted_features_reader(self, sample_data_dict: Dict[str, SampleData],
                                          data_package: SampleDataPackage) -> FeaturesReader:
        pass
