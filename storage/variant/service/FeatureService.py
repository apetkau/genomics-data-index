import abc
import logging
from pathlib import Path
from typing import List, Set, Any, Dict, Tuple, Optional

import pandas as pd

from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.model.db import Sample
from storage.variant.service import DatabaseConnection
from storage.variant.service import EntityExistsError
from storage.variant.service.SampleService import SampleService
from storage.variant.io.SampleFiles import SampleFiles

logger = logging.getLogger(__name__)

AUTO_SCOPE = '__AUTO_SCOPE__'


class FeatureService(abc.ABC):

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService):
        self._connection = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir

    @abc.abstractmethod
    def get_correct_reader(self) -> Any:
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

    def _create_feature_objects(self, features_df: pd.DataFrame, sample_names: Set[str]) -> List[Any]:
        sample_name_ids = self._sample_service.find_sample_name_ids(sample_names)

        features_df['_FEATURE_ID'] = features_df.apply(self._create_feature_identifier, axis='columns')
        features_df['_SAMPLE_ID'] = self._get_sample_id_series(features_df, sample_name_ids)

        index_df = features_df.groupby('_FEATURE_ID').agg(self.aggregate_feature_column()).reset_index()

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

    def _verify_correct_reader(self, features_reader: FeaturesReader) -> None:
        if not isinstance(features_reader, self.get_correct_reader()):
            raise Exception(f'features_reader=[{features_reader}] is not of type'
                            f' {self.get_correct_reader().__name__}')

    def _get_or_create_sample(self, sample_name: str) -> Sample:
        if self._sample_service.exists(sample_name):
            return self._sample_service.get_sample(sample_name)
        else:
            return Sample(name=sample_name)

    def _verify_correct_feature_scope(self, feature_scope_name: str) -> None:
        if feature_scope_name is None:
            raise Exception('feature_scope_name cannot be None')

    def progress_hook(self, number: int, print_every: int, total: int) -> None:
        if number % print_every == 0:
            logger.info(f'Proccessed {number/total*100:0.0f}% ({number}/{total}) samples')

    def insert(self, features_reader: FeaturesReader, feature_scope_name: str = AUTO_SCOPE) -> None:
        self._verify_correct_reader(features_reader=features_reader)
        self._verify_correct_feature_scope(feature_scope_name)

        sample_names = features_reader.samples_set()
        num_samples = len(sample_names)

        if self.check_samples_have_features(sample_names, feature_scope_name):
            raise EntityExistsError(f'Passed samples already have features for feature scope [{feature_scope_name}], '
                                    f'will not insert any new features')

        interval = max(1, int(num_samples / 50))
        processed_samples = 0
        persisted_sample_files_dict = {}
        for sample_files in features_reader.iter_sample_files():
            self.progress_hook(processed_samples, print_every=interval, total=num_samples)
            logger.debug(f'Loading sample {sample_files.sample_name}')

            sample = self._get_or_create_sample(sample_files.sample_name)
            persisted_sample_files = self._persist_sample_files(sample_files)
            sample_feature_object = self.build_sample_feature_object(sample=sample,
                                                                     sample_files=persisted_sample_files,
                                                                     features_reader=features_reader,
                                                                     feature_scope_name=feature_scope_name)

            persisted_sample_files_dict[persisted_sample_files.sample_name] = persisted_sample_files
            self._connection.get_session().add(sample_feature_object)
            processed_samples += 1
        self._connection.get_session().commit()
        logger.info(f'Finished processings {num_samples} samples')

        persisted_features_reader = self._create_persisted_features_reader(sample_files_dict=persisted_sample_files_dict,
                                                                           features_reader=features_reader)
        self.index_features(features_reader=persisted_features_reader, feature_scope_name=feature_scope_name)

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
    def build_sample_feature_object(self, sample: Sample,
                                    sample_files: SampleFiles,
                                    features_reader: FeaturesReader, feature_scope_name: str) -> Any:
        pass

    def _persist_sample_files(self, sample_files: SampleFiles) -> Optional[SampleFiles]:
        if sample_files is not None:
            return sample_files.persist(self._features_dir)
        else:
            return None

    @abc.abstractmethod
    def _create_persisted_features_reader(self, sample_files_dict: Dict[str, SampleFiles],
                                          features_reader: FeaturesReader) -> FeaturesReader:
        pass
