import abc
import logging
from pathlib import Path
from typing import List, Set, Any, Dict, Optional, Generator

import pandas as pd

from genomics_data_index.storage.io.FeaturesReader import FeaturesReader
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.model.db import Sample, FeatureSamples
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.util import TRACE_LEVEL
from genomics_data_index.storage.util.SamplesProgressLogger import SamplesProgressLogger

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
    def read_index(self, feature_ids: List[str]) -> Dict[str, FeatureSamples]:
        """
        Reads the passed list of feature IDs from the database and returns a dictionary containing those features that
        already exist. This dictionary maps {'id' => 'feature_object_in_database'}.
        :param feature_ids: The list of feature object ids to search through.
        :return: A dictionary containing any features that already existing in the database index.
        """
        pass

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
    def _create_feature_object(self, features_df: pd.DataFrame) -> FeatureSamples:
        pass

    @abc.abstractmethod
    def aggregate_feature_column(self) -> Dict[str, Any]:
        pass

    def _modify_df_types(self, features_df: pd.DataFrame) -> pd.DataFrame:
        return features_df

    def _create_feature_objects(self, features_df: pd.DataFrame, sample_names: Set[str]) -> List[FeatureSamples]:
        logger.debug(f'Aggregating samples by feature identifiers within a table of {len(features_df)} rows'
                     f' across {len(sample_names)} samples')

        sample_name_ids = self._sample_service.find_sample_name_ids(sample_names)

        features_df['_FEATURE_ID'] = features_df.apply(self._create_feature_identifier, axis='columns')
        features_df['_SAMPLE_ID'] = self._get_sample_id_series(features_df, sample_name_ids)

        index_df = features_df.groupby('_FEATURE_ID').agg(self.aggregate_feature_column())
        index_df = self._modify_df_types(index_df).reset_index()

        feature_objects_list = index_df.apply(self._create_feature_object, axis='columns').tolist()

        logger.debug(f'Finished aggregating {len(sample_names)} samples into a table of {len(index_df)} features '
                     f'from an original table size of {len(features_df)} rows.')
        return feature_objects_list

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
        progress_logger = SamplesProgressLogger(stage_name='Insert', stage_number=1, total_samples=num_samples)

        if self.check_samples_have_features(sample_names, feature_scope_name):
            raise EntityExistsError(f'Passed samples already have features for feature scope [{feature_scope_name}], '
                                    f'will not insert any new features')

        batch_size = self._set_batch_size(num_samples)
        logger.debug(f'Batch size {batch_size}')
        processed_samples = 0
        persisted_sample_data_dict = {}
        progress_logger.update_progress(processed_samples)
        for sample_data_batch in self._sample_data_iter_batch(data_package, batch_size=batch_size):
            persisted_sample_data = self._handle_batch(sample_data_batch, feature_scope_name)
            persisted_sample_data_dict.update(persisted_sample_data)
            processed_samples += len(sample_data_batch)
            progress_logger.update_progress(processed_samples)

        self._connection.get_session().commit()
        logger.info(f'Finished processing {num_samples} samples')

        persisted_features_reader = self._create_persisted_features_reader(sample_data_dict=persisted_sample_data_dict,
                                                                           data_package=data_package)
        self.index_features(features_reader=persisted_features_reader, feature_scope_name=feature_scope_name)

    def _handle_batch(self, sample_data_batch: Dict[str, SampleData], feature_scope_name: str) -> Dict[str, SampleData]:
        samples_dict = self._get_or_create_samples(list(sample_data_batch.keys()))
        persisted_sample_data_dict = {}
        for sample_name in samples_dict:
            logger.log(TRACE_LEVEL, f'Persisting sample {sample_name}')

            sample = samples_dict[sample_name]
            persisted_sample_data = self._persist_sample_data(sample_data_batch[sample.name])
            sample_feature_object = self.build_sample_feature_object(sample=sample, sample_data=persisted_sample_data,
                                                                     feature_scope_name=feature_scope_name)

            persisted_sample_data_dict[persisted_sample_data.sample_name] = persisted_sample_data
            self._connection.get_session().add(sample_feature_object)

        return persisted_sample_data_dict

    def _update_scope(self, features_df: pd.DataFrame, feature_scope_name: str) -> pd.DataFrame:
        return features_df

    def update_feature_objects_inplace(self, feature_objects_dict_dest: Dict[str, FeatureSamples],
                                       feature_objects_dict_source: Dict[str, FeatureSamples]) -> None:
        """
        Updates matching features from feature_objects_dict_source into feature_objects_dict_dest.
        This will modify objects in feature_objects_dict_dest. This will *not* insert any new objects from
        feature_objects_dict_source. This method is primarily used to update the existing objects in a database.
        :param feature_objects_dict_dest: The dictionary of feature ids to feature objects to merge into.
        :param feature_objects_dict_source: The dictionary of feature ids to feature objects to merge from.
        :return: Nothing. Modifies objects in feature_objects_dict_dest.
        """
        logger.debug(f'Scanning {len(feature_objects_dict_dest)} destination features for corresponding'
                     f' source features and updating if any matches are found.')
        count_updated_features = 0
        for feature_id in feature_objects_dict_dest:
            feature_object_dest = feature_objects_dict_dest[feature_id]
            if feature_id in feature_objects_dict_source:
                feature_object_source = feature_objects_dict_source[feature_id]
                feature_object_dest.update_sample_ids(feature_object_source)
                count_updated_features += 1

        logger.debug(f'Updated {count_updated_features}/{len(feature_objects_dict_dest)} destination feature objects '
                     f'inplace out of a total of {len(feature_objects_dict_source)} supplied source feature objects.')

    def index_features(self, features_reader: FeaturesReader, feature_scope_name: str) -> None:
        logger.info(f'Reading features from {len(features_reader.samples_list())} samples')
        features_df = features_reader.get_features_table()
        features_df = self._update_scope(features_df, feature_scope_name)
        sample_names = features_reader.samples_set()

        # Create new features and merge with existing features if they exist
        logger.info(f'Aggregating {len(features_df)} features found in {len(sample_names)} samples')
        created_feature_objects = self._create_feature_objects(features_df, sample_names)
        created_feature_objects_dict = {f.id: f for f in created_feature_objects}
        logger.info(f'A total of {len(created_feature_objects_dict)} unique features across {len(sample_names)} '
                    f'samples will be updated. Searching for existing features in the database.')
        db_feature_objects_dict = self.read_index(list(created_feature_objects_dict.keys()))
        logger.info(f'Found {len(db_feature_objects_dict)}/{len(created_feature_objects_dict)} features in database '
                    f'which already exist. These features will be updated to map to the new samples.')
        self.update_feature_objects_inplace(db_feature_objects_dict, created_feature_objects_dict)

        # Subtract out existing feature objects in database
        new_feature_object_ids = set(created_feature_objects_dict.keys()) - set(db_feature_objects_dict.keys())
        new_feature_objects = [created_feature_objects_dict[fid] for fid in new_feature_object_ids]

        # bulk update existing feature objects
        logger.info(
            f'Updating {len(db_feature_objects_dict)}/{len(created_feature_objects_dict)} existing features in database')
        self._connection.get_session().bulk_save_objects(db_feature_objects_dict.values())

        # bulk save new feature objects
        logger.info(
            f'Saving {len(new_feature_objects)}/{len(created_feature_objects_dict)} new features in database')
        self._connection.get_session().bulk_save_objects(new_feature_objects)

        self._connection.get_session().commit()
        logger.info(f'Finished indexing features from {len(sample_names)} samples')

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
