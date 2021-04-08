from typing import List, Set, Any
import abc
import logging
from pathlib import Path

import pandas as pd

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io import NucleotideFeaturesReader
from storage.variant.SampleSet import SampleSet
from storage.variant.model import Sample, SampleNucleotideVariation, NucleotideVariantsSamples
from storage.variant.service import DatabaseConnection
from storage.variant.service import EntityExistsError
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.service.SampleService import SampleService
from storage.variant.io.VariationFile import VariationFile

logger = logging.getLogger(__name__)


class FeatureService(abc.ABC):

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService):
        self._connection = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir

    @abc.abstractmethod
    def _create_feature_objects(self, features_df: pd.DataFrame) -> List[Any]:
        pass

    @abc.abstractmethod
    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        """
        Checks if any of the passed sample names already have features.
        :param sample_names: The set of sample names.
        :param feature_scope_name: A name defining the feature scope (i.e., reference genome name).
        :return: True if any of the samples have features, false otherwise.
        """
        pass

    @abc.abstractmethod
    def _verify_correct_reader(self, features_reader: FeaturesReader) -> None:
        pass

    def _get_or_create_sample(self, sample_name: str) -> Sample:
        if self._sample_service.exists(sample_name):
            return self._sample_service.get_sample(sample_name)
        else:
            return Sample(name=sample_name)

    def insert(self, features_reader: FeaturesReader, feature_scope_name: str = None) -> None:
        # reference = self._reference_service.find_reference_genome(reference_name)
        # sample_variant_files = variants_reader.sample_feature_files()
        # genomic_masked_regions = variants_reader.get_genomic_masked_regions()
        self._verify_correct_reader(features_reader=features_reader)

        sample_names = features_reader.samples_set()

        if self.check_samples_have_features(sample_names, feature_scope_name):
            raise EntityExistsError(f'Passed samples already have features for feature scope [{feature_scope_name}], '
                                    f'will not insert any new features')

        # TODO: keep track of saved feature files to index these ones
        saved_variation_files = {}
        saved_masked_regions = {}

        for sample_name in sample_names:
            sample = self._get_or_create_sample(sample_name)
            sample_feature_object = self.build_sample_feature_object(sample=sample,
                                                                     features_reader=features_reader,
                                                                     feature_scope_name=feature_scope_name)

            self._connection.get_session().add(sample_feature_object)
        self._connection.get_session().commit()

        self.index_features(features_reader=features_reader)

    def index_features(self, features_reader: FeaturesReader):
        features_df = features_reader.get_features_table()
        self._connection.get_session().bulk_save_objects(self._create_feature_objects(features_df))
        self._connection.get_session().commit()

    @abc.abstractmethod
    def build_sample_feature_object(self, sample: Sample,
                                    features_reader: FeaturesReader, feature_scope_name: str) -> Any:
        pass
