from typing import Dict, Set, List, Any
from pathlib import Path
import pandas as pd
import logging

from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.service import DatabaseConnection
from storage.variant.service.SampleService import SampleService
from storage.variant.service.FeatureService import FeatureService, AUTO_SCOPE
from storage.variant.model import Sample
from storage.variant.SampleSet import SampleSet
from storage.variant.model import MLSTScheme
from storage.variant.model import MLSTAllelesSamples
from storage.variant.model import SampleMLSTAlleles
from storage.variant.io.MLSTFeaturesReader import MLSTFeaturesReader
from storage.variant.service.QueryService import verify_columns_match


logger = logging.getLogger(__name__)


class MLSTService(FeatureService):

    def __init__(self, database_connection: DatabaseConnection, sample_service: SampleService, mlst_dir: Path):
        super().__init__(database_connection=database_connection,
                         features_dir=mlst_dir,
                         sample_service=sample_service)
        self._database = database_connection
        self._sample_service = sample_service

    def find_mlst_scheme(self, name: str) -> MLSTScheme:
        return self._database.get_session().query(MLSTScheme)\
            .filter(MLSTScheme.name == name)\
            .one()

    def exists_mlst_scheme(self, name: str):
        return self._connection.get_session().query(MLSTScheme.id).filter_by(name=name).scalar() is not None

    def get_or_create_mlst_scheme(self, name: str) -> MLSTScheme:
        if self.exists_mlst_scheme(name):
            return self.find_mlst_scheme(name)
        else:
            return MLSTScheme(name=name)

    def find_mlst_schemes(self, scheme_names: Set[str]) -> Dict[str, MLSTScheme]:
        schemes = {}
        for name in scheme_names:
            scheme = self.find_mlst_scheme(name)
            schemes[name] = scheme

        return schemes

    def _create_feature_identifier(self, features_df: pd.DataFrame) -> str:
        return MLSTAllelesSamples.to_sla(
            scheme_name=features_df['Scheme'],
            locus=features_df['Locus'],
            allele=features_df['Allele']
        )

    def aggregate_feature_column(self) -> Dict[str, Any]:
        return {'_SAMPLE_ID': SampleSet}

    def _get_sample_id_series(self, features_df: pd.DataFrame, sample_name_ids: Dict[str, int]) -> pd.Series:
        return features_df.apply(lambda x: sample_name_ids[x['Sample']], axis='columns')

    def _create_feature_object(self, features_df: pd.DataFrame):
        return MLSTAllelesSamples(sla=features_df['_FEATURE_ID'], sample_ids=features_df['_SAMPLE_ID'])

    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        logger.warning('TODO: implement check_samples_have_features')
        return False

    def get_correct_reader(self) -> Any:
        return MLSTFeaturesReader

    def _update_scope(self, features_df: pd.DataFrame, feature_scope_name: str) -> pd.DataFrame:
        if feature_scope_name != AUTO_SCOPE:
            features_df['Scheme'] = feature_scope_name
        return features_df

    def build_sample_feature_object(self, sample: Sample, features_reader: FeaturesReader,
                                    feature_scope_name: str) -> Any:
        self._verify_correct_reader(features_reader=features_reader)
        mlst_reader : MLSTFeaturesReader = features_reader

        mlst_scheme = self.get_or_create_mlst_scheme(feature_scope_name)
        sample_mlst_alleles = SampleMLSTAlleles(sample=sample, scheme=mlst_scheme)

        return sample_mlst_alleles
