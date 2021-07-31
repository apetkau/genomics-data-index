import logging
from pathlib import Path
from typing import Dict, Set, Any, Tuple, cast, List

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.io.FeaturesReader import FeaturesReader
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTSampleData import MLSTSampleData
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.model import MLST_UNKNOWN_ALLELE
from genomics_data_index.storage.model.db import MLSTScheme, SampleMLSTAlleles, MLSTAllelesSamples, Sample, \
    FeatureSamples
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.FeatureService import FeatureService, AUTO_SCOPE
from genomics_data_index.storage.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class MLSTService(FeatureService):

    def __init__(self, database_connection: DatabaseConnection, sample_service: SampleService, mlst_dir: Path):
        super().__init__(database_connection=database_connection,
                         features_dir=mlst_dir,
                         sample_service=sample_service)
        self._database = database_connection
        self._sample_service = sample_service

    def find_mlst_scheme(self, name: str) -> MLSTScheme:
        return self._database.get_session().query(MLSTScheme) \
            .filter(MLSTScheme.name == name) \
            .one()

    def get_mlst_schemes(self) -> List[MLSTScheme]:
        return self._database.get_session().query(MLSTScheme).all()

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

    def get_features(self, scheme: str = None, locus: str = None, include_present: bool = True,
                     include_unknown: bool = False) -> Dict[str, MLSTAllelesSamples]:
        query = self._database.get_session().query(MLSTAllelesSamples)
        if scheme is not None:
            query = query.filter(MLSTAllelesSamples.scheme == scheme)

        if not include_present:
            if include_unknown:
                query = query.filter(MLSTAllelesSamples.allele == MLST_UNKNOWN_ALLELE)
            else:
                # Here, I'm not including unknown or present, so don't query database
                return dict()
        elif not include_unknown:
            query = query.filter(MLSTAllelesSamples.allele != MLST_UNKNOWN_ALLELE)

        if locus is not None:
            query = query.filter(MLSTAllelesSamples.locus == locus)

        allele_samples = query.all()

        return {f.id: f for f in allele_samples}

    def get_all_alleles(self, scheme: str, locus: str) -> Set[str]:
        return {a for a, in self._database.get_session().query(MLSTAllelesSamples.allele) \
            .filter(MLSTAllelesSamples.scheme == scheme, MLSTAllelesSamples.locus == locus) \
            .all()}

    def get_all_loci_alleles(self, scheme: str) -> Set[Tuple[str, str]]:
        """
        Gets all (loci, allele) pairs from the database given a scheme.
        :param scheme: The scheme.
        :return: Gets a list of tuples of the form (loci, allele).
        """
        return {a for a in self._database.get_session().query(MLSTAllelesSamples.locus, MLSTAllelesSamples.allele) \
            .filter(MLSTAllelesSamples.scheme == scheme) \
            .all()}

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

    def _create_feature_object(self, features_df: pd.DataFrame) -> FeatureSamples:
        return MLSTAllelesSamples(sla=features_df['_FEATURE_ID'], sample_ids=features_df['_SAMPLE_ID'])

    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        samples_with_mlst = {sample.name for sample in
                             self._sample_service.get_samples_with_mlst_alleles(feature_scope_name)}
        return len(samples_with_mlst.intersection(sample_names)) != 0

    def get_correct_data_package(self) -> Any:
        return MLSTSampleDataPackage

    def get_correct_sample_data(self) -> Any:
        return MLSTSampleData

    def _update_scope(self, features_df: pd.DataFrame, feature_scope_name: str) -> pd.DataFrame:
        if feature_scope_name != AUTO_SCOPE:
            features_df['Scheme'] = feature_scope_name
        return features_df

    def build_sample_feature_object(self, sample: Sample, sample_data: SampleData, feature_scope_name: str) -> Any:
        self._verify_correct_sample_data(sample_data=sample_data)
        mlst_sample_data = cast(MLSTSampleData, sample_data)

        if feature_scope_name == AUTO_SCOPE:
            scheme_name = mlst_sample_data.get_scheme()
        else:
            scheme_name = feature_scope_name

        mlst_scheme = self.get_or_create_mlst_scheme(scheme_name)
        sample_mlst_alleles = SampleMLSTAlleles(sample=sample, scheme=mlst_scheme)

        return sample_mlst_alleles

    def _create_persisted_features_reader(self, sample_data_dict: Dict[str, SampleData],
                                          data_package: SampleDataPackage) -> FeaturesReader:
        self._verify_correct_data_package(data_package=data_package)
        mlst_data_package = cast(MLSTSampleDataPackage, data_package)
        return mlst_data_package.get_features_reader()

    def read_index(self, feature_ids: List[str]) -> Dict[str, FeatureSamples]:
        feature_samples = self._connection.get_session().query(MLSTAllelesSamples) \
            .filter(MLSTAllelesSamples._sla.in_(feature_ids)) \
            .all()

        return {f.id: f for f in feature_samples}
