from typing import Dict, Set, List
from pathlib import Path
import pandas as pd

from storage.variant.service import DatabaseConnection
from storage.variant.service.SampleService import SampleService
from storage.variant.model import Sample
from storage.variant.SampleSet import SampleSet
from storage.variant.model import MLSTScheme
from storage.variant.model import MLSTAllelesSamples
from storage.variant.service.QueryService import verify_columns_match


class MLSTService:

    def __init__(self, database_connection: DatabaseConnection, sample_service: SampleService):
        self._database = database_connection
        self._sample_service = sample_service

    def find_mlst_scheme(self, name: str) -> MLSTScheme:
        return self._database.get_session().query(MLSTScheme)\
            .filter(MLSTScheme.name == name)\
            .one()

    def find_mlst_schemes(self, scheme_names: Set[str]) -> Dict[str, MLSTScheme]:
        schemes = {}
        for name in scheme_names:
            scheme = self.find_mlst_scheme(name)
            schemes[name] = scheme

        return schemes

    def _create_mlst_alleles(self, mlst_df: pd.DataFrame) -> List[MLSTAllelesSamples]:
        samples_names = set(mlst_df['Sample'].tolist())
        sample_name_ids = self._sample_service.find_sample_name_ids(samples_names)

        mlst_df['SLA'] = mlst_df.apply(lambda x: MLSTAllelesSamples.to_sla(
            scheme=x['Scheme'],
            locu=x['Locus'],
            allele=x['Allele'],
        ), axis='columns')
        mlst_df['SAMPLE_ID'] = mlst_df.apply(lambda x: sample_name_ids[x['Sample']], axis='columns')

        index_df = mlst_df.groupby('SLA').agg({'TYPE': 'first', 'SAMPLE_ID': SampleSet}).reset_index()

        return index_df.apply(lambda x: MLSTAllelesSamples(
            sla=x['SLA'], sample_ids=x['SAMPLE_ID']),
            axis='columns').tolist()

    # def insert_mlst_results(self, mlst_reader) -> None:
    #     """
    #     Inserts MLST results into the database.
    #     :param mlst_df: The dataframe containing the MLST results.
    #     :return: None
    #     """
    #     if mlst_df is None:
    #         raise Exception(f'mlst_df is None')
    #     verify_columns_match({'Sample', 'Scheme', 'Locus', 'Allele'}, mlst_df)
    #
    #     schemes = self.find_mlst_schemes(set(mlst_df['Scheme'].tolist()))
    #     mlst_df['Scheme_obj'] = mlst_df['Scheme'].apply(lambda x: schemes[x])
    #
    #     mlst_df = mlst_reader.get_variants_table()
    #     self._connection.get_session().bulk_save_objects(self._create_mlst_alleles(mlst_df))
    #     self._connection.get_session().commit()
