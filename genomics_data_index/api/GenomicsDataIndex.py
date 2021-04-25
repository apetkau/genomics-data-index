from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import List

import pandas as pd

from genomics_data_index.api.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from genomics_data_index.api.query.impl.SamplesQueryIndex import SamplesQueryIndex
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.Project import Project
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater

logger = logging.getLogger(__name__)


class GenomicsDataIndex:
    QUERY_UNIVERSE = ['all', 'mutations', 'dataframe']
    MUTATION_ID_TYPES = ['spdi_ref', 'spdi']

    def __init__(self, connection: DataIndexConnection):
        self._connection = connection
        self._samples_count = self.count_samples()

    @property
    def connection(self) -> DataIndexConnection:
        return self._connection

    def count_samples(self) -> int:
        return self._connection.sample_service.count_samples()

    def sample_names(self) -> List[str]:
        return [s.name for s in self._connection.sample_service.get_samples()]

    def count_references(self) -> int:
        return self._connection.reference_service.count_reference_genomes()

    def reference_names(self) -> List[str]:
        return [r.name for r in self._connection.reference_service.get_reference_genomes()]

    def count_mutations(self, reference_genome: str, include_unknown: bool = False) -> int:
        return self._connection.variation_service.count_on_reference(reference_genome,
                                                                     include_unknown=include_unknown)

    def mutations_summary(self, reference_genome: str, id_type: str = 'spdi_ref',
                          include_unknown: bool = False) -> pd.DataFrame:
        rs = self._connection.reference_service
        if id_type not in self.MUTATION_ID_TYPES:
            raise Exception(f'id_type={id_type} must be one of {self.MUTATION_ID_TYPES}')

        vs = self._connection.variation_service
        mutation_counts = vs.mutation_counts_on_reference(reference_genome,
                                                          include_unknown=include_unknown)
        if id_type == 'spdi_ref':
            translated_ids = rs.translate_spdi(mutation_counts.keys(), to=id_type)
            mutation_counts = {translated_ids[m]: mutation_counts[m] for m in mutation_counts}
            convert_deletion = False
        else:
            convert_deletion = True

        data = []
        for mutation in mutation_counts:
            seq, pos, deletion, insertion = NucleotideMutationTranslater.from_spdi(mutation,
                                                                                   convert_deletion=convert_deletion)
            data.append([mutation, seq, pos, deletion, insertion, mutation_counts[mutation]])

        return pd.DataFrame(data,
                            columns=['Mutation', 'Sequence', 'Position',
                                     'Deletion', 'Insertion', 'Count']).set_index('Mutation')

    def db_size(self, unit: str = 'B') -> pd.DataFrame:
        return self._connection.db_size(unit)

    def samples_query(self, universe: str = 'all', **kwargs) -> SamplesQuery:
        if universe == 'all':
            return self._query_all_samples(self._connection)
        elif universe == 'mutations':
            return self._query_reference(connection=self._connection, **kwargs)
        elif universe == 'dataframe':
            return self._query_data_frame(connection=self._connection, **kwargs)
        else:
            raise Exception(f'Invalid universe=[{universe}]. Must be one of {self.QUERY_UNIVERSE}')

    def _query_all_samples(self, connection: DataIndexConnection):
        all_samples = connection.sample_service.get_all_sample_ids()
        return SamplesQueryIndex(connection=connection, sample_set=all_samples, universe_set=all_samples)

    def _query_reference(self, connection: DataIndexConnection, reference_name: str):
        reference_samples = connection.sample_service.get_samples_set_associated_with_reference(reference_name)
        reference_genome = connection.reference_service.find_reference_genome(reference_name)

        sample_query = SamplesQueryIndex(connection=connection, sample_set=reference_samples,
                                         universe_set=reference_samples)

        if reference_genome.has_tree():
            sample_query = TreeSamplesQuery(connection=connection, wrapped_query=sample_query,
                                            tree=reference_genome.tree,
                                            alignment_length=reference_genome.tree_alignment_length)

        return sample_query

    def _query_data_frame(self, connection: DataIndexConnection,
                          data_frame: pd.DataFrame = None,
                          sample_ids_column=None,
                          sample_names_column=None
                          ):
        if data_frame is None:
            raise Exception('data_frame must be set when querying with universe="dataframe"')
        if sample_ids_column is None and sample_names_column is None:
            raise Exception(
                'If querying with universe="dataframe", then one of sample_names_column or sample_ids_column '
                'must be set')
        elif sample_ids_column is not None:
            sample_query = self._query_all_samples(connection=connection)
            return DataFrameSamplesQuery.create_with_sample_ids_column(sample_ids_column,
                                                                       data_frame=data_frame,
                                                                       wrapped_query=sample_query,
                                                                       connection=connection)
        else:
            sample_query = self._query_all_samples(connection=connection)
            return DataFrameSamplesQuery.create_with_sample_names_column(sample_names_column,
                                                                         data_frame=data_frame,
                                                                         wrapped_query=sample_query,
                                                                         connection=connection)

    @classmethod
    def connect(cls, project_dir: Path = None, project: Project = None) -> GenomicsDataIndex:
        if project_dir is None and project is None:
            project_dir = os.getcwd()
            logger.warning(f'No project_dir or project specified. Assuming project is current dir [{project_dir}]')
            data_store_project = Project(project_dir)
        elif project_dir is not None and project is None:
            if not project_dir.exists():
                raise Exception(f'project_dir=[{project_dir}] does not exist')
            data_store_project = Project(root_dir=project_dir)
        else:
            data_store_project = project

        database_connection = data_store_project.create_connection()
        return GenomicsDataIndex(connection=database_connection)

    def __str__(self) -> str:
        samples_count = self._samples_count
        return f'<{self.__class__.__name__}(samples={samples_count})>'

    def __repr__(self) -> str:
        return str(self)
