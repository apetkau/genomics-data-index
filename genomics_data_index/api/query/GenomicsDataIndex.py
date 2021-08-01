from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import List, Union

import pandas as pd
from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from genomics_data_index.api.query.impl.SamplesQueryIndex import SamplesQueryIndex
from genomics_data_index.api.query.impl.TreeSamplesQueryFactory import TreeSamplesQueryFactory
from genomics_data_index.configuration.Project import Project
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.VariationService import VariationService

logger = logging.getLogger(__name__)


class GenomicsDataIndex:
    """
    The root class for the Python API. Used to connect to a genomics index and build query objects.
    """

    QUERY_UNIVERSE = ['all', 'mutations', 'mutations_experimental', 'dataframe']
    MUTATION_ID_TYPES = ['spdi_ref', 'spdi']
    FEAUTRE_KINDS = ['mutations', 'mlst']

    def __init__(self, connection: DataIndexConnection):
        """
        Builds a new GenomicsDataIndex.
        :param connection: The database connection to use.

        :return: A new GenomicsDataIndex.
        """
        self._connection = connection
        self._samples_count = self.count_samples()

    @property
    def connection(self) -> DataIndexConnection:
        """
        Gets the database connection object.

        :return: The database connection object.
        """
        return self._connection

    def count_samples(self) -> int:
        """
        Counts the samples stored in this genomics index.

        :return: The count of samples stored in this genomics index.
        """
        return self._connection.sample_service.count_samples()

    def sample_names(self) -> List[str]:
        """
        Gets a list of all sample names stored in this index.

        :return: A list of all sample names stored in this index.
        """
        return [s.name for s in self._connection.sample_service.get_samples()]

    def count_references(self) -> int:
        """
        Counts the number of reference genomes stored in this genomics index.

        :return: The count of reference genomes stored in this genomics index.
        """
        return self._connection.reference_service.count_reference_genomes()

    def reference_names(self) -> List[str]:
        """
        Gets a list of the names of the reference genomes used in this index.
        :return: A list of all reference genome names.
        """
        return [r.name for r in self._connection.reference_service.get_reference_genomes()]

    def mlst_schemes(self) -> List[str]:
        """
        Gets a list of the MLST schemes used in this index.
        :return: A list of all MLST scheme names.
        """
        return [r.name for r in self._connection.mlst_service.get_mlst_schemes()]

    def reference_tree(self, reference_name: str) -> Tree:
        """
        Gets a tree associated with a reference genome.
        :param reference_name: The name of the reference genome to get a tree for.
        :return: The tree (if any) associated with the reference genome.
        """
        reference_genome = self._connection.reference_service.find_reference_genome(reference_name)
        if reference_genome.has_tree():
            return reference_genome.tree
        else:
            raise EntityExistsError(f'reference_genome=[{reference_genome}] does not have a tree')

    def count_mutations(self, reference_genome: str, include_unknown: bool = False) -> int:
        """
        Counts all mutations indexed relative to the given reference genome.

        :param reference_genome: The reference genome.
        :param include_unknown: Whether or not unknown mutations should be included.

        :return: A count of all mutations relative to the given reference genome.
        """
        return self._connection.variation_service.count_on_reference(reference_genome,
                                                                     include_unknown=include_unknown)

    def mutations_summary(self, reference_name: str, id_type: str = 'spdi_ref', include_present: bool = True,
                          include_unknown: bool = False, ignore_annotations: bool = False) -> pd.DataFrame:
        """
        Summarizes all mutations stored in this index relative to a string for the passed scope.
        Shorthand for features_summary(kind='mutations', ...)

        :param reference_name: The reference genome name.
        :param id_type: The type of identifier to use.
        :param include_present: Whether or not mutation features present in this index (i.e., not unknown/missing)
                                should be included.
        :param include_unknown: Whether or not unknown mutations should be included.
        :param ignore_annotations: Whether or not mutation annotations should be ignored.
        :return: A summary of all mutations in this index as a DataFrame.
        """
        return self.features_summary(kind='mutations', scope=reference_name,
                                     include_present=include_present, include_unknown=include_unknown,
                                     id_type=id_type, ignore_annotations=ignore_annotations)

    def mlst_summary(self, scheme_name: str, locus: str = None, include_present: bool = True,
                     include_unknown: bool = False) -> pd.DataFrame:
        """
        Summarizes all MLST alleles stored in this index relative to the passed scheme name.
        Shorthand for features_summary(kind='mlst', ...)

        :param scheme_name: The MLST scheme to summarize.
        :param locus: The locus id to restrict summaries to. Defaults to all locus IDs.
        :param include_present: Whether or not MLST features present in this index (i.e., not unknown/missing)
                                should be included.
        :param include_unknown: Whether or not unknown MLST alleles should be included.
        :return: A summary of all MLST alleles in this index as a DataFrame.
        """
        return self.features_summary(kind='mlst', scope=scheme_name, include_present=include_present,
                                     include_unknown=include_unknown, locus=locus)

    def features_summary(self, kind: str = 'mutations', scope: str = None,
                         include_present: bool = True, include_unknown: bool = False, **kwargs) -> pd.DataFrame:
        """
        Summarizes all features stored in this index relative to a string for the passed scope.

        :param kind: The kind of feature (e.g., 'mutations' or 'mlst').
        :param scope: The scope (e.g., reference genome or MLST scheme).
        :param include_present: Whether or not features present in this index (i.e., not unknown/missing)
                                should be included.
        :param include_unknown: Whether or not unknown/missing features should be included.
        :return: A summary of all features in this index as a DataFrame.
        """
        if kind == 'mutations' or kind == 'mutation':
            return self._mutations_summary_internal(reference_name=scope, include_present=include_present,
                                                    include_unknown=include_unknown, **kwargs)
        elif kind == 'mlst':
            return self._mlst_summary_internal(scheme_name=scope, include_present=include_present,
                                               include_unknown=include_unknown, **kwargs)
        else:
            raise Exception(f'Unknown value for kind=[{kind}]. Must be one of {self.FEAUTRE_KINDS}.')

    def _mlst_summary_internal(self, scheme_name: str, locus: str = None, include_present: bool = True,
                               include_unknown: bool = False) -> pd.DataFrame:
        """
        Summarizes all MLST alleles stored in this index relative to the passed scheme name.

        :param scheme_name: The MLST scheme to summarize.
        :param locus: The locus id to restrict summaries to. Defaults to all locus IDs.
        :param include_present: Whether or not MLST features present in this index (i.e., not unknown/missing)
                                should be included.
        :param include_unknown: Whether or not unknown MLST alleles should be included.
        :return: A summary of all MLST alleles in this index as a DataFrame.
        """
        mlst_service: MLSTService = self._connection.mlst_service
        mlst_features = mlst_service.get_features(scheme_name, locus=locus,
                                                  include_present=include_present,
                                                  include_unknown=include_unknown)
        total_samples = self._connection.sample_service.count_samples_associated_with_mlst_scheme(scheme_name)

        data = []
        for mlst_feature_id in mlst_features:
            mlst_feature = mlst_features[mlst_feature_id]
            count = len(mlst_feature.sample_ids)
            data.append([mlst_feature.query_id, mlst_feature.scheme, mlst_feature.locus, mlst_feature.allele, count,
                         total_samples])

        features_df = pd.DataFrame(data,
                                   columns=['MLST Feature', 'Scheme', 'Locus',
                                            'Allele', 'Count', 'Total']).set_index('MLST Feature')

        features_df['Percent'] = 100 * (features_df['Count'] / features_df['Total'])

        return features_df

    def _mutations_summary_internal(self, reference_name: str, id_type: str = 'spdi_ref',
                                    include_present: bool = True,
                                    include_unknown: bool = False, ignore_annotations: bool = False) -> pd.DataFrame:
        """
        Summarizes all mutations stored in this index relative to the passed reference genome.

        :param reference_name: The reference genome.
        :param id_type: The type of identifier to use.
        :param include_present: Whether or not mutation features present in this index (i.e., not unknown/missing)
                                should be included.
        :param include_unknown: Whether or not unknown mutations should be included.
        :param ignore_annotations: Whether or not mutation annotations should be ignored.

        :return: A summary of all mutations in this index as a DataFrame.
        """
        rs = self._connection.reference_service
        if id_type not in self.MUTATION_ID_TYPES:
            raise Exception(f'id_type={id_type} must be one of {self.MUTATION_ID_TYPES}')

        vs: VariationService = self._connection.variation_service
        mutations = vs.get_variants_on_reference(reference_name, include_present=include_present,
                                                 include_unknown=include_unknown)

        if id_type == 'spdi_ref':
            translated_ids = rs.translate_spdi(mutations.keys(), to=id_type)
            mutations = {translated_ids[m]: mutations[m] for m in mutations}
            convert_deletion = False
        else:
            convert_deletion = True

        total_samples = self._connection.sample_service.count_samples_associated_with_reference(reference_name)

        data = []
        for mutation in mutations:
            seq, pos, deletion, insertion = NucleotideMutationTranslater.from_spdi(mutation,
                                                                                   convert_deletion=convert_deletion)
            count = len(mutations[mutation].sample_ids)
            data.append([mutation, seq, pos, deletion, insertion, count, total_samples])

        features_df = pd.DataFrame(data,
                                   columns=['Mutation', 'Sequence', 'Position',
                                            'Deletion', 'Insertion', 'Count', 'Total']).set_index('Mutation')

        features_df['Percent'] = 100 * (features_df['Count'] / features_df['Total'])

        if not ignore_annotations:
            features_df = vs.append_mutation_annotations(features_df)

        return features_df

    def db_size(self, unit: str = 'B') -> pd.DataFrame:
        """
        Gets the size of this genomics index broken apart into separate sections.

        :param unit: The units to use for reporting the sizes.

        :return: A DataFrame summarizing the size of this genomics index.
        """
        return self._connection.db_size(unit)

    def samples_query(self, universe: str = 'all', **kwargs) -> SamplesQuery:
        """
        Constructs a new SamplesQuery with respect to this index.

        :param universe: The universe of samples under consideration. Options include ['all', 'mutations', 'dataframe'].
        :param kwargs: Additional arguments depending on the universe of samples.
        :return: A SamplesQuery object which can be used to further refine the query.
        """
        if universe == 'all':
            return self._query_all_samples(self._connection)
        elif universe == 'mutations':
            return self._query_reference(kind=universe, connection=self._connection, **kwargs)
        elif universe == 'mutations_experimental':
            return self._query_reference(kind=universe, connection=self._connection, **kwargs)
        elif universe == 'dataframe':
            return self._query_data_frame(connection=self._connection, **kwargs)
        else:
            raise Exception(f'Invalid universe=[{universe}]. Must be one of {self.QUERY_UNIVERSE}')

    def _query_all_samples(self, connection: DataIndexConnection):
        all_samples = connection.sample_service.get_all_sample_ids()
        return SamplesQueryIndex(connection=connection, sample_set=all_samples, universe_set=all_samples,
                                 unknown_set=SampleSet.create_empty())

    def _query_reference(self, kind: str, connection: DataIndexConnection, reference_name: str):
        reference_samples = connection.sample_service.get_samples_set_associated_with_reference(reference_name)
        reference_genome = connection.reference_service.find_reference_genome(reference_name)

        sample_query = SamplesQueryIndex(connection=connection, sample_set=reference_samples,
                                         unknown_set=SampleSet.create_empty(),
                                         universe_set=reference_samples)

        if reference_genome.has_tree():
            sample_query = TreeSamplesQueryFactory.instance().create_from_reference_genome(kind=kind,
                                                                                           reference_genome=reference_genome,
                                                                                           connection=connection,
                                                                                           wrapped_query=sample_query,
                                                                                           include_reference=True)

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
    def connect(cls, project_dir: Union[Path, str] = None, project: Project = None) -> GenomicsDataIndex:
        """
        Connects to a new genomics index. One of either project_dir or project must be set (by default assumes
        project_dir is the current directory).

        :param project_dir: The project/index directory. Defaults to the current directory.
        :param project: The :py:class:`genomics_data_index.configuration.Project` to connect to.
        :return: A new GenomicsDataIndex connected to the given project.
        """
        if project_dir is None and project is None:
            project_dir = os.getcwd()
            logger.warning(f'No project_dir or project specified. Assuming project is current dir [{project_dir}]')
            data_store_project = Project(project_dir)
        elif project_dir is not None and project is None:
            if isinstance(project_dir, str):
                project_dir = Path(project_dir)

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
