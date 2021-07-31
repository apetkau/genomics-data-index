import logging
from pathlib import Path
from typing import List, Set, Any, Dict, cast, Union

import pandas as pd
import sqlalchemy

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.io.FeaturesReader import FeaturesReader
from genomics_data_index.storage.io.SampleData import SampleData
from genomics_data_index.storage.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.VariationFile import VariationFile
from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN_TYPE
from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import NucleotideVariantsSamples, SampleNucleotideVariation, Sample, \
    FeatureSamples
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.FeatureService import FeatureService
from genomics_data_index.storage.service.ReferenceService import ReferenceService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.util import TRACE_LEVEL
from genomics_data_index.storage.util.ListSliceIter import ListSliceIter
from genomics_data_index.storage.util.SamplesProgressLogger import SamplesProgressLogger

logger = logging.getLogger(__name__)


class VariationService(FeatureService):
    MUTATION_TYPES = ['snp', 'indel', 'all', 'other']
    MUTATION_ID_TYPES = ['spdi_ref', 'spdi']

    def __init__(self, database_connection: DatabaseConnection, variation_dir: Path,
                 reference_service: ReferenceService, sample_service: SampleService,
                 index_unknown_missing: bool, sql_select_limit: int):
        super().__init__(database_connection=database_connection,
                         features_dir=variation_dir,
                         sample_service=sample_service)
        self._reference_service = reference_service
        self._index_unknown_missing = index_unknown_missing
        self._sql_select_limit = sql_select_limit

    def _reference_sequence_names(self, reference_name: str) -> List[str]:
        return list(self._reference_service.get_reference_sequences(reference_name).keys())

    def _query_include_unknown(self, query, include_present: bool = True, include_unknown: bool = False):
        if not include_present:
            if include_unknown:
                return query.filter(NucleotideVariantsSamples.var_type == NUCLEOTIDE_UNKNOWN_TYPE)
            else:
                # Should return a query that gives no results
                return query.filter(sqlalchemy.sql.false())
        elif not include_unknown:
            return query.filter(NucleotideVariantsSamples.var_type != NUCLEOTIDE_UNKNOWN_TYPE)
        else:
            return query

    def count_on_reference(self, reference_name: str, include_unknown: bool = False) -> int:
        reference_sequence_names = self._reference_sequence_names(reference_name)

        query = self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence.in_(reference_sequence_names))

        return self._query_include_unknown(query, include_unknown=include_unknown).count()

    def mutation_counts_on_reference(self, reference_name: str, include_unknown: bool = False) -> Dict[str, int]:
        reference_sequence_names = self._reference_sequence_names(reference_name)

        query = self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence.in_(reference_sequence_names))
        mutations = self._query_include_unknown(query, include_unknown=include_unknown).all()

        return {m.spdi: len(m.sample_ids) for m in mutations}

    def get_features(self, include_present: bool = True,
                     include_unknown: bool = False,
                     id_type: str = 'spdi') -> Dict[str, NucleotideVariantsSamples]:
        query = self._connection.get_session().query(NucleotideVariantsSamples)

        if not include_present:
            if include_unknown:
                query = query.filter(NucleotideVariantsSamples.var_type == NUCLEOTIDE_UNKNOWN_TYPE)
            else:
                # Here, I'm not including unknown or present, so don't query database
                return dict()
        elif not include_unknown:
            query = query.filter(NucleotideVariantsSamples.var_type != NUCLEOTIDE_UNKNOWN_TYPE)

        features = query.all()
        if id_type == 'spdi_ref':
            feature_ids_db = {f.id: f for f in features}
            feature_ids_translated = self._reference_service.translate_spdi(feature_ids_db.keys(), to=id_type)
            return {feature_ids_translated[fid]: feature_ids_db[fid] for fid in feature_ids_db}
        elif id_type == 'spdi':
            return {f.id: f for f in features}
        else:
            raise Exception(f'Unknown value for id_type={id_type}. Must be one of {self.MUTATION_ID_TYPES}')

    def get_variants_on_reference(self, reference_name: str, include_present: bool = True,
                                  include_unknown: bool = False) -> Dict[str, NucleotideVariantsSamples]:
        reference_sequence_names = self._reference_sequence_names(reference_name)

        query = self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence.in_(reference_sequence_names))

        mutations = self._query_include_unknown(query, include_present=include_present,
                                                include_unknown=include_unknown).all()

        return {m.spdi: m for m in mutations}

    def count_mutations_in_sample_ids_dataframe(self, sample_ids: Union[SampleSet, List[int]],
                                                ncores: int = 1,
                                                batch_size: int = 50,
                                                mutation_type: str = 'all',
                                                include_unknown: bool = False) -> pd.DataFrame:
        if include_unknown:
            raise Exception(f'support for include_unknown is not implemented')

        if mutation_type == 'all':
            include_expression = None
        elif mutation_type == 'snp':
            include_expression = 'TYPE="SNP"'
        elif mutation_type == 'indel':
            include_expression = 'TYPE="INDEL"'
        elif mutation_type == 'other':
            include_expression = 'TYPE="OTHER"'
        else:
            raise Exception(f'Unsupported option mutation_type=[{mutation_type}]. Must be one of {self.MUTATION_TYPES}')

        if isinstance(sample_ids, SampleSet):
            sample_ids = list(sample_ids)

        sample_nucleotide_variation = self._connection.get_session().query(SampleNucleotideVariation) \
            .filter(SampleNucleotideVariation.sample_id.in_(sample_ids)) \
            .all()

        nucleotide_variants_files = [snv.nucleotide_variants_file for snv in sample_nucleotide_variation]
        mutation_df = VariationFile.union_all_files(nucleotide_variants_files,
                                                    include_expression=include_expression,
                                                    batch_size=batch_size,
                                                    ncores=ncores)

        mutation_df = mutation_df.rename(columns={
            'ID': 'Mutation',
            'CHROM': 'Sequence',
            'POS': 'Position',
            'REF': 'Deletion',
            'ALT': 'Insertion',
            'COUNT': 'Count',
        })

        return mutation_df.set_index('Mutation')

    def append_mutation_annotations(self, features_df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds annotations to the mutations stored within the passed dataframe.
        :param features_df: The dataframe to add annotations.
                            Assumes the index of the dataframe contains the SPDI identifier.
        :return: A new dataframe with mutation annotations.
        """
        if features_df.index.name != 'Mutation':
            raise Exception(f'Index name not equal to "Mutation" for features_df={features_df}')

        mutation_ids = features_df.index.tolist()
        query_features = [QueryFeatureMutationSPDI(i) for i in mutation_ids]
        id_to_nucleotide_variants_samples: Dict[str, NucleotideVariantsSamples] = \
            self._sample_service.get_variants_samples_by_variation_features(query_features)

        annotation_data = []
        for mutation_id in id_to_nucleotide_variants_samples:
            variants_samples = id_to_nucleotide_variants_samples[mutation_id]

            if variants_samples.id_hgvs_c is not None:
                id_hgvs_gn_c = QueryFeatureHGVSGN.create(sequence_name=variants_samples.sequence,
                                                         gene=variants_samples.annotation_gene_name,
                                                         variant=variants_samples.annotation_hgvs_c).id
            else:
                id_hgvs_gn_c = pd.NA

            if variants_samples.id_hgvs_p is not None:
                id_hgvs_gn_p = QueryFeatureHGVSGN.create(sequence_name=variants_samples.sequence,
                                                         gene=variants_samples.annotation_gene_name,
                                                         variant=variants_samples.annotation_hgvs_p).id
            else:
                id_hgvs_gn_p = pd.NA

            annotation_data.append([mutation_id,
                                    variants_samples.annotation if variants_samples.annotation is not None else pd.NA,
                                    variants_samples.annotation_impact if variants_samples.annotation_impact is not None else pd.NA,
                                    variants_samples.annotation_gene_name if variants_samples.annotation_gene_name is not None else pd.NA,
                                    variants_samples.annotation_gene_id if variants_samples.annotation_gene_id is not None else pd.NA,
                                    variants_samples.annotation_feature_type if variants_samples.annotation_feature_type is not None else pd.NA,
                                    variants_samples.annotation_transcript_biotype if variants_samples.annotation_transcript_biotype is not None else pd.NA,
                                    variants_samples.annotation_hgvs_c if variants_samples.annotation_hgvs_c is not None else pd.NA,
                                    variants_samples.annotation_hgvs_p if variants_samples.annotation_hgvs_p is not None else pd.NA,
                                    variants_samples.id_hgvs_c if variants_samples.id_hgvs_c is not None else pd.NA,
                                    variants_samples.id_hgvs_p if variants_samples.id_hgvs_p is not None else pd.NA,
                                    id_hgvs_gn_c,
                                    id_hgvs_gn_p,
                                    ])

        annotation_df = pd.DataFrame(data=annotation_data,
                                     columns=['Mutation',
                                              'Annotation',
                                              'Annotation_Impact',
                                              'Gene_Name',
                                              'Gene_ID',
                                              'Feature_Type',
                                              'Transcript_BioType',
                                              'HGVS.c',
                                              'HGVS.p',
                                              'ID_HGVS.c',
                                              'ID_HGVS.p',
                                              'ID_HGVS_GN.c',
                                              'ID_HGVS_GN.p',
                                              ]).set_index('Mutation')

        return features_df.merge(annotation_df, how='left', left_index=True, right_index=True)

    def get_variants_ordered(self, sequence_name: str, type: str = 'SNP') -> List[NucleotideVariantsSamples]:
        return self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence == sequence_name) \
            .filter(NucleotideVariantsSamples.var_type == type) \
            .order_by(NucleotideVariantsSamples.position) \
            .all()

    def get_sample_nucleotide_variation(self, sample_names: List[str]) -> List[SampleNucleotideVariation]:
        return self._connection.get_session().query(SampleNucleotideVariation) \
            .join(SampleNucleotideVariation.sample) \
            .filter(Sample.name.in_(sample_names)) \
            .all()

    def _create_feature_identifier(self, features_df: pd.DataFrame) -> str:
        return NucleotideVariantsSamples.to_spdi(
            sequence_name=features_df['CHROM'],
            position=features_df['POS'],
            ref=features_df['REF'],
            alt=features_df['ALT']
        )

    def _get_sample_id_series(self, features_df: pd.DataFrame, sample_name_ids: Dict[str, int]) -> pd.Series:
        return features_df.apply(lambda x: sample_name_ids[x['SAMPLE']], axis='columns')

    def _create_feature_object(self, features_df: pd.DataFrame) -> FeatureSamples:
        return NucleotideVariantsSamples(spdi=features_df['_FEATURE_ID'],
                                         var_type=features_df['TYPE'],
                                         sample_ids=features_df['_SAMPLE_ID'],
                                         annotation=features_df['ANN.Annotation'],
                                         annotation_impact=features_df['ANN.Annotation_Impact'],
                                         annotation_gene_name=features_df['ANN.Gene_Name'],
                                         annotation_gene_id=features_df['ANN.Gene_ID'],
                                         annotation_feature_type=features_df['ANN.Feature_Type'],
                                         annotation_transcript_biotype=features_df['ANN.Transcript_BioType'],
                                         annotation_hgvs_c=features_df['ANN.HGVS.c'],
                                         annotation_hgvs_p=features_df['ANN.HGVS.p'])

    def get_correct_data_package(self) -> Any:
        return NucleotideSampleDataPackage

    def get_correct_sample_data(self) -> Any:
        return NucleotideSampleData

    def aggregate_feature_column(self) -> Dict[str, Any]:
        return {'TYPE': 'first', '_SAMPLE_ID': SampleSet,
                'ANN.Annotation': 'first', 'ANN.Annotation_Impact': 'first',
                'ANN.Gene_Name': 'first', 'ANN.Gene_ID': 'first', 'ANN.Feature_Type': 'first',
                'ANN.Transcript_BioType': 'first', 'ANN.HGVS.c': 'first', 'ANN.HGVS.p': 'first'}

    def _modify_df_types(self, features_df: pd.DataFrame) -> pd.DataFrame:
        features_df = features_df.copy()
        columns = ['ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                   'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

        # Convert NA to None to properly save in database
        features_df[columns] = features_df[columns].apply(lambda x: x.astype('object'))
        features_df[columns] = features_df[columns].where(pd.notnull(features_df[columns]), None)

        return features_df

    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        samples_with_variants = {sample.name for sample in
                                 self._sample_service.get_samples_with_variants(feature_scope_name)}
        return len(samples_with_variants.intersection(sample_names)) != 0

    def _verify_correct_feature_scope(self, feature_scope_name: str) -> None:
        if feature_scope_name is None:
            raise Exception('feature_scope_name must not be None')

    def build_sample_feature_object(self, sample: Sample, sample_data: SampleData, feature_scope_name: str) -> Any:
        nucleotide_sample_files = cast(NucleotideSampleData, sample_data)

        reference = self._reference_service.find_reference_genome(feature_scope_name)

        sample_nucleotide_variation = SampleNucleotideVariation(reference=reference)
        vcf_file, vcf_index = nucleotide_sample_files.get_vcf_file()
        sample_nucleotide_variation.nucleotide_variants_file = vcf_file
        sample_nucleotide_variation.sample = sample
        sample_nucleotide_variation.masked_regions_file = nucleotide_sample_files.get_mask_file()

        return sample_nucleotide_variation

    def _create_persisted_features_reader(self, sample_data_dict: Dict[str, SampleData],
                                          data_package: SampleDataPackage) -> FeaturesReader:
        sample_data_dict = cast(Dict[str, NucleotideSampleData], sample_data_dict)
        data_package = cast(NucleotideSampleDataPackage, data_package)
        index_unknown_missing = self._index_unknown_missing and data_package.index_unknown_missing()
        if not index_unknown_missing:
            logger.debug(f'index_unknown_missing={index_unknown_missing} so will not '
                         'index missing/unknown positions')
        else:
            logger.debug(f'index_unknown_missing={index_unknown_missing}')
        progress_logger = SamplesProgressLogger(stage_name='Index', stage_number=2, total_samples=len(sample_data_dict))
        return VcfVariantsReader(sample_data_dict, include_masked_regions=index_unknown_missing,
                                 progress_logger=progress_logger,
                                 variants_processor_factory=data_package.get_variants_processor_factory())

    def read_index(self, feature_ids: Union[List[str], Set[str]]) -> Dict[str, FeatureSamples]:
        if isinstance(feature_ids, set):
            feature_ids = list(feature_ids)

        feature_ids_slicer = ListSliceIter(feature_ids, slice_size=self._sql_select_limit)

        feature_ids_to_feature_samples_dict = {}
        logger.debug(f'Reading nucleotide/variation index for {len(feature_ids)} feature ids '
                     f'dividing up into a maximum of {self._sql_select_limit} feature ids per SQL query')
        slice_number = 0
        for feature_ids_slice in feature_ids_slicer.islice():
            logger.log(TRACE_LEVEL, f'Reading feature ids slice={slice_number}')
            feature_samples = self._connection.get_session().query(NucleotideVariantsSamples) \
                .filter(NucleotideVariantsSamples._spdi.in_(feature_ids_slice)) \
                .all()
            feature_ids_samples_subdict = {f.id: f for f in feature_samples}
            feature_ids_to_feature_samples_dict.update(feature_ids_samples_subdict)

            slice_number = slice_number + 1
        logger.debug(f'Finished reading nucleotide/variation index for {len(feature_ids)} feature ids,'
                     f' found a total of {len(feature_ids_to_feature_samples_dict)} already in the database')

        return feature_ids_to_feature_samples_dict
