import logging
from pathlib import Path
from typing import List, Set, Any, Dict, cast, Union

import pandas as pd

from genomics_data_index.variant.SampleSet import SampleSet
from genomics_data_index.variant.io.FeaturesReader import FeaturesReader
from genomics_data_index.variant.io.SampleData import SampleData
from genomics_data_index.variant.io.SampleDataPackage import SampleDataPackage
from genomics_data_index.variant.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.variant.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.variant.io.mutation.VariationFile import VariationFile
from genomics_data_index.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.variant.model.db import NucleotideVariantsSamples, SampleNucleotideVariation, Sample
from genomics_data_index.variant.service import DatabaseConnection
from genomics_data_index.variant.service.FeatureService import FeatureService
from genomics_data_index.variant.service.ReferenceService import ReferenceService
from genomics_data_index.variant.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class VariationService(FeatureService):
    MUTATION_TYPES = ['snp', 'indel', 'all', 'other']

    def __init__(self, database_connection: DatabaseConnection, variation_dir: Path,
                 reference_service: ReferenceService, sample_service: SampleService):
        super().__init__(database_connection=database_connection,
                         features_dir=variation_dir,
                         sample_service=sample_service)
        self._reference_service = reference_service

    def _reference_sequence_names(self, reference_name: str) -> List[str]:
        return list(self._reference_service.get_reference_sequences(reference_name).keys())

    def count_on_reference(self, reference_name: str, include_unknown: bool) -> int:
        reference_sequence_names = self._reference_sequence_names(reference_name)
        return self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence.in_(reference_sequence_names)) \
            .count()

    def mutation_counts_on_reference(self, reference_name: str, include_unknown: bool) -> Dict[str, int]:
        reference_sequence_names = self._reference_sequence_names(reference_name)
        mutations = self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence.in_(reference_sequence_names)) \
            .all()

        return {m.spdi: len(m.sample_ids) for m in mutations}

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

    def _create_feature_object(self, features_df: pd.DataFrame):
        return NucleotideVariantsSamples(spdi=features_df['_FEATURE_ID'], var_type=features_df['TYPE'],
                                         sample_ids=features_df['_SAMPLE_ID'])

    def get_correct_data_package(self) -> Any:
        return NucleotideSampleDataPackage

    def get_correct_sample_data(self) -> Any:
        return NucleotideSampleData

    def aggregate_feature_column(self) -> Dict[str, Any]:
        return {'TYPE': 'first', '_SAMPLE_ID': SampleSet}

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
        return VcfVariantsReader(sample_data_dict)
