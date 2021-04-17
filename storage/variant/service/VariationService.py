import logging
from pathlib import Path
from typing import List, Set, Any, Dict, cast

import pandas as pd

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.SampleSet import SampleSet
from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.io.mutation.NucleotideFeaturesReader import NucleotideFeaturesReader
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.mutation.VariationFile import VariationFile
from storage.variant.model.db import NucleotideVariantsSamples, SampleNucleotideVariation, Sample
from storage.variant.service import DatabaseConnection
from storage.variant.service.FeatureService import FeatureService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.io.SampleData import SampleData
from storage.variant.io.mutation.NucleotideSampleData import NucleotideSampleData
from storage.variant.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor

logger = logging.getLogger(__name__)


class VariationService(FeatureService):

    def __init__(self, database_connection: DatabaseConnection, variation_dir: Path,
                 reference_service: ReferenceService, sample_service: SampleService):
        super().__init__(database_connection=database_connection,
                         features_dir=variation_dir,
                         sample_service=sample_service)
        self._reference_service = reference_service

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

    def get_correct_reader(self) -> Any:
        return NucleotideFeaturesReader

    def aggregate_feature_column(self) -> Dict[str, Any]:
        return {'TYPE': 'first', '_SAMPLE_ID': SampleSet}

    def check_samples_have_features(self, sample_names: Set[str], feature_scope_name: str) -> bool:
        samples_with_variants = {sample.name for sample in
                                 self._sample_service.get_samples_with_variants(feature_scope_name)}
        return len(samples_with_variants.intersection(sample_names)) != 0

    def _verify_correct_feature_scope(self, feature_scope_name: str) -> None:
        if feature_scope_name is None:
            raise Exception('feature_scope_name must not be None')

    def build_sample_feature_object(self, sample: Sample,
                                    sample_files: SampleData,
                                    feature_scope_name: str) -> Any:
        nucleotide_sample_files = cast(NucleotideSampleData, sample_files)

        reference = self._reference_service.find_reference_genome(feature_scope_name)

        sample_nucleotide_variation = SampleNucleotideVariation(reference=reference)
        vcf_file, vcf_index = nucleotide_sample_files.get_vcf_file()
        sample_nucleotide_variation.nucleotide_variants_file = vcf_file
        sample_nucleotide_variation.sample = sample
        sample_nucleotide_variation.masked_regions_file = nucleotide_sample_files.get_mask_file()

        return sample_nucleotide_variation

    def _create_persisted_features_reader(self, sample_files_dict: Dict[str, SampleData]) -> FeaturesReader:
        file_processor = NullSampleFilesProcessor()
        for sample in sample_files_dict:
            file_processor.add(sample_files_dict[sample])
        return VcfVariantsReader(sample_files_dict)
