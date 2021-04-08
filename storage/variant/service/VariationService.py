import logging
from pathlib import Path
from typing import List, Set, Any, Dict

import pandas as pd

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.SampleSet import SampleSet
from storage.variant.io.FeaturesReader import FeaturesReader
from storage.variant.io.NucleotideFeaturesReader import NucleotideFeaturesReader
from storage.variant.io.VariationFile import VariationFile
from storage.variant.model import Sample, SampleNucleotideVariation, NucleotideVariantsSamples
from storage.variant.service import DatabaseConnection
from storage.variant.service.FeatureService import FeatureService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class VariationService(FeatureService):

    def __init__(self, database_connection: DatabaseConnection, variation_dir: Path,
                 reference_service: ReferenceService, sample_service: SampleService):
        super().__init__(database_connection=database_connection,
                         features_dir=variation_dir,
                         sample_service=sample_service)
        self._reference_service = reference_service

    def get_variants_ordered(self, sequence_name: str, type: str = 'snp') -> List[NucleotideVariantsSamples]:
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
                                    features_reader: FeaturesReader, feature_scope_name: str) -> Any:
        self._verify_correct_reader(features_reader=features_reader)
        variants_reader: NucleotideFeaturesReader = features_reader

        reference = self._reference_service.find_reference_genome(feature_scope_name)

        feature_file = variants_reader.get_or_create_feature_file(sample.name)
        sample_nucleotide_variation = SampleNucleotideVariation(reference=reference)
        sample_nucleotide_variation.nucleotide_variants_file = self._save_variation_file(feature_file, sample)
        sample_nucleotide_variation.sample = sample

        genomic_masked_regions = variants_reader.get_genomic_masked_regions()

        if sample.name in genomic_masked_regions:
            masked_regions = genomic_masked_regions[sample.name]
        else:
            masked_regions = MaskedGenomicRegions.empty_mask()

        sample_nucleotide_variation.masked_regions_file = self._save_masked_regions_file(masked_regions, sample)

        return sample_nucleotide_variation

    def _save_variation_file(self, original_file: Path, sample: Sample) -> Path:
        new_file = self._features_dir / f'{sample.name}.vcf.gz'
        if new_file.exists():
            raise Exception(f'File {new_file} already exists')

        return VariationFile(original_file).write(new_file)

    def _save_masked_regions_file(self, masked_regions, sample: Sample):
        new_file = self._features_dir / f'{sample.name}.bed.gz'
        if new_file.exists():
            raise Exception(f'File {new_file} already exists')

        masked_regions.write(new_file)
        return new_file
