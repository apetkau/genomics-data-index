from typing import Dict, List, Any, Set
import logging
from pathlib import Path
import shutil

import pandas as pd

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.io import VariantsReader
from storage.variant.SampleSet import SampleSet
from storage.variant.model import ReferenceSequence, Sample, SampleNucleotideVariation, NucleotideVariantsSamples
from storage.variant.service import DatabaseConnection
from storage.variant.service import EntityExistsError
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class VariationService:

    def __init__(self, database_connection: DatabaseConnection, variation_dir: Path,
                 reference_service: ReferenceService, sample_service: SampleService):
        self._connection = database_connection
        self._reference_service = reference_service
        self._sample_service = sample_service
        self._variation_dir = variation_dir

    def get_variants_ordered(self, sequence_name: str, type: str = 'snp') -> List[NucleotideVariantsSamples]:
        return self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples.sequence == sequence_name) \
            .filter(NucleotideVariantsSamples.var_type == type) \
            .order_by(NucleotideVariantsSamples.position) \
            .all()

    def _create_nucleotide_variants(self, var_df: pd.DataFrame) -> List[NucleotideVariantsSamples]:
        samples_names = set(var_df['SAMPLE'].tolist())
        sample_name_ids = self._sample_service.find_sample_name_ids(samples_names)

        var_df['SPDI'] = var_df.apply(lambda x: NucleotideVariantsSamples.to_spdi(
            sequence_name=x['CHROM'],
            position=x['POS'],
            ref=x['REF'],
            alt=x['ALT']
        ), axis='columns')
        var_df['SAMPLE_ID'] = var_df.apply(lambda x: sample_name_ids[x['SAMPLE']], axis='columns')

        index_df = var_df.groupby('SPDI').agg({'TYPE': 'first', 'SAMPLE_ID': SampleSet}).reset_index()

        return index_df.apply(lambda x: NucleotideVariantsSamples(
            spdi=x['SPDI'], var_type=x['TYPE'], sample_ids=x['SAMPLE_ID']),
            axis='columns').tolist()

    def check_samples_have_variants(self, sample_names: Set[str], reference_name: str) -> bool:
        """
        Checks if any of the passed sample names already have variants.
        :param sample_names: The dataframe of variants
        :param reference_name: The reference genome name.
        :return: True if any of the samples have variants, false otherwise.
        """
        samples_with_variants = {sample.name for sample in
                                 self._sample_service.get_samples_with_variants(reference_name)}
        return len(samples_with_variants.intersection(sample_names)) != 0

    def insert_variants(self, reference_name: str, variants_reader: VariantsReader) -> None:
        reference = self._reference_service.find_reference_genome(reference_name)
        sample_variant_files = variants_reader.sample_variant_files()

        if self.check_samples_have_variants(set(sample_variant_files.keys()), reference_name):
            raise EntityExistsError(f'Passed samples already have variants for reference genome [{reference_name}], '
                                    f'will not insert any new variants')

        for sample_name in sample_variant_files:
            variant_file = sample_variant_files[sample_name]
            sample = Sample(name=sample_name)
            sample_nucleotide_variation = SampleNucleotideVariation(reference=reference)
            sample_nucleotide_variation.nucleotide_variants_file = self._save_variation_file(variant_file, sample)
            sample_nucleotide_variation.sample = sample
            self._connection.get_session().add(sample_nucleotide_variation)
        self._connection.get_session().commit()

        self.index_variants(variants_reader=variants_reader)

    def index_variants(self, variants_reader: VariantsReader):
        variants_df = variants_reader.get_variants_table()
        self._connection.get_session().bulk_save_objects(self._create_nucleotide_variants(variants_df))
        self._connection.get_session().commit()

    def _save_variation_file(self, original_file: Path, sample: Sample) -> Path:
        new_file = self._variation_dir / f'{sample.name}.vcf.gz'
        if new_file.exists():
            raise Exception(f'File {new_file} already exists')

        shutil.copyfile(original_file, new_file)
        return new_file

    def pairwise_distance(self, samples: List[str], var_type='all', distance_type='jaccard') -> pd.DataFrame:
        sample_objs = self._connection.get_session().query(Sample).filter(Sample.name.in_(samples)).all()

        if var_type == 'all':
            sample_variants = {s.name: {v.to_spdi() for v in s.variants} for s in sample_objs}
        else:
            sample_variants = {s.name: {v.to_spdi() for v in s.variants if v.var_type == var_type} for s in sample_objs}

        names = list(sample_variants.keys())
        distances = []
        for name1 in names:
            row = []
            for name2 in names:
                if name1 == name2:
                    row.append(0)
                else:
                    if distance_type == 'jaccard':
                        logger.debug(f'variants1=[{sample_variants[name1]}]')
                        logger.debug(f'variants2=[{sample_variants[name2]}]')
                        intersection = sample_variants[name1].intersection(sample_variants[name2])
                        union = sample_variants[name1].union(sample_variants[name2])

                        row.append(1 - (len(intersection) / len(union)))
                    else:
                        raise Exception(f'Unsupported distance_type=[{distance_type}]')
            distances.append(row)

        return pd.DataFrame(distances, columns=names, index=names)
