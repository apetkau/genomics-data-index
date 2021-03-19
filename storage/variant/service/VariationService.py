from typing import Dict, List, Any
import logging
from pathlib import Path

import pandas as pd

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.io import check_variants_table_columns
from storage.variant.model import ReferenceSequence, VariationAllele, Sample, SampleSequence
from storage.variant.service import DatabaseConnection
from storage.variant.service import EntityExistsError
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class VariationService:

    def __init__(self, database_connection: DatabaseConnection, reference_service: ReferenceService,
                 sample_service: SampleService):
        self._connection = database_connection
        self._reference_service = reference_service
        self._sample_service = sample_service

    def get_variants(self, sequence_name: str, type: str = 'snp') -> Dict[int, Dict[str, VariationAllele]]:
        variants = self._connection.get_session().query(VariationAllele) \
            .join(ReferenceSequence) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .filter(VariationAllele.var_type == type) \
            .order_by(VariationAllele.position) \
            .all()

        variants_dict = {}

        for variant in variants:
            if variant.position not in variants_dict:
                variants_dict[variant.position] = {}
            for sample in variant.samples:
                variants_dict[variant.position][sample.name] = variant

        return variants_dict

    def _create_file_variants(self, var_df: pd.DataFrame, ref_contigs: Dict[str, ReferenceSequence]) -> Dict[
        str, List[VariationAllele]]:
        variant_table = {}
        file_variants = {}
        for row in var_df.iterrows():
            sample_name = row[1]['SAMPLE']

            ref_contig = ref_contigs[row[1]['CHROM']]
            variant_id = VariationAllele.spdi(sequence_name=ref_contig.sequence_name,
                                              position=row[1]['POS'],
                                              ref=row[1]['REF'],
                                              alt=row[1]['ALT']
                                              )

            if variant_id not in variant_table:
                variant = VariationAllele(sequence=ref_contig, position=row[1]['POS'],
                                          ref=row[1]['REF'], alt=row[1]['ALT'], var_type=row[1]['TYPE'])
                variant_table[variant.id] = variant
            else:
                variant = variant_table[variant_id]

            if sample_name not in file_variants:
                file_variants[sample_name] = []

            file_variants[sample_name].append(variant)

        return file_variants

    def check_samples_have_variants(self, var_df: pd.DataFrame, reference_name: str) -> bool:
        """
        Checks if any of the samples loaded in the passed dataframe already have variants.
        :param var_df: The dataframe of variants
        :param reference_name: The reference genome name.
        :return: True if any of the samples have variants, false otherwise.
        """
        samples_with_variants = {sample.name for sample in
                                 self._sample_service.get_samples_with_variants(reference_name)}
        samples_from_df = set(var_df['SAMPLE'].tolist())
        return len(samples_with_variants.intersection(samples_from_df)) != 0

    # def insert_variant_files(self, variant_data: Dict[str, Dict[str, Any]]):
    #     """
    #     Inserts nucleotide variation files from the given data stucture.
    #     :param reference_name: The name of the reference genome.
    #     :param variant_data: A dictionary mapping a sample name to the file path and core bitmask.
    #     :return: None
    #     """
    #     for sample_name in variant_data:
    #
    #
    #
    #     ref_contigs = self._reference_service.get_reference_contigs(reference_name)
    #     file_variants = self._create_file_variants(var_df, ref_contigs)
    #
    #     for s in file_variants:
    #         ref_objects = {ref_contigs[v.sequence.sequence_name] for v in file_variants[s]}
    #         sample_core_masks = core_masks[s]
    #         sample_sequences = []
    #         for r in ref_objects:
    #             sample_sequence = SampleSequence(sequence=r)
    #             sample_sequence.core_mask = sample_core_masks[r.sequence_name]
    #             sample_sequences.append(sample_sequence)
    #         sample = Sample(name=s, variants=file_variants[s], sample_sequences=sample_sequences)
    #         self._connection.get_session().add(sample)
    #
    #     self._connection.get_session().commit()

    def insert_variants(self, var_df: pd.DataFrame, reference_name: str,
                        core_masks: Dict[str, Dict[str, CoreBitMask]]) -> None:
        check_variants_table_columns(var_df)
        if self.check_samples_have_variants(var_df, reference_name):
            raise EntityExistsError(f'Passed samples already have variants for reference genome [{reference_name}], '
                                    f'will not insert any new variants')

        ref_contigs = self._reference_service.get_reference_contigs(reference_name)
        file_variants = self._create_file_variants(var_df, ref_contigs)

        for s in file_variants:
            ref_objects = {ref_contigs[v.sequence.sequence_name] for v in file_variants[s]}
            sample_core_masks = core_masks[s]
            sample_sequences = []
            for r in ref_objects:
                sample_sequence = SampleSequence(sequence=r)
                sample_sequence.core_mask = sample_core_masks[r.sequence_name]
                sample_sequences.append(sample_sequence)
            sample = Sample(name=s, variants=file_variants[s], sample_sequences=sample_sequences)
            self._connection.get_session().add(sample)

        self._connection.get_session().commit()

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
