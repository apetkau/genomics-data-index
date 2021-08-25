import copy
import logging
import tempfile
import time
from pathlib import Path
from typing import List, Dict, Generator

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.VariationFile import VariationFile
from genomics_data_index.storage.model.db import SampleNucleotideVariation
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.ReferenceService import ReferenceService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.VariationService import VariationService

logger = logging.getLogger(__name__)


class CoreAlignmentService:
    ALIGN_TYPES = ['core', 'full']
    INCLUDE_VARIANT_SUBSTITUTION = ['SNP', 'MNP']
    INCLUDE_VARIANT_DELETION = ['DELETION', 'DELETION_OTHER']
    INCLUDE_VARIANT_TYPES = INCLUDE_VARIANT_SUBSTITUTION + INCLUDE_VARIANT_DELETION
    INCLUDE_VARIANT_DEFAULT = ['SNP', 'MNP', 'DELETION']

    # Mask generated sequence with this character (which should not appear anywhere else) so I can remove
    # all positions with this character later (when generating core alignment)
    CORE_MASK_CHAR = '?'

    def __init__(self, database: DatabaseConnection, reference_service: ReferenceService,
                 sample_service: SampleService, variation_service: VariationService):
        self._database = database
        self._reference_service = reference_service
        self._variation_service = variation_service
        self._sample_service = sample_service

    def _all_sample_names(self, reference_name: str) -> List[str]:
        samples = self._sample_service.get_samples_with_variants(reference_name)
        return [s.name for s in samples]

    def _create_core_mask(self, sample_variations: List[SampleNucleotideVariation]) -> MaskedGenomicRegions:
        masked_regions = [v.masked_regions for v in sample_variations]
        return MaskedGenomicRegions.union_all(masked_regions)

    def _get_core_positions(self, sample_variations: List[SampleNucleotideVariation],
                            core_mask: MaskedGenomicRegions, include_expression: str) -> Dict[str, List[int]]:
        variation_files = [v.nucleotide_variants_file for v in sample_variations]
        union_df = VariationFile.union_all_files(variation_files, include_expression=include_expression)
        union_df = union_df.sort_values(['CHROM', 'POS'])
        core_positions = {}
        for index, row in union_df.iterrows():
            seq_name = row['CHROM']
            position = row['POS']
            if not core_mask.contains(seq_name, position, start_position_index='1'):
                if seq_name not in core_positions:
                    core_positions[seq_name] = [position]
                else:
                    core_positions[seq_name].append(position)

        return core_positions

    def _extract_sequence_from_positions(self, input_seq: Seq, positions: List[int]) -> Seq:
        sequence_string = ''

        for position in positions:
            sequence_string += str(input_seq[position - 1:position])

        return Seq(sequence_string)

    def _core_alignment_sequence_generator(self, reference_file: Path, core_positions: Dict[str, List[int]],
                                           sample_variations: List[SampleNucleotideVariation],
                                           include_expression: str) -> Generator[
        Dict[str, SeqRecord], None, None]:
        for sample_variation in sample_variations:
            seq_records = {}

            consensus_records = VariationFile(sample_variation.nucleotide_variants_file).consensus(
                reference_file, include_expression=include_expression)
            for record in consensus_records:
                sequence_name = record.id
                record.id = sample_variation.sample.name
                record.description = 'generated automatically'
                record.seq = self._extract_sequence_from_positions(input_seq=record.seq,
                                                                   positions=core_positions[sequence_name])

                seq_records[sequence_name] = record
            yield seq_records

    def _full_alignment_sequence_generator(self, reference_file: Path,
                                           sample_variations: List[SampleNucleotideVariation],
                                           include_expression: str) -> Generator[
        Dict[str, SeqRecord], None, None]:
        for sample_variation in sample_variations:
            seq_records = {}

            consensus_records = VariationFile(sample_variation.nucleotide_variants_file).consensus(
                reference_file=reference_file,
                include_expression=include_expression,
                mask_file=sample_variation.masked_regions_file
            )
            for record in consensus_records:
                seq_records[record.id] = record
                record.id = sample_variation.sample.name
                record.description = 'generated automatically'
            yield seq_records

    def construct_alignment(self, reference_name: str, samples: List[str] = None,
                            include_reference: bool = True, align_type: str = 'core',
                            include_variants: List[str] = None) -> MultipleSeqAlignment:
        if samples is None or len(samples) == 0:
            samples = self._all_sample_names(reference_name)

        if include_variants is None:
            include_variants = self.INCLUDE_VARIANT_DEFAULT

        if align_type == 'core' and include_variants != ['SNP']:
            raise Exception(f'align_type={align_type} and include_variants={include_variants}. Currently'
                            f' align_type=core only works with include_variants=["SNP"]')
        else:
            subtitution_types = []
            deletion_types = []
            for variant_type in include_variants:
                if variant_type not in self.INCLUDE_VARIANT_TYPES:
                    raise Exception(f'variant_type={variant_type} found in include_variants={include_variants}. '
                                    f'Only {self.INCLUDE_VARIANT_TYPES} are supported')
                elif variant_type in self.INCLUDE_VARIANT_SUBSTITUTION:
                    subtitution_types.append(variant_type)
                elif variant_type in self.INCLUDE_VARIANT_DELETION:
                    deletion_types.append(variant_type)

        include_expression = ''
        if len(deletion_types) == 1:
            deletion_type = deletion_types[0]
            # To make sure I only include deletions (out of all INDELs) I use (ILEN<0) which
            # represents INDELs that are deletions (see <http://samtools.github.io/bcftools/bcftools.html#expressions>).
            if deletion_type == 'DELETION':
                include_expression = '((INFO/TYPE=="INDEL") & (ILEN<0))'
            elif deletion_type == 'DELETION_OTHER':
                include_expression = '(ILEN<0)'
            else:
                raise Exception(f'deletion_type={deletion_type} is invalid. Must be one of '
                                f'{self.INCLUDE_VARIANT_DELETION}')
        elif len(deletion_types) > 1:
            raise Exception(f'Can only set one of {self.INCLUDE_VARIANT_DELETION}, got {deletion_types}')

        if len(subtitution_types) > 0 and len(include_expression) > 0:
            include_expression = include_expression + ' | '

        include_expression = include_expression + ' | '.join(
            f'(INFO/TYPE=="{variant}")' for variant in subtitution_types)

        sample_nucleotide_variants = self._variation_service.get_sample_nucleotide_variation(samples)

        alignment_seqs = {}
        alignments = {}

        start_time = time.time()
        logger.info(f'Started building alignment for {len(samples)} samples with include_variants={include_variants}')
        logger.debug(f'Alignment include_expression=\'{include_expression}\'')

        with tempfile.TemporaryDirectory() as tmp_dir:
            # Write reference genome to file for creating consensus sequences
            reference_file = Path(tmp_dir) / 'reference.fasta'
            reference_records = self._reference_service.get_reference_genome_records(reference_name)
            with open(reference_file, 'w') as h:
                SeqIO.write(reference_records, h, 'fasta')

            if align_type == 'core':
                core_mask = self._create_core_mask(sample_nucleotide_variants)
                core_positions = self._get_core_positions(sample_variations=sample_nucleotide_variants,
                                                          core_mask=core_mask,
                                                          include_expression=include_expression)

                alignment_generator = self._core_alignment_sequence_generator(
                    reference_file=reference_file,
                    sample_variations=sample_nucleotide_variants,
                    core_positions=core_positions,
                    include_expression=include_expression
                )

                if include_reference:
                    for record in reference_records:
                        sequence_name = record.id
                        core_seq = self._extract_sequence_from_positions(input_seq=record.seq,
                                                                         positions=core_positions[sequence_name])
                        new_record = SeqRecord(id=reference_name,
                                               seq=core_seq,
                                               description='[reference genome]')
                        alignment_seqs[sequence_name] = [new_record]

            elif align_type == 'full':
                alignment_generator = self._full_alignment_sequence_generator(
                    reference_file=reference_file,
                    sample_variations=sample_nucleotide_variants,
                    include_expression=include_expression
                )

                if include_reference:
                    # Add the reference sequence in
                    reference_records_align = copy.deepcopy(reference_records)
                    for record in reference_records_align:
                        alignment_seqs[record.id] = [record]
                        record.id = reference_name
                        record.description = '[reference genome]'
            else:
                raise Exception(f'Unknown value for align_type=[{align_type}]. Must be one of {self.ALIGN_TYPES}')

            # Actually generate the data
            for sequence_record in alignment_generator:
                for sequence in sequence_record:
                    if sequence == reference_name:
                        raise Exception(
                            f'Data contains a sequence with name [{sequence}], which is the same as the reference genome name')
                    elif sequence not in alignment_seqs:
                        alignment_seqs[sequence] = [sequence_record[sequence]]
                    else:
                        alignment_seqs[sequence].append(sequence_record[sequence])

        for sequence_name in alignment_seqs:
            alignments[sequence_name] = MultipleSeqAlignment(
                alignment_seqs[sequence_name])
            alignments[sequence_name].sort()

        sequence_names = sorted(alignments.keys())
        seq1 = sequence_names.pop()
        snv_align = alignments[seq1]

        for sequence_name in sequence_names:
            snv_align += alignments[sequence_name]

        end_time = time.time()
        logger.info(f'Finished building alignment for {len(samples)} samples. '
                    f'Took {end_time - start_time:0.2f} seconds')

        return snv_align
