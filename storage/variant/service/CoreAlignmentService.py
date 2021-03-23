import copy
import logging
import time
from typing import List, Dict

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.model import Sample, Reference, ReferenceSequence
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.VariationService import VariationService

logger = logging.getLogger(__name__)


class CoreAlignmentService:
    ALIGN_TYPES = ['core', 'full']

    def __init__(self, database: DatabaseConnection, reference_service: ReferenceService,
                 sample_service: SampleService, variation_service: VariationService):
        self._database = database
        self._reference_service = reference_service
        self._variation_service = variation_service
        self._sample_service = sample_service

    def _all_sample_names(self, reference_name: str) -> List[str]:
        samples = self._sample_service.get_samples_with_variants(reference_name)
        return [s.name for s in samples]

    def _create_core_mask(self, sequences: List[SampleSequence]) -> MaskedGenomicRegions:
        if sequences is None or len(sequences) == 0:
            raise Exception('Cannot create bitmask of empty sequences')
        else:
            start_time = time.time()
            logger.debug(f'Started creating core mask for {len(sequences)} sequences')

            core_mask = sequences[0].core_mask
            for i in range(1, len(sequences)):
                sequence = sequences[i]
                core_mask = core_mask.append_bitmask(sequence.core_mask)

            end_time = time.time()
            logger.debug(f'Finished creating core mask for {len(sequences)} sequences. '
                         f'Took {end_time - start_time:0.2f} seconds')

            return core_mask

    def _build_core_alignment_sequence(self, ref_sequence: SeqRecord,
                                       variants_dict: Dict[int, Dict[str, VariationAllele]],
                                       core_mask: MaskedGenomicRegions,
                                       samples: List[str],
                                       include_reference: bool) -> Dict[str, SeqRecord]:
        start_time = time.time()
        logger.debug(f'Started building core alignment for {ref_sequence.id} using {len(variants_dict)} '
                     'variant positions')

        seq_records = {}
        for position in variants_dict:
            ref = ref_sequence[position - 1:position].seq

            # if in core
            if position in core_mask:
                variant_samples = variants_dict[position]
                if len(set(samples).intersection(set(variant_samples.keys()))) == 0:
                    continue

                for sample in samples:
                    if sample not in seq_records:
                        seq_records[sample] = SeqRecord(
                            seq=Seq(''), id=sample, description='generated automatically')

                    if sample in variant_samples:
                        seq_records[sample] += variant_samples[sample].alt
                    else:
                        seq_records[sample] += ref

                if include_reference:
                    # Add the reference sequence in
                    if 'reference' not in seq_records:
                        seq_records['reference'] = SeqRecord(
                            seq=Seq(''), id='reference', description='generated automatically')

                    seq_records['reference'] += ref

        end_time = time.time()
        logger.debug(f'Finished building core alignment for {ref_sequence.id}. '
                     f'Took {end_time - start_time:0.2f} seconds')

        return seq_records

    def _build_full_alignment_sequence(self, ref_sequence: SeqRecord,
                                       variants_dict: Dict[int, Dict[str, VariationAllele]],
                                       samples: List[str],
                                       sample_sequences: List[SampleSequence],
                                       include_reference: bool) -> Dict[str, SeqRecord]:
        seq_records = {}

        start_time = time.time()
        logger.debug(f'Started building full alignment for {ref_sequence.id}'
                     f' using {len(variants_dict)} variant positions')

        for sample in samples:
            seq_records[sample] = copy.deepcopy(ref_sequence)
            seq_records[sample].id = sample
            seq_records[sample].description = 'generated automatically'
            if isinstance(seq_records[sample].seq, str):
                seq_records[sample].seq = MutableSeq(seq_records[sample].seq)
            else:
                seq_records[sample].seq = seq_records[sample].seq.tomutable()

        # Add all variants
        for position in variants_dict:
            variant_samples = variants_dict[position]
            for sample in samples:
                if sample not in seq_records:
                    raise Exception(f'Alignment for sample {sample} is not valid')

                if sample in variant_samples:
                    seq_records[sample].seq[position - 1] = variant_samples[sample].alt

        # Mask positions
        for sample in samples:
            seq = seq_records[sample].seq

            # Search for correct sample sequence
            sample_sequence = None
            for s in sample_sequences:
                if s.sample.name == sample:
                    sample_sequence = s
                    break

            core_mask = sample_sequence.core_mask

            for missing_pos in core_mask.iter_missing_positions():
                seq[missing_pos - 1] = 'N'

        # Change back to immutable sequences
        for sample in samples:
            seq_records[sample].seq = seq_records[sample].seq.toseq()

        if include_reference and 'reference' in seq_records:
            raise Exception('Error, [reference] is a sample name so cannot add "reference" sequence to alignment')
        elif include_reference:
            seq_records['reference'] = copy.deepcopy(ref_sequence)
            seq_records['reference'].id = 'reference'
            seq_records['reference'].description = 'generated automatically'
            if isinstance(seq_records['reference'].seq, str):
                seq_records['reference'].seq = Seq(seq_records['reference'].seq)

        end_time = time.time()
        logger.debug(f'Finished building full alignment for {ref_sequence.id}. '
                     f'Took {end_time - start_time:0.2f} seconds')

        return seq_records

    def construct_alignment(self, reference_name: str, samples: List[str] = None,
                            include_reference: bool = True, align_type: str = 'core') -> MultipleSeqAlignment:
        if samples is None or len(samples) == 0:
            samples = self._all_sample_names(reference_name)

        reference_sequences = self._sample_sequence_service.get_sample_sequences(reference_name, samples)

        alignment_seqs = {}
        alignments = {}

        for sequence_name in sample_sequences:
            seq = self._reference_service.get_sequence(sequence_name)
            variants_dict = self._variation_service.get_variants(sequence_name, type='snp')

            if align_type == 'core':
                core_mask = self._create_core_mask(sample_sequences[sequence_name])
                alignment_seqs[sequence_name] = self._build_core_alignment_sequence(
                    ref_sequence=seq,
                    variants_dict=variants_dict,
                    core_mask=core_mask,
                    samples=samples,
                    include_reference=include_reference
                )
            elif align_type == 'full':
                alignment_seqs[sequence_name] = self._build_full_alignment_sequence(
                    ref_sequence=seq,
                    variants_dict=variants_dict,
                    samples=samples,
                    sample_sequences=sample_sequences[sequence_name],
                    include_reference=include_reference
                )
            else:
                raise Exception(f'Unknown value for align_type=[{align_type}]. Must be one of {self.ALIGN_TYPES}')

        for sequence_name in alignment_seqs:
            alignments[sequence_name] = MultipleSeqAlignment(
                alignment_seqs[sequence_name].values())
            alignments[sequence_name].sort()

        sequence_names = sorted(alignments.keys())
        seq1 = sequence_names.pop()
        core_snv_align = alignments[seq1]

        for sequence_name in sequence_names:
            core_snv_align += alignments[sequence_name]

        return core_snv_align
