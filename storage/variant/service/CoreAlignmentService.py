from typing import List, Dict
import copy

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from bitarray import bitarray

from storage.variant.CoreBitMask import CoreBitMask
from storage.variant.model import Sample, SampleSequence, Reference, ReferenceSequence, VariationAllele
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService


class CoreAlignmentService:
    ALIGN_TYPES = ['core', 'full']

    def __init__(self, database: DatabaseConnection, reference_service: ReferenceService):
        self._database = database
        self._reference_service = reference_service

    def _all_sample_names(self, reference_name: str) -> List[str]:
        samples = self._database.get_session().query(Sample) \
            .join(SampleSequence) \
            .join(ReferenceSequence) \
            .join(Reference) \
            .filter(Reference.name == reference_name) \
            .all()

        return [s.name for s in samples]

    def _sample_sequence(self, reference_name: str, samples: List[str]) -> Dict[str, List[SampleSequence]]:
        sample_sequences = self._database.get_session().query(SampleSequence) \
            .join(Sample) \
            .join(ReferenceSequence) \
            .join(Reference) \
            .filter(Sample.name.in_(samples)) \
            .filter(Reference.name == reference_name) \
            .all()

        sample_sequences_map = {}
        for ss in sample_sequences:
            if ss.sequence.sequence_name in sample_sequences_map:
                sample_sequences_map[ss.sequence.sequence_name].append(ss)
            else:
                sample_sequences_map[ss.sequence.sequence_name] = [ss]

        return sample_sequences_map

    def _get_variants(self, sequence_name: str) -> Dict[int, Dict[str, VariationAllele]]:
        variants = self._database.get_session().query(VariationAllele) \
            .join(ReferenceSequence) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .filter(VariationAllele.var_type == 'snp') \
            .order_by(VariationAllele.position) \
            .all()

        variants_dict = {}

        for variant in variants:
            if variant.position not in variants_dict:
                variants_dict[variant.position] = {}
            for sample in variant.samples:
                variants_dict[variant.position][sample.name] = variant

        return variants_dict

    def _create_core_mask(self, sequences: List[SampleSequence]) -> CoreBitMask:
        if sequences is None or len(sequences) == 0:
            raise Exception('Cannot create bitmask of empty sequences')
        else:
            core_mask = sequences[0].get_core_mask()
            for i in range(1, len(sequences)):
                sequence = sequences[i]
                core_mask = core_mask.append_bitmask(sequence.get_core_mask())

            return core_mask

    def _build_core_alignment_sequence(self, ref_sequence: SeqRecord,
                                       variants_dict: Dict[int, Dict[str, VariationAllele]],
                                       core_mask: CoreBitMask,
                                       samples: List[str],
                                       include_reference: bool) -> Dict[str, SeqRecord]:
        seq_records = {}
        for position in variants_dict:
            ref = ref_sequence[position - 1:position].seq

            # if in core
            if core_mask[position - 1]:
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

        return seq_records

    def _build_full_alignment_sequence(self, ref_sequence: SeqRecord,
                                       variants_dict: Dict[int, Dict[str, VariationAllele]],
                                       samples: List[str],
                                       sample_sequences: List[SampleSequence],
                                       include_reference: bool) -> Dict[str, SeqRecord]:
        seq_records = {}

        for sample in samples:
            seq_records[sample] = copy.deepcopy(ref_sequence)
            seq_records[sample].id = sample
            seq_records[sample].description = 'generated automatically'
            if isinstance(seq_records[sample].seq, str):
                seq_records[sample].seq = MutableSeq(seq_records[sample].seq)
            else:
                seq_records[sample].seq = seq_records[sample].seq.tomutable()

        if include_reference:
            seq_records['reference'] = ref_sequence

        # Add all variants
        for position in variants_dict:
            variant_samples = variants_dict[position]
            for sample in samples:
                if sample not in seq_records:
                    raise Exception(f'Alignment for sample {sample} is not valid')

                if sample in variant_samples:
                    seq_records[sample].seq[position-1] = variant_samples[sample].alt

        # Mask positions
        for sample in samples:
            seq = seq_records[sample].seq

            # Search for correct sample sequence
            sample_sequence = None
            for s in sample_sequences:
                if s.sample.name == sample:
                    sample_sequence = s
                    break

            core_mask = sample_sequence.get_core_mask()

            for missing_pos in core_mask._core_bitmask.itersearch(bitarray('0')):
                seq[missing_pos] = 'N'

        # Change back to immutable sequences
        for sample in samples:
            seq_records[sample].seq = seq_records[sample].seq.toseq()

        return seq_records

    def construct_alignment(self, reference_name: str, samples: List[str] = None,
                            include_reference: bool = True, align_type: str = 'core') -> MultipleSeqAlignment:
        if samples is None or len(samples) == 0:
            samples = self._all_sample_names(reference_name)

        sample_sequences = self._sample_sequence(reference_name, samples)

        alignment_seqs = {}
        alignments = {}

        for sequence_name in sample_sequences:
            seq = self._reference_service.get_sequence(sequence_name)
            variants_dict = self._get_variants(sequence_name)

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
