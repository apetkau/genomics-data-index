import copy
import logging
import time
from typing import List, Dict, Generator
from pathlib import Path
import tempfile

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.VariationFile import VariationFile
from storage.variant.model import Sample, Reference, ReferenceSequence, NucleotideVariantsSamples, SampleNucleotideVariation
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

    def _create_core_mask(self, sample_variations: List[SampleNucleotideVariation]) -> MaskedGenomicRegions:
        masked_regions = [v.masked_regions for v in sample_variations]
        return MaskedGenomicRegions.union_all(masked_regions)

    def _build_core_alignment_sequence(self, ref_sequence: SeqRecord,
                                       variants_ordered: List[NucleotideVariantsSamples],
                                       core_mask: MaskedGenomicRegions,
                                       samples: List[Sample],
                                       include_reference: bool) -> Dict[str, SeqRecord]:
        start_time = time.time()
        logger.debug(f'Started building core alignment for {ref_sequence.id} using {len(variants_ordered)} '
                     'variant positions')

        seq_records = {}
        selected_sample_ids = {s.id for s in samples}
        for variant in variants_ordered:
            position = variant.position
            ref = ref_sequence[position - 1:position].seq

            # if in core
            if position in core_mask:
                if len(selected_sample_ids.intersection(variant.sample_ids)) == 0:
                    continue

                for sample in samples:
                    if sample.name not in seq_records:
                        seq_records[sample.name] = SeqRecord(
                            seq=Seq(''), id=sample.name, description='generated automatically')

                    if sample.id in variant.sample_ids:
                        seq_records[sample] += variant.insertion
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

    def _full_alignment_sequence_generator(self, reference_file: Path,
                                           sample_variations: List[SampleNucleotideVariation]) -> Generator[Dict[str, SeqRecord], None, None]:
        for sample_variation in sample_variations:
            seq_records = {}

            consensus_records = VariationFile(sample_variation.nucleotide_variants_file).consensus(
                reference_file=reference_file,
                mask_file=sample_variation.masked_regions_file
            )
            for record in consensus_records:
                seq_records[record.id] = record
                record.id = sample_variation.sample.name
                record.description = 'generated automatically'
            yield seq_records

    def construct_alignment(self, reference_name: str, samples: List[str] = None,
                            include_reference: bool = True, align_type: str = 'core') -> MultipleSeqAlignment:
        if samples is None or len(samples) == 0:
            samples = self._all_sample_names(reference_name)

        sample_nucleotide_variants = self._variation_service.get_sample_nucleotide_variation(samples)

        alignment_seqs = {}
        alignments = {}

        start_time = time.time()
        logger.debug(f'Started building alignment for {len(samples)} samples')

        with tempfile.TemporaryDirectory() as tmp_dir:
            # Write reference genome to file for creating consensus sequences
            reference_file = Path(tmp_dir) / 'reference.fasta'
            reference_records = self._reference_service.get_reference_genome_records(reference_name)
            with open(reference_file, 'w') as h:
                SeqIO.write(reference_records, h, 'fasta')

            if align_type == 'core':
                raise Exception('align_type=core not implemented')
            elif align_type == 'full':
                alignment_generator = self._full_alignment_sequence_generator(
                    reference_file=reference_file,
                    sample_variations=sample_nucleotide_variants
                )

                if include_reference:
                    # Add the reference sequence in
                    reference_records_align = copy.deepcopy(reference_records)
                    for record in reference_records_align:
                        alignment_seqs[record.id] = [record]
                        record.id = 'reference'
                        record.description = 'generated automatically'

                for sequence_record in alignment_generator:
                    for sequence in sequence_record:
                        if sequence not in alignment_seqs:
                            alignment_seqs[sequence] = [sequence_record[sequence]]
                        else:
                            alignment_seqs[sequence].append(sequence_record[sequence])
            else:
                raise Exception(f'Unknown value for align_type=[{align_type}]. Must be one of {self.ALIGN_TYPES}')

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
        logger.debug(f'Finished building alignment for {len(samples)} samples. '
                     f'Took {end_time - start_time:0.2f} seconds')

        return snv_align
