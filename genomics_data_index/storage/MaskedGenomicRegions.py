from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List, Set, Dict, Generator, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pybedtools import BedTool

from genomics_data_index.storage.model import NUCLEOTIDE_UNKNOWN, NUCLEOTIDE_UNKNOWN_TYPE


class MaskedGenomicRegions:

    def __init__(self, mask: BedTool):
        self._mask = mask.sort().merge()

    def intersect(self, other: MaskedGenomicRegions) -> MaskedGenomicRegions:
        return MaskedGenomicRegions(self._mask.intersect(other._mask))

    def union(self, other: MaskedGenomicRegions) -> MaskedGenomicRegions:
        union = self._mask.cat(other._mask, postmerge=True, force_truncate=True)
        return MaskedGenomicRegions(union)

    def mask_to_features(self) -> pd.DataFrame:
        mask_features = []
        ref = 1
        alt = NUCLEOTIDE_UNKNOWN
        nuc_type = NUCLEOTIDE_UNKNOWN_TYPE
        for sequence_name, position in self.positions_iter(start_position_index='1'):
            variant_id = f'{sequence_name}:{position}:{ref}:{alt}'
            mask_features.append([sequence_name, position, ref, alt, nuc_type, variant_id])

        return pd.DataFrame(mask_features, columns=['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'VARIANT_ID'])

    def mask_genome(self, genome_file: Path, mask_char: str = '?', remove: bool = True) -> Dict[str, SeqRecord]:
        """
        Gets a SeqRecord with all those regions on the passed genome that are in the masked regions removed
        (or masked with mask_char).
        :param genome_file: The genome file to mask.
        :param mask_char: The character to mask with.
        :param remove: Whether or not to remove masked sequence data.
        :return: A Dictionary mapping a sequence name to a SeqRecord containing all those regions on the sequence
                 within the masked regions removed (or masked with mask_char)
        """
        with tempfile.TemporaryDirectory() as out_f:
            seq_records = {}
            output_fasta = Path(out_f) / 'masked.fasta'
            self._mask.mask_fasta(fi=str(genome_file), fo=str(output_fasta), mc=mask_char)
            for record in SeqIO.parse(output_fasta, 'fasta'):
                if remove:
                    record.seq = record.seq.ungap(mask_char)
                seq_records[record.id] = record
            return seq_records

    def write(self, file: Path):
        self._mask.saveas(str(file), compressed=True)

    @classmethod
    def union_all(cls, masked_regions: List[MaskedGenomicRegions]):
        if len(masked_regions) == 0:
            raise Exception('Cannot merge empty list')
        elif len(masked_regions) == 1:
            return masked_regions[0]
        else:
            start_mask = masked_regions.pop()
            union = start_mask._mask.cat(*[o._mask for o in masked_regions], postmerge=True, force_truncate=True)
            return MaskedGenomicRegions(union)

    @classmethod
    def from_sequences(cls, sequences: List[SeqRecord]) -> MaskedGenomicRegions:
        def is_missing(char):
            return char.upper() == 'N' or char == '-'

        # pybedtools internally stores as 0-based BED file intervals
        # https://daler.github.io/pybedtools/intervals.html#bed-is-0-based-others-are-1-based
        mask_intervals = []

        for record in sequences:
            start = 0
            in_mask = False
            for idx, char in enumerate(record.seq):
                if in_mask:
                    if not is_missing(char):
                        in_mask = False
                        # pybedtools stop position is not included in interval
                        stop = idx
                        mask_intervals.append((record.id, start, stop))
                else:
                    if is_missing(char):
                        in_mask = True
                        start = idx

            # Finish recording last interval if it exists (e.g., if last bit of sequence was like 'NNNN')
            if in_mask:
                stop = len(record)
                mask_intervals.append((record.id, start, stop))

        bedtool_intervals = BedTool(mask_intervals)
        return MaskedGenomicRegions(bedtool_intervals)

    @classmethod
    def from_file(cls, file: Path) -> MaskedGenomicRegions:
        bed_file_data = BedTool(str(file))
        return MaskedGenomicRegions(bed_file_data)

    @classmethod
    def empty_mask(cls):
        return MaskedGenomicRegions(BedTool('', from_string=True))

    def is_empty(self):
        return len(self) == 0

    def sequence_names(self) -> Set[str]:
        """
        Gets a set of sequence names from this genomic regions mask.
        :return: A set of sequence names.
        """
        return {x.chrom for x in self._mask}

    def contains(self, sequence: str, position: int, start_position_index: str = '0') -> bool:
        if start_position_index != '0' and start_position_index != '1':
            raise Exception((f'Unknown value start_position_index=[{start_position_index}].'
                             'Should be "0" or "1" to indicate which is the starting base position'))
        elif start_position_index == '1':
            position = position - 1

        for i in self._mask:
            if i.chrom == sequence and i.start <= position < i.end:
                return True
        return False

    def _validate_start_position_index(self, start_position_index: str) -> None:
        if start_position_index not in ['0', '1']:
            raise Exception((f'Unknown value start_position_index=[{start_position_index}].'
                             'Should be "0" or "1" to indicate which is the starting base position'))

    def overlaps_range(self, sequence: str, start: int, stop: int, start_position_index: str = '0') -> bool:
        self._validate_start_position_index(start_position_index)

        if start_position_index == '1':
            start = start - 1
            stop = stop - 1

        if stop <= start:
            raise Exception(f'start=[{start}] is less than stop=[{stop}]')

        for i in self._mask:
            if i.chrom == sequence:
                if i.start <= start and i.end > start:
                    return True
                elif start < i.end and stop > i.end:
                    return True
        return False

    def positions_iter(self, start_position_index: str = '0') -> Generator[Tuple[str, int], None, None]:
        """
        Creates an iterator to iterate over all ('sequence', 'position') in this mask.
        :param start_position_index: Whether positions should be in 0-base coordinates or 1-base coordinates.
                                     See https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format
                                     for a description of the differences in coordinates.
        :return: An iterator which will return tuples like ('sequence', 'position') for every
                 position in this mask.
        """
        self._validate_start_position_index(start_position_index)

        for sequence in self._mask:
            sequence_name = sequence.chrom
            start = sequence.start
            end = sequence.end

            if start_position_index == '1':
                start = start + 1
                end = end + 1

            for pos in range(start, end):
                yield sequence_name, pos

    def __len__(self) -> int:
        """
        Calculates length of underlying masked intervals. Assumes the intervals have been merged beforehand.
        :return: The length of the masked intervals.
        """
        total = 0
        for i in self._mask:
            total += len(i)
        return total
