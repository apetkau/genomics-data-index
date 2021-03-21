from __future__ import annotations

from typing import List

from pybedtools import BedTool
from Bio.SeqRecord import SeqRecord


class MaskedGenomicRegions:

    def __init__(self, mask: BedTool):
        self._mask = mask.sort().merge()

    def intersect(self, other: MaskedGenomicRegions) -> MaskedGenomicRegions:
        return MaskedGenomicRegions(self._mask.intersect(other._mask))

    def union(self, other: MaskedGenomicRegions) -> MaskedGenomicRegions:
        union = self._mask.cat(other._mask, postmerge=True, force_truncate=True)
        return MaskedGenomicRegions(union)

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
    def empty_mask(cls):
        return MaskedGenomicRegions(BedTool('', from_string=True))

    def is_empty(self):
        return len(self) == 0

    def contains(self, sequence: str, position: int) -> bool:
        for i in self._mask:
            if i.chrom == sequence and i.start <= position < i.end:
                return True
        return False

    def __len__(self) -> int:
        """
        Calculates length of underlying masked intervals. Assumes the intervals have been merged beforehand.
        :return: The length of the masked intervals.
        """
        total = 0
        for i in self._mask:
            total += len(i)
        return total
