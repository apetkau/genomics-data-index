from pathlib import Path
from typing import Tuple, List

from Bio.SeqRecord import SeqRecord

from genomics_data_index.storage.util import parse_sequence_file


class SequenceFile:

    def __init__(self, file: Path):
        self._file = file

    def parse_sequence_file(self) -> Tuple[str, List[SeqRecord]]:
        return parse_sequence_file(self._file)
