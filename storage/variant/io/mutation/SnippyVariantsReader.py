from __future__ import annotations

import logging
from typing import Dict

from storage.variant.io.mutation.NucleotideSampleData import NucleotideSampleData
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader

logger = logging.getLogger(__name__)


class SnippyVariantsReader(VcfVariantsReader):

    def __init__(self, sample_files_map: Dict[str, NucleotideSampleData]):
        super().__init__(sample_files_map=sample_files_map)
