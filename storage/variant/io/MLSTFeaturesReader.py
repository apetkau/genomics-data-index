from typing import Set

import pandas as pd

from storage.variant.io.FeaturesReader import FeaturesReader


class MLSTFeaturesReader(FeaturesReader):

    def __init__(self):
        super().__init__()

    def _minimal_expected_columns(self) -> Set[str]:
        return {'Sample', 'Scheme', 'Locus', 'Allele'}
