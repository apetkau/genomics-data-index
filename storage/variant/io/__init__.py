from typing import Dict
import abc
import pandas as pd

from storage.variant.CoreBitMask import CoreBitMask


class VariantsReader(abc.ABC):

    def __init(self):
        pass

    @abc.abstractmethod
    def get_variants_table(self) -> pd.DataFrame:
        pass

    @abc.abstractmethod
    def get_core_masks(self) -> Dict[str, Dict[str, CoreBitMask]]:
        pass