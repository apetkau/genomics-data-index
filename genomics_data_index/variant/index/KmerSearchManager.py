import logging
import tempfile
import time
from pathlib import Path
from typing import List

import pandas as pd

from genomics_data_index.variant.util import execute_commands

logger = logging.getLogger(__name__)


class KmerSearchManagerSourmash:

    def __init__(self):
        pass

    def search(self, kmer_size: int, query_file: Path, search_files: List[Path]) -> pd.DataFrame:
        start_time = time.time()
        logger.debug(f'Start search for matches to [{query_file}].')

        with tempfile.TemporaryDirectory() as tmp:
            out_file = Path(tmp, 'output.csv')
            command = ['sourmash', 'search', '-k', str(kmer_size), '-o', str(out_file),
                       str(query_file)]
            command.extend([str(f) for f in search_files])

            execute_commands([command])
            end_time = time.time()
            logger.debug(f'Finished search for matches. Took {end_time - start_time:0.2f} seconds')

            return pd.read_csv(out_file)
