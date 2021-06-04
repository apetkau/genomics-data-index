import abc
import logging
import os
import time
from pathlib import Path
from typing import List, Union, Tuple, Dict

from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)


class KmerIndexer(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def index(self, index_name: str, index_path: Path, files: List[Path]):
        pass

    @abc.abstractmethod
    def is_compress(self):
        pass


class KmerIndexerSourmash(KmerIndexer):

    def __init__(self, k: Union[int, List[int]] = 31, scaled: int = None, num: int = None,
                 abund: bool = False, compress: bool = True):
        super().__init__()

        if k is None:
            raise Exception('Cannot index with k=None')
        elif not isinstance(k, list) and not isinstance(k, int):
            raise Exception(f'Invalid parameter k={k}')

        if scaled is None and num is None:
            raise Exception('Must set one of "scaled" or "num" to index')
        elif scaled is not None and num is not None:
            raise Exception(f'Cannot set both scaled={scaled} and num={num}')

        self._k = k
        self._scaled = scaled
        self._num = num
        self._abund = abund

        self._params = ''
        if isinstance(k, list):
            self._params += ','.join([f'k={v}' for v in k])
        else:
            self._params += f'k={k}'

        if scaled is not None:
            self._params += f',scaled={scaled}'
        else:
            self._params += f',num={num}'

        if abund:
            self._params += ',abund'

        self._compress = compress

    def is_compress(self):
        return self._compress

    def index(self, index_name: str, index_path: Path, files: List[Path]) -> Path:
        if self._compress:
            if index_path.suffix != '.gz':
                logger.warning((f'index_path=[{index_path}] does not end in ".gz" but compress={self._compress}. '
                                'Adding ".gz" to index_path.'))
                index_path = Path(str(index_path) + '.gz')
            index_out = Path(str(index_path) + '.tmp')
        else:
            index_out = index_path

        command = ['sourmash', 'sketch', 'dna', '-p', self._params,
                   '--name', index_name, '--output', str(index_out)]
        command.extend([str(f) for f in files])
        if self._compress:
            execute_commands([command, ['gzip', str(index_out)]])
            index_gzip_out = Path(str(index_out) + '.gz')
            os.rename(index_gzip_out, index_path)
        else:
            execute_commands([command])

        return index_path


class KmerIndexManager:

    def __init__(self, kmer_index_directory: Path, kmer_indexer: KmerIndexer):
        self._kmer_index_directory = kmer_index_directory
        self._kmer_indexer = kmer_indexer

    def index_single_genome(self, index_name: str, files: List[Path]) -> Path:
        """
        Indexes reads and constructs an index file with the given name from the given files.
        :param index_name: The name of the index.
        :param files: The list of files to index.
        :return: The path to the index.
        """
        if self._kmer_indexer.is_compress():
            index_out = self._kmer_index_directory / f'{index_name}.sig.gz'
        else:
            index_out = self._kmer_index_directory / f'{index_name}.sig'

        if index_out.exists():
            raise Exception(f'Index output file {index_out} already exists')

        indexed_path = self._kmer_indexer.index(index_name=index_name,
                                                index_path=index_out,
                                                files=files)

        return indexed_path

    def index_all_genomes(self, genomes_files: List[Tuple[str, List[Path]]]) -> Dict[str, Path]:
        start_time = time.time()
        logger.debug('Start building kmer indexes')

        indexed_genomes = {}

        for genome_files in genomes_files:
            index_name = genome_files[0]
            index_files = genome_files[1]
            logger.debug(f'Working on {index_name}')
            indexed_path = self.index_single_genome(index_name=index_name, files=index_files)
            logger.debug(f'Finished creating index {indexed_path}')

            indexed_genomes[index_name] = indexed_path

        end_time = time.time()
        logger.debug(f'Finished building kmer indexes. Took {end_time - start_time:0.2f} seconds')

        return indexed_genomes
