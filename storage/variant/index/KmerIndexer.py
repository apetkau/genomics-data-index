from typing import List, Union
import abc
from pathlib import Path
import logging
import subprocess


logger = logging.getLogger(__name__)


class KmerIndexer(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def index(self, index_name: str, index_path: Path, files: List[Path]):
        pass


class KmerIndexerSourmash(KmerIndexer):

    def __init__(self, k: Union[int, List[int]] = 31, scaled: int = None, num: int = None, abund: bool = False):
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

    def index(self, index_name: str, index_path: Path, files: List[Path]):
        command = ['sourmash', 'sketch', 'dna', '-p', self._params,
                   '--name', index_name, '--output', str(index_path)]
        command.extend([str(f) for f in files])
        try:
            logger.debug(f'Running: {" ".join(command)}')
            subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       check=True, text=True)
        except subprocess.CalledProcessError as e:
            err_msg = str(e.stderr.strip())
            raise Exception(f'Could not run sourmash on alignment: error {err_msg}')


class KmerIndexManager:

    def __init__(self, kmer_index_directory: Path):
        self._kmer_index_directory = kmer_index_directory

    def index_genome_files(self, index_name: str, files: List[Path], kmer_indexer: KmerIndexer) -> Path:
        """
        Indexes reads and constructs an index file with the given name from the given files.
        :param index_name: The name of the index.
        :param files: The list of files to index.
        :return: The path to the index.
        """
        index_out = self._kmer_index_directory / f'{index_name}.sig'
        if index_out.exists():
            raise Exception(f'Index output file {index_out} already exists')

        kmer_indexer.index(index_name=index_name,
                           index_path=index_out,
                           files=files)

        return index_out
