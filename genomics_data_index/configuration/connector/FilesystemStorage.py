import logging
import os
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class FilesystemStorage:
    ALL_SUBDIRECTORIES = ['reference', 'kmer', 'variation', 'mlst']

    def __init__(self, root_dir: Path):
        self._root_dir = root_dir

        if not root_dir.exists():
            os.mkdir(root_dir)

    def _check_make_dir(self, name) -> Path:
        d = self._root_dir / name
        if not d.exists():
            os.mkdir(d)
        return d

    def get_storage_size(self) -> pd.DataFrame:
        '''
        Gets the total size of all files on the disk in bytes.
        :return: A Dataframe of the size of all files on the disk (in bytes).
        '''

        # There is likely a Python-specific way of doing this but I haven't found a simple solution
        # This thread (https://stackoverflow.com/questions/1392413/calculating-a-directorys-size-using-python/1392549)
        # has a lot of different implementations but they all have complexities with symlinks, subdirectories,
        # broken links, etc.
        logger.warning(
            'A reminder to myself to look for a Python solution for directory sizes (instead of running `du`)')

        sizes_list = []
        for dname in self.ALL_SUBDIRECTORIES:
            dpath = self._check_make_dir(dname)
            # Line from https://stackoverflow.com/a/47930319
            number_of_files = file_count = sum(len(files) for _, _, files in os.walk(dpath))
            size = subprocess.check_output(['du', '-s', '--block-size=1', dpath]).split()[0].decode('utf-8')
            sizes_list.append(['Filesystem',
                               self._root_dir.name,
                               dpath.name,
                               int(size),
                               number_of_files])

        return pd.DataFrame(sizes_list, columns=['Type', 'Name', 'Division', 'Data Size', 'Number of Items'])

    @property
    def root_dir(self):
        return self._root_dir

    @property
    def reference_dir(self):
        return self._check_make_dir('reference')

    @property
    def kmer_dir(self):
        return self._check_make_dir('kmer')

    @property
    def variation_dir(self):
        return self._check_make_dir('variation')

    @property
    def mlst_dir(self):
        return self._check_make_dir('mlst')
