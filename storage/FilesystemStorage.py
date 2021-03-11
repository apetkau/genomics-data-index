from pathlib import Path
import os

class FilesystemStorage:

    def __init__(self, root_dir: Path):
        self._root_dir = root_dir

    @property
    def reference_dir(self):
        ref_dir = self._root_dir / 'reference'
        if not ref_dir.exists():
            os.mkdir(ref_dir)
        return ref_dir

    @property
    def kmer_dir(self):
        k_dir = self._root_dir / 'kmer'
        if not k_dir.exists():
            os.mkdir(k_dir)
        return k_dir
