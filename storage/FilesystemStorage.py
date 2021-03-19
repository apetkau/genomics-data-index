from pathlib import Path
import os

class FilesystemStorage:

    def __init__(self, root_dir: Path):
        self._root_dir = root_dir

        if not root_dir.exists():
            os.mkdir(root_dir)

    def _check_make_dir(self, name) -> Path:
        d = self._root_dir / name
        if not d.exists():
            os.mkdir(d)
        return d

    @property
    def reference_dir(self):
        return self._check_make_dir('reference')

    @property
    def kmer_dir(self):
        return self._check_make_dir('kmer')

    @property
    def variation_dir(self):
        return self._check_make_dir('variation')
