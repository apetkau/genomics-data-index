from pathlib import Path


class FilesystemStorage:

    def __init__(self, root_dir: Path):
        self._root_dir = root_dir

    @property
    def reference_dir(self):
        return self._root_dir / 'reference'

    @property
    def kmer_dir(self):
        return self._root_dir / 'kmer'
