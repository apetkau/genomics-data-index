from pathlib import Path
from typing import Dict


class ExecutorResults:

    def __init__(self, files: Dict[str, Path]):
        self._files = files

    def get_file(self, key: str) -> Path:
        return self._files[key]
