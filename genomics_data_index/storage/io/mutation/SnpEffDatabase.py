from pathlib import Path


class SnpEffDatabase:

    def __init__(self, snpeff_config: Path, database_dir: Path, genome_name: str):
        self._snpeff_config = snpeff_config
        self._database_dir = database_dir
        self._genome_name = genome_name

    @property
    def config(self) -> Path:
        return self._snpeff_config

    @property
    def genome_name(self) -> str:
        return self._genome_name
