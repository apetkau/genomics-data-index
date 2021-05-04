import abc
from typing import List
from pathlib import Path


class PipelineExecutor(abc.ABC):
    """
    This class defines an interface for all analysis pipelines which need to execute.
    """

    def __init__(self):
        pass

    @abc.abstractmethod
    def execute(self, input_files: List[Path], reference_file: Path,
                working_directory: Path, ncores: int = 1) -> Path:
        """
        Executes a pipeline and returns a file which can be used as input for GDI (a file of file names).
        """
        pass
