import abc
from pathlib import Path
from typing import List

from genomics_data_index.pipelines.ExecutorResults import ExecutorResults


class PipelineExecutor(abc.ABC):
    """
    This class defines an interface for all analysis pipelines which need to execute.
    """

    def __init__(self):
        pass

    @abc.abstractmethod
    def execute(self, input_files: List[Path], reference_file: Path, ncores: int = 1) -> ExecutorResults:
        """
        Executes a pipeline and returns a file which can be used as input for GDI (a file of file names).
        """
        pass
