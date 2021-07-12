import abc
from pathlib import Path
from typing import List
import pandas as pd

from genomics_data_index.pipelines.ExecutorResults import ExecutorResults
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile


class PipelineExecutor(abc.ABC):
    """
    This class defines an interface for all analysis pipelines which need to execute.
    """

    INPUT_SAMPLE_FILE_COLUMNS = ['Sample', 'Assemblies', 'Reads1', 'Reads2']

    def __init__(self):
        pass

    @abc.abstractmethod
    def execute(self, input_files: List[Path], reference_file: Path, ncores: int = 1) -> ExecutorResults:
        """
        Executes a pipeline and returns a file which can be used as input for GDI (a file of file names).
        """
        pass

    def validate_input_sample_files(self, input_sample_files: pd.DataFrame) -> None:
        if input_sample_files.columns.tolist() != self.INPUT_SAMPLE_FILE_COLUMNS:
            raise Exception(f'input_sample_files={input_sample_files} has invalid columns. '
                            f'Expected {self.INPUT_SAMPLE_FILE_COLUMNS}')

    def _sample_name_from_file(self, sample_file: Path) -> str:
        return SequenceFile(sample_file).get_genome_name()

    def create_input_sample_files(self, input_files: List[Path]) -> pd.DataFrame:
        sample_files = []
        for file in input_files:
            sample_files.append([self._sample_name_from_file(file), str(file.absolute()), pd.NA, pd.NA])

        sample_files_df = pd.DataFrame(sample_files, columns=self.INPUT_SAMPLE_FILE_COLUMNS)

        return sample_files_df

    def write_input_sample_files(self, input_sample_files: pd.DataFrame, output_file: Path) -> None:
        input_sample_files.to_csv(output_file, sep='\t', index=False)
