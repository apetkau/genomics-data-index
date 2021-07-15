import abc
import io
import logging
import sys
from pathlib import Path
from typing import List, Union

import pandas as pd

from genomics_data_index.pipelines.ExecutorResults import ExecutorResults
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile

logger = logging.getLogger(__name__)


class PipelineExecutor(abc.ABC):
    """
    This class defines an interface for all analysis pipelines which need to execute.
    """

    INPUT_SAMPLE_FILE_COLUMNS = ['Sample', 'Assemblies', 'Reads1', 'Reads2']

    def __init__(self):
        pass

    @abc.abstractmethod
    def execute(self, sample_files: pd.DataFrame, reference_file: Path, ncores: int = 1) -> ExecutorResults:
        """
        Executes a pipeline and returns a file which can be used as input for GDI (a file of file names).
        """
        pass

    def validate_input_sample_files(self, input_sample_files: pd.DataFrame) -> None:
        if input_sample_files.columns.tolist() != self.INPUT_SAMPLE_FILE_COLUMNS:
            raise Exception(f'input_sample_files={input_sample_files} has invalid columns. '
                            f'Expected {self.INPUT_SAMPLE_FILE_COLUMNS}')

        for col in ['Assemblies', 'Reads1', 'Reads2']:
            are_paths = set(input_sample_files[col].apply(lambda x: isinstance(x, Path) or pd.isna(x)).tolist())
            if not all(are_paths):
                raise Exception(
                    f'column=[{col}] in input_sample_files={input_sample_files} does not contain Path or NA')

    def create_input_sample_files(self, input_files: List[Path]) -> pd.DataFrame:
        """

        :rtype: object
        """
        assemblies = {}
        reads = {}
        sample_names = set()
        data = []

        # Initial pass of files to break up into assemblies/reads
        for file in input_files:
            sf = SequenceFile(file)
            sample_name = sf.get_genome_name(exclude_paired_end_indicators=True)
            if sf.is_assembly():
                if sample_name in sample_names:
                    if sample_name in assemblies:
                        previous_files = [assemblies[sample_name]]
                    else:
                        previous_files = reads[sample_name]
                    raise Exception(f'Duplicate sample with name [{sample_name}]. current_file=[{file}], '
                                    f'previous_file(s)={previous_files}')
                else:
                    sample_names.add(sample_name)
                    assemblies[sample_name] = file
            elif sf.is_reads():
                if sample_name in assemblies:
                    previous_files = assemblies[sample_name]
                    raise Exception(f'Duplicate sample with name [{sample_name}]. current_file=[{file}], '
                                    f'previous_file(s)={previous_files}')
                elif sample_name in reads:
                    if len(reads[sample_name]) != 1:
                        raise Exception(f'Invalid number of files for sample with name [{sample_name}]. '
                                        f'current_file=[{file}], previous_files={reads[sample_name]}')
                    else:
                        reads[sample_name].append(file)
                else:
                    reads[sample_name] = [file]

                sample_names.add(sample_name)
            else:
                logger.warning(f'Input file [{file}] with unknown file type (not assembly or reads). Ignoring.')

        # Now we iterate over samples to insert into an array to create the final dataframe
        for sample in assemblies:
            data.append([sample, assemblies[sample], pd.NA, pd.NA])

        # Iterate over reads to insert into array for final dataframe
        for sample in reads:
            if len(reads[sample]) == 1:
                data.append([sample, pd.NA, reads[sample][0], pd.NA])
            elif len(reads[sample]) == 2:
                file1 = SequenceFile(reads[sample][0])
                file2 = SequenceFile(reads[sample][1])

                file1_differences = file1.name_differences(file2)
                file2_differences = file2.name_differences(file1)

                if len(file1_differences) != 1 or len(file2_differences) != 1:
                    raise Exception(
                        f'Files [{reads[sample]}] do not have exactly one difference between names, cannot determine'
                        f' paired structure.')
                else:
                    f1d = file1_differences[0].lower()
                    f2d = file2_differences[0].lower()

                    if f1d == '1' and f2d == '2':
                        forward = file1
                        reverse = file2
                    elif f1d == 'f' and f2d == 'r':
                        forward = file1
                        reverse = file2
                    elif f2d == '1' and f1d == '2':
                        reverse = file1
                        forward = file2
                    elif f1d == 'r' and f2d == 'f':
                        reverse = file1
                        forward = file2
                    else:
                        raise Exception(f'Cannot determine pair structure for files [{reads[sample]}]')

                    data.append([sample, pd.NA, forward.file, reverse.file])
            else:
                raise Exception(f'Invalid number of files for sample [{sample}], files={reads[sample]}')

        return pd.DataFrame(data, columns=self.INPUT_SAMPLE_FILE_COLUMNS)

    def write_input_sample_files(self, input_sample_files: pd.DataFrame,
                                 output_file: Union[Path, io.TextIOWrapper] = sys.stdout,
                                 abolute_paths: bool = False) -> None:
        input_sample_files = input_sample_files.copy()
        # Convert paths to strings
        for col in ['Assemblies', 'Reads1', 'Reads2']:
            if abolute_paths:
                input_sample_files[col] = input_sample_files[col].apply(
                    lambda x: str(x.absolute()) if not pd.isna(x) else pd.NA)
            else:
                input_sample_files[col] = input_sample_files[col].apply(lambda x: str(x) if not pd.isna(x) else pd.NA)
        input_sample_files.to_csv(output_file, sep='\t', index=False)

    def read_input_sample_files(self, input_file: Path) -> pd.DataFrame:
        df = pd.read_csv(input_file, sep='\t')
        for col in ['Assemblies', 'Reads1', 'Reads2']:
            df[col] = df[col].apply(lambda x: Path(x) if not pd.isna(x) else pd.NA)

        return df
