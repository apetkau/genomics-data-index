import abc
import gzip
import io
import logging
import random
import sys
from pathlib import Path
from typing import List, Union, Tuple, Optional, Set

import pandas as pd
from Bio import SeqIO
from pathvalidate import sanitize_filename

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

    def fix_sample_names(self, sample_files: pd.DataFrame) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        fixed_sample_files = sample_files.copy()

        fixed_sample_files['Sample_fixed'] = fixed_sample_files['Sample'].apply(
            lambda x: sanitize_filename(x, replacement_text='__'))

        samples_changed = (fixed_sample_files['Sample'] != fixed_sample_files['Sample_fixed'])
        if not samples_changed.any():
            return sample_files, None
        else:
            logger.debug(f'Modified {samples_changed.sum()} sample names so they are valid file names for '
                         f'analysis. These will be restored when the analysis is complete.')
            agg_samples_fixed = fixed_sample_files.groupby('Sample_fixed').agg({
                'Sample_fixed': 'count',
                'Sample': 'first'
            })
            duplicate_names_fixed = agg_samples_fixed[agg_samples_fixed['Sample_fixed'] > 1]
            if len(duplicate_names_fixed) > 0:
                duplicate_samples = duplicate_names_fixed['Sample'].tolist()
                raise Exception(f'Sanitizing sample names to use as file names leads to duplicates for the '
                                f'following samples: {duplicate_samples}. Please rename these samples and try'
                                f' executing again.')
            else:
                fixed_sample_files = fixed_sample_files.rename({'Sample': 'Sample_original'}, axis='columns')
                sample_sample_fixed_df = fixed_sample_files[['Sample_original', 'Sample_fixed']]
                fixed_sample_files = fixed_sample_files.rename({'Sample_fixed': 'Sample'}, axis='columns')

                return fixed_sample_files[['Sample', 'Assemblies', 'Reads1', 'Reads2']], \
                       sample_sample_fixed_df

    def restore_sample_names(self, data: pd.DataFrame, sample_column: str,
                             samples_original_fixed: pd.DataFrame) -> pd.DataFrame:
        data_columns_to_keep = data.columns.tolist()

        # Rename column so I know which is the original set of sample names after merging
        if sample_column == 'Sample_original':
            original_sample_column = 'Sample_original_renamed'
            samples_original_fixed = samples_original_fixed.rename({'Sample_original': original_sample_column},
                                                                   axis='columns')
        else:
            original_sample_column = 'Sample_original'

        restored_data = data.merge(samples_original_fixed, left_on=sample_column, right_on='Sample_fixed')
        restored_data = restored_data.drop(sample_column, axis='columns')  # Remove sample names I don't want
        restored_data = restored_data.rename({original_sample_column: sample_column}, axis='columns')
        restored_data = restored_data[data_columns_to_keep]
        return restored_data

    def create_input_sample_files(self, input_files: List[Path], skip_samples: Set[str] = None) -> pd.DataFrame:
        """
        Create a dataframe which associates the given files with a sample entry.
        :param input_files: The list of files.
        :param skip_samples: An optional set of sample names to skip.
        :return: A pandas.DataFrame which has one sample per row associated with the input files.
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

        samples_df = pd.DataFrame(data, columns=self.INPUT_SAMPLE_FILE_COLUMNS)
        samples_df = self.skip_samples_from_input_files(samples_df, skip_samples=skip_samples)

        return samples_df

    def skip_samples_from_input_files(self, samples_df: pd.DataFrame, skip_samples: Set[str] = None):
        if not (skip_samples is None or len(skip_samples) == 0):
            full_length = len(samples_df)
            samples_df = samples_df.loc[~samples_df['Sample'].isin(skip_samples)]
            reduced_length = len(samples_df)
            logger.info(f'Skipping {full_length - reduced_length}/{full_length} samples '
                         f'since they are already indexed')
        return samples_df

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

    def select_random_samples(self, input_sequence_files: List[Path],
                              number_samples: float, random_seed: int = None) -> Set[str]:
        """
        Selects a random set of samples from the passed sequence files (assumes one sample is one sequence record).
        :param input_sequence_files: The list of sequence files to select samples from.
        :param number_samples: If >= 1, then represents the number of samples to select (e.g., 5 for 5 samples).
                               If < 1, then represents the proportion of samples in all files to select
                               (e.g., 0.5 means 50% of samples across all files).
        :param random_seed: The random seed. If None uses the default Python random seed
                            (system time or some other method, see Python documentation on random.seed()).
        :return: A randomly selected set of sample names from the passed sequence files (no check is performed
                 for duplicate sequences across files).
        """
        if number_samples < 0:
            raise Exception(f'number_samples={number_samples} must be non-negative')

        samples = []
        for input_sequence_file in input_sequence_files:
            sequence_file = SequenceFile(input_sequence_file)
            for record in sequence_file.records():
                samples.append(record.id)

        if number_samples < 1:
            number_samples = round(number_samples * len(samples))
        else:
            number_samples = round(number_samples)

        random.seed(random_seed)
        return set(random.sample(samples, k=number_samples))

    def split_input_sequence_files(self, input_sequence_files: List[Path], output_dir: Path,
                                   samples: Set[str] = None) -> pd.DataFrame:
        if samples is not None and isinstance(samples, list):
            samples = set(samples)

        sample_data = []
        print_on = 200
        for input_sequence_file in input_sequence_files:
            count = 0
            sequence_file = SequenceFile(input_sequence_file)

            if sequence_file.is_fasta():
                extension = '.fasta.gz'
                filetype = 'fasta'
            elif sequence_file.is_genbank():
                extension = '.gbk.gz'
                filetype = 'genbank'
            else:
                raise Exception(f'Unknown filetype for {input_sequence_file}, must be either fasta or genbank')

            for record in sequence_file.records():
                if count % print_on == 0:
                    logger.debug(f'Processed {count} sequences from {input_sequence_file}')

                sample_name = record.id

                if samples is None or sample_name in samples:
                    sample_filename = sanitize_filename(sample_name, replacement_text='__')
                    sample_path = output_dir / (sample_filename + extension)
                    with gzip.open(sample_path, "wt") as oh:
                        SeqIO.write(record, oh, filetype)

                    sample_data.append([sample_name, sample_path, pd.NA, pd.NA])

                count += 1
        return pd.DataFrame(sample_data, columns=self.INPUT_SAMPLE_FILE_COLUMNS)

    def read_input_sample_files(self, input_file: Path) -> pd.DataFrame:
        df = pd.read_csv(input_file, sep='\t')
        for col in ['Assemblies', 'Reads1', 'Reads2']:
            df[col] = df[col].apply(lambda x: Path(x) if not pd.isna(x) else pd.NA)

        return df
