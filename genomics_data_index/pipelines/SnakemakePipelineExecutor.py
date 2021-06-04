import logging
import math
from os import path
from pathlib import Path
from typing import List

import pandas as pd
import yaml

from genomics_data_index.pipelines.ExecutorResults import ExecutorResults
from genomics_data_index.pipelines.PipelineExecutor import PipelineExecutor
from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)

snakemake_file = Path(path.dirname(__file__), 'assembly_input', 'workflow', 'Snakefile')


class SnakemakePipelineExecutor(PipelineExecutor):

    def __init__(self, working_directory: Path, use_conda: bool = True,
                 include_kmer: bool = True, include_mlst: bool = True,
                 kmer_sizes: List[int] = None, kmer_scaled: int = 1000,
                 snakemake_batch_size: int = 5000):
        super().__init__()
        if kmer_sizes is None:
            kmer_sizes = [31]

        self._working_directory = working_directory
        self._use_conda = use_conda
        self._include_kmer = include_kmer
        self._include_mlst = include_mlst
        self._snakemake_batch_size = snakemake_batch_size
        self._sourmash_params = self._prepare_sourmash_params(kmer_sizes=kmer_sizes, kmer_scaled=kmer_scaled)

    def _prepare_sourmash_params(self, kmer_sizes: List[int], kmer_scaled: int) -> str:
        params = ','.join([f'k={v}' for v in kmer_sizes])
        params += f',scaled={kmer_scaled}'

        return params

    def _sample_name_from_file(self, sample_file: Path) -> str:
        return sample_file.stem

    def _prepare_working_directory(self, reference_file: Path,
                                   input_files: List[Path]) -> Path:
        config_dir = self._working_directory / 'config'
        config_dir.mkdir()

        config_file = config_dir / 'config.yaml'
        samples_file = config_dir / 'samples.tsv'

        logger.debug(f'Writing snakemake config [{config_file}]')
        with open(config_file, 'w') as fh:
            config = {
                'reference': str(reference_file.absolute()),
                'samples': str(samples_file.absolute()),
                'include_mlst': self._include_mlst,
                'include_kmer': self._include_kmer,
                'sourmash_params': self._sourmash_params,
            }
            yaml.dump(config, fh)

        logger.debug(f'Writing samples list [{samples_file}]')

        sample_names = []
        for file in input_files:
            sample_names.append([self._sample_name_from_file(file), str(file.absolute())])

        sample_names_df = pd.DataFrame(sample_names, columns=['Sample', 'File'])
        sample_names_df.to_csv(samples_file, sep='\t', index=False)

        return config_file

    def _apply_use_conda(self, command: List[str]) -> List[str]:
        if self._use_conda:
            logger.debug('Enabling --use-conda for snakemake pipeline')
            return command + ['--use-conda']
        else:
            return command

    def _get_number_batches(self, number_input_files: int) -> int:
        return int(math.ceil(number_input_files / self._snakemake_batch_size))

    def execute(self, input_files: List[Path], reference_file: Path, ncores: int = 1) -> ExecutorResults:
        working_directory = self._working_directory
        logger.debug(f'Preparing working directory [{working_directory}] for snakemake')
        config_file = self._prepare_working_directory(reference_file=reference_file,
                                                      input_files=input_files)

        logger.debug(f'Executing snakemake on {len(input_files)} files with reference_file=[{reference_file}]'
                     f' using {ncores} cores in [{working_directory}]')
        snakemake_output_fofn = working_directory / 'gdi-input.fofn'
        snakemake_output_mlst = working_directory / 'mlst.tsv'

        command = ['snakemake', '--configfile', str(config_file)]
        command = self._apply_use_conda(command)
        command = command + ['-j', str(ncores), '--directory', str(working_directory),
                             '--snakefile', str(snakemake_file)]

        # I found with very large numbers of input files snakemake would slow down significantly
        # So I'm executing batches of files in snakemake if you pass too many input files.
        number_batches = self._get_number_batches(len(input_files))
        if number_batches > 1:
            for batch_number in range(1, number_batches + 1):
                logger.info(f'Running snakemake batch: {batch_number}/{number_batches}')
                batch_command_gdi = command + ['--batch', f'gdi_input_fofn={batch_number}/{number_batches}']
                execute_commands([batch_command_gdi])

                # Also batch MLST results
                if self._include_mlst:
                    batch_command_mlst = command + ['--batch', f'mlst_full_table={batch_number}/{number_batches}']
                    execute_commands([batch_command_mlst])

        logger.info('Running Snakemake for rule all')
        execute_commands([command])
        logger.info('Finished running snakemake.')
        logger.debug(f'Output file [{snakemake_output_fofn}]. '
                     f'MLST file [{snakemake_output_mlst}]')

        return ExecutorResults({'gdi-fofn': snakemake_output_fofn, 'mlst': snakemake_output_mlst})
