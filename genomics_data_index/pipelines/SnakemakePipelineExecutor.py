from typing import List
from pathlib import Path
from os import path
import logging

from genomics_data_index.pipelines.PipelineExecutor import PipelineExecutor
from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)


snakemake_file = Path(path.dirname(__file__), 'assembly_input', 'workflow', 'Snakefile')


class SnakemakePipelineExecutor(PipelineExecutor):

    def __init__(self, snakemake_directory: Path):
        super().__init__()
        self._snakemake_directory = snakemake_directory

    def execute(self, input_files: List[Path], reference_file: Path,
                working_directory: Path, ncores: int = 1) -> Path:
        logger.debug(f'Executing snakemake on {len(input_files)} files with reference_file=[{reference_file}]'
                     f' using {ncores} cores in [{working_directory}]')
        snakemake_output = working_directory / 'gdi-input.fofn'
        command = ['snakemake', '--use-conda', '-j', str(ncores), '--directory', str(working_directory),
                   '--snakefile', str(snakemake_file)]
        execute_commands([command])
        logger.debug(f'Finished executing snakemake. Output file [{snakemake_output}]')

        return snakemake_output
