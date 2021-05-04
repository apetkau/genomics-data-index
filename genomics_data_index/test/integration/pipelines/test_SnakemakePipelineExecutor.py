from pathlib import Path
from tempfile import TemporaryDirectory

from genomics_data_index.test.integration.pipelines import assemblies_samples, assemblies_reference
from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor


def test_create_fofn_file():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir)

        input_fofn = pipeline_executor.execute(input_files=assemblies_samples,
                                               reference_file=assemblies_reference,
                                               ncores=1)

        assert input_fofn.exists()
