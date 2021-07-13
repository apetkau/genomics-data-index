from pathlib import Path
import pandas as pd


from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor


def test_create_input_sample_files_assemblies_only():
    files = [Path('SampleA.fasta'),
             Path('SampleB.fasta'),
             Path('SampleC_xyz.fa'),
             Path('SampleD.something.fna.gz'),
             Path('/tmp', 'SampleE.fasta.gz')]

    executor = SnakemakePipelineExecutor(Path('/tmp'))

    df = executor.create_input_sample_files(files)
    assert 5 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', 'SampleA.fasta', pd.NA, pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', 'SampleB.fasta', pd.NA, pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC_xyz', 'SampleC_xyz.fa', pd.NA, pd.NA] == df.loc['SampleC_xyz'].tolist()
    assert ['SampleD.something', 'SampleD.something.fna.gz', pd.NA, pd.NA] == df.loc['SampleD.something'].tolist()
    assert ['SampleE', '/tmp/SampleE.fasta.gz', pd.NA, pd.NA] == df.loc['SampleE'].tolist()
