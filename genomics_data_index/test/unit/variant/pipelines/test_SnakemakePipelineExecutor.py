import pytest
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


def test_create_input_sample_files_reads_only():
    files = [Path('SampleA.fastq'),
             Path('SampleB.fastq.gz'),
             Path('SampleC_1.fastq.gz'),
             Path('SampleC_2.fastq.gz'),
             Path('/tmp', 'SampleD_xyz_R1_001.fastq'),
             Path('/tmp', 'SampleD_xyz_R2_001.fastq')]

    executor = SnakemakePipelineExecutor(Path('/tmp'))

    df = executor.create_input_sample_files(files)
    assert 4 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', pd.NA, 'SampleA.fastq', pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', pd.NA, 'SampleB.fastq.gz', pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC', pd.NA, 'SampleC_1.fastq.gz', 'SampleC_2.fastq.gz'] == df.loc['SampleC'].tolist()
    assert ['SampleD_xyz', pd.NA, '/tmp/SampleD_xyz_R1_001.fastq',
            '/tmp/SampleD_xyz_R2_001.fastq'] == df.loc['SampleD_xyz'].tolist()


def test_create_input_sample_files_reads_and_assemblies():
    files = [Path('SampleA.fastq'),
             Path('SampleB.fastq.gz'),
             Path('SampleC_1.fastq.gz'),
             Path('SampleC_2.fastq.gz'),
             Path('/tmp', 'SampleD_xyz_R1_001.fastq'),
             Path('/tmp', 'SampleD_xyz_R2_001.fastq'),
             Path('..', 'SampleE_xyz.fasta'),
             Path('SampleF.fna.gz')]

    executor = SnakemakePipelineExecutor(Path('/tmp'))

    df = executor.create_input_sample_files(files)
    assert 6 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', pd.NA, 'SampleA.fastq', pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', pd.NA, 'SampleB.fastq.gz', pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC', pd.NA, 'SampleC_1.fastq.gz', 'SampleC_2.fastq.gz'] == df.loc['SampleC'].tolist()
    assert ['SampleD_xyz', pd.NA, '/tmp/SampleD_xyz_R1_001.fastq',
            '/tmp/SampleD_xyz_R2_001.fastq'] == df.loc['SampleD_xyz'].tolist()
    assert ['SampleE_xyz', '../SampleE_xyz.fasta', pd.NA, pd.NA] == df.loc['SampleE_xyz'].tolist()
    assert ['SampleF', 'SampleF.fna.gz', pd.NA, pd.NA] == df.loc['SampleF'].tolist()


def test_create_input_sample_files_reads_and_assemblies_duplicate_sample_names():
    files = [Path('SampleA.fasta'),
             Path('SampleB.fasta'),
             Path('SampleA_1.fastq.gz'),
             Path('SampleA_2.fastq.gz')]

    executor = SnakemakePipelineExecutor(Path('/tmp'))

    with pytest.raises(Exception) as execinfo:
        executor.create_input_sample_files(files)

    assert 'Duplicate sample with name [SampleA]' in str(execinfo.value)


def test_create_input_sample_files_reads_and_assemblies_duplicate_sample_names_2():
    files = [Path('/tmp', 'SampleA_1.fastq.gz'),
             Path('/tmp', 'SampleA_2.fastq.gz'),
             Path('SampleA.fasta'),
             Path('SampleB.fasta')]

    executor = SnakemakePipelineExecutor(Path('/tmp'))

    with pytest.raises(Exception) as execinfo:
        executor.create_input_sample_files(files)

    assert 'Duplicate sample with name [SampleA]' in str(execinfo.value)
