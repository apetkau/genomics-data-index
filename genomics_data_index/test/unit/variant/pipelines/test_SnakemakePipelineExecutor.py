import tempfile
from pathlib import Path

import pandas as pd
import pytest

from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor


def test_create_input_sample_files_assemblies_only():
    files = [Path('SampleA.fasta'),
             Path('SampleB.fasta'),
             Path('SampleC_xyz.fa'),
             Path('SampleD.something.fna.gz'),
             Path('/tmp', 'SampleE.fasta.gz')]

    executor = SnakemakePipelineExecutor()

    df = executor.create_input_sample_files(files)
    assert 5 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', Path('SampleA.fasta'), pd.NA, pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', Path('SampleB.fasta'), pd.NA, pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC_xyz', Path('SampleC_xyz.fa'), pd.NA, pd.NA] == df.loc['SampleC_xyz'].tolist()
    assert ['SampleD.something', Path('SampleD.something.fna.gz'), pd.NA, pd.NA] == df.loc['SampleD.something'].tolist()
    assert ['SampleE', Path('/tmp/SampleE.fasta.gz'), pd.NA, pd.NA] == df.loc['SampleE'].tolist()


def test_create_input_sample_files_reads_only():
    files = [Path('SampleA.fastq'),
             Path('SampleB.fastq.gz'),
             Path('SampleC_1.fastq.gz'),
             Path('SampleC_2.fastq.gz'),
             Path('/tmp', 'SampleD_xyz_R1_001.fastq'),
             Path('/tmp', 'SampleD_xyz_R2_001.fastq')]

    executor = SnakemakePipelineExecutor()

    df = executor.create_input_sample_files(files)
    assert 4 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', pd.NA, Path('SampleA.fastq'), pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', pd.NA, Path('SampleB.fastq.gz'), pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC', pd.NA, Path('SampleC_1.fastq.gz'), Path('SampleC_2.fastq.gz')] == df.loc['SampleC'].tolist()
    assert ['SampleD_xyz', pd.NA, Path('/tmp/SampleD_xyz_R1_001.fastq'),
            Path('/tmp/SampleD_xyz_R2_001.fastq')] == df.loc['SampleD_xyz'].tolist()


def test_create_input_sample_files_reads_and_assemblies():
    files = [Path('SampleA.fastq'),
             Path('SampleB.fastq.gz'),
             Path('SampleC_1.fastq.gz'),
             Path('SampleC_2.fastq.gz'),
             Path('/tmp', 'SampleD_xyz_R1_001.fastq'),
             Path('/tmp', 'SampleD_xyz_R2_001.fastq'),
             Path('..', 'SampleE_xyz.fasta'),
             Path('SampleF.fna.gz')]

    executor = SnakemakePipelineExecutor()

    df = executor.create_input_sample_files(files)
    assert 6 == len(df)
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()

    df = df.set_index('Sample', drop=False)
    assert ['SampleA', pd.NA, Path('SampleA.fastq'), pd.NA] == df.loc['SampleA'].tolist()
    assert ['SampleB', pd.NA, Path('SampleB.fastq.gz'), pd.NA] == df.loc['SampleB'].tolist()
    assert ['SampleC', pd.NA, Path('SampleC_1.fastq.gz'), Path('SampleC_2.fastq.gz')] == df.loc['SampleC'].tolist()
    assert ['SampleD_xyz', pd.NA, Path('/tmp/SampleD_xyz_R1_001.fastq'),
            Path('/tmp/SampleD_xyz_R2_001.fastq')] == df.loc['SampleD_xyz'].tolist()
    assert ['SampleE_xyz', Path('../SampleE_xyz.fasta'), pd.NA, pd.NA] == df.loc['SampleE_xyz'].tolist()
    assert ['SampleF', Path('SampleF.fna.gz'), pd.NA, pd.NA] == df.loc['SampleF'].tolist()


def test_create_input_sample_files_reads_and_assemblies_duplicate_sample_names():
    files = [Path('SampleA.fasta'),
             Path('SampleB.fasta'),
             Path('SampleA_1.fastq.gz'),
             Path('SampleA_2.fastq.gz')]

    executor = SnakemakePipelineExecutor()

    with pytest.raises(Exception) as execinfo:
        executor.create_input_sample_files(files)

    assert 'Duplicate sample with name [SampleA]' in str(execinfo.value)


def test_create_input_sample_files_reads_and_assemblies_duplicate_sample_names_2():
    files = [Path('/tmp', 'SampleA_1.fastq.gz'),
             Path('/tmp', 'SampleA_2.fastq.gz'),
             Path('SampleA.fasta'),
             Path('SampleB.fasta')]

    executor = SnakemakePipelineExecutor()

    with pytest.raises(Exception) as execinfo:
        executor.create_input_sample_files(files)

    assert 'Duplicate sample with name [SampleA]' in str(execinfo.value)


def test_validate_input_sample_files():
    input_samples = pd.DataFrame([
        ['A', Path('file.fasta'), pd.NA, pd.NA],
        ['B', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    executor.validate_input_sample_files(input_samples)


def test_validate_input_sample_files_fail():
    input_samples = pd.DataFrame([
        ['A', 'file.fasta', pd.NA, pd.NA],
        ['B', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    with pytest.raises(Exception) as execinfo:
        executor.validate_input_sample_files(input_samples)

    assert 'column=[Assemblies]' in str(execinfo.value)
    assert 'does not contain Path or NA' in str(execinfo.value)


def test_validate_input_sample_files_fail2():
    input_samples = pd.DataFrame([
        ['A', pd.NA, pd.NA, 'file'],
        ['B', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    with pytest.raises(Exception) as execinfo:
        executor.validate_input_sample_files(input_samples)

    assert 'column=[Reads2]' in str(execinfo.value)
    assert 'does not contain Path or NA' in str(execinfo.value)


def test_write_input_sample_files():
    input_samples = pd.DataFrame([
        ['A', Path('file.fasta'), pd.NA, pd.NA],
        ['B', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir)
        out_file = out_dir / 'output.tsv'

        executor.write_input_sample_files(input_samples, output_file=out_file)

        df = pd.read_csv(out_file, sep='\t').fillna('NA')

        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 2 == len(df)
        assert ['A', 'file.fasta', 'NA', 'NA'] == df.iloc[0].tolist()
        assert ['B', 'file2.fasta', 'NA', 'NA'] == df.iloc[1].tolist()


def test_write_read_input_sample_files():
    input_samples = pd.DataFrame([
        ['A', Path('file.fasta'), pd.NA, pd.NA],
        ['B', pd.NA, Path('file.fastq'), pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir)
        out_file = out_dir / 'output.tsv'

        executor.write_input_sample_files(input_samples, output_file=out_file)

        df = executor.read_input_sample_files(out_file)
        executor.validate_input_sample_files(df)

        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 2 == len(df)
        assert ['A', Path('file.fasta'), pd.NA, pd.NA] == df.iloc[0].tolist()
        assert ['B', pd.NA, Path('file.fastq'), pd.NA] == df.iloc[1].tolist()
