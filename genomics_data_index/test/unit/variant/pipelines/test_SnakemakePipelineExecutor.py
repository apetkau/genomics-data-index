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


def test_skip_samples_from_input_files():
    input_samples = pd.DataFrame([
        ['A', Path('file.fasta'), pd.NA, pd.NA],
        ['B', Path('file2.fasta'), pd.NA, pd.NA],
        ['C', pd.NA, Path('file_1.fastq'), Path('file_2.fastq')]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples=None)
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 3 == len(df)
    assert ['A', Path('file.fasta'), 'NA', 'NA'] == df.iloc[0].tolist()
    assert ['B', Path('file2.fasta'), 'NA', 'NA'] == df.iloc[1].tolist()
    assert ['C', 'NA', Path('file_1.fastq'), Path('file_2.fastq')] == df.iloc[2].tolist()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples=set())
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 3 == len(df)
    assert ['A', Path('file.fasta'), 'NA', 'NA'] == df.iloc[0].tolist()
    assert ['B', Path('file2.fasta'), 'NA', 'NA'] == df.iloc[1].tolist()
    assert ['C', 'NA', Path('file_1.fastq'), Path('file_2.fastq')] == df.iloc[2].tolist()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples={'A'})
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 2 == len(df)
    assert ['B', Path('file2.fasta'), 'NA', 'NA'] == df.iloc[0].tolist()
    assert ['C', 'NA', Path('file_1.fastq'), Path('file_2.fastq')] == df.iloc[1].tolist()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples=set('B'))
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 2 == len(df)
    assert ['A', Path('file.fasta'), 'NA', 'NA'] == df.iloc[0].tolist()
    assert ['C', 'NA', Path('file_1.fastq'), Path('file_2.fastq')] == df.iloc[1].tolist()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples={'A', 'C', 'invalid'})
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 1 == len(df)
    assert ['B', Path('file2.fasta'), 'NA', 'NA'] == df.iloc[0].tolist()

    df = executor.skip_samples_from_input_files(input_samples, skip_samples={'A', 'B', 'C'})
    df = df.fillna('NA')
    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
    assert 0 == len(df)


def test_skip_missing_sample_files():
    with tempfile.NamedTemporaryFile() as f1_tmp:
        f1 = Path(f1_tmp.name)

        executor = SnakemakePipelineExecutor()

        input_samples = pd.DataFrame([
            ['A', f1, pd.NA, pd.NA],
            ['B', f1, pd.NA, pd.NA],
            ['C', pd.NA, f1, f1]
        ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

        df = executor.skip_missing_sample_files(input_samples)
        df = df.fillna('NA')
        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 3 == len(df)
        assert ['A', f1, 'NA', 'NA'] == df.iloc[0].tolist()
        assert ['B', f1, 'NA', 'NA'] == df.iloc[1].tolist()
        assert ['C', 'NA', f1, f1] == df.iloc[2].tolist()


        input_samples = pd.DataFrame([
            ['A', f1, pd.NA, pd.NA],
            ['B', Path('file2.fasta'), pd.NA, pd.NA],
            ['C', pd.NA, Path('file_1.fastq'), Path('file_2.fastq')]
        ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

        df = executor.skip_missing_sample_files(input_samples)
        df = df.fillna('NA')
        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 1 == len(df)
        assert ['A', f1, 'NA', 'NA'] == df.iloc[0].tolist()


        input_samples = pd.DataFrame([
            ['A', f1, pd.NA, pd.NA],
            ['B', f1, pd.NA, pd.NA],
            ['C', pd.NA, Path('file_1.fastq'), Path('file_2.fastq')]
        ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

        df = executor.skip_missing_sample_files(input_samples)
        df = df.fillna('NA')
        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 2 == len(df)
        assert ['A', f1, 'NA', 'NA'] == df.iloc[0].tolist()
        assert ['B', f1, 'NA', 'NA'] == df.iloc[1].tolist()


        input_samples = pd.DataFrame([
            ['A', f1, pd.NA, pd.NA],
            ['B', f1, pd.NA, pd.NA],
            ['C', pd.NA, f1, Path('file_2.fastq')]
        ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

        df = executor.skip_missing_sample_files(input_samples)
        df = df.fillna('NA')
        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 2 == len(df)
        assert ['A', f1, 'NA', 'NA'] == df.iloc[0].tolist()
        assert ['B', f1, 'NA', 'NA'] == df.iloc[1].tolist()


        input_samples = pd.DataFrame([
            ['A', f1, pd.NA, pd.NA],
            ['B', Path('file2.fasta'), pd.NA, pd.NA],
            ['C', pd.NA, f1, f1]
        ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

        df = executor.skip_missing_sample_files(input_samples)
        df = df.fillna('NA')
        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == df.columns.tolist()
        assert 2 == len(df)
        assert ['A', f1, 'NA', 'NA'] == df.iloc[0].tolist()
        assert ['C', 'NA', f1, f1] == df.iloc[1].tolist()


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


def test_fix_sample_names_no_change():
    input_samples = pd.DataFrame([
        ['A', Path('file1.fasta'), pd.NA, pd.NA],
        ['B', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    input_samples_fixed, original_fixed_names = executor.fix_sample_names(input_samples)
    assert original_fixed_names is None

    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == input_samples_fixed.columns.tolist()
    assert 2 == len(input_samples_fixed)
    assert ['A', 'B'] == input_samples_fixed['Sample'].tolist()
    assert [Path('file1.fasta'), Path('file2.fasta')] == input_samples_fixed['Assemblies'].tolist()


def test_fix_sample_names_with_change():
    input_samples = pd.DataFrame([
        ['A/1', Path('file1.fasta'), pd.NA, pd.NA],
        ['B2', Path('file2.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    input_samples_fixed, original_fixed_names = executor.fix_sample_names(input_samples)
    assert original_fixed_names is not None

    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == input_samples_fixed.columns.tolist()
    assert 2 == len(input_samples_fixed)
    assert ['A__1', 'B2'] == input_samples_fixed['Sample'].tolist()
    assert [Path('file1.fasta'), Path('file2.fasta')] == input_samples_fixed['Assemblies'].tolist()

    assert ['Sample_original', 'Sample_fixed'] == original_fixed_names.columns.tolist()
    assert 2 == len(original_fixed_names)
    assert ['A/1', 'B2'] == original_fixed_names['Sample_original'].tolist()
    assert ['A__1', 'B2'] == original_fixed_names['Sample_fixed'].tolist()


def test_fix_sample_names_with_change2():
    input_samples = pd.DataFrame([
        ['new/A/1/', Path('file1.fasta'), pd.NA, pd.NA],
        ['new/B/2', Path('file2.fasta'), pd.NA, pd.NA],
        ['new/C//2', Path('file3.fasta'), pd.NA, pd.NA]
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    input_samples_fixed, original_fixed_names = executor.fix_sample_names(input_samples)

    assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == input_samples_fixed.columns.tolist()
    assert 3 == len(input_samples_fixed)
    assert ['new__A__1__', 'new__B__2', 'new__C____2'] == input_samples_fixed['Sample'].tolist()
    assert [Path('file1.fasta'), Path('file2.fasta'), Path('file3.fasta')] == input_samples_fixed['Assemblies'].tolist()

    assert ['Sample_original', 'Sample_fixed'] == original_fixed_names.columns.tolist()
    assert 3 == len(original_fixed_names)
    assert ['new/A/1/', 'new/B/2', 'new/C//2'] == original_fixed_names['Sample_original'].tolist()
    assert ['new__A__1__', 'new__B__2', 'new__C____2'] == original_fixed_names['Sample_fixed'].tolist()


def test_fix_sample_names_with_duplicates():
    input_samples = pd.DataFrame([
        ['A/1', Path('file1.fasta'), pd.NA, pd.NA],
        ['A__1', Path('file2.fasta'), pd.NA, pd.NA],
    ], columns=['Sample', 'Assemblies', 'Reads1', 'Reads2'])

    executor = SnakemakePipelineExecutor()

    with pytest.raises(Exception) as execinfo:
        executor.fix_sample_names(input_samples)

    assert "Sanitizing sample names to use as file names leads to duplicates" in str(execinfo.value)


def test_restore_sample_names():
    original_new_names = pd.DataFrame([
        ['new/A/1/', 'new__A__1__'],
        ['new/B/2', 'new__B__2'],
        ['new/C//2', 'new__C____2'],
    ], columns=['Sample_original', 'Sample_fixed'])

    data = pd.DataFrame([
        ['new__A__1__', Path('file1.vcf'), Path('file1.fasta')],
        ['new__B__2', Path('file2.vcf'), Path('file2.fasta')],
        ['new__C____2', Path('file3.vcf'), Path('file3.fasta')]
    ], columns=['Sample', 'VCF', 'Mask'])

    executor = SnakemakePipelineExecutor()

    actual_df = executor.restore_sample_names(data, sample_column='Sample', samples_original_fixed=original_new_names)

    assert ['Sample', 'VCF', 'Mask'] == actual_df.columns.tolist()
    assert 3 == len(actual_df)
    assert ['new/A/1/', 'new/B/2', 'new/C//2'] == actual_df['Sample'].tolist()
    assert [Path('file1.vcf'), Path('file2.vcf'), Path('file3.vcf')] == actual_df['VCF'].tolist()
    assert [Path('file1.fasta'), Path('file2.fasta'), Path('file3.fasta')] == actual_df['Mask'].tolist()
