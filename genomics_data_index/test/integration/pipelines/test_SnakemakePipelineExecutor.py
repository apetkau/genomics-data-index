from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

import pandas as pd
import vcf
from Bio.SeqRecord import SeqRecord

from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration.pipelines import assemblies_samples, assemblies_reference, expected_mutations
from genomics_data_index.test.integration.pipelines import snpeff_input_sampleA, snpeff_reference_genome
from genomics_data_index.test.integration.pipelines import snpeff_reads_paired, snpeff_reads_single


def vcf_to_mutations_list(vcf_file: Path) -> List[str]:
    reader = vcf.Reader(filename=str(vcf_file))
    df = pd.DataFrame([vars(r) for r in reader])
    df['ALT'] = df['ALT'].apply(lambda x: x[0])
    mutations = df.apply(lambda x: f'{x["CHROM"]}:{x["POS"]}:{x["REF"]}:{x["ALT"]}', axis='columns')

    return mutations.tolist()


def read_expected_mutations(file: Path) -> List[str]:
    with open(file, 'r') as fh:
        mutations = [l.strip() for l in fh.readlines()]
        return mutations


def get_consensus_sequences(file: Path) -> List[SeqRecord]:
    ref_name, sequences = SequenceFile(file).parse_sequence_file()
    return sequences


def assert_vcf(vcf_file: Path, expected_mutations_file: Path):
    sample_expected_mutations = read_expected_mutations(expected_mutations_file)
    actual_mutations = vcf_to_mutations_list(vcf_file)
    assert len(actual_mutations) == len(sample_expected_mutations)
    assert actual_mutations == sample_expected_mutations


def assert_consensus(consensus_file: Path, expected_length: int, expected_Ns: int, expected_gaps: int):
    actual_consensus_records = get_consensus_sequences(consensus_file)
    assert 1 == len(actual_consensus_records)

    actual_consensus_record = actual_consensus_records[0]
    assert expected_length == len(actual_consensus_record)
    assert 'reference' == actual_consensus_record.id
    assert expected_Ns == actual_consensus_record.upper().seq.count('N')
    assert expected_gaps == actual_consensus_record.seq.count('-')


def test_create_fofn_file_single_sample():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'assemblies' / 'variant' / 'SampleA.vcf.gz'
        actual_consensus_file = tmp_dir / 'assemblies' / 'consensus' / 'SampleA.fasta.gz'
        unexpected_mlst_file = tmp_dir / 'mlst.tsv'
        input_samples = [assemblies_samples['SampleA']]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_mutations_file.exists()
        assert actual_consensus_file.exists()
        assert not unexpected_mlst_file.exists()  # Make sure MLST file does not exist in this case

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df

        assert_vcf(actual_mutations_file, expected_mutations['SampleA'])
        assert_consensus(actual_consensus_file, expected_length=5180, expected_Ns=0, expected_gaps=0)


def test_create_fofn_file_multiple_samples():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples = ['SampleA', 'SampleB', 'SampleC']

        input_samples = [assemblies_samples[s] for s in samples]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 3 == len(fofn_df)
        assert ['SampleA', 'SampleB', 'SampleC'] == fofn_df['Sample'].tolist()

        actual_mutations_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_A, expected_mutations['SampleA'])
        assert_consensus(actual_consensus_A, expected_length=5180, expected_Ns=0, expected_gaps=0)

        actual_mutations_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['VCF'].tolist()[0])
        actual_consensus_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_B, expected_mutations['SampleB'])
        assert_consensus(actual_consensus_B, expected_length=5180, expected_Ns=0, expected_gaps=0)

        actual_mutations_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['VCF'].tolist()[0])
        actual_consensus_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_C, expected_mutations['SampleC'])
        assert_consensus(actual_consensus_C, expected_length=5180, expected_Ns=0, expected_gaps=0)


def test_create_fofn_file_multiple_samples_batching():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples = ['SampleA', 'SampleB', 'SampleC']

        input_samples = [assemblies_samples[s] for s in samples]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False, snakemake_input_batch_size=2)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 3 == len(fofn_df)
        assert ['SampleA', 'SampleB', 'SampleC'] == fofn_df['Sample'].tolist()

        actual_mutations_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_A, expected_mutations['SampleA'])
        assert_consensus(actual_consensus_A, expected_length=5180, expected_Ns=0, expected_gaps=0)

        actual_mutations_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['VCF'].tolist()[0])
        actual_consensus_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_B, expected_mutations['SampleB'])
        assert_consensus(actual_consensus_B, expected_length=5180, expected_Ns=0, expected_gaps=0)

        actual_mutations_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['VCF'].tolist()[0])
        actual_consensus_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_C, expected_mutations['SampleC'])
        assert_consensus(actual_consensus_C, expected_length=5180, expected_Ns=0, expected_gaps=0)


def test_create_fofn_file_multiple_samples_with_ns():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples = ['SampleD', 'SampleE', 'SampleF']

        input_samples = [assemblies_samples[s] for s in samples]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 3 == len(fofn_df)
        assert ['SampleD', 'SampleE', 'SampleF'] == fofn_df['Sample'].tolist()

        actual_mutations_D = Path(fofn_df[fofn_df['Sample'] == 'SampleD']['VCF'].tolist()[0])
        actual_consensus_D = Path(fofn_df[fofn_df['Sample'] == 'SampleD']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_D, expected_mutations['SampleD'])
        assert_consensus(actual_consensus_D, expected_length=5180, expected_Ns=480, expected_gaps=0)

        actual_mutations_E = Path(fofn_df[fofn_df['Sample'] == 'SampleE']['VCF'].tolist()[0])
        actual_consensus_E = Path(fofn_df[fofn_df['Sample'] == 'SampleE']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_E, expected_mutations['SampleE'])
        assert_consensus(actual_consensus_E, expected_length=5180, expected_Ns=960, expected_gaps=0)

        actual_mutations_F = Path(fofn_df[fofn_df['Sample'] == 'SampleF']['VCF'].tolist()[0])
        actual_consensus_F = Path(fofn_df[fofn_df['Sample'] == 'SampleF']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_F, expected_mutations['SampleF'])
        assert_consensus(actual_consensus_F, expected_length=5180, expected_Ns=740, expected_gaps=0)


def test_create_fofn_file_multiple_samples_multiple_cores_and_use_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples = ['SampleA', 'SampleB', 'SampleC', 'SampleD', 'SampleE', 'SampleF']

        input_samples = [assemblies_samples[s] for s in samples]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=True,
                                                      include_mlst=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=2)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 6 == len(fofn_df)
        assert ['SampleA', 'SampleB', 'SampleC', 'SampleD', 'SampleE', 'SampleF'] == fofn_df['Sample'].tolist()


def test_create_fofn_file_single_sketch_mlst():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'assemblies' / 'variant' / 'SampleA.vcf.gz'
        actual_consensus_file = tmp_dir / 'assemblies' / 'consensus' / 'SampleA.fasta.gz'
        actual_sketch_file = tmp_dir / 'assemblies' / 'sketch' / 'SampleA.sig.gz'
        actual_mlst_file = tmp_dir / 'mlst.tsv'
        input_samples = [assemblies_samples['SampleA']]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=True)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')
        mlst_file = results.get_file('mlst')

        assert input_fofn.exists()
        assert actual_mutations_file.exists()
        assert actual_consensus_file.exists()
        assert actual_sketch_file.exists()
        assert actual_mlst_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        actual_sketch_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Sketch File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df
        assert actual_sketch_file == actual_sketch_file_from_df
        assert actual_mlst_file == mlst_file


def test_create_fofn_file_single_no_sketch_with_mlst():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'assemblies' / 'variant' / 'SampleA.vcf.gz'
        actual_consensus_file = tmp_dir / 'assemblies' / 'consensus' / 'SampleA.fasta.gz'
        unexpected_sketch_file = tmp_dir / 'assemblies' / 'sketch' / 'SampleA.sig.gz'
        actual_mlst_file = tmp_dir / 'mlst.tsv'
        input_samples = [assemblies_samples['SampleA']]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=True, include_kmer=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=assemblies_reference,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')
        mlst_file = results.get_file('mlst')

        assert input_fofn.exists()
        assert actual_mutations_file.exists()
        assert actual_consensus_file.exists()
        assert not unexpected_sketch_file.exists()
        assert actual_mlst_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df
        assert actual_mlst_file == mlst_file


def test_create_fofn_file_snpeff_no_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file = tmp_dir / 'assemblies' / 'variant-snpeff' / 'SampleA.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_consensus_file = tmp_dir / 'assemblies' / 'consensus' / 'SampleA.fasta.gz'
        input_samples = [snpeff_input_sampleA]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False, include_kmer=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file.exists()
        assert actual_consensus_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df

        # snpeff annotations should be added in headers
        reader = vcf.Reader(filename=str(actual_mutations_snpeff_file))
        assert 'ANN' in reader.infos


def test_create_fofn_file_snpeff_with_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file = tmp_dir / 'assemblies' / 'variant-snpeff' / 'SampleA.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_consensus_file = tmp_dir / 'assemblies' / 'consensus' / 'SampleA.fasta.gz'
        input_samples = [snpeff_input_sampleA]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=True,
                                                      include_mlst=False, include_kmer=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file.exists()
        assert actual_consensus_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df

        # snpeff annotations should be added in headers
        reader = vcf.Reader(filename=str(actual_mutations_snpeff_file))
        assert 'ANN' in reader.infos


def test_create_fofn_file_snpeff_reads_with_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file_paired = tmp_dir / 'reads' / 'paired' / 'variant-snpeff' / 'SampleA-snpeff.vcf.gz'
        actual_mutations_snpeff_file_single = tmp_dir / 'reads' / 'single' / 'variant-snpeff' / 'SampleA-single-snpeff.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_consensus_file_paired = tmp_dir / 'reads' / 'paired' / 'consensus' / 'SampleA-snpeff.fasta.gz'
        actual_consensus_file_single = tmp_dir / 'reads' / 'single' / 'consensus' / 'SampleA-single-snpeff.fasta.gz'
        input_samples = snpeff_reads_paired + snpeff_reads_single

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=True,
                                                      include_mlst=False, include_kmer=False,
                                                      reads_mincov=5)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file_paired.exists()
        assert actual_mutations_snpeff_file_single.exists()
        assert actual_consensus_file_paired.exists()
        assert actual_consensus_file_single.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 2 == len(fofn_df)
        assert {'SampleA-snpeff', 'SampleA-single-snpeff'} == set(fofn_df['Sample'].tolist())
        actual_mutations_file_paired_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA-snpeff']['VCF'].tolist()[0])
        actual_mutations_file_single_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-single-snpeff']['VCF'].tolist()[0])
        actual_consensus_file_paired_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-snpeff']['Mask File'].tolist()[0])
        actual_consensus_file_single_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-single-snpeff']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file_paired == actual_mutations_file_paired_from_df
        assert actual_mutations_snpeff_file_single == actual_mutations_file_single_from_df
        assert actual_consensus_file_paired == actual_consensus_file_paired_from_df
        assert actual_consensus_file_single == actual_consensus_file_single_from_df

        # snpeff annotations should be added in headers
        reader = vcf.Reader(filename=str(actual_mutations_snpeff_file_paired))
        assert 'ANN' in reader.infos

        reader = vcf.Reader(filename=str(actual_mutations_snpeff_file_single))
        assert 'ANN' in reader.infos
