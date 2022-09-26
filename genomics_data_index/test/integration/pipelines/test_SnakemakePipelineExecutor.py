from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Union

import pandas as pd
import vcfpy
from pybedtools import BedTool

from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.test.integration import reference_file, reference_file_5000_snpeff_2
from genomics_data_index.test.integration.pipelines import assemblies_samples, assemblies_reference, expected_mutations
from genomics_data_index.test.integration.pipelines import empty_bed_file, expected_beds
from genomics_data_index.test.integration.pipelines import snpeff_input_sampleA, snpeff_reference_genome
from genomics_data_index.test.integration.pipelines import snpeff_reads_paired, snpeff_reads_single


def vcf_to_mutations_list(vcf_file: Path) -> List[str]:
    reader = vcfpy.Reader.from_path(path=str(vcf_file))
    df = pd.DataFrame([vars(r) for r in reader])
    df['ALT'] = df['ALT'].apply(lambda x: x[0].value)
    mutations = df.apply(lambda x: f'{x["CHROM"]}:{x["POS"]}:{x["REF"]}:{x["ALT"]}', axis='columns')

    return mutations.tolist()


def read_expected_mutations(file: Path) -> List[str]:
    with open(file, 'r') as fh:
        mutations = [l.strip() for l in fh.readlines()]
        return mutations


def assert_vcf(vcf_file: Path, expected_mutations_file: Path):
    sample_expected_mutations = read_expected_mutations(expected_mutations_file)
    actual_mutations = vcf_to_mutations_list(vcf_file)
    assert len(actual_mutations) == len(sample_expected_mutations)
    assert actual_mutations == sample_expected_mutations


def assert_masks_equal(actual_mask_file: Path, expected_bed: Union[Path, BedTool]):
    actual_mask = BedTool(str(actual_mask_file))

    if isinstance(expected_bed, Path):
        expected_bed = BedTool(str(expected_bed))

    assert expected_bed == actual_mask, f'expected=\n{expected_bed} != actual=\n{actual_mask}'


def test_create_fofn_file_single_sample():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'assemblies' / 'variant' / 'SampleA.vcf.gz'
        actual_mask_file = tmp_dir / 'assemblies' / 'mask' / 'SampleA.bed.gz'
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
        assert actual_mask_file.exists()
        assert not unexpected_mlst_file.exists()  # Make sure MLST file does not exist in this case

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_mask_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_mask_file == actual_mask_file_from_df

        assert_vcf(actual_mutations_file, expected_mutations['SampleA'])
        assert_masks_equal(actual_mask_file, empty_bed_file)


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
        actual_mask_file_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_A, expected_mutations['SampleA'])
        assert_masks_equal(actual_mask_file_A, empty_bed_file)

        actual_mutations_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['VCF'].tolist()[0])
        actual_mask_file_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_B, expected_mutations['SampleB'])
        assert_masks_equal(actual_mask_file_B, empty_bed_file)

        actual_mutations_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['VCF'].tolist()[0])
        actual_mask_file_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_C, expected_mutations['SampleC'])
        assert_masks_equal(actual_mask_file_C, empty_bed_file)


def test_create_fofn_file_multiple_samples_with_invalid_charcaters_in_name():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples_original = ['SampleA', 'SampleB', 'SampleC']
        input_samples = [assemblies_samples[s] for s in samples_original]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)

        # Modify sample names so they contain extra characters that can't be used for a file name
        # (in this case, a slash '/')
        sample_files['Sample'] = sample_files['Sample'] + '/extra|1'
        expected_mutations_local = {s + '/extra|1': expected_mutations[s] for s in expected_mutations}

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
        assert ['SampleA/extra|1', 'SampleB/extra|1', 'SampleC/extra|1'] == fofn_df['Sample'].tolist()

        actual_mutations_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA/extra|1']['VCF'].tolist()[0])
        actual_mask_file_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA/extra|1']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_A, expected_mutations_local['SampleA/extra|1'])
        assert_masks_equal(actual_mask_file_A, empty_bed_file)

        actual_mutations_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB/extra|1']['VCF'].tolist()[0])
        actual_mask_file_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB/extra|1']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_B, expected_mutations_local['SampleB/extra|1'])
        assert_masks_equal(actual_mask_file_B, empty_bed_file)

        actual_mutations_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC/extra|1']['VCF'].tolist()[0])
        actual_mask_file_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC/extra|1']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_C, expected_mutations_local['SampleC/extra|1'])
        assert_masks_equal(actual_mask_file_C, empty_bed_file)


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
        actual_mask_file_A = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_A, expected_mutations['SampleA'])
        assert_masks_equal(actual_mask_file_A, empty_bed_file)

        actual_mutations_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['VCF'].tolist()[0])
        actual_mask_file_B = Path(fofn_df[fofn_df['Sample'] == 'SampleB']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_B, expected_mutations['SampleB'])
        assert_masks_equal(actual_mask_file_B, empty_bed_file)

        actual_mutations_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['VCF'].tolist()[0])
        actual_mask_file_C = Path(fofn_df[fofn_df['Sample'] == 'SampleC']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_C, expected_mutations['SampleC'])
        assert_masks_equal(actual_mask_file_C, empty_bed_file)


def test_create_fofn_file_multiple_samples_with_ns_deletions():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        samples = ['SampleD', 'SampleE', 'SampleF', 'SampleH']

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

        assert 4 == len(fofn_df)
        assert ['SampleD', 'SampleE', 'SampleF', 'SampleH'] == fofn_df['Sample'].tolist()

        actual_mutations_D = Path(fofn_df[fofn_df['Sample'] == 'SampleD']['VCF'].tolist()[0])
        actual_mask_D = Path(fofn_df[fofn_df['Sample'] == 'SampleD']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_D, expected_mutations['SampleD'])
        assert_masks_equal(actual_mask_D, expected_beds['SampleD'])

        actual_mutations_E = Path(fofn_df[fofn_df['Sample'] == 'SampleE']['VCF'].tolist()[0])
        actual_mask_E = Path(fofn_df[fofn_df['Sample'] == 'SampleE']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_E, expected_mutations['SampleE'])
        assert_masks_equal(actual_mask_E, expected_beds['SampleE'])

        actual_mutations_F = Path(fofn_df[fofn_df['Sample'] == 'SampleF']['VCF'].tolist()[0])
        actual_mask_F = Path(fofn_df[fofn_df['Sample'] == 'SampleF']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_F, expected_mutations['SampleF'])
        assert_masks_equal(actual_mask_F, expected_beds['SampleF'])

        actual_mutations_H = Path(fofn_df[fofn_df['Sample'] == 'SampleH']['VCF'].tolist()[0])
        actual_mask_H = Path(fofn_df[fofn_df['Sample'] == 'SampleH']['Mask File'].tolist()[0])
        assert_vcf(actual_mutations_H, expected_mutations['SampleH'])
        assert_masks_equal(actual_mask_H, expected_beds['SampleH'])


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
        actual_mask_file = tmp_dir / 'assemblies' / 'mask' / 'SampleA.bed.gz'
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
        assert actual_mask_file.exists()
        assert actual_sketch_file.exists()
        assert actual_mlst_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File', 'Sketch File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_mask_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])
        actual_sketch_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Sketch File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_mask_file == actual_mask_file_from_df
        assert actual_sketch_file == actual_sketch_file_from_df
        assert actual_mlst_file == mlst_file


def test_create_fofn_file_single_no_sketch_with_mlst():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'assemblies' / 'variant' / 'SampleA.vcf.gz'
        actual_mask_file = tmp_dir / 'assemblies' / 'mask' / 'SampleA.bed.gz'
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
        assert actual_mask_file.exists()
        assert not unexpected_sketch_file.exists()
        assert actual_mlst_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_mask_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_mask_file == actual_mask_file_from_df
        assert actual_mlst_file == mlst_file


def test_create_fofn_file_snpeff_no_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file = tmp_dir / 'assemblies' / 'variant-snpeff' / 'SampleA.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_mask_file = tmp_dir / 'assemblies' / 'mask' / 'SampleA.bed.gz'
        input_samples = [snpeff_input_sampleA]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False,
                                                      include_mlst=False, include_kmer=False,
                                                      snpeff_no_check=True)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file.exists()
        assert actual_mask_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_mask_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file == actual_mutations_file_from_df
        assert actual_mask_file == actual_mask_file_from_df

        # snpeff annotations should be added in headers
        reader = vcfpy.Reader.from_path(path=str(actual_mutations_snpeff_file))
        assert 'ANN' in reader.header.info_ids()


def test_create_fofn_file_snpeff_with_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file = tmp_dir / 'assemblies' / 'variant-snpeff' / 'SampleA.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_mask_file = tmp_dir / 'assemblies' / 'mask' / 'SampleA.bed.gz'
        input_samples = [snpeff_input_sampleA]

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=True,
                                                      include_mlst=False, include_kmer=False,
                                                      snpeff_no_check=True)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file.exists()
        assert actual_mask_file.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_mask_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file == actual_mutations_file_from_df
        assert actual_mask_file == actual_mask_file_from_df

        # snpeff annotations should be added in headers
        reader = vcfpy.Reader.from_path(path=str(actual_mutations_snpeff_file))
        assert 'ANN' in reader.header.info_ids()


def test_create_fofn_file_snpeff_reads_with_conda():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_snpeff_file_paired = tmp_dir / 'reads' / 'paired' / 'variant-snpeff' / 'SampleA-snpeff.vcf.gz'
        actual_mutations_snpeff_file_single = tmp_dir / 'reads' / 'single' / 'variant-snpeff' / 'SampleA-single-snpeff.vcf.gz'
        actual_snpeff_config = tmp_dir / 'snpeff_db' / 'snpEff.config'
        actual_mask_file_paired = tmp_dir / 'reads' / 'paired' / 'mask' / 'SampleA-snpeff.bed.gz'
        actual_mask_file_single = tmp_dir / 'reads' / 'single' / 'mask' / 'SampleA-single-snpeff.bed.gz'
        input_samples = snpeff_reads_paired + snpeff_reads_single

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=True,
                                                      include_mlst=False, include_kmer=False,
                                                      reads_mincov=5, snpeff_no_check=True)

        sample_files = pipeline_executor.create_input_sample_files(input_samples)
        results = pipeline_executor.execute(sample_files=sample_files,
                                            reference_file=snpeff_reference_genome,
                                            ncores=1)

        input_fofn = results.get_file('gdi-fofn')

        assert input_fofn.exists()
        assert actual_snpeff_config.exists()
        assert actual_mutations_snpeff_file_paired.exists()
        assert actual_mutations_snpeff_file_single.exists()
        assert actual_mask_file_paired.exists()
        assert actual_mask_file_single.exists()

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 2 == len(fofn_df)
        assert {'SampleA-snpeff', 'SampleA-single-snpeff'} == set(fofn_df['Sample'].tolist())
        actual_mutations_file_paired_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA-snpeff']['VCF'].tolist()[0])
        actual_mutations_file_single_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-single-snpeff']['VCF'].tolist()[0])
        actual_mask_file_paired_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-snpeff']['Mask File'].tolist()[0])
        actual_mask_file_single_from_df = Path(
            fofn_df[fofn_df['Sample'] == 'SampleA-single-snpeff']['Mask File'].tolist()[0])

        assert actual_mutations_snpeff_file_paired == actual_mutations_file_paired_from_df
        assert actual_mutations_snpeff_file_single == actual_mutations_file_single_from_df
        assert actual_mask_file_paired == actual_mask_file_paired_from_df
        assert actual_mask_file_single == actual_mask_file_single_from_df

        # snpeff annotations should be added in headers
        reader = vcfpy.Reader.from_path(path=str(actual_mutations_snpeff_file_paired))
        assert 'ANN' in reader.header.info_ids()

        reader = vcfpy.Reader.from_path(path=str(actual_mutations_snpeff_file_single))
        assert 'ANN' in reader.header.info_ids()


def test_split_input_sequence_files_single_record():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir)
        sample_data = pipeline_executor.split_input_sequence_files([reference_file], output_dir=tmp_dir)

        expected_out1 = tmp_dir / 'reference.fasta.gz'

        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == sample_data.columns.tolist()
        assert 1 == len(sample_data)
        assert ['reference'] == sample_data['Sample'].tolist()
        assert [expected_out1] == sample_data['Assemblies'].tolist()
        assert sample_data['Reads1'].isna().all()
        assert sample_data['Reads2'].isna().all()

        assert expected_out1.exists()

        sf = SequenceFile(expected_out1)
        name, records = sf.parse_sequence_file()
        assert 5180 == len(records[0])


def test_split_input_sequence_files_multiple_records():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir)
        sample_data = pipeline_executor.split_input_sequence_files([reference_file,
                                                                    reference_file_5000_snpeff_2], output_dir=tmp_dir)
        sample_data = sample_data.sort_values('Sample')

        expected_out1 = tmp_dir / 'CP001602.2.gbk.gz'
        expected_out2 = tmp_dir / 'NC_011083.1.gbk.gz'
        expected_out3 = tmp_dir / 'reference.fasta.gz'

        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == sample_data.columns.tolist()
        assert 3 == len(sample_data)
        assert ['CP001602.2', 'NC_011083.1', 'reference'] == sample_data['Sample'].tolist()
        assert [expected_out1, expected_out2, expected_out3] == sample_data['Assemblies'].tolist()
        assert sample_data['Reads1'].isna().all()
        assert sample_data['Reads2'].isna().all()

        assert expected_out1.exists()
        assert expected_out2.exists()
        assert expected_out3.exists()

        sf1 = SequenceFile(expected_out1)
        name, records = sf1.parse_sequence_file()
        assert 5000 == len(records[0])

        sf2 = SequenceFile(expected_out1)
        name, records = sf2.parse_sequence_file()
        assert 5000 == len(records[0])

        sf3 = SequenceFile(expected_out3)
        name, records = sf3.parse_sequence_file()
        assert 5180 == len(records[0])


def test_split_input_sequence_files_multiple_records_subsample():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        samples = {'CP001602.2', 'reference'}

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir)
        sample_data = pipeline_executor.split_input_sequence_files([reference_file,
                                                                    reference_file_5000_snpeff_2],
                                                                   output_dir=tmp_dir,
                                                                   samples=samples)
        sample_data = sample_data.sort_values('Sample')

        expected_out1 = tmp_dir / 'CP001602.2.gbk.gz'
        unexpected_out2 = tmp_dir / 'NC_011083.1.gbk.gz'
        expected_out3 = tmp_dir / 'reference.fasta.gz'

        assert ['Sample', 'Assemblies', 'Reads1', 'Reads2'] == sample_data.columns.tolist()
        assert 2 == len(sample_data)
        assert ['CP001602.2', 'reference'] == sample_data['Sample'].tolist()
        assert [expected_out1, expected_out3] == sample_data['Assemblies'].tolist()
        assert sample_data['Reads1'].isna().all()
        assert sample_data['Reads2'].isna().all()

        assert expected_out1.exists()
        assert not unexpected_out2.exists()
        assert expected_out3.exists()

        sf1 = SequenceFile(expected_out1)
        name, records = sf1.parse_sequence_file()
        assert 5000 == len(records[0])

        sf3 = SequenceFile(expected_out3)
        name, records = sf3.parse_sequence_file()
        assert 5180 == len(records[0])


def test_select_random_samples():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        samples = {'CP001602.2', 'reference', 'NC_011083.1'}

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir)

        # Test select 1 sample
        subsamples = pipeline_executor.select_random_samples([reference_file, reference_file_5000_snpeff_2],
                                                             number_samples=1)
        assert 1 == len(subsamples)
        assert subsamples.issubset(samples)

        # Test select 2 samples
        subsamples = pipeline_executor.select_random_samples([reference_file, reference_file_5000_snpeff_2],
                                                             number_samples=2)
        assert 2 == len(subsamples)
        assert subsamples.issubset(samples)

        # Test select 3 samples
        subsamples = pipeline_executor.select_random_samples([reference_file, reference_file_5000_snpeff_2],
                                                             number_samples=3)
        assert 3 == len(subsamples)
        assert subsamples.issubset(samples)

        # Test select 33% of samples
        subsamples = pipeline_executor.select_random_samples([reference_file, reference_file_5000_snpeff_2],
                                                             number_samples=0.33)
        assert 1 == len(subsamples)
        assert subsamples.issubset(samples)

        # Test select 66% of samples
        subsamples = pipeline_executor.select_random_samples([reference_file, reference_file_5000_snpeff_2],
                                                             number_samples=0.66)
        assert 2 == len(subsamples)
        assert subsamples.issubset(samples)
