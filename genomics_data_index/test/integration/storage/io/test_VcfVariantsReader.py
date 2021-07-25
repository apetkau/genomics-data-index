import tempfile
from os import path, listdir
from pathlib import Path
from typing import List, Dict, cast

import pytest

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration import data_dir_empty
from genomics_data_index.test.integration import snpeff_sample_vcfs


@pytest.fixture
def sample_dirs() -> List[Path]:
    return [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


@pytest.fixture
def sample_dirs_empty() -> List[Path]:
    return [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


def vcf_and_mask_files(sample_dirs) -> Dict[str, Dict[str, Path]]:
    sample_vcf_map = {}
    sample_genomic_files_mask = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        genomic_file_mask = Path(d, 'snps.aligned.fa')

        sample_vcf_map[sample_name] = vcf_file
        sample_genomic_files_mask[sample_name] = genomic_file_mask

    return {
        'vcfs': sample_vcf_map,
        'masks': sample_genomic_files_mask
    }


def variants_reader_internal(sample_dirs, include_masked_regions: bool) -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    vcf_masks = vcf_and_mask_files(sample_dirs)
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                          masked_genomic_files_map=vcf_masks[
                                                                              'masks'],
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=include_masked_regions)


@pytest.fixture
def variants_reader(sample_dirs) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs, include_masked_regions=False)


@pytest.fixture
def variants_reader_empty(sample_dirs_empty) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs_empty, include_masked_regions=False)


@pytest.fixture
def variants_reader_snpeff_annotations_single_sample() -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)

    vcfs_map = {
        'SH10-014': snpeff_sample_vcfs['SH10-014']
    }

    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcfs_map,
                                                                          masked_genomic_files_map=None,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=False)


@pytest.fixture
def variants_reader_snpeff_annotations_multiple_samples() -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                          masked_genomic_files_map=None,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=False)


@pytest.fixture
def variants_reader_empty_masks(sample_dirs) -> VcfVariantsReader:
    sample_vcf_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')

        sample_vcf_map[sample_name] = vcf_file

    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcf_map,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=False)


@pytest.fixture
def variants_reader_default_no_data() -> VcfVariantsReader:
    return VcfVariantsReader(sample_files_map={})


def variants_reader_from_snippy_internal(sample_dirs, include_masked_regions: bool,
                                         ncores: int) -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs=sample_dirs,
                                                                  sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=include_masked_regions,
                                    ncores=ncores)


@pytest.fixture
def variants_reader_from_snippy(sample_dirs) -> VcfVariantsReader:
    return variants_reader_from_snippy_internal(sample_dirs, include_masked_regions=False, ncores=1)


@pytest.fixture
def variants_reader_from_snippy_masked(sample_dirs) -> VcfVariantsReader:
    return variants_reader_from_snippy_internal(sample_dirs, include_masked_regions=True, ncores=1)


@pytest.fixture
def variants_reader_from_snippy_masked_multicore(sample_dirs) -> VcfVariantsReader:
    return variants_reader_from_snippy_internal(sample_dirs, include_masked_regions=True, ncores=2)


@pytest.fixture
def variants_reader_snpeff() -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files, include_masked_regions=False)


def test_get_variants_table(variants_reader):
    df = variants_reader.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    assert ['SNP'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 293), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 302), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 324), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 374), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 461), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleB') & (df['POS'] == 1325), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleC') & (df['POS'] == 1984), 'TYPE'].tolist()


def test_get_genomic_masks_and_mask_to_features(variants_reader):
    mask = variants_reader.get_genomic_masked_region('SampleA')
    assert 437 == len(mask)
    assert {'reference'} == mask.sequence_names()
    assert 437 == len(variants_reader.mask_to_features(mask))

    mask = variants_reader.get_genomic_masked_region('SampleB')
    assert 276 == len(mask)
    assert {'reference'} == mask.sequence_names()
    assert 276 == len(variants_reader.mask_to_features(mask))

    mask = variants_reader.get_genomic_masked_region('SampleC')
    assert 329 == len(mask)
    assert {'reference'} == mask.sequence_names()
    assert 329 == len(variants_reader.mask_to_features(mask))


def test_get_genomic_masks_empty(variants_reader_empty_masks):
    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleA')
    assert mask.is_empty()
    assert 0 == len(variants_reader_empty_masks.mask_to_features(mask))

    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleB')
    assert mask.is_empty()
    assert 0 == len(variants_reader_empty_masks.mask_to_features(mask))

    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleC')
    assert mask.is_empty()
    assert 0 == len(variants_reader_empty_masks.mask_to_features(mask))


@pytest.mark.skip()
def test_combine_vcf_mask(variants_reader, sample_dirs):
    expected_columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                        'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                        'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    vcf_mask_files_dict = vcf_and_mask_files(sample_dirs)

    sampleA_vcf = vcf_mask_files_dict['vcfs']['SampleA']
    sampleA_vcf_df = variants_reader.read_vcf(sampleA_vcf, 'SampleA')
    assert 46 == len(sampleA_vcf_df)
    assert 461 in sampleA_vcf_df['POS'].tolist()

    mask = variants_reader.get_genomic_masked_region('SampleA')
    sampleA_mask_df = variants_reader.mask_to_features(mask)
    assert 437 == len(sampleA_mask_df)
    assert 461 in sampleA_mask_df['POS'].tolist()

    combined_df = variants_reader.combine_vcf_mask(vcf_frame=sampleA_vcf_df, mask_frame=sampleA_mask_df)

    assert expected_columns == combined_df.columns.tolist()
    assert 482 == len(combined_df)
    assert 436 == len(combined_df[combined_df['TYPE'] == 'UNKNOWN_MISSING'])
    assert 46 == len(combined_df[combined_df['TYPE'] != 'UNKNOWN_MISSING'])

    # Make sure it's the SNV/SNP from the VCF file that I keep and not the missing position in this file
    assert combined_df[combined_df['POS'] == 461]['TYPE'] == 'INDEL'


def test_get_samples_list(variants_reader):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader.samples_list())


def test_get_variants_table_empty(variants_reader_empty):
    df = variants_reader_empty.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'


def test_get_or_create_feature_file(variants_reader):
    file = variants_reader.get_or_create_feature_file('SampleA')
    assert file.exists()


def test_snippy_read_vcf(variants_reader_from_snippy):
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'

    df = variants_reader_from_snippy.read_vcf(vcf_file, 'SampleA')

    assert 46 == len(df), 'Data fram has incorrect length'

    assert {'snps.vcf.gz'} == set(df['FILE'].tolist()), 'Incorrect filename'
    assert {'SampleA'} == set(df['SAMPLE'].tolist()), 'Incorrect sample name'

    v = df[df['POS'] == 461]
    assert 'AAAT' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1048]
    assert 'C' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1253]
    assert 'T' == v['REF'].values[0], 'Incorrect reference'
    assert 'TAA' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 3656]
    assert 'CATT' == v['REF'].values[0], 'Incorrect reference'
    assert 'C' == v['ALT'].values[0], 'Incorrect alt'


def test_snippy_get_variants_table(variants_reader_from_snippy):
    df = variants_reader_from_snippy.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'


def test_snippy_get_variants_table_masked(variants_reader_from_snippy_masked):
    df = variants_reader_from_snippy_masked.get_features_table()

    assert 1170 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    # Missing/unknown
    assert 437 == len(df[(df['SAMPLE'] == 'SampleA') & (df['TYPE'] == 'UNKNOWN_MISSING')])
    assert 276 == len(df[(df['SAMPLE'] == 'SampleB') & (df['TYPE'] == 'UNKNOWN_MISSING')])
    assert 329 == len(df[(df['SAMPLE'] == 'SampleC') & (df['TYPE'] == 'UNKNOWN_MISSING')])

    # Variants
    assert 45 == len(df[(df['SAMPLE'] == 'SampleA') & (df['TYPE'] != 'UNKNOWN_MISSING')])
    assert 50 == len(df[(df['SAMPLE'] == 'SampleB') & (df['TYPE'] != 'UNKNOWN_MISSING')])
    assert 33 == len(df[(df['SAMPLE'] == 'SampleC') & (df['TYPE'] != 'UNKNOWN_MISSING')])


def test_snippy_get_variants_table_masked_multicore(variants_reader_from_snippy_masked_multicore):
    df = variants_reader_from_snippy_masked_multicore.get_features_table()

    assert 1170 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    # Missing/unknown
    assert 437 == len(df[(df['SAMPLE'] == 'SampleA') & (df['TYPE'] == 'UNKNOWN_MISSING')])
    assert 276 == len(df[(df['SAMPLE'] == 'SampleB') & (df['TYPE'] == 'UNKNOWN_MISSING')])
    assert 329 == len(df[(df['SAMPLE'] == 'SampleC') & (df['TYPE'] == 'UNKNOWN_MISSING')])

    # Variants
    assert 45 == len(df[(df['SAMPLE'] == 'SampleA') & (df['TYPE'] != 'UNKNOWN_MISSING')])
    assert 50 == len(df[(df['SAMPLE'] == 'SampleB') & (df['TYPE'] != 'UNKNOWN_MISSING')])
    assert 33 == len(df[(df['SAMPLE'] == 'SampleC') & (df['TYPE'] != 'UNKNOWN_MISSING')])


def test_snippy_get_genomic_masks(variants_reader_from_snippy):
    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleA')
    assert 437 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleB')
    assert 276 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleC')
    assert 329 == len(mask)
    assert {'reference'} == mask.sequence_names()


def test_snippy_get_samples_list(variants_reader_from_snippy):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader_from_snippy.samples_list())


def test_snippy_get_samples_list_two_files():
    sample_dirs = [data_dir / 'SampleA', data_dir / 'SampleB']
    reader = variants_reader_from_snippy_internal(sample_dirs, include_masked_regions=False,
                                                  ncores=1)

    assert {'SampleA', 'SampleB'} == set(reader.samples_list())


def test_snippy_read_empty_vcf(sample_dirs_empty):
    reader = variants_reader_from_snippy_internal(sample_dirs_empty, include_masked_regions=False,
                                                  ncores=1)
    df = reader.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'


def test_snippy_read_empty_vcf_include_masked_regions(sample_dirs_empty):
    reader = variants_reader_from_snippy_internal(sample_dirs_empty, include_masked_regions=True,
                                                  ncores=1)
    df = reader.get_features_table()

    assert 437 == len(df), 'Data has incorrect length'
    assert {'UNKNOWN_MISSING'} == set(df['TYPE'])
    assert {'?'} == set(df['ALT'])


def test_snippy_read_empty_vcf_multicore_include_masked_regions(sample_dirs_empty):
    reader = variants_reader_from_snippy_internal(sample_dirs_empty, include_masked_regions=True,
                                                  ncores=2)
    df = reader.get_features_table()

    assert 437 == len(df), 'Data has incorrect length'
    assert {'UNKNOWN_MISSING'} == set(df['TYPE'])
    assert {'?'} == set(df['ALT'])


def test_read_snpeff(variants_reader_default_no_data: VcfVariantsReader):
    vr = variants_reader_default_no_data
    sample_10_014 = vr.read_vcf(file=snpeff_sample_vcfs['SH10-014'], sample_name='SH10-014').sort_values('POS')
    sample_14_001 = vr.read_vcf(file=snpeff_sample_vcfs['SH14-001'], sample_name='SH14-001').sort_values('POS')
    sample_14_014 = vr.read_vcf(file=snpeff_sample_vcfs['SH14-014'], sample_name='SH14-014').sort_values('POS')

    assert 1434 == len(sample_10_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_10_014.columns)

    # snv/snp
    sample_10_014_varA = sample_10_014[sample_10_014['POS'] == 140658]
    assert 9 == len(sample_10_014_varA)
    assert ['SH10-014', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH10-014.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_10_014_varA[
               sample_10_014_varA['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    # del
    sample_10_014_varB = sample_10_014[sample_10_014['POS'] == 1125996]
    assert 14 == len(sample_10_014_varB)
    assert ['SH10-014', 'NC_011083', 1125996, 'CG', 'C', 'del', 'SH10-014.vcf.gz', 'NC_011083:1125996:CG:C',
            'C', 'frameshift_variant', 'HIGH', 'SEHA_RS05995', 'SEHA_RS05995', 'transcript', 'protein_coding',
            'c.418delG', 'p.Glu140fs'] == sample_10_014_varB[
               sample_10_014_varB['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # ins
    sample_10_014_varC = sample_10_014[sample_10_014['POS'] == 1246085]
    assert 11 == len(sample_10_014_varC)
    assert ['SH10-014', 'NC_011083', 1246085, 'C', 'CG', 'ins', 'SH10-014.vcf.gz', 'NC_011083:1246085:C:CG',
            'CG', 'frameshift_variant', 'HIGH', 'mdtG', 'SEHA_RS06605', 'transcript', 'protein_coding',
            'c.722dupC', 'p.Leu242fs'] == sample_10_014_varC[
               sample_10_014_varC['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # complex
    sample_10_014_varD = sample_10_014[sample_10_014['POS'] == 3535121]
    assert 10 == len(sample_10_014_varD)
    assert ['SH10-014', 'NC_011083', 3535121, 'CGCGA', 'TGTGG', 'complex', 'SH10-014.vcf.gz',
            'NC_011083:3535121:CGCGA:TGTGG',
            'TGTGG', 'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1119_1123delTCGCGinsCCACA', 'p.ArgAla374HisThr'] == sample_10_014_varD[
               sample_10_014_varD['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 1218 == len(sample_14_001)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_001.columns)
    sample_14_001_var = sample_14_001[sample_14_001['POS'] == 140658]
    assert 9 == len(sample_14_001_var)
    assert ['SH14-001', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH14-001.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_14_001_var[
               sample_14_001_var['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 1128 == len(sample_14_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_014.columns)
    sample_14_014_var = sample_14_014[sample_14_014['POS'] == 298472]
    assert 6 == len(sample_14_014_var)
    assert ['SH14-014', 'NC_011083', 298472, 'A', 'C', 'snp', 'SH14-014.vcf.gz', 'NC_011083:298472:A:C',
            'C', 'intergenic_region', 'MODIFIER', 'SEHA_RS01880-SEHA_RS01885', 'SEHA_RS01880-SEHA_RS01885',
            'intergenic_region', 'n.298472A>C'] == sample_14_014_var[
               sample_14_014_var['ANN.Annotation'] == 'intergenic_region'].drop(
        ['ANN.Transcript_BioType', 'ANN.HGVS.p'], axis='columns').iloc[0].tolist()
    assert {True} == set(sample_14_014_var[sample_14_014_var['ANN.Annotation'] == 'intergenic_region'] \
                             [['ANN.Transcript_BioType', 'ANN.HGVS.p']].iloc[0].isna().tolist())


def test_get_variants_table_snpeff_annotations_single_sample(
        variants_reader_snpeff_annotations_single_sample: VcfVariantsReader):
    vr = variants_reader_snpeff_annotations_single_sample
    vcf_df = vr.get_features_table()

    assert 139 == len(vcf_df)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'VARIANT_ID'] == list(
        vcf_df.columns)

    # Get rid of 'FILE' since these are processed files
    vcf_df = vcf_df.drop(columns='FILE')

    # # missense variant
    vcf_df_varA = vcf_df[vcf_df['POS'] == 140658]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 140658, 'C', 'A', 'SNP',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu', 'NC_011083:140658:C:A'] == vcf_df_varA.iloc[0].tolist()

    # Synonymous variant
    vcf_df_varA = vcf_df[vcf_df['POS'] == 723772]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 723772, 'T', 'C', 'SNP',
            'C', 'synonymous_variant', 'LOW', 'SEHA_RS03990', 'SEHA_RS03990', 'transcript', 'protein_coding',
            'c.1638T>C', 'p.Cys546Cys', 'NC_011083:723772:T:C'] == vcf_df_varA.iloc[0].tolist()

    # Intergenic variant (snp)
    vcf_df_varA = vcf_df[vcf_df['POS'] == 1031571]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 1031571, 'T', 'C', 'SNP',
            'C', 'intergenic_region', 'MODIFIER', 'SEHA_RS05505-SEHA_RS05515', 'SEHA_RS05505-SEHA_RS05515',
            'intergenic_region', 'NA',
            'n.1031571T>C', 'NA', 'NC_011083:1031571:T:C'] == vcf_df_varA.iloc[0].fillna('NA').tolist()

    # Intergenic variant (del)
    vcf_df_varA = vcf_df[vcf_df['POS'] == 4555461]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 4555461, 'T', 'TC', 'INDEL',
            'TC', 'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA', 'NC_011083:4555461:T:TC'] == vcf_df_varA.iloc[0].fillna('NA').tolist()

    # frameshift variant
    vcf_df_varA = vcf_df[vcf_df['POS'] == 1770751]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 1770751, 'CA', 'C', 'INDEL',
            'C', 'frameshift_variant', 'HIGH', 'SEHA_RS09270', 'SEHA_RS09270',
            'transcript', 'protein_coding',
            'c.335delT', 'p.Leu112fs', 'NC_011083:1770751:CA:C'] == vcf_df_varA.iloc[0].fillna('NA').tolist()

    # stop gained variant
    vcf_df_varA = vcf_df[vcf_df['POS'] == 4882099]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 4882099, 'C', 'T', 'SNP',
            'T', 'stop_gained', 'HIGH', 'SEHA_RS24155', 'SEHA_RS24155',
            'transcript', 'protein_coding',
            'c.1287G>A', 'p.Trp429*', 'NC_011083:4882099:C:T'] == vcf_df_varA.iloc[0].fillna('NA').tolist()

    # stop lost variant
    vcf_df_varA = vcf_df[vcf_df['POS'] == 4824790]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 4824790, 'T', 'C', 'SNP',
            'C', 'stop_lost', 'HIGH', 'SEHA_RS23860', 'SEHA_RS23860',
            'transcript', 'protein_coding',
            'c.1378T>C', 'p.Ter460Glnext*?', 'NC_011083:4824790:T:C'] == vcf_df_varA.iloc[0].fillna('NA').tolist()

    # inframe deletion
    vcf_df_varA = vcf_df[vcf_df['POS'] == 4465400]
    assert 1 == len(vcf_df_varA)
    assert ['SH10-014', 'NC_011083', 4465400, 'GGCCGAA', 'G', 'INDEL',
            'G', 'conservative_inframe_deletion', 'MODERATE', 'tyrB', 'SEHA_RS22180',
            'transcript', 'protein_coding',
            'c.157_162delGAAGCC', 'p.Glu53_Ala54del', 'NC_011083:4465400:GGCCGAA:G'] == vcf_df_varA.iloc[0].fillna(
        'NA').tolist()


def test_get_variants_table_snpeff_annotations_multiple_samples(
        variants_reader_snpeff_annotations_multiple_samples: VcfVariantsReader):
    vr = variants_reader_snpeff_annotations_multiple_samples
    vcf_df = vr.get_features_table()

    assert 361 == len(vcf_df)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'VARIANT_ID'] == list(
        vcf_df.columns)

    # Get rid of 'FILE' since these are processed files
    vcf_df = vcf_df.drop(columns='FILE')

    # disruptive inframe deletion
    vcf_df_varA = vcf_df[vcf_df['POS'] == 3167187]
    assert 3 == len(vcf_df_varA)
    vcf_df_varA = vcf_df_varA.set_index('SAMPLE')
    assert ['NC_011083', 3167187, 'AACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCACG', 'A', 'INDEL',
            'A', 'disruptive_inframe_deletion', 'MODERATE', 'SEHA_RS15905', 'SEHA_RS15905', 'transcript',
            'protein_coding',
            'c.417_464delCGACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCA', 'p.Asp140_His155del',
            'NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCACG:A'] == vcf_df_varA.loc[
               'SH10-014'].tolist()
    assert ['NC_011083', 3167187, 'AACCACGACCACGACCACGACCACGACCACGACCACG', 'A', 'INDEL',
            'A', 'disruptive_inframe_deletion', 'MODERATE', 'SEHA_RS15905', 'SEHA_RS15905', 'transcript',
            'protein_coding',
            'c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA', 'p.Asp144_His155del',
            'NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACG:A'] == vcf_df_varA.loc['SH14-001'].tolist()
    assert ['NC_011083', 3167187, 'AACCACGACCACGACCACGACCACGACCACGACCACG', 'A', 'INDEL',
            'A', 'disruptive_inframe_deletion', 'MODERATE', 'SEHA_RS15905', 'SEHA_RS15905', 'transcript',
            'protein_coding',
            'c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA', 'p.Asp144_His155del',
            'NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACG:A'] == vcf_df_varA.loc['SH14-014'].tolist()
