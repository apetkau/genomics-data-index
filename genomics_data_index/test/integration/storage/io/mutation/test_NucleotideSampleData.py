import tempfile

import pytest
from pathlib import Path
import logging

from pybedtools import BedTool

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.test.integration import snippy_sample_vcfs_dict, snippy_sample_mask_sequences_dict
from genomics_data_index.test.integration import snpeff_sample_vcfs


logger = logging.getLogger(__name__)


@pytest.fixture
def sample_dataA():
    tmp_file = Path(tempfile.mktemp())
    seq_name, seq_file = SequenceFile(snippy_sample_mask_sequences_dict['SampleA']).parse_sequence_file()
    mask_region = MaskedGenomicRegions.from_sequences(seq_file)
    mask_region.write(tmp_file)
    sample_data = NucleotideSampleData(sample_name='SampleA',
                                       vcf_file=snippy_sample_vcfs_dict['SampleA'],
                                       vcf_file_index=None,
                                       mask_bed_file=tmp_file,
                                       preprocessed=True)
    return sample_data


@pytest.fixture
def sample_data_sh10_014():
    tmp_file = Path(tempfile.mktemp())
    mask_region = MaskedGenomicRegions(BedTool([('NC_011083', 0, 100), ('NC_011083', 199, 250)]))
    mask_region.write(tmp_file)
    sample_data = NucleotideSampleData(sample_name='SH10-014',
                                       vcf_file=snpeff_sample_vcfs['SH10-014'],
                                       vcf_file_index=None,
                                       mask_bed_file=tmp_file,
                                       preprocessed=True)
    return sample_data


def test_read_sample_data_features_no_unknown(sample_dataA):
    df = sample_dataA.read_sample_data_features(include_masked_regions=False).sort_values('POS')

    assert 46 == len(df)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(df.columns)
    assert {'SampleA'} == set(df['SAMPLE'].tolist())

    # Variant position
    v = df[df['POS'] == 461]
    assert ['reference'] == v['CHROM'].tolist()
    assert ['AAAT'] == v['REF'].tolist()
    assert ['G'] == v['ALT'].tolist()
    assert ['reference:461:AAAT:G'] == v['VARIANT_ID'].tolist()

    # Variant position
    v = df[df['POS'] == 1048]
    assert ['reference'] == v['CHROM'].tolist()
    assert ['C'] == v['REF'].tolist()
    assert ['G'] == v['ALT'].tolist()
    assert ['reference:1048:C:G'] == v['VARIANT_ID'].tolist()

    # No snpeff annotations
    assert all(df['ANN.Annotation'].isna())

    # No missing/unknown
    assert not any(df['ALT'] == '?')


def test_read_sample_data_features_with_unknown(sample_dataA):
    df = sample_dataA.read_sample_data_features(include_masked_regions=True).sort_values('POS')

    assert 482 == len(df)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(df.columns)
    assert {'SampleA'} == set(df['SAMPLE'].tolist())

    # Number actual variants
    assert 45 == len(df[df['ALT'] != '?'])

    # Number missing
    assert 437 == len(df[df['ALT'] == '?'])

    # 461... is missing
    assert ['reference:460:1:?', 'reference:461:1:?',
            'reference:462:1:?'] == df[
        (df['POS'] >= 460) & (df['POS'] <= 465)]['VARIANT_ID'].tolist()
    logger.warning('I need to verify that mask coordinates are correct')

    # Variant position
    v = df[df['POS'] == 1048]
    assert ['reference'] == v['CHROM'].tolist()
    assert ['C'] == v['REF'].tolist()
    assert ['G'] == v['ALT'].tolist()
    assert ['reference:1048:C:G'] == v['VARIANT_ID'].tolist()

    # No snpeff annotations
    assert all(df['ANN.Annotation'].isna())

    # One of the missing
    assert ['reference:1:1:?'] == df[df['POS'] == 1]['VARIANT_ID'].tolist()


def test_read_sample_data_features_snpeff_with_unknown(sample_data_sh10_014):
    df = sample_data_sh10_014.read_sample_data_features(include_masked_regions=True).sort_values('POS')

    num_variants = 139
    num_masked = 151
    assert num_variants + num_masked == len(df)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(df.columns)
    assert {'SH10-014'} == set(df['SAMPLE'].tolist())

    # Number actual variants
    assert num_variants == len(df[df['ALT'] != '?'])

    # Number missing
    assert num_masked == len(df[df['ALT'] == '?'])

    variant_ids = set(df['VARIANT_ID'].tolist())

    # Test to verify some are mossing
    assert 'NC_011083:1:1:?' in variant_ids
    assert 'NC_011083:100:1:?' in variant_ids
    assert 'NC_011083:101:1:?' not in variant_ids
    assert 'NC_011083:199:1:?' not in variant_ids
    assert 'NC_011083:200:1:?' in variant_ids
    assert 'NC_011083:250:1:?' in variant_ids
    assert 'NC_011083:251:1:?' not in variant_ids

    # Variant position
    v = df[df['POS'] == 140658]
    assert ['NC_011083'] == v['CHROM'].tolist()
    assert ['C'] == v['REF'].tolist()
    assert ['A'] == v['ALT'].tolist()
    assert ['NC_011083:140658:C:A'] == v['VARIANT_ID'].tolist()
    assert ['missense_variant'] == v['ANN.Annotation'].tolist()
    assert ['murF'] == v['ANN.Gene_Name'].tolist()
    assert ['c.497C>A'] == v['ANN.HGVS.c'].tolist()
    assert ['p.Ala166Glu'] == v['ANN.HGVS.p'].tolist()
