import pytest
import pandas as pd

from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser, \
    InvalidSnpEffVcfError


@pytest.fixture
def vcf_snpeff_annotation_parser() -> VcfSnpEffAnnotationParser:
    return VcfSnpEffAnnotationParser()


@pytest.fixture
def mock_snpeff_infos():
    class MockAnn():
        def __init__(self):
            pass

        desc = ("Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID"
                " | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p"
                " | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance"
                " | ERRORS / WARNINGS / INFO'")

    return {
        'ANN': MockAnn()
    }


@pytest.fixture
def mock_vcf_df_with_ann_single() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A',
         {'ANN': ('A|missense_variant|MODERATE|murF|SEHA_RS01180|transcript|SEHA_RS01180|'
                  'protein_coding|1/1|c.497C>A|p.Ala166Glu|497/1359|497/1359|166/452||')}
         ],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_vcf_df_with_ann_multiple_entries_single_sample() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A',
         {'ANN': [('A|missense_variant|MODERATE|murF|SEHA_RS01180|transcript|SEHA_RS01180|'
                  'protein_coding|1/1|c.497C>A|p.Ala166Glu|497/1359|497/1359|166/452||'),
                  ('A|upstream_gene_variant|MODIFIER|mraY|SEHA_RS01185|transcript|SEHA_RS01185|'
                  'protein_coding||c.-856C>A|||||856|'),
                  ('A|upstream_gene_variant|MODIFIER|murD|SEHA_RS01190|transcript|SEHA_RS01190|'
                  'protein_coding||c.-1941C>A|||||1941|')]}
         ],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_snpeff_infos_invalid():
    class MockAnn():
        def __init__(self):
            pass

        desc = 'invalid'

    return {
        'ANN': MockAnn()
    }


def test_parse_annotation_headers(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser, mock_snpeff_infos):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_infos)

    assert ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID',
            'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length',
            'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO'] == headers_list


def test_parse_annotation_headers_invalid(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                          mock_snpeff_infos_invalid):
    with pytest.raises(InvalidSnpEffVcfError) as execinfo:
        vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_infos_invalid)
    assert "Found 'ANN' in VCF information but description" in str(execinfo.value)


def test_parse_annotation_entries_single(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                         mock_snpeff_infos, mock_vcf_df_with_ann_single: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_infos)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_ann_single)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert [0] == list(ann_entries_df.index)
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])


def test_parse_annotation_entries_multiple_entries_single_sample(
        vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
        mock_snpeff_infos, mock_vcf_df_with_ann_multiple_entries_single_sample: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_infos)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_ann_multiple_entries_single_sample)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(ann_entries_df.columns)
    print(ann_entries_df)
    assert 3 == len(ann_entries_df)
    assert [0, 0, 0] == list(ann_entries_df.index)
    print(ann_entries_df['ANN.HGVS.p'].isna())
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
            'c.-856C>A', pd.NA] == list(ann_entries_df.iloc[1])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
            'c.-1941C>A', pd.NA] == list(ann_entries_df.iloc[2])
