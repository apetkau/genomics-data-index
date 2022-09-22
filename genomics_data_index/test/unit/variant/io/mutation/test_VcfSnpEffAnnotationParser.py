import pandas as pd
import pytest

from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser, \
    InvalidSnpEffVcfError


@pytest.fixture
def vcf_snpeff_annotation_parser() -> VcfSnpEffAnnotationParser:
    return VcfSnpEffAnnotationParser()


@pytest.fixture
def mock_snpeff_headers():
    class MockAnn():
        def __init__(self):
            pass

        def info_ids(self):
            return ['ANN']

        def get_info_field_info(self, field: str):
            if field == 'ANN':
                obj = type('', (), {})()
                obj.description = ("Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID"
                    " | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p"
                    " | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance"
                    " | ERRORS / WARNINGS / INFO'")
                return obj
            else:
                return None

    return MockAnn()


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
def mock_vcf_df_without_ann_single() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A', {}],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_vcf_df_without_ann_multiple() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A', {}],
        ['NC_011083', 203200, 'C', 'T', {}],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_vcf_df_with_and_without_ann() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A',
         {'ANN': [('A|missense_variant|MODERATE|murF|SEHA_RS01180|transcript|SEHA_RS01180|'
                   'protein_coding|1/1|c.497C>A|p.Ala166Glu|497/1359|497/1359|166/452||'),
                  ('A|upstream_gene_variant|MODIFIER|mraY|SEHA_RS01185|transcript|SEHA_RS01185|'
                   'protein_coding||c.-856C>A|||||856|'),
                  ('A|upstream_gene_variant|MODIFIER|murD|SEHA_RS01190|transcript|SEHA_RS01190|'
                   'protein_coding||c.-1941C>A|||||1941|')]}
         ],
        ['NC_011083', 203200, 'C', 'T', {'ANN': []}],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_vcf_df_empty() -> pd.DataFrame:
    return pd.DataFrame(columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_snpeff_headers_empty():
    class MockAnn():
        def __init__(self):
            pass

        def info_ids(self):
            return []

        def get_info_field_info(self, field: str):
            return None

    return MockAnn()


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
def mock_vcf_df_single_sample_annotations() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
         'c.-856C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
         'c.497C>A', 'p.Ala166Glu', '497/1359'],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
         'c.-1941C>A', pd.NA, pd.NA]
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_single_sample_annotations_some_na_values() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
         'c.-856C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', pd.NA,
         'c.497C>A', 'p.Ala166Glu', '497/1359'],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
         'c.-1941C>A', pd.NA, pd.NA]
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_multiple_sample_annotations() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
         'c.-856C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
         'c.497C>A', 'p.Ala166Glu', '497/1359'],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
         'c.-1941C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T',
         'T', 'missense_variant', 'MODERATE', 'SEHA_RS01460', 'SEHA_RS01460', 'transcript',
         'protein_coding', 'c.602C>T', 'p.Thr201Met', '602/927'],
        ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T',
         'T', 'upstream_gene_variant', 'MODIFIER', 'SEHA_RS01455', 'SEHA_RS01455',
         'transcript', 'protein_coding', 'c.-2172G>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T',
         'T', 'upstream_gene_variant', 'MODIFIER', 'SEHA_RS01455', 'SEHA_RS01455',
         'transcript', 'protein_coding', 'c.-710G>A', pd.NA, pd.NA],
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_single_sample_annotations_sars_cov_2() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.3683A>G', 'p.D1228G', '3683/21291'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.3683A>G', 'p.D1228G', '3683/13218'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.1229A>G', 'p.D410G', '1229/5835'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.1229A>G', 'p.D410G', '1229/5835'],
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_single_sample_annotations_rv_sars_cov_2() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.1229A>G', 'p.D410G', '1229/5835'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.1229A>G', 'p.D410G', '1229/5835'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.3683A>G', 'p.D1228G', '3683/13218'],
        ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
         'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
         'c.3683A>G', 'p.D1228G', '3683/21291'],
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_multiple_sample_one_empty() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
         'c.-856C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
         'c.497C>A', 'p.Ala166Glu', '497/1359'],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
         'c.-1941C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T'] + [pd.NA] * 10,
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length',
                ])


@pytest.fixture
def mock_vcf_df_multiple_sample_one_invalid() -> pd.DataFrame:
    return pd.DataFrame([
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
         'c.-856C>A', pd.NA, pd.NA],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
         'c.497C>A', 'p.Ala166Glu', '497/1359'],
        ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
         'A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
         'c.-1941C>A', pd.NA, pd.NA],

        # This one is invalid because snpeff 'Allele' refers to some other variant (indicating it's a compound variant)
        # I want to ignore all compound variants because it would be very difficult to index them in my software.
        ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T',
         'T:123456_A>T', 'missense_variant', 'MODERATE', 'SEHA_RS01460', 'SEHA_RS01460', 'transcript',
         'protein_coding', 'c.602C>T', 'p.Thr201Met', '602/927'],
    ], columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p',
                'ANN.cDNA.pos / cDNA.length'
                ])


@pytest.fixture
def mock_vcf_df_with_ann_multiple_entries_multiple_samples() -> pd.DataFrame:
    return pd.DataFrame([
        ['NC_011083', 140658, 'C', 'A',
         {'ANN': [('A|missense_variant|MODERATE|murF|SEHA_RS01180|transcript|SEHA_RS01180|'
                   'protein_coding|1/1|c.497C>A|p.Ala166Glu|497/1359|497/1359|166/452||'),
                  ('A|upstream_gene_variant|MODIFIER|mraY|SEHA_RS01185|transcript|SEHA_RS01185|'
                   'protein_coding||c.-856C>A|||||856|'),
                  ('A|upstream_gene_variant|MODIFIER|murD|SEHA_RS01190|transcript|SEHA_RS01190|'
                   'protein_coding||c.-1941C>A|||||1941|')]}
         ],
        ['NC_011083', 203200, 'C', 'T',
         {'ANN': [('T|missense_variant|MODERATE|SEHA_RS01460|SEHA_RS01460|transcript|SEHA_RS01460|'
                   'protein_coding|1/1|c.602C>T|p.Thr201Met|602/927|602/927|201/308||'),
                  ('T|upstream_gene_variant|MODIFIER|SEHA_RS01445|SEHA_RS01445|transcript|SEHA_RS01445|'
                   'protein_coding||c.-2172G>A|||||2172|'),
                  ('T|upstream_gene_variant|MODIFIER|can|SEHA_RS01455|transcript|SEHA_RS01455|'
                   'protein_coding||c.-710G>A|||||710|')]}
         ],
    ], columns=[
        'CHROM', 'POS', 'REF', 'ALT', 'INFO',
    ])


@pytest.fixture
def mock_snpeff_headers_invalid():
    class MockAnn():
        def __init__(self):
            pass

        def info_ids(self):
            return ['ANN']

        def get_info_field_info(self, field: str):
            if field == 'ANN':
                obj = type('', (), {})()
                obj.description = "invalid"
                return obj
            else:
                return None

    return MockAnn()


def test_parse_annotation_headers(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser, mock_snpeff_headers):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers)

    assert ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID',
            'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length',
            'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO'] == headers_list


def test_parse_annotation_headers_invalid(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                          mock_snpeff_headers_invalid):
    with pytest.raises(InvalidSnpEffVcfError) as execinfo:
        vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers_invalid)
    assert "Found 'ANN' in VCF information but description" in str(execinfo.value)


def test_parse_annotation_headers_no_annotation(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                mock_snpeff_headers_empty):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers_empty)

    assert [] == headers_list


def test_parse_annotation_entries_single(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                         mock_snpeff_headers, mock_vcf_df_with_ann_single: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_ann_single)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert [0] == list(ann_entries_df.index)
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu', '497/1359', 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[0])


def test_parse_annotation_entries_multiple_entries_single_sample(
        vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
        mock_snpeff_headers, mock_vcf_df_with_ann_multiple_entries_single_sample: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_ann_multiple_entries_single_sample)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 3 == len(ann_entries_df)
    ann_entries_df = ann_entries_df.sort_values(['original_index', 'ANN.Gene_ID'])
    assert [0, 0, 0] == list(ann_entries_df.index)
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu', '497/1359', 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[0])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
            'c.-856C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[1])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
            'c.-1941C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[2])


def test_parse_annotation_entries_multiple_entries_multiple_samples(
        vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
        mock_snpeff_headers, mock_vcf_df_with_ann_multiple_entries_multiple_samples: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_ann_multiple_entries_multiple_samples)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 6 == len(ann_entries_df)
    ann_entries_df = ann_entries_df.sort_values(['original_index', 'ANN.Gene_ID'])
    assert [0, 0, 0, 1, 1, 1] == list(ann_entries_df.index)
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu', '497/1359', 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[0])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
            'c.-856C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[1])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
            'c.-1941C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[2])

    assert ['T', 'upstream_gene_variant', 'MODIFIER', 'SEHA_RS01445', 'SEHA_RS01445', 'transcript', 'protein_coding',
            'c.-2172G>A', pd.NA, pd.NA, 'NC_011083:203200:C:T'] == list(ann_entries_df.iloc[3])
    assert ['T', 'upstream_gene_variant', 'MODIFIER', 'can', 'SEHA_RS01455', 'transcript', 'protein_coding',
            'c.-710G>A', pd.NA, pd.NA, 'NC_011083:203200:C:T'] == list(ann_entries_df.iloc[4])
    assert ['T', 'missense_variant', 'MODERATE', 'SEHA_RS01460', 'SEHA_RS01460', 'transcript', 'protein_coding',
            'c.602C>T', 'p.Thr201Met', '602/927', 'NC_011083:203200:C:T'] == list(ann_entries_df.iloc[5])


def test_parse_annotation_entries_no_annotation_single_sample(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                              mock_vcf_df_without_ann_single: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=[],
                                                                           vcf_df=mock_vcf_df_without_ann_single)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert [0] == list(ann_entries_df.index)
    assert {True} == set(ann_entries_df.drop('VARIANT_ID', axis='columns').iloc[0].isna())
    assert 'NC_011083:140658:C:A' == ann_entries_df['VARIANT_ID'].iloc[0]


def test_parse_annotation_entries_no_annotation_multiple_sample(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                                mock_vcf_df_without_ann_multiple: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=[],
                                                                           vcf_df=mock_vcf_df_without_ann_multiple)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 2 == len(ann_entries_df)
    assert [0, 1] == list(ann_entries_df.index)
    assert {True} == set(ann_entries_df.drop('VARIANT_ID', axis='columns').iloc[0].isna())
    assert 'NC_011083:140658:C:A' == ann_entries_df['VARIANT_ID'].iloc[0]
    assert {True} == set(ann_entries_df.drop('VARIANT_ID', axis='columns').iloc[1].isna())
    assert 'NC_011083:203200:C:T' == ann_entries_df['VARIANT_ID'].iloc[1]


def test_parse_annotation_entries_some_with_some_without_annotations(
        vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
        mock_snpeff_headers,
        mock_vcf_df_with_and_without_ann: pd.DataFrame):
    headers_list = vcf_snpeff_annotation_parser.parse_annotation_headers(mock_snpeff_headers)
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=headers_list,
                                                                           vcf_df=mock_vcf_df_with_and_without_ann)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 4 == len(ann_entries_df)
    ann_entries_df = ann_entries_df.sort_values(['original_index', 'ANN.Gene_ID'])
    assert [0, 0, 0, 1] == list(ann_entries_df.index)
    assert ['A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu', '497/1359', 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[0])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'mraY', 'SEHA_RS01185', 'transcript', 'protein_coding',
            'c.-856C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[1])
    assert ['A', 'upstream_gene_variant', 'MODIFIER', 'murD', 'SEHA_RS01190', 'transcript', 'protein_coding',
            'c.-1941C>A', pd.NA, pd.NA, 'NC_011083:140658:C:A'] == list(ann_entries_df.iloc[2])

    assert {True} == set(ann_entries_df.drop('VARIANT_ID', axis='columns').iloc[3].isna())
    assert 'NC_011083:203200:C:T' == ann_entries_df['VARIANT_ID'].iloc[3]


def test_parse_annotation_entries_empty(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                        mock_vcf_df_empty: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.parse_annotation_entries(vcf_ann_headers=[],
                                                                           vcf_df=mock_vcf_df_empty)

    assert ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p', 'ANN.cDNA.pos / cDNA.length',
            'VARIANT_ID'] == list(
        ann_entries_df.columns)
    assert 0 == len(ann_entries_df)


def test_select_variant_annotations_single_sample(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                  mock_vcf_df_single_sample_annotations: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(mock_vcf_df_single_sample_annotations)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])


def test_select_variant_annotations_some_na_values(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                   mock_vcf_df_single_sample_annotations_some_na_values: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(
        mock_vcf_df_single_sample_annotations_some_na_values)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'NA',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0].fillna('NA'))


def test_select_variant_annotations_single_sample_sars_cov_2(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                             mock_vcf_df_single_sample_annotations_sars_cov_2: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(
        mock_vcf_df_single_sample_annotations_sars_cov_2)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
            'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
            'c.3683A>G', 'p.D1228G'] == list(ann_entries_df.iloc[0])


def test_select_variant_annotations_rv_single_sample_sars_cov_2(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                                mock_vcf_df_single_sample_annotations_rv_sars_cov_2: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(
        mock_vcf_df_single_sample_annotations_rv_sars_cov_2)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 1 == len(ann_entries_df)
    assert ['SampleA', 'NC_045512.2', 3948, 'A', 'G', 'SNP', 'SampleA.vcf', 'NC_045512.2:3948:A:G',
            'G', 'missense_variant', 'MODERATE', 'ORF1ab', 'GU280_gp01', 'transcript', 'protein_coding',
            'c.3683A>G', 'p.D1228G'] == list(ann_entries_df.iloc[0])


def test_select_variant_annotations_multiple_sample(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                    mock_vcf_df_multiple_sample_annotations: pd.DataFrame):
    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(mock_vcf_df_multiple_sample_annotations)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 2 == len(ann_entries_df)
    assert ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])
    assert ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf', 'NC_011083:203200:C:T',
            'T', 'missense_variant', 'MODERATE', 'SEHA_RS01460', 'SEHA_RS01460', 'transcript', 'protein_coding',
            'c.602C>T', 'p.Thr201Met'] == list(ann_entries_df.iloc[1])


def test_select_variant_annotations_one_no_annotation(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                      mock_vcf_df_multiple_sample_one_empty: pd.DataFrame):
    annotation_columns = ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                          'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(mock_vcf_df_multiple_sample_one_empty)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 2 == len(ann_entries_df)
    assert ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])
    assert ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf',
            'NC_011083:203200:C:T'] == list(ann_entries_df.drop(annotation_columns, axis='columns').iloc[1])
    assert {True} == set(ann_entries_df[annotation_columns].iloc[1].isna().tolist())


def test_select_variant_annotations_one_invalid_annotation(vcf_snpeff_annotation_parser: VcfSnpEffAnnotationParser,
                                                           mock_vcf_df_multiple_sample_one_invalid: pd.DataFrame):
    annotation_columns = ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                          'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    ann_entries_df = vcf_snpeff_annotation_parser.select_variant_annotations(mock_vcf_df_multiple_sample_one_invalid)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
        ann_entries_df.columns)
    assert 2 == len(ann_entries_df)
    assert ['SampleA', 'NC_011083', 140658, 'C', 'A', 'snp', 'SampleA.vcf', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == list(ann_entries_df.iloc[0])
    assert ['SampleA', 'NC_011083', 203200, 'C', 'T', 'snp', 'SampleA.vcf',
            'NC_011083:203200:C:T'] == list(ann_entries_df.drop(annotation_columns, axis='columns').iloc[1])
    assert {True} == set(ann_entries_df[annotation_columns].iloc[1].isna().tolist())
