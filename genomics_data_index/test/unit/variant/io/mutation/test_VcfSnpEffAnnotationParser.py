import pytest

from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser, InvalidSnpEffVcfError


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
