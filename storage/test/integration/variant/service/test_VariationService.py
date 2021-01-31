from os import path, listdir
from pathlib import Path
from typing import Dict, Any

import pytest

from storage.variant.VariantsReader import SnippyVariantsReader
from storage.variant.model import VariationAllele
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService

data_dir = Path(path.dirname(__file__), '..', '..', 'data', 'snippy')
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
reference_file = data_dir / 'genome.fasta.gz'


@pytest.fixture
def setup() -> Dict[str, Any]:
    database = DatabaseConnection('sqlite:///:memory:')

    return {
        'database': database,
        'reader': SnippyVariantsReader(sample_dirs),
        'reference_service': ReferenceService(database),
        'variation_service': VariationService(database)
    }


def test_insert_variants(setup):
    setup['reference_service'].create_reference_genome(reference_file)
    reference = setup['reference_service'].find_reference_genome('genome')

    ref_contigs = {s.sequence_name: s for s in reference.sequences}

    core_masks = setup['reader'].get_core_masks()
    var_df = setup['reader'].get_variants_table()

    session = setup['database'].get_session()

    setup['variation_service'].insert_variants(var_df=var_df,
                                               ref_contigs=ref_contigs,
                                               core_masks=core_masks)

    assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'

    # Check to make sure some variation alleles are stored
    v = session.query(VariationAllele).get('reference:1048:C:G')
    assert v is not None, 'Particular variant does not exist'
    assert 1048 == v.position, 'Position is incorrect'
    assert 'C' == v.ref, 'Reference is incorrect'
    assert 'G' == v.alt, 'Alt is incorrect'
    assert 'reference' == v.sequence.sequence_name, 'Sequence name is incorrect'
    assert 5180 == v.sequence.sequence_length, 'Sequence length is incorrect'
    assert 'snp' == v.var_type, 'Type is incorrect'

    v = session.query(VariationAllele).get('reference:1135:CCT:C')
    assert v is not None, 'Particular variant does not exist'
    assert 1135 == v.position, 'Position is incorrect'
    assert 'CCT' == v.ref, 'Reference is incorrect'
    assert 'C' == v.alt, 'Alt is incorrect'
    assert 'reference' == v.sequence.sequence_name, 'Sequence name is incorrect'
    assert 5180 == v.sequence.sequence_length, 'Sequence length is incorrect'
    assert 'del' == v.var_type, 'Type is incorrect'
