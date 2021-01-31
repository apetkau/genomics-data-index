import math

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
    reference_service = ReferenceService(database)
    reference_service.create_reference_genome(reference_file)
    reference = reference_service.find_reference_genome('genome')

    ref_contigs = {s.sequence_name: s for s in reference.sequences}

    return {
        'database': database,
        'reader': SnippyVariantsReader(sample_dirs),
        'ref_contigs': ref_contigs,
        'variation_service': VariationService(database)
    }


def test_insert_variants(setup):
    core_masks = setup['reader'].get_core_masks()
    var_df = setup['reader'].get_variants_table()
    ref_contigs = setup['ref_contigs']
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


def test_pariwise_distance(setup):
    core_masks = setup['reader'].get_core_masks()
    var_df = setup['reader'].get_variants_table()
    ref_contigs = setup['ref_contigs']
    variation_service = setup['variation_service']

    variation_service.insert_variants(var_df=var_df,
                                      ref_contigs=ref_contigs,
                                      core_masks=core_masks)

    distances = variation_service.pairwise_distance(['SampleA', 'SampleC', 'SampleB'])

    expected_distAB = 1
    expected_distAC = 1
    expected_distBC = (1 - 17 / 66)

    assert math.isclose(expected_distAB, distances.loc['SampleA', 'SampleB']), 'Incorrect pairwise distance'
    assert math.isclose(expected_distAC, distances.loc['SampleA', 'SampleC']), 'Incorrect pairwise distance'
    assert math.isclose(expected_distBC, distances.loc['SampleB', 'SampleC']), 'Incorrect pairwise distance'
