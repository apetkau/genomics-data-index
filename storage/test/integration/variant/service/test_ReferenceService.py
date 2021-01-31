from os import path
from pathlib import Path
from typing import Dict, Any

import pytest

from storage.variant.model import Reference
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService

data_dir = Path(path.dirname(__file__), '..', '..', 'data', 'snippy')
reference_file = data_dir / 'genome.fasta.gz'


@pytest.fixture
def setup() -> Dict[str, Any]:
    database = DatabaseConnection('sqlite:///:memory:')
    reference_service = ReferenceService(database)

    return dict(database=database, reference_service=reference_service)


def test_insert_reference_genome(setup):
    assert 0 == setup['database'].get_session().query(Reference).count(), 'Database should be empty initially'
    setup['reference_service'].create_reference_genome(reference_file)
    assert 1 == setup['database'].get_session().query(Reference).count(), 'Database should have one entry'
    assert 'genome' == setup['database'].get_session().query(Reference).all()[0].name, 'Name should match'


def test_find_reference_genome(setup):
    assert 0 == setup['database'].get_session().query(Reference).count(), 'Database should be empty initially'
    setup['reference_service'].create_reference_genome(reference_file)

    reference = setup['reference_service'].find_reference_genome('genome')
    assert 'genome' == reference.name, 'Reference name should match'
    assert 5180 == reference.length, 'Reference length should match'
