import pytest
import gzip
from pathlib import Path
from os import path

from Bio import SeqIO

from genomics_data_index.api.query.kind.isa.typing.ExperimentalSARSCov2ConstellationsTyper import ExperimentalSARSCov2ConstellationsTyper

root_data_dir = Path(path.dirname(__file__)) / '..' / '..' / '..' / 'data'
sequence_path = root_data_dir / 'NC_045512.2.gb.gz'


@pytest.fixture
def sequence_record():
    with gzip.open(sequence_path, mode='rt') as fh:
        sequences = list(SeqIO.parse(fh, 'genbank'))
        return sequences[0]


def test_mutation_to_identifier(sequence_record):
    typer = ExperimentalSARSCov2ConstellationsTyper(constellation_files=[],
                                                    sequence=sequence_record)

    assert 'hgvs_gn:NC_045512.2:S:p.H69_V70del' == typer.mutation_to_identifier('s:HV69-')
    assert 'hgvs_gn:NC_045512.2:ORF1ab:p.S3675_F3677del' == typer.mutation_to_identifier('1ab:SGF3675-')
    assert 'NC_045512.2:5986:C:T' == typer.mutation_to_identifier('nuc:C5986T')
    assert 'hgvs_gn:NC_045512.2:ORF8:p.Q27*' == typer.mutation_to_identifier('8:Q27*')
    assert 'hgvs_gn:NC_045512.2:S:p.Y144del' == typer.mutation_to_identifier('s:Y144-')

    assert 'NC_045512.2:3:AA:A' == typer.mutation_to_identifier('del:4:1')
    assert 'NC_045512.2:3:AAA:A' == typer.mutation_to_identifier('del:4:2')
    assert 'NC_045512.2:3:AAAG:A' == typer.mutation_to_identifier('del:4:3')
    assert 'NC_045512.2:20:CAGGTAACAAA:C' == typer.mutation_to_identifier('del:21:10')
