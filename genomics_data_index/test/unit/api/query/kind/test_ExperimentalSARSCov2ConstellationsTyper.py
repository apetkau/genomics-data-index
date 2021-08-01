from genomics_data_index.api.query.kind.isa.typing.ExperimentalSARSCov2ConstellationsTyper import ExperimentalSARSCov2ConstellationsTyper


def test_mutation_to_identifier():
    typer = ExperimentalSARSCov2ConstellationsTyper(constellation_files=[],
                                                    sequence_name='MN908947.3')

    assert 'hgvs_gn:MN908947.3:S:p.H69_V70del' == typer.mutation_to_identifier('s:HV69-')
    assert 'hgvs_gn:MN908947.3:orf1ab:p.S3675_F3677del' == typer.mutation_to_identifier('1ab:SGF3675-')
    assert 'MN908947.3:5986:C:T' == typer.mutation_to_identifier('nuc:C5986T')
    assert 'hgvs_gn:MN908947.3:ORF8:p.Q27*' == typer.mutation_to_identifier('8:Q27*')
    assert 'hgvs_gn:MN908947.3:S:p.Y144del' == typer.mutation_to_identifier('s:Y144-')
