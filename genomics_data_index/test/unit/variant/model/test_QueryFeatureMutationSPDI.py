from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI


def test_to_unknown():
    f = QueryFeatureMutationSPDI('ref:1:A:T')
    assert f.to_unknown().id == 'ref:1:A:?'

    f = QueryFeatureMutationSPDI('ref:1:1:T')
    assert f.to_unknown().id == 'ref:1:1:?'

    f = QueryFeatureMutationSPDI('ref:10:AT:C')
    assert f.to_unknown().id == 'ref:10:AT:?'

    f = QueryFeatureMutationSPDI('ref:1:A:TT')
    assert f.to_unknown().id == 'ref:1:A:?'

    f = QueryFeatureMutationSPDI('ref:1:AT:TCG')
    assert f.to_unknown().id == 'ref:1:AT:?'

    f = QueryFeatureMutationSPDI('ref:1:2:TCG')
    assert f.to_unknown().id == 'ref:1:2:?'
