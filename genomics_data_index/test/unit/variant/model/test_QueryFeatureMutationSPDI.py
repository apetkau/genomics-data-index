from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI


def test_create():
    f = QueryFeatureMutationSPDI('ref:1:A:T')
    assert 'ref:1:A:T' == f.id
    assert 'ref:1:A:T' == f.id_no_prefix


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


def test_to_unknown_explode():
    f = QueryFeatureMutationSPDI('ref:1:A:T')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:1:A:?']

    f = QueryFeatureMutationSPDI('ref:1:1:T')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:1:1:?']

    f = QueryFeatureMutationSPDI('ref:10:AT:T')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:10:A:?', 'ref:11:T:?']

    f = QueryFeatureMutationSPDI('ref:10:T:TC')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:10:T:?']

    f = QueryFeatureMutationSPDI('ref:10:ATC:TCC')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:10:A:?', 'ref:11:T:?', 'ref:12:C:?']

    f = QueryFeatureMutationSPDI('ref:10:3:TCC')
    assert [u.id for u in f.to_unknown_explode()] == ['ref:10:1:?', 'ref:11:1:?', 'ref:12:1:?']


def test_has_deletion_sequence():
    assert QueryFeatureMutationSPDI('ref:1:A:T').has_deletion_sequence()
    assert QueryFeatureMutationSPDI('ref:10:ATT:T').has_deletion_sequence()
    assert not QueryFeatureMutationSPDI('ref:10:1:T').has_deletion_sequence()
    assert not QueryFeatureMutationSPDI('ref:10:3:T').has_deletion_sequence()


def test_deletion_length():
    assert 1 == QueryFeatureMutationSPDI('ref:10:A:T').deletion_length()
    assert 3 == QueryFeatureMutationSPDI('ref:10:ATT:T').deletion_length()
    assert 1 == QueryFeatureMutationSPDI('ref:10:1:T').deletion_length()
    assert 3 == QueryFeatureMutationSPDI('ref:10:3:T').deletion_length()
