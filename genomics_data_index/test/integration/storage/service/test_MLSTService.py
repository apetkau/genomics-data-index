import pytest

from genomics_data_index.storage.model.db import SampleMLSTAlleles, MLSTAllelesSamples, Sample
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.MLSTService import MLSTService


def test_insert_mlst_results(database, mlst_data_package_single_scheme, sample_service, filesystem_storage):
    num_loci = 7

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_single_scheme)

    samples = session.query(Sample).all()

    assert 2 == len(samples)
    assert {'CFSAN002349', 'CFSAN023463'} == {s.name for s in samples}

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 2 == len(sample_mlst_alleles)

    mlst_alleles_samples = session.query(MLSTAllelesSamples).all()
    mlst_alleles_samples_id_allele = {allele.sla: allele for allele in mlst_alleles_samples}

    assert num_loci == len(mlst_alleles_samples)
    assert {'lmonocytogenes:abcZ:1', 'lmonocytogenes:bglA:51', 'lmonocytogenes:cat:11',
            'lmonocytogenes:dapE:13', 'lmonocytogenes:dat:2', 'lmonocytogenes:ldh:5',
            'lmonocytogenes:lhkA:5'} == set(mlst_alleles_samples_id_allele.keys())

    assert 2 == len(mlst_alleles_samples_id_allele['lmonocytogenes:abcZ:1'].sample_ids)


def test_insert_mlst_results_multiple_schemes(database, mlst_data_package_basic, sample_service, filesystem_storage):
    num_loci = 7
    num_schemes = 3

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(data_package=mlst_data_package_basic)

    samples = session.query(Sample).all()

    assert 6 == len(samples)
    assert {'CFSAN002349', 'CFSAN023463', '2014C-3598', '2014C-3599',
            '2014D-0067', '2014D-0068'} == {s.name for s in samples}

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 6 == len(sample_mlst_alleles)

    mlst_alleles_all = session.query(MLSTAllelesSamples).all()
    mlst_alleles_all_id = {allele.sla: allele for allele in mlst_alleles_all}
    assert num_loci * num_schemes + 1 == len(mlst_alleles_all)

    assert {'lmonocytogenes', 'ecoli', 'campylobacter'} == {allele.scheme for allele in mlst_alleles_all}

    assert 2 == len(mlst_alleles_all_id['lmonocytogenes:abcZ:1'].sample_ids)
    assert 2 == len(mlst_alleles_all_id['ecoli:fumC:23'].sample_ids)

    mlst_alleles_lmono = session.query(MLSTAllelesSamples).filter(MLSTAllelesSamples.scheme == 'lmonocytogenes').all()
    mlst_alleles_id_lmono = {allele.sla for allele in mlst_alleles_lmono}

    assert {'lmonocytogenes:abcZ:1', 'lmonocytogenes:bglA:51', 'lmonocytogenes:cat:11',
            'lmonocytogenes:dapE:13', 'lmonocytogenes:dat:2', 'lmonocytogenes:ldh:5',
            'lmonocytogenes:lhkA:4', 'lmonocytogenes:lhkA:5'} == mlst_alleles_id_lmono


def test_insert_mlst_results_multiple_schemes_override_scheme(database, mlst_data_package_basic,
                                                              sample_service, filesystem_storage):
    num_loci = 7
    num_alt_schemes = 3

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_basic)

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 6 == len(sample_mlst_alleles)

    mlst_alleles_samples = session.query(MLSTAllelesSamples).all()
    mlst_alleles_samples_id_allele = {allele.sla: allele for allele in mlst_alleles_samples}

    assert num_loci * num_alt_schemes + 1 == len(mlst_alleles_samples)
    assert 2 == len(mlst_alleles_samples_id_allele['lmonocytogenes:abcZ:1'].sample_ids)
    assert 2 == len(mlst_alleles_samples_id_allele['lmonocytogenes:fumC:23'].sample_ids)


def test_double_insert_mlst(database, mlst_data_package_single_scheme, sample_service, filesystem_storage):
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_single_scheme)

    samples = session.query(Sample).all()

    assert 2 == len(samples)
    assert {'CFSAN002349', 'CFSAN023463'} == {s.name for s in samples}

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 2 == len(sample_mlst_alleles)

    with pytest.raises(EntityExistsError) as execinfo:
        mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_single_scheme)

    assert 'Passed samples already have features for feature scope [lmonocytogenes]' in str(execinfo.value)

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 2 == len(sample_mlst_alleles)


def test_multiple_insert_different_samples_same_scheme(database, mlst_data_package_single_scheme,
                                                       mlst_data_package_single_scheme2,
                                                       sample_service, filesystem_storage):
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    # Insert first set of data
    mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_single_scheme)

    sample_CFSAN002349 = session.query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_CFSAN023463 = session.query(Sample).filter(Sample.name == 'CFSAN023463').one()

    ## Test feature where all new samples get added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'cat',
        'allele': '11',
    })
    assert m.id == 'lmonocytogenes:cat:11'
    assert m.sla == 'lmonocytogenes:cat:11'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id}

    ## Test feature where only one of the new samples gets added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'abcZ',
        'allele': '1',
    })
    assert m.id == 'lmonocytogenes:abcZ:1'
    assert m.sla == 'lmonocytogenes:abcZ:1'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id}

    ## Test feature where none of the new samples gets added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'dapE',
        'allele': '13',
    })
    assert m.id == 'lmonocytogenes:dapE:13'
    assert m.sla == 'lmonocytogenes:dapE:13'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id}

    ## Test feature that only exists after the 2nd addition
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'dapE',
        'allele': '12',
    })
    assert m is None

    # Insert second set of data
    mlst_service.insert(feature_scope_name='lmonocytogenes', data_package=mlst_data_package_single_scheme2)

    sample_CFSAN002349_2 = session.query(Sample).filter(Sample.name == 'CFSAN002349-2').one()
    sample_CFSAN023463_2 = session.query(Sample).filter(Sample.name == 'CFSAN023463-2').one()

    ## Test feature where all new samples get added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'cat',
        'allele': '11',
    })
    assert m.id == 'lmonocytogenes:cat:11'
    assert m.sla == 'lmonocytogenes:cat:11'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id,
                                 sample_CFSAN002349_2.id, sample_CFSAN023463_2.id}

    ## Test feature where only one of the new samples gets added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'abcZ',
        'allele': '1',
    })
    assert m.id == 'lmonocytogenes:abcZ:1'
    assert m.sla == 'lmonocytogenes:abcZ:1'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id,
                                 sample_CFSAN002349_2.id}

    ## Test feature where none of the new samples gets added
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'dapE',
        'allele': '13',
    })
    assert m.id == 'lmonocytogenes:dapE:13'
    assert m.sla == 'lmonocytogenes:dapE:13'
    assert set(m.sample_ids) == {sample_CFSAN002349.id, sample_CFSAN023463.id}

    ## Test feature that only exists after the 2nd addition
    m: MLSTAllelesSamples = session.query(MLSTAllelesSamples).get({
        'scheme': 'lmonocytogenes',
        'locus': 'dapE',
        'allele': '12',
    })
    assert m.id == 'lmonocytogenes:dapE:12'
    assert m.sla == 'lmonocytogenes:dapE:12'
    assert set(m.sample_ids) == {sample_CFSAN002349_2.id, sample_CFSAN023463_2.id}


def test_get_all_alleles(mlst_service_loaded: MLSTService):
    assert {'1'} == mlst_service_loaded.get_all_alleles('lmonocytogenes', 'abcZ')


def test_get_all_alleles_multiple_alleles(mlst_service_loaded: MLSTService):
    assert {'4', '5'} == mlst_service_loaded.get_all_alleles('lmonocytogenes', 'lhkA')


def test_get_all_alleles_unknown_alleles(mlst_service_loaded_unknown: MLSTService):
    assert {'1', '?'} == mlst_service_loaded_unknown.get_all_alleles('lmonocytogenes', 'abcZ')


def test_get_all_alleles_unknown_alleles2(mlst_service_loaded_unknown: MLSTService):
    assert {'?'} == mlst_service_loaded_unknown.get_all_alleles('lmonocytogenes', 'bglA')


def test_get_all_loci_alleles(mlst_service_loaded: MLSTService):
    assert {('abcZ', '1'), ('bglA', '51'), ('cat', '11'),
            ('dapE', '13'), ('dat', '2'), ('ldh', '5'),
            ('lhkA', '4'), ('lhkA', '5')} == mlst_service_loaded.get_all_loci_alleles('lmonocytogenes')


def test_get_all_loci_alleles_unknown(mlst_service_loaded_unknown: MLSTService):
    assert {('abcZ', '1'), ('abcZ', '?'), ('bglA', '?'), ('cat', '11'),
            ('dapE', '13'), ('dat', '2'), ('ldh', '5'),
            ('lhkA', '5')} == mlst_service_loaded_unknown.get_all_loci_alleles('lmonocytogenes')
