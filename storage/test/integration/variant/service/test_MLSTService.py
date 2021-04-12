import pytest

from storage.variant.model import MLSTAllelesSamples
from storage.variant.model import Sample
from storage.variant.model import SampleMLSTAlleles
from storage.variant.service.MLSTService import MLSTService
from storage.variant.service import EntityExistsError


def test_insert_mlst_results(database, mlst_reader_single_scheme, sample_service, filesystem_storage):
    num_loci = 7

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', features_reader=mlst_reader_single_scheme)

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


def test_insert_mlst_results_multiple_schemes(database, mlst_reader_basic, sample_service, filesystem_storage):
    num_loci = 7
    num_schemes = 3

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(features_reader=mlst_reader_basic)

    samples = session.query(Sample).all()

    assert 6 == len(samples)
    assert {'CFSAN002349', 'CFSAN023463', '2014C-3598', '2014C-3599',
            '2014D-0067', '2014D-0068'} == {s.name for s in samples}

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 6 == len(sample_mlst_alleles)

    mlst_alleles_all = session.query(MLSTAllelesSamples).all()
    mlst_alleles_all_id = {allele.sla: allele for allele in mlst_alleles_all}
    assert num_loci * num_schemes == len(mlst_alleles_all)

    assert {'lmonocytogenes', 'ecoli', 'campylobacter'} == {allele.scheme for allele in mlst_alleles_all}

    assert 2 == len(mlst_alleles_all_id['lmonocytogenes:abcZ:1'].sample_ids)
    assert 2 == len(mlst_alleles_all_id['ecoli:fumC:23'].sample_ids)

    mlst_alleles_lmono = session.query(MLSTAllelesSamples).filter(MLSTAllelesSamples.scheme == 'lmonocytogenes').all()
    mlst_alleles_id_lmono = {allele.sla: allele for allele in mlst_alleles_lmono}

    assert {'lmonocytogenes:abcZ:1', 'lmonocytogenes:bglA:51', 'lmonocytogenes:cat:11',
            'lmonocytogenes:dapE:13', 'lmonocytogenes:dat:2', 'lmonocytogenes:ldh:5',
            'lmonocytogenes:lhkA:5'} == set(mlst_alleles_id_lmono.keys())


def test_insert_mlst_results_multiple_schemes_override_scheme(database, mlst_reader_basic,
                                                              sample_service, filesystem_storage):
    num_loci = 7
    num_alt_schemes = 3

    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', features_reader=mlst_reader_basic)

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 6 == len(sample_mlst_alleles)

    mlst_alleles_samples = session.query(MLSTAllelesSamples).all()
    mlst_alleles_samples_id_allele = {allele.sla: allele for allele in mlst_alleles_samples}

    assert num_loci * num_alt_schemes == len(mlst_alleles_samples)
    assert 2 == len(mlst_alleles_samples_id_allele['lmonocytogenes:abcZ:1'].sample_ids)
    assert 2 == len(mlst_alleles_samples_id_allele['lmonocytogenes:fumC:23'].sample_ids)


def test_double_insert_mlst(database, mlst_reader_single_scheme, sample_service, filesystem_storage):
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)

    session = database.get_session()

    mlst_service.insert(feature_scope_name='lmonocytogenes', features_reader=mlst_reader_single_scheme)

    samples = session.query(Sample).all()

    assert 2 == len(samples)
    assert {'CFSAN002349', 'CFSAN023463'} == {s.name for s in samples}

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 2 == len(sample_mlst_alleles)

    with pytest.raises(EntityExistsError) as execinfo:
        mlst_service.insert(feature_scope_name='lmonocytogenes', features_reader=mlst_reader_single_scheme)

    assert 'Passed samples already have features for feature scope [lmonocytogenes]' in str(execinfo.value)

    sample_mlst_alleles = session.query(SampleMLSTAlleles).all()
    assert 2 == len(sample_mlst_alleles)
