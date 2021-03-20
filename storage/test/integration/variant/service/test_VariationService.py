import math
import pytest

from storage.variant.model import Sample, NucleotideVariantsSamples, SampleNucleotideVariation
from storage.variant.service import EntityExistsError
from storage.variant.service.VariationService import VariationService
from storage.variant.io.SnippyVariantsReader import SnippyVariantsReader

from storage.test.integration.variant import data_dir


def test_insert_variants_tmp(database, snippy_variants_reader, reference_service_with_data,
                         sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    session = database.get_session()

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    samples = session.query(Sample).all()

    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}

    variation_files = {v.nucleotide_variants_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(variation_files)


def test_insert_variants(database, snippy_variants_reader, reference_service_with_data,
                         sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    sample_name_ids = sample_service.find_sample_name_ids({'SampleA', 'SampleB', 'SampleC'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'

    # Check to make sure some variants are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1048,
        'deletion': 'C',
        'insertion': 'G'
    })
    assert v is not None, 'Particular variant does not exist'
    assert 'snp' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleA']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1135,
        'deletion': 'CCT',
        'insertion': 'C'
    })
    assert 'del' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB'], sample_name_ids['SampleC']} == set(v.sample_ids)


def test_insert_variants_duplicates(database, snippy_variants_reader, reference_service_with_data,
                                    sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    assert 'Passed samples already have variants for reference genome [genome]' in str(execinfo.value)
    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'


def test_insert_variants_duplicates_subset(database, snippy_variants_reader, reference_service_with_data,
                                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    # Select a subset of samples
    sample_dirs_subset = [data_dir / 'SampleA']
    snippy_variants_reader_subset = SnippyVariantsReader(sample_dirs_subset)

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader_subset)

    assert 'Passed samples already have variants for reference genome [genome]' in str(execinfo.value)
    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'


@pytest.mark.skip()
def test_get_variants(database, snippy_variants_reader, reference_service_with_data,
                      sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    variants = variation_service.get_variants(sequence_name='reference')

    assert 60 == len(variants.keys()), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    assert 'reference:4265:G:C' == variants[4265]['SampleA'].id, 'Incorrect variant returned'
    assert 'SampleB' not in variants[4265], 'Incorrect variant returned'
    assert 'SampleC' not in variants[4265], 'Incorrect variant returned'

    assert 'reference:839:C:G' == variants[839]['SampleB'].id, 'Incorrect variant returned'
    assert 'reference:839:C:G' == variants[839]['SampleC'].id, 'Incorrect variant returned'
    assert 'SampleA' not in variants[839], 'Incorrect variant returned'


@pytest.mark.skip()
def test_pariwise_distance(database, snippy_variants_reader, reference_service_with_data,
                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    distances = variation_service.pairwise_distance(['SampleA', 'SampleC', 'SampleB'])

    expected_distAB = 1
    expected_distAC = 1
    expected_distBC = (1 - 17 / 66)

    assert math.isclose(expected_distAB, distances.loc['SampleA', 'SampleB']), 'Incorrect pairwise distance'
    assert math.isclose(expected_distAC, distances.loc['SampleA', 'SampleC']), 'Incorrect pairwise distance'
    assert math.isclose(expected_distBC, distances.loc['SampleB', 'SampleC']), 'Incorrect pairwise distance'
