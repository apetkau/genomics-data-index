from typing import List

import pytest

from storage.test.integration.variant import data_dir
from storage.variant.io.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.model import Sample, NucleotideVariantsSamples, SampleNucleotideVariation
from storage.variant.service import EntityExistsError
from storage.variant.service.VariationService import VariationService


def test_insert_variants_saved_files(database, snippy_variants_reader, reference_service_with_data,
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

    genomic_mask_files = {v.masked_regions_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(genomic_mask_files)


def test_insert_variants_masked_regions(database, snippy_variants_reader, reference_service_with_data,
                                        sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    session = database.get_session()

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    samples = session.query(Sample).all()

    genomic_masks = {s.name: v.masked_regions for s in samples for v in s.sample_nucleotide_variation}
    assert 437 == len(genomic_masks['SampleA'])
    assert 276 == len(genomic_masks['SampleB'])
    assert 329 == len(genomic_masks['SampleC'])


def test_insert_variants_examine_variation(database, snippy_variants_reader, reference_service_with_data,
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


def test_get_variants_ordered(database, snippy_variants_reader, reference_service_with_data,
                              sample_service, filesystem_storage):
    def find_variant_by_position(variants: List[NucleotideVariantsSamples], position: int) -> NucleotideVariantsSamples:
        for v in variants:
            if v.position == position:
                return v

        raise Exception(f'Could not find variant with position {position}')

    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)

    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    variants = variation_service.get_variants_ordered(sequence_name='reference')

    assert 60 == len(variants), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    v1 = find_variant_by_position(variants, 4265)
    assert 'reference:4265:G:C' == v1.spdi, 'Incorrect variant returned'
    assert {sampleA.id} == set(v1.sample_ids)

    v2 = find_variant_by_position(variants, 839)
    assert 'reference:839:C:G' == v2.spdi, 'Incorrect variant returned'
    assert {sampleB.id, sampleC.id} == set(v2.sample_ids)


def test_get_sample_nucleotide_variation_one_sample(variation_service):
    sample_variations = variation_service.get_sample_nucleotide_variation(['SampleA'])

    assert 1 == len(sample_variations)
    assert isinstance(sample_variations[0], SampleNucleotideVariation)
    assert 'SampleA' == sample_variations[0].sample.name


def test_get_sample_nucleotide_variation_all_samples(variation_service):
    sample_variations = variation_service.get_sample_nucleotide_variation(['SampleA', 'SampleB', 'SampleC'])

    assert 3 == len(sample_variations)
    assert {'SampleA', 'SampleB', 'SampleC'} == {v.sample.name for v in sample_variations}
