import math

import pytest

from storage.variant.model import Sample
from storage.variant.service import EntityExistsError
from storage.variant.service.VariationService import VariationService


def test_insert_variants(database, snippy_variants_reader, reference_service_with_data,
                         sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    variation_service.insert_variants(reference_name='genome', variants_reader=snippy_variants_reader)
    session = database.get_session()
    samples = session.query(Sample).all()

    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}

    variation_files = {v.nucleotide_variants_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(variation_files)


# def test_insert_variants(database, snippy_variants_reader, reference_service_with_data, sample_service):
#     variation_service = VariationService(database, reference_service_with_data, sample_service)
#
#     core_masks = snippy_variants_reader.get_core_masks()
#     var_df = snippy_variants_reader.get_variants_table()
#     session = database.get_session()
#
#     variation_service.insert_variants(var_df=var_df,
#                                       reference_name='genome',
#                                       core_masks=core_masks)
#
#     assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 129 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#     # Check to make sure some variation alleles are stored
#     v = session.query(VariationAllele).get('reference:1048:C:G')
#     assert v is not None, 'Particular variant does not exist'
#     assert 1048 == v.position, 'Position is incorrect'
#     assert 'C' == v.ref, 'Reference is incorrect'
#     assert 'G' == v.alt, 'Alt is incorrect'
#     assert 'reference' == v.sequence.sequence_name, 'Sequence name is incorrect'
#     assert 5180 == v.sequence.sequence_length, 'Sequence length is incorrect'
#     assert 'snp' == v.var_type, 'Type is incorrect'
#
#     v = session.query(VariationAllele).get('reference:1135:CCT:C')
#     assert v is not None, 'Particular variant does not exist'
#     assert 1135 == v.position, 'Position is incorrect'
#     assert 'CCT' == v.ref, 'Reference is incorrect'
#     assert 'C' == v.alt, 'Alt is incorrect'
#     assert 'reference' == v.sequence.sequence_name, 'Sequence name is incorrect'
#     assert 5180 == v.sequence.sequence_length, 'Sequence length is incorrect'
#     assert 'del' == v.var_type, 'Type is incorrect'
#
#
# def test_insert_variants_duplicates(database, snippy_variants_reader, reference_service_with_data, sample_service):
#     variation_service = VariationService(database, reference_service_with_data, sample_service)
#
#     core_masks = snippy_variants_reader.get_core_masks()
#     var_df = snippy_variants_reader.get_variants_table()
#     session = database.get_session()
#
#     assert 0 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 0 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 0 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#     variation_service.insert_variants(var_df=var_df,
#                                       reference_name='genome',
#                                       core_masks=core_masks)
#
#     assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 3 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 129 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#     with pytest.raises(EntityExistsError) as execinfo:
#         variation_service.insert_variants(var_df=var_df,
#                                           reference_name='genome',
#                                           core_masks=core_masks)
#     assert 'Passed samples already have variants for reference genome [genome]' in str(execinfo.value)
#     assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 3 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 129 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#
# def test_insert_variants_duplicates_subset(database, snippy_variants_reader, reference_service_with_data,
#                                            sample_service):
#     variation_service = VariationService(database, reference_service_with_data, sample_service)
#
#     core_masks = snippy_variants_reader.get_core_masks()
#     var_df = snippy_variants_reader.get_variants_table()
#     session = database.get_session()
#
#     assert 0 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 0 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 0 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#     variation_service.insert_variants(var_df=var_df,
#                                       reference_name='genome',
#                                       core_masks=core_masks)
#
#     assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 3 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 129 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#     # Select a subset of samples
#     var_df_2 = var_df[var_df['SAMPLE'].isin(['SampleA', 'SampleB'])]
#
#     with pytest.raises(EntityExistsError) as execinfo:
#         variation_service.insert_variants(var_df=var_df_2,
#                                           reference_name='genome',
#                                           core_masks=core_masks)
#     assert 'Passed samples already have variants for reference genome [genome]' in str(execinfo.value)
#     assert 112 == session.query(VariationAllele).count(), 'Incorrect number of variant entries'
#     assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
#     assert 3 == session.query(SampleSequence).count(), 'Incorrect number of SampleSequences'
#     assert 129 == session.query(sample_variation_association).count(), 'Incorrect number of sample variants'
#
#
# def test_get_variants(database, snippy_variants_reader, reference_service_with_data, sample_service):
#     variation_service = VariationService(database, reference_service_with_data, sample_service)
#
#     core_masks = snippy_variants_reader.get_core_masks()
#     var_df = snippy_variants_reader.get_variants_table()
#
#     variation_service.insert_variants(var_df=var_df,
#                                       reference_name='genome',
#                                       core_masks=core_masks)
#
#     variants = variation_service.get_variants(sequence_name='reference')
#
#     assert 60 == len(variants.keys()), 'Incorrect number of variants returned (counting only SNV/SNPs)'
#
#     assert 'reference:4265:G:C' == variants[4265]['SampleA'].id, 'Incorrect variant returned'
#     assert 'SampleB' not in variants[4265], 'Incorrect variant returned'
#     assert 'SampleC' not in variants[4265], 'Incorrect variant returned'
#
#     assert 'reference:839:C:G' == variants[839]['SampleB'].id, 'Incorrect variant returned'
#     assert 'reference:839:C:G' == variants[839]['SampleC'].id, 'Incorrect variant returned'
#     assert 'SampleA' not in variants[839], 'Incorrect variant returned'
#
#
# def test_pariwise_distance(database, snippy_variants_reader, reference_service_with_data, sample_service):
#     variation_service = VariationService(database, reference_service_with_data, sample_service)
#
#     core_masks = snippy_variants_reader.get_core_masks()
#     var_df = snippy_variants_reader.get_variants_table()
#
#     variation_service.insert_variants(var_df=var_df,
#                                       reference_name='genome',
#                                       core_masks=core_masks)
#
#     distances = variation_service.pairwise_distance(['SampleA', 'SampleC', 'SampleB'])
#
#     expected_distAB = 1
#     expected_distAC = 1
#     expected_distBC = (1 - 17 / 66)
#
#     assert math.isclose(expected_distAB, distances.loc['SampleA', 'SampleB']), 'Incorrect pairwise distance'
#     assert math.isclose(expected_distAC, distances.loc['SampleA', 'SampleC']), 'Incorrect pairwise distance'
#     assert math.isclose(expected_distBC, distances.loc['SampleB', 'SampleC']), 'Incorrect pairwise distance'
