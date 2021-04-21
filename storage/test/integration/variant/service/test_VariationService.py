from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

import pytest
import pandas as pd
from sqlalchemy.orm.exc import NoResultFound

from storage.test.integration import data_dir, snippy_snps_dataframes, snippy_all_dataframes
from storage.variant.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from storage.variant.model.db import NucleotideVariantsSamples, SampleNucleotideVariation, Sample
from storage.variant.service import EntityExistsError
from storage.variant.service.VariationService import VariationService


def test_insert_variants_saved_files(database, snippy_nucleotide_data_package, reference_service_with_data,
                                     sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    samples = session.query(Sample).all()

    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}

    variation_files = {v.nucleotide_variants_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(variation_files)

    genomic_mask_files = {v.masked_regions_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(genomic_mask_files)


def test_insert_variants_masked_regions(database, snippy_nucleotide_data_package, reference_service_with_data,
                                        sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)

    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    samples = session.query(Sample).all()

    genomic_masks = {s.name: v.masked_regions for s in samples for v in s.sample_nucleotide_variation}
    assert 437 == len(genomic_masks['SampleA'])
    assert 276 == len(genomic_masks['SampleB'])
    assert 329 == len(genomic_masks['SampleC'])


def test_summarize_variants(database, snippy_nucleotide_data_package, reference_service_with_data,
                                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == variation_service.count_on_reference('genome', include_unknown=False)

    with pytest.raises(NoResultFound) as execinfo:
        variation_service.count_on_reference('no_exists', include_unknown=False)
    assert 'No row was found' in str(execinfo.value)

    mutation_counts = variation_service.mutation_counts_on_reference('genome', include_unknown=False)
    assert 112 == len(mutation_counts)
    assert 2 == mutation_counts['reference:839:1:G']
    assert 1 == mutation_counts['reference:866:9:G']
    assert 1 == mutation_counts['reference:1048:1:G']
    assert 2 == mutation_counts['reference:3897:5:G']


def test_insert_variants_examine_variation(database, snippy_nucleotide_data_package, reference_service_with_data,
                                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    sample_name_ids = sample_service.find_sample_name_ids({'SampleA', 'SampleB', 'SampleC'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'

    # Check to make sure some variants are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1048,
        'deletion': len('C'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular variant does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleA']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1135,
        'deletion': len('CCT'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB'], sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 883,
        'deletion': len('CACATG'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1984,
        'deletion': len('GTGATTG'),
        'insertion': 'TTGA'
    })
    assert 'OTHER' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)


def test_insert_variants_regular_vcf_reader_examine_variation(database, regular_nucleotide_data_package,
                                                              reference_service_with_data,
                                                              sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=regular_nucleotide_data_package)

    sample_name_ids = sample_service.find_sample_name_ids({'SampleA', 'SampleB', 'SampleC'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'

    # Check to make sure some variants are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1048,
        'deletion': len('C'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular variant does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleA']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1135,
        'deletion': len('CCT'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB'], sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 883,
        'deletion': len('CACATG'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1984,
        'deletion': len('GTGATTG'),
        'insertion': 'TTGA'
    })
    assert 'OTHER' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)


def test_insert_variants_duplicates(database, snippy_nucleotide_data_package, reference_service_with_data,
                                    sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 'Passed samples already have features for feature scope [genome]' in str(execinfo.value)
    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'


def test_insert_variants_duplicates_subset(database, snippy_nucleotide_data_package, reference_service_with_data,
                                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    # Select a subset of samples
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        sample_dirs_subset = [data_dir / 'SampleA']
        subset_snippy_data_package = NucleotideSampleDataPackage.create_from_snippy(
            sample_dirs_subset,
            sample_files_processor=SerialSampleFilesProcessor(tmp_dir)
        )

        with pytest.raises(EntityExistsError) as execinfo:
            variation_service.insert(feature_scope_name='genome', data_package=subset_snippy_data_package)

        assert 'Passed samples already have features for feature scope [genome]' in str(execinfo.value)
        assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of variant entries'
        assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
        assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'


def test_get_variants_ordered(database, snippy_nucleotide_data_package, reference_service_with_data,
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

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    variants = variation_service.get_variants_ordered(sequence_name='reference')

    assert 60 == len(variants), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    v1 = find_variant_by_position(variants, 4265)
    assert 'reference:4265:1:C' == v1.spdi, 'Incorrect variant returned'
    assert {sampleA.id} == set(v1.sample_ids)

    v2 = find_variant_by_position(variants, 839)
    assert 'reference:839:1:G' == v2.spdi, 'Incorrect variant returned'
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


def test_count_mutations_in_sample_ids_one_sample(database, variation_service: VariationService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()

    expected_df = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t').groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    print('Expected')
    print(expected_df)

    mutations_df = variation_service.count_mutations_in_sample_ids_dataframe([sampleA.id])
    mutations_df = mutations_df.sort_index()
    print('Actual')
    print(mutations_df)

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    print(list(mutations_df.index))
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Position']) == list(mutations_df['Position'])
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Insertion']) == list(mutations_df['Insertion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])


def test_count_mutations_in_sample_ids_three_samples_only_snps(database, variation_service: VariationService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    dfA = pd.read_csv(snippy_snps_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_snps_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_snps_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfA, dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()

    mutations_df = variation_service.count_mutations_in_sample_ids_dataframe([sampleA.id, sampleB.id, sampleC.id],
                                                                             mutation_type='snp')
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    print(list(mutations_df.index))
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Position']) == list(mutations_df['Position'])
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Insertion']) == list(mutations_df['Insertion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])


def test_count_mutations_in_sample_ids_three_samples_all_mutations(database, variation_service: VariationService):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfA, dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()

    mutations_df = variation_service.count_mutations_in_sample_ids_dataframe([sampleA.id, sampleB.id, sampleC.id],
                                                                             )
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    print(list(mutations_df.index))
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Position']) == list(mutations_df['Position'])
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Insertion']) == list(mutations_df['Insertion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])
