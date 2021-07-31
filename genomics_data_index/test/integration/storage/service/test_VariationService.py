from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

import pandas as pd
import pytest

from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.model.db import NucleotideVariantsSamples, SampleNucleotideVariation, Sample
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.VariationService import VariationService
from genomics_data_index.test.integration import data_dir, snippy_snps_dataframes, snippy_all_dataframes
from genomics_data_index.test.integration import reference_file_snpeff


def test_insert_variants_saved_files(database, snippy_nucleotide_data_package, reference_service_with_data,
                                     sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=50)

    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    samples = session.query(Sample).all()

    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}

    variation_files = {v.nucleotide_variants_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(variation_files)

    genomic_mask_files = {v.masked_regions_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(genomic_mask_files)


def test_insert_variants_saved_files_set_no_index_package(database, snippy_nucleotide_data_package_no_index_missing,
                                                          reference_service_with_data,
                                                          sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)

    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package_no_index_missing)

    samples = session.query(Sample).all()

    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}

    variation_files = {v.nucleotide_variants_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(variation_files)

    genomic_mask_files = {v.masked_regions_file for s in samples for v in s.sample_nucleotide_variation}
    assert 3 == len(genomic_mask_files)


def test_read_index(database, snippy_nucleotide_data_package, reference_service_with_data,
                    sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    session = database.get_session()
    sampleB = session.query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = session.query(Sample).filter(Sample.name == 'SampleC').one()

    # Read 2 that do exist
    feature_samples = variation_service.read_index(['reference:839:1:G', 'reference:866:9:G'])
    assert len(feature_samples) == 2
    assert {'reference:839:1:G', 'reference:866:9:G'} == set(feature_samples.keys())
    assert isinstance(feature_samples['reference:839:1:G'], NucleotideVariantsSamples)
    assert isinstance(feature_samples['reference:866:9:G'], NucleotideVariantsSamples)
    assert {sampleB.id, sampleC.id} == set(feature_samples['reference:839:1:G'].sample_ids)
    assert {sampleC.id} == set(feature_samples['reference:866:9:G'].sample_ids)

    # Read empty list
    feature_samples = variation_service.read_index([])
    assert len(feature_samples) == 0

    # Read list with additional non-existing IDs
    feature_samples = variation_service.read_index(['reference:839:1:G', 'reference:866:9:G', 'not_exist:100:A:T'])
    assert len(feature_samples) == 2
    assert {'reference:839:1:G', 'reference:866:9:G'} == set(feature_samples.keys())


def test_count_on_reference_no_index_unknowns(database, snippy_nucleotide_data_package, reference_service_with_data,
                                              sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == variation_service.count_on_reference(reference_name='genome', include_unknown=False)
    assert 112 == variation_service.count_on_reference(reference_name='genome', include_unknown=True)

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.count_on_reference('no_exists')
    assert 'No reference genome with name=[no_exists]' in str(execinfo.value)


def test_count_on_reference_with_index_unknowns_low_slice_limit(database, snippy_nucleotide_data_package,
                                                                reference_service_with_data,
                                                                sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=5)

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 111 == variation_service.count_on_reference(reference_name='genome', include_unknown=False)
    assert 632 == variation_service.count_on_reference(reference_name='genome', include_unknown=True)

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.count_on_reference('no_exists')
    assert 'No reference genome with name=[no_exists]' in str(execinfo.value)


def test_count_on_reference_with_index_unknowns_high_slice_limit(database, snippy_nucleotide_data_package,
                                                                 reference_service_with_data,
                                                                 sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=1000)

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 111 == variation_service.count_on_reference(reference_name='genome', include_unknown=False)
    assert 632 == variation_service.count_on_reference(reference_name='genome', include_unknown=True)

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.count_on_reference('no_exists')
    assert 'No reference genome with name=[no_exists]' in str(execinfo.value)


def test_insert_variants_masked_regions(database, snippy_nucleotide_data_package, reference_service_with_data,
                                        sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)

    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    samples = session.query(Sample).all()

    genomic_masks = {s.name: v.masked_regions for s in samples for v in s.sample_nucleotide_variation}
    assert 437 == len(genomic_masks['SampleA'])
    assert 276 == len(genomic_masks['SampleB'])
    assert 329 == len(genomic_masks['SampleC'])

    assert 632 == session.query(NucleotideVariantsSamples).count()
    assert 521 == session.query(NucleotideVariantsSamples) \
        .filter(NucleotideVariantsSamples.var_type == 'UNKNOWN_MISSING') \
        .count()
    assert 111 == session.query(NucleotideVariantsSamples) \
        .filter(NucleotideVariantsSamples.var_type != 'UNKNOWN_MISSING') \
        .count()


def test_mutation_counts_on_reference(database, snippy_nucleotide_data_package, reference_service_with_data,
                                      sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    mutation_counts = variation_service.mutation_counts_on_reference('genome')
    assert 111 == len(mutation_counts)
    assert 2 == mutation_counts['reference:839:1:G']
    assert 1 == mutation_counts['reference:866:9:G']
    assert 1 == mutation_counts['reference:1048:1:G']
    assert 2 == mutation_counts['reference:3897:5:G']

    mutation_counts = variation_service.mutation_counts_on_reference('genome', include_unknown=True)
    assert 632 == len(mutation_counts)
    assert 2 == mutation_counts['reference:839:1:G']
    assert 1 == mutation_counts['reference:866:9:G']
    assert 1 == mutation_counts['reference:1048:1:G']
    assert 2 == mutation_counts['reference:3897:5:G']
    assert 2 == mutation_counts['reference:3897:5:G']

    assert 3 == mutation_counts['reference:9:1:?']
    assert 2 == mutation_counts['reference:5100:1:?']
    assert 1 == mutation_counts['reference:888:1:?']


def test_get_variants_on_reference_no_index_unknowns(database, snippy_nucleotide_data_package,
                                                     reference_service_with_data,
                                                     sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    mutations = variation_service.get_variants_on_reference('genome', include_unknown=False)
    assert 112 == len(mutations)
    assert 2 == len(mutations['reference:839:1:G'].sample_ids)
    assert 839 == mutations['reference:839:1:G'].position

    assert 1 == len(mutations['reference:866:9:G'].sample_ids)
    assert 866 == mutations['reference:866:9:G'].position

    assert 1 == len(mutations['reference:1048:1:G'].sample_ids)
    assert 1048 == mutations['reference:1048:1:G'].position

    assert 2 == len(mutations['reference:3897:5:G'].sample_ids)
    assert 3897 == mutations['reference:3897:5:G'].position

    # Test should still have the same number of mutations when including unknowns since none were indexed
    mutations = variation_service.get_variants_on_reference('genome', include_unknown=True)
    assert 112 == len(mutations)
    assert 'reference:839:1:G' in mutations

    # Test setting include_present to False, which should give no results in this case
    mutations = variation_service.get_variants_on_reference('genome', include_present=False, include_unknown=True)
    assert 0 == len(mutations)


def test_get_variants_on_reference_index_unknowns(database, snippy_nucleotide_data_package, reference_service_with_data,
                                                  sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    mutations = variation_service.get_variants_on_reference('genome', include_unknown=True)
    assert 632 == len(mutations)
    assert 2 == len(mutations['reference:839:1:G'].sample_ids)
    assert 839 == mutations['reference:839:1:G'].position

    assert 1 == len(mutations['reference:866:9:G'].sample_ids)
    assert 866 == mutations['reference:866:9:G'].position

    assert 1 == len(mutations['reference:1048:1:G'].sample_ids)
    assert 1048 == mutations['reference:1048:1:G'].position

    assert 2 == len(mutations['reference:3897:5:G'].sample_ids)
    assert 3897 == mutations['reference:3897:5:G'].position

    assert 3 == len(mutations['reference:87:1:?'].sample_ids)
    assert 87 == mutations['reference:87:1:?'].position

    mutations = variation_service.get_variants_on_reference('genome', include_unknown=False)
    assert 111 == len(mutations)
    assert 2 == len(mutations['reference:839:1:G'].sample_ids)
    assert 839 == mutations['reference:839:1:G'].position

    assert 1 == len(mutations['reference:866:9:G'].sample_ids)
    assert 866 == mutations['reference:866:9:G'].position

    assert 1 == len(mutations['reference:1048:1:G'].sample_ids)
    assert 1048 == mutations['reference:1048:1:G'].position

    assert 2 == len(mutations['reference:3897:5:G'].sample_ids)
    assert 3897 == mutations['reference:3897:5:G'].position

    assert 'reference:87:1:?' not in mutations

    # Test only including unknowns
    mutations = variation_service.get_variants_on_reference('genome', include_present=False, include_unknown=True)
    assert 521 == len(mutations)
    assert 'reference:87:1:?' in mutations
    assert 'reference:839:1:G' not in mutations

    # Test inluding neigther present nor unknowns
    mutations = variation_service.get_variants_on_reference('genome', include_present=False, include_unknown=False)
    assert 0 == len(mutations)


def test_get_features(database, snippy_nucleotide_data_package, reference_service_with_data,
                      sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    # Test include unknown
    mutations = variation_service.get_features(include_unknown=True)
    assert 632 == len(mutations)
    assert 2 == len(mutations['reference:839:1:G'].sample_ids)
    assert 1 == len(mutations['reference:866:9:G'].sample_ids)
    assert 1 == len(mutations['reference:1048:1:G'].sample_ids)
    assert 2 == len(mutations['reference:3897:5:G'].sample_ids)
    assert 3 == len(mutations['reference:87:1:?'].sample_ids)

    # Test include unknown ids translated
    mutations = variation_service.get_features(include_unknown=True, id_type='spdi_ref')
    assert 632 == len(mutations)
    assert 2 == len(mutations['reference:839:C:G'].sample_ids)
    assert 1 == len(mutations['reference:866:GCCAGATCC:G'].sample_ids)
    assert 1 == len(mutations['reference:1048:C:G'].sample_ids)
    assert 2 == len(mutations['reference:3897:GCGCA:G'].sample_ids)
    assert 3 == len(mutations['reference:87:G:?'].sample_ids)

    # Test no include unknown
    mutations = variation_service.get_features(include_unknown=False)
    assert 111 == len(mutations)
    assert 2 == len(mutations['reference:839:1:G'].sample_ids)
    assert 1 == len(mutations['reference:866:9:G'].sample_ids)
    assert 1 == len(mutations['reference:1048:1:G'].sample_ids)
    assert 2 == len(mutations['reference:3897:5:G'].sample_ids)
    assert 'reference:87:1:?' not in mutations

    # Test no include unknown ids translated
    mutations = variation_service.get_features(include_unknown=False,
                                               id_type='spdi_ref')
    assert 111 == len(mutations)
    assert 2 == len(mutations['reference:839:C:G'].sample_ids)
    assert 1 == len(mutations['reference:866:GCCAGATCC:G'].sample_ids)
    assert 1 == len(mutations['reference:1048:C:G'].sample_ids)
    assert 2 == len(mutations['reference:3897:GCGCA:G'].sample_ids)
    assert 'reference:87:G:?' not in mutations
    assert 'reference:87:1:?' not in mutations

    # Test only including unknowns
    mutations = variation_service.get_features(include_present=False, include_unknown=True)
    assert 521 == len(mutations)
    assert 3 == len(mutations['reference:87:1:?'].sample_ids)
    assert 'reference:839:1:G' not in mutations

    # Test only including unknowns, id translated
    mutations = variation_service.get_features(include_present=False, include_unknown=True,
                                               id_type='spdi_ref')
    assert 521 == len(mutations)
    assert 3 == len(mutations['reference:87:G:?'].sample_ids)
    assert 'reference:839:1:G' not in mutations

    # Test including neither present nor unknowns
    mutations = variation_service.get_features(include_present=False, include_unknown=False)
    assert 0 == len(mutations)


def test_insert_variants_examine_variation_with_unknown(database, snippy_nucleotide_data_package,
                                                        reference_service_with_data,
                                                        sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    sample_name_ids = sample_service.find_sample_name_ids({'SampleA', 'SampleB', 'SampleC'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 632 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'

    # Check to make sure some variants are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1048,
        'deletion': len('C'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
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
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleB']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1984,
        'deletion': len('GTGATTG'),
        'insertion': 'TTGA'
    })
    assert 'OTHER' == v.var_type, 'Type is incorrect'
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)

    # Check to make sure unknown/missing positions are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 9,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleA'], sample_name_ids['SampleB'],
            sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 5100,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleA'], sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 873,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)


def do_test_insert_variants_examine_variation_annotations(variation_service, sample_service,
                                                          database, snpeff_nucleotide_data_package):
    session = database.get_session()

    variation_service.insert(feature_scope_name='NC_011083', data_package=snpeff_nucleotide_data_package)

    sample_name_ids = sample_service.find_sample_name_ids({'SH10-014', 'SH14-001', 'SH14-014'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 177 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'

    # SNV/SNP
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'NC_011083',
        'position': 602110,
        'deletion': len('T'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SH10-014'], sample_name_ids['SH14-001'],
            sample_name_ids['SH14-014']} == set(v.sample_ids)

    assert 'missense_variant' == v.annotation
    assert 'MODERATE' == v.annotation_impact
    assert 'SEHA_RS03410' == v.annotation_gene_id
    assert 'SEHA_RS03410' == v.annotation_gene_name
    assert 'transcript' == v.annotation_feature_type
    assert 'protein_coding' == v.annotation_transcript_biotype
    assert 'c.1229T>G' == v.annotation_hgvs_c
    assert 'p.Leu410Arg' == v.annotation_hgvs_p
    assert 'hgvs:NC_011083:SEHA_RS03410:c.1229T>G' == v.id_hgvs_c
    assert 'hgvs:NC_011083:SEHA_RS03410:p.Leu410Arg' == v.id_hgvs_p

    # Deletion
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'NC_011083',
        'position': 833114,
        'deletion': len('GT'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SH14-001'],
            sample_name_ids['SH14-014']} == set(v.sample_ids)

    assert 'frameshift_variant' == v.annotation
    assert 'HIGH' == v.annotation_impact
    assert 'SEHA_RS04540' == v.annotation_gene_id
    assert 'glf' == v.annotation_gene_name
    assert 'transcript' == v.annotation_feature_type
    assert 'protein_coding' == v.annotation_transcript_biotype
    assert 'c.696delT' == v.annotation_hgvs_c
    assert 'p.Phe232fs' == v.annotation_hgvs_p
    assert 'hgvs:NC_011083:SEHA_RS04540:c.696delT' == v.id_hgvs_c
    assert 'hgvs:NC_011083:SEHA_RS04540:p.Phe232fs' == v.id_hgvs_p

    # Insertion
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'NC_011083',
        'position': 2380144,
        'deletion': len('A'),
        'insertion': 'AATTTTAT'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SH14-001'],
            sample_name_ids['SH14-014']} == set(v.sample_ids)

    assert 'frameshift_variant' == v.annotation
    assert 'HIGH' == v.annotation_impact
    assert 'SEHA_RS12315' == v.annotation_gene_id
    assert 'SEHA_RS12315' == v.annotation_gene_name
    assert 'transcript' == v.annotation_feature_type
    assert 'protein_coding' == v.annotation_transcript_biotype
    assert 'c.431_437dupTTTATAT' == v.annotation_hgvs_c
    assert 'p.Ter149fs' == v.annotation_hgvs_p
    assert 'hgvs:NC_011083:SEHA_RS12315:c.431_437dupTTTATAT' == v.id_hgvs_c
    assert 'hgvs:NC_011083:SEHA_RS12315:p.Ter149fs' == v.id_hgvs_p

    # MNP
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'NC_011083',
        'position': 3535656,
        'deletion': len('CTCAAAA'),
        'insertion': 'GTCATAG'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'MNP' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SH14-001'],
            sample_name_ids['SH10-014']} == set(v.sample_ids)

    assert 'missense_variant' == v.annotation
    assert 'MODERATE' == v.annotation_impact
    assert 'SEHA_RS17780' == v.annotation_gene_id
    assert 'oadA' == v.annotation_gene_name
    assert 'transcript' == v.annotation_feature_type
    assert 'protein_coding' == v.annotation_transcript_biotype
    assert 'c.582_588delTTTTGAGinsCTATGAC' == v.annotation_hgvs_c
    assert 'p.PheGlu195TyrAsp' == v.annotation_hgvs_p
    assert 'hgvs:NC_011083:SEHA_RS17780:c.582_588delTTTTGAGinsCTATGAC' == v.id_hgvs_c
    assert 'hgvs:NC_011083:SEHA_RS17780:p.PheGlu195TyrAsp' == v.id_hgvs_p


def test_insert_variants_examine_variation_annotations(database, snpeff_nucleotide_data_package,
                                                       reference_service_with_snpeff_data,
                                                       sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_snpeff_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    do_test_insert_variants_examine_variation_annotations(variation_service=variation_service,
                                                          sample_service=sample_service,
                                                          database=database,
                                                          snpeff_nucleotide_data_package=snpeff_nucleotide_data_package)


def test_insert_variants_examine_variation_annotations_parallel_variants(database,
                                                                         snpeff_nucleotide_data_package_parallel_variants,
                                                                         reference_service_with_snpeff_data,
                                                                         sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_snpeff_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    do_test_insert_variants_examine_variation_annotations(variation_service=variation_service,
                                                          sample_service=sample_service,
                                                          database=database,
                                                          snpeff_nucleotide_data_package=snpeff_nucleotide_data_package_parallel_variants)


def test_insert_variants_regular_vcf_reader_examine_variation(database, regular_nucleotide_data_package,
                                                              reference_service_with_data,
                                                              sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=regular_nucleotide_data_package)

    sample_name_ids = sample_service.find_sample_name_ids({'SampleA', 'SampleB', 'SampleC'})
    assert 3 == len(sample_name_ids.values())

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    assert 632 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'

    # Check to make sure some variants are stored
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1048,
        'deletion': len('C'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleA']} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1135,
        'deletion': len('CCT'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB'], sample_name_ids['SampleC']} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 883,
        'deletion': len('CACATG'),
        'insertion': 'C'
    })
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleB']} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1984,
        'deletion': len('GTGATTG'),
        'insertion': 'TTGA'
    })
    assert 'OTHER' == v.var_type, 'Type is incorrect'
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    # Check to make sure unknown/missing positions are stored
    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 9,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleA'], sample_name_ids['SampleB'],
            sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 5100,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleA'], sample_name_ids['SampleC']} == set(v.sample_ids)

    v = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 873,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None
    assert 'UNKNOWN_MISSING' == v.var_type
    assert v.id_hgvs_c is None
    assert v.id_hgvs_p is None
    assert {sample_name_ids['SampleC']} == set(v.sample_ids)


def test_insert_variants_regular_vcf_reader_examine_variation_no_unknown(database, regular_nucleotide_data_package,
                                                                         reference_service_with_data,
                                                                         sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)
    session = database.get_session()

    variation_service.insert(feature_scope_name='genome', data_package=regular_nucleotide_data_package)

    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'
    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'


def test_insert_variants_duplicates(database, snippy_nucleotide_data_package, reference_service_with_data,
                                    sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    with pytest.raises(EntityExistsError) as execinfo:
        variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 'Passed samples already have features for feature scope [genome]' in str(execinfo.value)
    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'


def test_insert_variants_duplicates_subset(database, snippy_nucleotide_data_package, reference_service_with_data,
                                           sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
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
        assert 112 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
        assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
        assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleSequences'


def test_multiple_insert_variants_different_samples(database, snippy_nucleotide_data_package_AB,
                                                    snippy_nucleotide_data_package_C,
                                                    reference_service_with_data,
                                                    sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=50)
    session = database.get_session()

    assert 0 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 0 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 0 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    # Insert samples A and B
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package_AB)

    assert 597 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 2 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 2 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    sampleA = session.query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = session.query(Sample).filter(Sample.name == 'SampleB').one()

    ## Variant which is present only in Sample A
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1708,
        'deletion': len('ATGCTGTTCAATAC'),
        'insertion': 'A'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sampleA.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Variant which is present only in Samples B
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4539,
        'deletion': len('T'),
        'insertion': 'A'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Variant which is present only in Samples B and C
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4585,
        'deletion': len('T'),
        'insertion': 'C'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample A
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 5071,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleA.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample B
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4247,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample B and C
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 350,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    # Insert sample C
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package_C)

    assert 632 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 3 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 3 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    sampleC = session.query(Sample).filter(Sample.name == 'SampleC').one()

    ## Variant which is present only in Sample A
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1708,
        'deletion': len('ATGCTGTTCAATAC'),
        'insertion': 'A'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sampleA.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Variant which is present only in Samples B
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4539,
        'deletion': len('T'),
        'insertion': 'A'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Variant which is present only in Samples B and C
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4585,
        'deletion': len('T'),
        'insertion': 'C'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sampleB.id, sampleC.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample A
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 5071,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleA.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample B
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 4247,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleB.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample B and C
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 350,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleB.id, sampleC.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    ## Missing only in Sample C
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1989,
        'deletion': 1,
        'insertion': '?'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'UNKNOWN_MISSING' == v.var_type, 'Type is incorrect'
    assert {sampleC.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])


def test_multiple_insert_variants_different_reference_genomes(database, snippy_nucleotide_data_package_AB,
                                                              snpeff_nucleotide_data_package,
                                                              reference_service_with_data,
                                                              sample_service, filesystem_storage):
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=True,
                                         sql_select_limit=500)
    session = database.get_session()

    # Add snpeff reference genome
    reference_service_with_data.add_reference_genome(reference_file_snpeff)

    # Insert samples A and B, and snpeff samples
    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package_AB)
    variation_service.insert(feature_scope_name='NC_011083', data_package=snpeff_nucleotide_data_package)

    assert 774 == session.query(NucleotideVariantsSamples).count(), 'Incorrect number of storage entries'
    assert 5 == session.query(Sample).count(), 'Incorrect number of Samples'
    assert 5 == session.query(SampleNucleotideVariation).count(), 'Incorrect number of SampleNucleotideVariation'

    sampleA = session.query(Sample).filter(Sample.name == 'SampleA').one()
    sample_SH10_014 = session.query(Sample).filter(Sample.name == 'SH10-014').one()
    sample_SH14_001 = session.query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_SH14_014 = session.query(Sample).filter(Sample.name == 'SH14-014').one()

    ## Variant which is present only in Sample A
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'reference',
        'position': 1708,
        'deletion': len('ATGCTGTTCAATAC'),
        'insertion': 'A'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'INDEL' == v.var_type, 'Type is incorrect'
    assert {sampleA.id} == set(v.sample_ids)
    assert all([x is None for x in [v.annotation, v.annotation_impact, v.annotation_gene_name, v.annotation_gene_id,
                                    v.annotation_feature_type, v.annotation_transcript_biotype, v.annotation_hgvs_c,
                                    v.annotation_hgvs_p]])

    # SNV/SNP for snpeff data
    v: NucleotideVariantsSamples = session.query(NucleotideVariantsSamples).get({
        'sequence': 'NC_011083',
        'position': 602110,
        'deletion': len('T'),
        'insertion': 'G'
    })
    assert v is not None, 'Particular storage does not exist'
    assert 'SNP' == v.var_type, 'Type is incorrect'
    assert {sample_SH10_014.id, sample_SH14_001.id, sample_SH14_014.id} == set(v.sample_ids)

    assert 'missense_variant' == v.annotation
    assert 'MODERATE' == v.annotation_impact
    assert 'SEHA_RS03410' == v.annotation_gene_id
    assert 'SEHA_RS03410' == v.annotation_gene_name
    assert 'transcript' == v.annotation_feature_type
    assert 'protein_coding' == v.annotation_transcript_biotype
    assert 'c.1229T>G' == v.annotation_hgvs_c
    assert 'p.Leu410Arg' == v.annotation_hgvs_p
    assert 'hgvs:NC_011083:SEHA_RS03410:c.1229T>G' == v.id_hgvs_c
    assert 'hgvs:NC_011083:SEHA_RS03410:p.Leu410Arg' == v.id_hgvs_p


def test_get_variants_ordered(database, snippy_nucleotide_data_package, reference_service_with_data,
                              sample_service, filesystem_storage):
    def find_variant_by_position(variants: List[NucleotideVariantsSamples], position: int) -> NucleotideVariantsSamples:
        for v in variants:
            if v.position == position:
                return v

        raise Exception(f'Could not find storage with position {position}')

    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service_with_data,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=500)

    variation_service.insert(feature_scope_name='genome', data_package=snippy_nucleotide_data_package)

    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    variants = variation_service.get_variants_ordered(sequence_name='reference')

    assert 60 == len(variants), 'Incorrect number of variants returned (counting only SNV/SNPs)'

    v1 = find_variant_by_position(variants, 4265)
    assert 'reference:4265:1:C' == v1.spdi, 'Incorrect storage returned'
    assert {sampleA.id} == set(v1.sample_ids)

    v2 = find_variant_by_position(variants, 839)
    assert 'reference:839:1:G' == v2.spdi, 'Incorrect storage returned'
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
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Position']) == list(mutations_df['Position'])
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Insertion']) == list(mutations_df['Insertion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])
