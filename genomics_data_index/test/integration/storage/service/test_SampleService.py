from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.storage.service.SampleService import SampleService


def count_samples(sample_service, variation_service):
    assert 3 == sample_service.count_samples()


def test_samples_with_variants(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_with_variants_on_sequence(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants_on_sequence('reference')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_associated_with_reference(sample_service: SampleService, variation_service):
    samples_reference = sample_service.get_samples_associated_with_reference('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_reference}


def test_samples_set_associated_with_reference(database, sample_service, variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    expected_sample_set = SampleSet([sampleA.id, sampleB.id, sampleC.id])

    samples_set_reference = sample_service.get_samples_set_associated_with_reference('genome')

    assert isinstance(samples_set_reference, SampleSet)
    assert set(samples_set_reference) == set(expected_sample_set)


def test_find_sample_name_ids(database, sample_service, variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    sample_ids_list = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])

    assert {sampleA.name: sampleA.id, sampleB.name: sampleB.id, sampleC.name: sampleC.id} == sample_ids_list


def test_count_samples_associated_with_reference(sample_service, variation_service):
    assert 3 == sample_service.count_samples_associated_with_reference('genome')


def test_count_samples_associated_with_reference_empty(sample_service, variation_service):
    assert 0 == sample_service.count_samples_associated_with_reference('no_exist')


def test_samples_associated_with_reference_empty(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_associated_with_reference('no_exist')
    assert set() == {sample.name for sample in samples_with_variants}


def test_samples_with_variants_empty(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('does_not_exist')
    assert set() == set(samples_with_variants)


def test_samples_which_exists(sample_service, variation_service):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(sample_service.which_exists(['SampleA', 'SampleB',
                                                                                 'SampleC', 'not_exist']))


def test_samples_which_exists_only_one(sample_service, variation_service):
    assert {'SampleA'} == set(sample_service.which_exists(['SampleA']))


def test_samples_which_exists_none(sample_service, variation_service):
    assert [] == sample_service.which_exists(['not_exist'])


def test_samples_which_exists_none2(sample_service, variation_service):
    assert [] == sample_service.which_exists([])


def test_get_samples(sample_service, variation_service):
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in sample_service.get_samples()}


def test_get_samples_empty(sample_service):
    assert [] == sample_service.get_samples()


def test_get_existing_samples_by_names_all(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names(['SampleA', 'SampleB', 'SampleC'])
    assert 3 == len(samples)
    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}


def test_get_existing_samples_by_names_two(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names(['SampleA', 'SampleC'])
    assert 2 == len(samples)
    assert {'SampleA', 'SampleC'} == {s.name for s in samples}


def test_get_existing_samples_by_names_one(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names(['SampleC'])
    assert 1 == len(samples)
    assert {'SampleC'} == {s.name for s in samples}


def test_get_existing_samples_by_names_none(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names([])
    assert 0 == len(samples)
    assert set() == {s.name for s in samples}


def test_get_existing_samples_by_names_sample_not_exist(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names(['SampleA', 'NotFound', 'SampleC'])
    assert 2 == len(samples)
    assert {'SampleA', 'SampleC'} == {s.name for s in samples}


def test_get_existing_samples_by_names_sample_not_exis2t(sample_service, variation_service):
    samples = sample_service.get_existing_samples_by_names(['SampleA', 'NotFound', 'SampleC', 'NotFound2'])
    assert 2 == len(samples)
    assert {'SampleA', 'SampleC'} == {s.name for s in samples}


def test_find_samples_by_ids(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA'} == {s.name for s in samples}


def test_find_samples_by_ids_2_samples(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA'], sample_name_ids['SampleC']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA', 'SampleC'} == {s.name for s in samples}


def test_find_samples_by_ids_3_samples(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA'], sample_name_ids['SampleB'],
                                       sample_name_ids['SampleC']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}


def test_find_samples_by_features_variations(database, sample_service, variation_service):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    variant_samples = sample_service.find_samples_by_features([QueryFeatureMutationSPDI('reference:5061:G:A')])

    assert f'reference:5061:G:A' in variant_samples
    assert {'SampleB'} == {s.name for s in variant_samples[f'reference:5061:G:A']}
    assert {sampleB.id} == {s.id for s in variant_samples[f'reference:5061:G:A']}


def test_find_sample_sets_by_features_variations(database, sample_service, variation_service):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    sample_sets = sample_service.find_sample_sets_by_features([QueryFeatureMutationSPDI('reference:5061:G:A')])

    assert f'reference:5061:G:A' in sample_sets
    assert {sampleB.id} == set(sample_sets[f'reference:5061:G:A'])


def test_find_sample_sets_by_features_variations_two_samples(database, sample_service, variation_service):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    variant_samples = sample_service.find_sample_sets_by_features([QueryFeatureMutationSPDI('reference:4975:T:C')])

    assert {'reference:4975:T:C'} == set(variant_samples.keys())
    assert {sampleC.id, sampleB.id} == set(variant_samples[f'reference:4975:T:C'])
    assert 2 == len(variant_samples[f'reference:4975:T:C'])


def test_find_samples_by_features_variations_numeric_deletion(database, sample_service, variation_service):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    variant_samples = sample_service.find_samples_by_features([QueryFeatureMutationSPDI('reference:5061:1:A')])

    assert f'reference:5061:1:A' in variant_samples
    assert {'SampleB'} == {s.name for s in variant_samples[f'reference:5061:1:A']}
    assert {sampleB.id} == {s.id for s in variant_samples[f'reference:5061:1:A']}


def test_count_samples_by_variation_features_single_feature(sample_service, variation_service):
    features = [QueryFeatureMutationSPDI('reference:5061:G:A')]

    variant_counts = sample_service.count_samples_by_features(features)

    assert 1 == len(variant_counts)
    assert 1 == variant_counts[f'reference:5061:G:A']


def test_count_samples_by_variation_features_multiple_features(sample_service, variation_service):
    features = [QueryFeatureMutationSPDI('reference:5061:G:A'), QueryFeatureMutationSPDI('reference:3063:A:ATGCAGC')]

    variant_counts = sample_service.count_samples_by_features(features)

    assert 2 == len(variant_counts)
    assert 1 == variant_counts[f'reference:5061:G:A']
    assert 2 == variant_counts[f'reference:3063:A:ATGCAGC']


def test_count_samples_by_variation_features_multiple_features_numeric_deletion(sample_service, variation_service):
    features = [QueryFeatureMutationSPDI('reference:5061:1:A'), QueryFeatureMutationSPDI('reference:3063:1:ATGCAGC')]

    variant_counts = sample_service.count_samples_by_features(features)

    assert 2 == len(variant_counts)
    assert 1 == variant_counts[f'reference:5061:1:A']
    assert 2 == variant_counts[f'reference:3063:1:ATGCAGC']


def test_get_samples_with_mlst_alleles(sample_service, mlst_service_loaded):
    samples = sample_service.get_samples_with_mlst_alleles('lmonocytogenes')

    assert {'CFSAN002349', 'CFSAN023463'} == {sample.name for sample in samples}


def test_get_samples_with_mlst_alleles2(sample_service, mlst_service_loaded):
    samples = sample_service.get_samples_with_mlst_alleles('ecoli')

    assert {'2014C-3598', '2014C-3599'} == {sample.name for sample in samples}


def test_find_samples_by_features_mlst(database, sample_service, mlst_service_loaded):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    mlst_samples = sample_service.find_samples_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1')])

    assert {'lmonocytogenes:abcZ:1'} == set(mlst_samples.keys())
    assert {'CFSAN002349', 'CFSAN023463'} == {s.name for s in mlst_samples['lmonocytogenes:abcZ:1']}
    assert {sample1.id, sample2.id} == {s.id for s in mlst_samples['lmonocytogenes:abcZ:1']}


def test_find_sample_sets_by_features_mlst(database, sample_service: SampleService, mlst_service_loaded):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    mlst_sample_sets = sample_service.find_sample_sets_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1')])

    assert {'lmonocytogenes:abcZ:1'} == set(mlst_sample_sets.keys())
    assert {sample1.id, sample2.id} == set(mlst_sample_sets['lmonocytogenes:abcZ:1'])


def test_find_samples_by_features_mlst_two(database, sample_service, mlst_service_loaded):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()
    sample3 = database.get_session().query(Sample).filter(Sample.name == '2014C-3599').one()
    sample4 = database.get_session().query(Sample).filter(Sample.name == '2014C-3598').one()

    mlst_samples = sample_service.find_samples_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1'),
                                                            QueryFeatureMLST('ecoli:adk:100')])

    assert {'lmonocytogenes:abcZ:1', 'ecoli:adk:100'} == set(mlst_samples.keys())

    assert {'CFSAN002349', 'CFSAN023463'} == {s.name for s in mlst_samples['lmonocytogenes:abcZ:1']}
    assert {sample1.id, sample2.id} == {s.id for s in mlst_samples['lmonocytogenes:abcZ:1']}

    assert {'2014C-3599', '2014C-3598'} == {s.name for s in mlst_samples['ecoli:adk:100']}
    assert {sample3.id, sample4.id} == {s.id for s in mlst_samples['ecoli:adk:100']}


def test_find_sample_sets_by_features_mlst_two(database, sample_service, mlst_service_loaded):
    sample1 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = database.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()
    sample3 = database.get_session().query(Sample).filter(Sample.name == '2014C-3599').one()
    sample4 = database.get_session().query(Sample).filter(Sample.name == '2014C-3598').one()

    mlst_sample_sets = sample_service.find_sample_sets_by_features([QueryFeatureMLST('lmonocytogenes:abcZ:1'),
                                                                    QueryFeatureMLST('ecoli:adk:100')])

    assert {'lmonocytogenes:abcZ:1', 'ecoli:adk:100'} == set(mlst_sample_sets.keys())

    assert {sample1.id, sample2.id} == set(mlst_sample_sets['lmonocytogenes:abcZ:1'])
    assert {sample3.id, sample4.id} == set(mlst_sample_sets['ecoli:adk:100'])


def test_count_samples_by_mlst_features_single_feature(sample_service, mlst_service_loaded):
    features = [QueryFeatureMLST('lmonocytogenes:abcZ:1')]

    mlst_counts = sample_service.count_samples_by_features(features)

    assert 1 == len(mlst_counts)
    assert 2 == mlst_counts['lmonocytogenes:abcZ:1']


def test_count_samples_by_mlst_features_multiple_features(sample_service, mlst_service_loaded):
    features = [QueryFeatureMLST('lmonocytogenes:abcZ:1'), QueryFeatureMLST('ecoli:adk:100')]

    mlst_counts = sample_service.count_samples_by_features(features)

    assert 2 == len(mlst_counts)
    assert 2 == mlst_counts['lmonocytogenes:abcZ:1']
    assert 2 == mlst_counts['ecoli:adk:100']


def test_create_dataframe_from_sample_set(database, sample_service: SampleService, variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    sample_set = SampleSet([sampleA.id, sampleB.id, sampleC.id])
    queries_expression = ''

    df = sample_service.create_dataframe_from_sample_set(sample_set=sample_set,
                                                         universe_set=sample_set,
                                                         exclude_absent=True,
                                                         queries_expression=queries_expression)
    assert len(df) == 3
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleA.id, sampleB.id, sampleC.id] == df['Sample ID'].tolist()


def test_create_dataframe_from_sample_set_subset_samples(database, sample_service: SampleService, variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    sample_set = SampleSet([sampleA.id, sampleC.id])
    queries_expression = ''

    df = sample_service.create_dataframe_from_sample_set(sample_set,
                                                         universe_set=sample_set,
                                                         exclude_absent=True,
                                                         queries_expression=queries_expression)
    assert len(df) == 2
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleA.id, sampleC.id] == df['Sample ID'].tolist()


def test_create_dataframe_from_sample_set_subset_samples_include_all(database,
                                                                     sample_service: SampleService,
                                                                     variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    universe_set = SampleSet([sampleA.id, sampleB.id, sampleC.id])
    sample_set = SampleSet([sampleA.id, sampleC.id])
    queries_expression = ''

    df = sample_service.create_dataframe_from_sample_set(sample_set,
                                                         universe_set=universe_set,
                                                         exclude_absent=False,
                                                         queries_expression=queries_expression)
    assert len(df) == 3
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleA.id, sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['Present', 'Absent', 'Present'] == df['Status'].tolist()


def test_create_dataframe_from_sample_set_subset_samples_include_all_matching_universe(database,
                                                                                       sample_service: SampleService,
                                                                                       variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    sample_set = SampleSet([sampleA.id, sampleC.id])
    queries_expression = ''

    df = sample_service.create_dataframe_from_sample_set(sample_set,
                                                         universe_set=sample_set,
                                                         exclude_absent=False,
                                                         queries_expression=queries_expression)
    assert len(df) == 2
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleA.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['Present', 'Present'] == df['Status'].tolist()


def test_create_dataframe_from_sample_set_empty(sample_service: SampleService, variation_service):
    sample_set = SampleSet.create_empty()
    queries_expression = ''

    df = sample_service.create_dataframe_from_sample_set(sample_set=sample_set,
                                                         universe_set=sample_set,
                                                         exclude_absent=True,
                                                         queries_expression=queries_expression)
    assert len(df) == 0
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()


def test_create_dataframe_from_sample_set_with_query_expression(database, sample_service: SampleService,
                                                                variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    sample_set = SampleSet([sampleA.id, sampleB.id])
    queries_expression = 'lmonocytogenes:abc:1'

    df = sample_service.create_dataframe_from_sample_set(sample_set,
                                                         universe_set=sample_set,
                                                         exclude_absent=True,
                                                         queries_expression=queries_expression)
    assert {'lmonocytogenes:abc:1'} == set(df['Query'].tolist())


def test_get_all_sample_ids(database, sample_service, variation_service):
    sampleA = database.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = database.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    assert {sampleA.id, sampleB.id, sampleC.id} == set(sample_service.get_all_sample_ids())


def test_get_variants_samples_by_variation_features_only_spdi(database, sample_service, variation_service):
    sampleB = database.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    feature_id_nucleotide_samples = sample_service.get_variants_samples_by_variation_features([QueryFeatureMutationSPDI('reference:5061:G:A')])

    assert 1 == len(feature_id_nucleotide_samples)
    assert f'reference:5061:G:A' in feature_id_nucleotide_samples
    assert {sampleB.id} == set(feature_id_nucleotide_samples[f'reference:5061:G:A'].sample_ids)


def test_get_variants_samples_by_variation_features_only_hgvs_c(database, sample_service_snpeff_annotations):
    sample_sh10_014 = database.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    features = [QueryFeatureHGVS.create_from_id('hgvs:NC_011083:n.882634G>A')]

    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)

    assert 1 == len(feature_id_nucleotide_samples)
    assert 'hgvs:NC_011083:n.882634G>A' in feature_id_nucleotide_samples
    assert {sample_sh10_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:n.882634G>A'].sample_ids)


def test_get_variants_samples_by_variation_features_only_hgvs_cp_spdi(database, sample_service_snpeff_annotations):
    sample_sh14_001 = database.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = database.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()

    # hgvs c (nucleotide)
    features = [QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:c.670dupA')]
    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)
    assert 1 == len(feature_id_nucleotide_samples)
    assert 'hgvs:NC_011083:SEHA_RS04550:c.670dupA' in feature_id_nucleotide_samples
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:SEHA_RS04550:c.670dupA'].sample_ids)


    # hgvs p (protein)
    features = [QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')]
    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)
    assert 1 == len(feature_id_nucleotide_samples)
    assert 'hgvs:NC_011083:SEHA_RS04550:p.Ile224fs' in feature_id_nucleotide_samples
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:SEHA_RS04550:p.Ile224fs'].sample_ids)


    # spdi
    features = [QueryFeatureMutationSPDI('NC_011083:835147:C:CA')]
    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)
    assert 1 == len(feature_id_nucleotide_samples)
    assert 'NC_011083:835147:C:CA' in feature_id_nucleotide_samples
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['NC_011083:835147:C:CA'].sample_ids)


def test_get_variants_samples_by_variation_features_both_hgvs_spdi(database, sample_service_snpeff_annotations):
    sample_sh14_001 = database.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = database.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = database.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()

    features = [QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:c.670dupA'),
                QueryFeatureHGVS.create_from_id('hgvs:NC_011083:n.882634G>A'),
                QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs'),
                QueryFeatureMutationSPDI('NC_011083:835147:C:CA')]
    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)

    assert 4 == len(feature_id_nucleotide_samples)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:SEHA_RS04550:c.670dupA'].sample_ids)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:SEHA_RS04550:p.Ile224fs'].sample_ids)
    assert {sample_sh10_014.id} == set(feature_id_nucleotide_samples['hgvs:NC_011083:n.882634G>A'].sample_ids)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(feature_id_nucleotide_samples['NC_011083:835147:C:CA'].sample_ids)


def test_get_variants_samples_by_variation_features_no_matches(database, sample_service_snpeff_annotations):
    features = [QueryFeatureHGVS.create_from_id('hgvs:NC_011083:n.unknown')]

    feature_id_nucleotide_samples = sample_service_snpeff_annotations.get_variants_samples_by_variation_features(features)

    assert 0 == len(feature_id_nucleotide_samples)
