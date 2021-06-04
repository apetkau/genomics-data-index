import shutil
import tempfile
from pathlib import Path
from tempfile import TemporaryDirectory

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.configuration.Project import Project
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.test.integration import sample_dirs, reference_file


def test_summaries_loaded_data(loaded_database_genomic_data_store: GenomicsDataIndex):
    gds = loaded_database_genomic_data_store

    # Samples
    assert 9 == gds.count_samples()
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'} == set(gds.sample_names())

    # References
    assert 1 == gds.count_references()
    assert ['genome'] == gds.reference_names()

    # Mutations
    assert 112 == gds.count_mutations('genome')
    ms = gds.mutations_summary('genome', id_type='spdi')
    print(ms)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count'] == list(ms.columns)
    assert ['reference', 839, 1, 'G', 2] == ms.loc['reference:839:1:G'].values.tolist()
    assert ['reference', 866, 9, 'G', 1] == ms.loc['reference:866:9:G'].values.tolist()
    assert ['reference', 1048, 1, 'G', 1] == ms.loc['reference:1048:1:G'].values.tolist()
    assert ['reference', 3897, 5, 'G', 2] == ms.loc['reference:3897:5:G'].values.tolist()

    ms = gds.mutations_summary('genome', id_type='spdi_ref')
    print(ms)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count'] == list(ms.columns)
    assert ['reference', 839, 'C', 'G', 2] == ms.loc['reference:839:C:G'].values.tolist()
    assert ['reference', 866, 'GCCAGATCC', 'G', 1] == ms.loc['reference:866:GCCAGATCC:G'].values.tolist()
    assert ['reference', 1048, 'C', 'G', 1] == ms.loc['reference:1048:C:G'].values.tolist()
    assert ['reference', 3897, 'GCGCA', 'G', 2] == ms.loc['reference:3897:GCGCA:G'].values.tolist()


def test_connect_to_project_from_dir():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        Project.initialize_project(project_dir)

        ds = GenomicsDataIndex.connect(project_dir=project_dir)
        assert ds is not None
        assert ds.connection.reference_service is not None
        assert ds.connection.filesystem_storage.variation_dir.parent == project_dir / '.gdi-data'


def test_connect_to_project_from_project():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        project = Project.initialize_project(project_dir)

        ds = GenomicsDataIndex.connect(project=project)
        assert ds is not None
        assert ds.connection.reference_service is not None
        assert ds.connection.filesystem_storage.variation_dir.parent == project_dir / '.gdi-data'


def test_connect_and_tree_after_moving_project(loaded_data_store_from_project_dir):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        Project.initialize_project(project_dir)
        ds = GenomicsDataIndex.connect(project_dir=project_dir)
        database_connection = ds.connection

        # Load Nucleotide variation
        database_connection.reference_service.add_reference_genome(reference_file)
        snippy_tmp_dir = Path(tempfile.mkdtemp())
        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                      SerialSampleFilesProcessor(snippy_tmp_dir))
        database_connection.variation_service.insert(data_package, feature_scope_name='genome')

        # I should be able to build tree initially
        tq1 = ds.samples_query().build_tree(kind='mutation', scope='genome')
        assert tq1.tree is not None

        # Move project
        project_dir_2 = tmp_file / 'project2'
        shutil.move(project_dir, project_dir_2)

        ds2 = GenomicsDataIndex.connect(project_dir=project_dir_2)

        # I should still be able to build tree even after project has been moved
        tq2 = ds2.samples_query().build_tree(kind='mutation', scope='genome')
        assert tq2.tree is not None
