import tempfile
from pathlib import Path

import vcf

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.io.mutation.SnpEffDatabase import SnpEffDatabase
from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.test.integration import reference_file_5000_snpeff, snpeff_vcf_file


def test_config_io():
    with tempfile.TemporaryDirectory() as out_dir:
        tmp_dir = Path(out_dir)
        snpeff_config = tmp_dir / 'snpEff.config'

        snpeff_database = SnpEffDatabase(
            snpeff_config=snpeff_config,
            database_dir=tmp_dir,
            genome_name='genome',
        )

        config_file = tmp_dir / 'config.ini'
        snpeff_database.write_database_config(config_file)

        assert config_file.exists()

        snpeff_database2 = SnpEffDatabase.create_from_config(config_file)

        assert snpeff_database2.database_dir == snpeff_database.database_dir
        assert snpeff_database2.config == snpeff_database.config
        assert snpeff_database2.genome_name == snpeff_database.genome_name


def test_annotate_vcf_file():
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        output_vcf_file = database_dir / 'output.vcf.gz'

        snpeff_database = SequenceFile(reference_file_5000_snpeff).create_snpeff_database(database_dir)

        returned_output = snpeff_database.annotate(input_vcf_file=snpeff_vcf_file, output_vcf_file=output_vcf_file)

        assert output_vcf_file == returned_output

        # snpeff annotations should be added in headers
        reader = vcf.Reader(filename=str(output_vcf_file))
        assert 'ANN' in reader.infos

        # Verify VCF annotation contents
        variants_reader = VcfVariantsReader(sample_files_map={})
        vcf_df = variants_reader.read_vcf(output_vcf_file, 'SampleA')
        assert 2 == len(vcf_df)

        sub_vcf_df = vcf_df[['POS', 'ANN.Annotation', 'ANN.Gene_ID', 'ANN.Gene_Name', 'ANN.HGVS.c', 'ANN.HGVS.p']]
        assert [195, 'missense_variant', 'SEHA_RS00560', 'thrL', 'c.6C>G', 'p.N2K'] == \
               sub_vcf_df[sub_vcf_df['POS'] == 195].values[0].tolist()
        assert [207, 'synonymous_variant', 'SEHA_RS00560', 'thrL', 'c.18C>G', 'p.T6T'] == \
               sub_vcf_df[sub_vcf_df['POS'] == 207].values[
                   0].tolist()
