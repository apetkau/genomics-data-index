import tempfile
from pathlib import Path
from typing import Optional

from pybedtools import BedTool

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataBedMask import NucleotideSampleDataBedMask
from genomics_data_index.test.integration import variation_dir

sampleA_mask_file = variation_dir / 'SampleA.bed.gz'
sampleA_vcf = variation_dir / 'SampleA.vcf.gz'


def create_sample_data(name: str, vcf: Path, mask: Optional[Path]) -> NucleotideSampleData:
    return NucleotideSampleDataBedMask.create(sample_name=name, vcf_file=vcf,
                                              mask_bed_file=mask)


def test_preprocess_mask():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        sample_data = NucleotideSampleDataBedMask.create(sample_name='SampleA', vcf_file=sampleA_vcf,
                                                         mask_bed_file=sampleA_mask_file)
        assert not sample_data.is_preprocessed()

        data_preprocessed = sample_data.preprocess(tmp_dir)

        assert data_preprocessed.is_preprocessed()
        vcf_file, vcf_file_index = data_preprocessed.get_vcf_file()
        assert tmp_dir == data_preprocessed.get_mask_file().parent
        assert tmp_dir == vcf_file.parent
        assert tmp_dir == vcf_file_index.parent
        assert 70 == len(data_preprocessed.get_mask())


def test_preprocess_mask_empty():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        sample_data = NucleotideSampleDataBedMask.create(sample_name='SampleA', vcf_file=sampleA_vcf,
                                                         mask_bed_file=None)
        assert not sample_data.is_preprocessed()

        data_preprocessed = sample_data.preprocess(tmp_dir)

        assert data_preprocessed.is_preprocessed()
        vcf_file, vcf_file_index = data_preprocessed.get_vcf_file()
        assert tmp_dir == data_preprocessed.get_mask_file().parent
        assert tmp_dir == vcf_file.parent
        assert tmp_dir == vcf_file_index.parent

        assert data_preprocessed.get_mask().is_empty()


def test_preprocess_other_mask():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        initial_mask = tmp_dir / 'initial.bed.gz'
        BedTool([('NC_011083', 0, 100)]).saveas(initial_mask)
        sample_data = NucleotideSampleDataBedMask.create(sample_name='SampleA', vcf_file=sampleA_vcf,
                                                         mask_bed_file=initial_mask)
        assert not sample_data.is_preprocessed()

        data_preprocessed = sample_data.preprocess(tmp_dir)

        assert data_preprocessed.is_preprocessed()
        vcf_file, vcf_file_index = data_preprocessed.get_vcf_file()
        assert tmp_dir == data_preprocessed.get_mask_file().parent
        assert tmp_dir == vcf_file.parent
        assert tmp_dir == vcf_file_index.parent

        assert 100 == len(data_preprocessed.get_mask())
