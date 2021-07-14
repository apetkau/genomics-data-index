from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile


# Method used to shorten creating sequence files with a given name
def sf_name(name: str) -> SequenceFile:
    return SequenceFile(Path(name))


def test_get_genome_name():
    assert 'file' == sf_name('file.gbk').get_genome_name()
    assert 'file_1' == sf_name('file_1.gbk').get_genome_name()
    assert 'file' == sf_name('file.gbk.gz').get_genome_name()
    assert 'file' == sf_name('file.fasta.gz').get_genome_name()

    # Should only exclude if it's reads
    assert 'file_1' == sf_name('file_1.gbk').get_genome_name(exclude_paired_end_indicators=True)

    assert 'file' == sf_name('file.fastq.gz').get_genome_name()
    assert 'file_A' == sf_name('file_A.fastq').get_genome_name()
    assert 'file_A_1' == sf_name('file_A_1.fastq').get_genome_name()
    assert 'file_A' == sf_name('file_A_1.fastq').get_genome_name(exclude_paired_end_indicators=True)
    assert 'file_A_R1' == sf_name('file_A_R1.fastq').get_genome_name()
    assert 'file_A' == sf_name('file_A_R1.fastq').get_genome_name(exclude_paired_end_indicators=True)
    assert 'file_A_R2' == sf_name('file_A_R2.fastq').get_genome_name()
    assert 'file_A' == sf_name('file_A_R2.fastq').get_genome_name(exclude_paired_end_indicators=True)
    assert 'file_A_R1_xyz' == sf_name('file_A_R1_xyz.fastq').get_genome_name()
    assert 'file_A' == sf_name('file_A_R1_xyz.fastq').get_genome_name(exclude_paired_end_indicators=True)
    assert 'file_A_R2_xyz' == sf_name('file_A_R2_xyz.fq.gz').get_genome_name()
    assert 'file_A' == sf_name('file_A_R2_xyz.fq.gz').get_genome_name(exclude_paired_end_indicators=True)


def test_is_genbank():
    assert SequenceFile(Path('file.gbk')).is_genbank()
    assert SequenceFile(Path('file.gb')).is_genbank()
    assert SequenceFile(Path('file.gbk.gz')).is_genbank()
    assert SequenceFile(Path('file.gb.gz')).is_genbank()

    assert not SequenceFile(Path('file.fasta')).is_genbank()
    assert not SequenceFile(Path('file.fasta.gz')).is_genbank()


def test_is_fasta():
    assert SequenceFile(Path('file.fasta')).is_fasta()
    assert SequenceFile(Path('file.fna')).is_fasta()
    assert SequenceFile(Path('file.fa')).is_fasta()
    assert SequenceFile(Path('file.fasta.gz')).is_fasta()
    assert SequenceFile(Path('file.fna.gz')).is_fasta()
    assert SequenceFile(Path('file.fa.gz')).is_fasta()

    assert not SequenceFile(Path('file.gb')).is_fasta()
    assert not SequenceFile(Path('file.fastq.gz')).is_fasta()


def test_is_assembly():
    assert SequenceFile(Path('file.fasta')).is_assembly()
    assert SequenceFile(Path('file.fna')).is_assembly()
    assert SequenceFile(Path('file.gb')).is_assembly()
    assert SequenceFile(Path('file.gbk.gz')).is_assembly()
    assert SequenceFile(Path('file.fna.gz')).is_assembly()
    assert SequenceFile(Path('file.fasta.gz')).is_assembly()

    assert not SequenceFile(Path('file.fastq')).is_assembly()
    assert not SequenceFile(Path('file.fastq.gz')).is_assembly()


def test_is_reads():
    assert SequenceFile(Path('file.fastq')).is_reads()
    assert SequenceFile(Path('file.fq')).is_reads()
    assert SequenceFile(Path('file.fastq.gz')).is_reads()
    assert SequenceFile(Path('file.fq.gz')).is_reads()

    assert not SequenceFile(Path('file.fasta')).is_reads()
    assert not SequenceFile(Path('file.gbk.gz')).is_reads()


def test_name_differences():
    assert [] == sf_name('file.fastq').name_differences(sf_name('file.fastq'))
    assert [] == sf_name('').name_differences(sf_name(''))
    assert [] == sf_name('.').name_differences(sf_name('.'))

    assert ['1'] == sf_name('file1.fastq').name_differences(sf_name('file2.fastq'))
    assert ['2'] == sf_name('file2.fastq').name_differences(sf_name('file1.fastq'))

    assert ['1'] == sf_name('file_1.fastq').name_differences(sf_name('file_2.fastq'))
    assert ['2'] == sf_name('file_2.fastq').name_differences(sf_name('file_1.fastq'))

    assert ['1'] == sf_name('file_xyz_R1.fastq').name_differences(sf_name('file_xyz_R2.fastq'))
    assert ['2'] == sf_name('file_xyz_R2.fastq').name_differences(sf_name('file_xyz_R1.fastq'))

    assert ['a', '1'] == sf_name('file_a_xyz_R1.fastq').name_differences(sf_name('file_b_xyz_R2.fastq'))
    assert ['b', '2'] == sf_name('file_b_xyz_R2.fastq').name_differences(sf_name('file_a_xyz_R1.fastq'))

    # Difference at end
    assert ['gz'] == sf_name('file.fastq.gz').name_differences(sf_name('file.fastq.by'))

    # Difference at beginning
    assert ['p.'] == sf_name('p.file.fastq.gz').name_differences(sf_name('x_file.fastq.gz'))

    # Many differences in beginning, middle, and end
    assert ['p.', '1', 'gz'] == sf_name('p.file1.fastq.gz').name_differences(sf_name('x_file2.fastq.by'))
    assert ['x_', '2', 'by'] == sf_name('x_file2.fastq.by').name_differences(sf_name('p.file1.fastq.gz'))
