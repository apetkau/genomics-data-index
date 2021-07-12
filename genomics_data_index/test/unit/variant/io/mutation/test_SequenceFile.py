from pathlib import Path

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile


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
