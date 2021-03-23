import tempfile
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pybedtools import BedTool

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions


def test_create_from_sequence():
    sequences = [SeqRecord(seq=Seq('ATCG-NN'), id='record1')]
    mask = MaskedGenomicRegions.from_sequences(sequences=sequences)

    assert 3 == len(mask), 'Invalid length'
    assert not mask.is_empty()
    assert not mask.contains('record1', 0)
    assert not mask.contains('record1', 3)
    assert mask.contains('record1', 4)
    assert mask.contains('record1', 6)
    assert not mask.contains('record1', 7)


def test_create_from_two_sequences():
    sequences = [
        SeqRecord(seq=Seq('ATCG-NN'), id='record1'),
        SeqRecord(seq=Seq('NN-GAT'), id='record2')
    ]
    mask = MaskedGenomicRegions.from_sequences(sequences=sequences)

    assert 6 == len(mask), 'Invalid length'
    assert not mask.is_empty()
    assert not mask.contains('record1', 0)
    assert not mask.contains('record1', 3)
    assert mask.contains('record1', 4)
    assert mask.contains('record1', 6)
    assert not mask.contains('record1', 7)

    assert mask.contains('record2', 0)
    assert mask.contains('record2', 2)
    assert not mask.contains('record2', 3)
    assert not mask.contains('record2', 5)


def test_sequence_names():
    sequences = [
        SeqRecord(seq=Seq('ATCG-NN'), id='record1'),
        SeqRecord(seq=Seq('NN-GAT'), id='record2')
    ]
    mask = MaskedGenomicRegions.from_sequences(sequences=sequences)
    assert {'record1', 'record2'} == mask.sequence_names()


def test_create_from_sequence_all_masked():
    sequences = [SeqRecord(seq=Seq('---------'), id='record1')]
    mask = MaskedGenomicRegions.from_sequences(sequences=sequences)

    assert 9 == len(mask), 'Invalid length'
    assert mask.contains('record1', 0)
    assert mask.contains('record1', 8)


def test_create_empty():
    mask = MaskedGenomicRegions.empty_mask()

    assert 0 == len(mask), 'Invalid length'
    assert mask.is_empty()
    assert not mask.contains('record1', 0)
    assert not mask.contains('record1', 3)


def test_intersection():
    sequences1 = [SeqRecord(seq=Seq('ATCG-NN'), id='record1')]
    sequences2 = [SeqRecord(seq=Seq('ATCGNNTCC'), id='record1')]
    mask1 = MaskedGenomicRegions.from_sequences(sequences=sequences1)
    mask2 = MaskedGenomicRegions.from_sequences(sequences=sequences2)

    mask_intersect = mask1.intersect(mask2)

    assert 2 == len(mask_intersect), 'Invalid length'
    assert not mask_intersect.contains('record1', 3)
    assert mask_intersect.contains('record1', 4)
    assert mask_intersect.contains('record1', 5)
    assert not mask_intersect.contains('record1', 6)


def test_union_all():
    sequences1 = [SeqRecord(seq=Seq('ATCG-NN'), id='record1')]
    sequences2 = [SeqRecord(seq=Seq('ATCGNNTCC'), id='record1')]
    mask1 = MaskedGenomicRegions.from_sequences(sequences=sequences1)
    mask2 = MaskedGenomicRegions.from_sequences(sequences=sequences2)

    mask_union = MaskedGenomicRegions.union_all([mask1, mask2])

    assert 3 == len(mask_union), 'Invalid length'
    assert not mask_union.contains('record1', 3)
    assert mask_union.contains('record1', 4)
    assert mask_union.contains('record1', 5)
    assert mask_union.contains('record1', 6)
    assert not mask_union.contains('record1', 7)


def test_union():
    sequences1 = [SeqRecord(seq=Seq('ATCG-NN'), id='record1')]
    sequences2 = [SeqRecord(seq=Seq('ATCGNNTCC'), id='record1')]
    mask1 = MaskedGenomicRegions.from_sequences(sequences=sequences1)
    mask2 = MaskedGenomicRegions.from_sequences(sequences=sequences2)

    mask_union = mask1.union(mask2)

    assert 3 == len(mask_union), 'Invalid length'
    assert not mask_union.contains('record1', 3)
    assert mask_union.contains('record1', 4)
    assert mask_union.contains('record1', 5)
    assert mask_union.contains('record1', 6)
    assert not mask_union.contains('record1', 7)


def test_write():
    sequences = [SeqRecord(seq=Seq('ATCG-NN'), id='record1')]
    mask = MaskedGenomicRegions.from_sequences(sequences=sequences)

    assert 3 == len(mask), 'Invalid length'
    assert not mask.contains('record1', 3)
    assert mask.contains('record1', 4)
    assert mask.contains('record1', 6)
    assert not mask.contains('record1', 7)

    with tempfile.NamedTemporaryFile() as f:
        file = Path(f.name)
        mask.write(file)

        mask2 = MaskedGenomicRegions.from_file(file)
        assert 3 == len(mask2), 'Invalid length'
        assert not mask2.contains('record1', 3)
        assert mask2.contains('record1', 4)
        assert mask2.contains('record1', 6)
        assert not mask2.contains('record1', 7)


def test_mask_genome():
    with tempfile.TemporaryDirectory() as tmp_dir:
        genome_file = Path(tmp_dir) / 'genome.fasta'
        genome_record = SeqRecord(id='reference', seq=Seq('ATCGAATT'))
        with open(genome_file, 'w') as f:
            SeqIO.write(genome_record, f, 'fasta')

        mask = MaskedGenomicRegions(BedTool('reference 0 4', from_string=True))
        masked_records = mask.mask_genome(genome_file, mask_char='?', remove=True)

        assert 1 == len(masked_records)
        assert 'reference' == masked_records['reference'].id
        assert Seq('AATT') == masked_records['reference'].seq


def test_mask_genome_no_remove():
    with tempfile.TemporaryDirectory() as tmp_dir:
        genome_file = Path(tmp_dir) / 'genome.fasta'
        genome_record = SeqRecord(id='reference', seq=Seq('ATCGAATT'))
        with open(genome_file, 'w') as f:
            SeqIO.write(genome_record, f, 'fasta')

        mask = MaskedGenomicRegions(BedTool('reference 0 4', from_string=True))
        masked_records = mask.mask_genome(genome_file, mask_char='?', remove=False)

        assert 1 == len(masked_records)
        assert 'reference' == masked_records['reference'].id
        assert Seq('????AATT') == masked_records['reference'].seq


def test_mask_genome_multiple_sequence():
    with tempfile.TemporaryDirectory() as tmp_dir:
        genome_file = Path(tmp_dir) / 'genome.fasta'
        genome_records = [
            SeqRecord(id='r1', seq=Seq('ATCGAATT')),
            SeqRecord(id='r2', seq=Seq('GGGGCCCCAAAA')),
        ]
        with open(genome_file, 'w') as f:
            SeqIO.write(genome_records, f, 'fasta')

        mask = MaskedGenomicRegions(BedTool('r1 2 4\nr2 4 8', from_string=True))
        masked_records = mask.mask_genome(genome_file, mask_char='?', remove=True)

        assert 2 == len(masked_records)
        assert 'r1' == masked_records['r1'].id
        assert Seq('ATAATT') == masked_records['r1'].seq
        assert 'r2' == masked_records['r2'].id
        assert Seq('GGGGAAAA') == masked_records['r2'].seq
