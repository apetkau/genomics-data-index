from typing import Tuple, Union
from pathlib import Path
from ete3 import Tree
from sqlalchemy import Column, Integer, String, ForeignKey, Table, LargeBinary, UnicodeText
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.SampleSet import SampleSet

Base = declarative_base()


class NucleotideVariantsSamples(Base):
    __tablename__ = 'nucleotide_variants_samples'
    sequence = Column(String(255), primary_key=True)
    position = Column(Integer, primary_key=True)
    deletion = Column(String(255), primary_key=True)
    insertion = Column(String(255), primary_key=True)
    _spdi = Column('spdi', String(255))
    var_type = Column(String(255))
    _sample_ids = Column(LargeBinary(length=500 * 10 ** 6))  # Max of 500 million

    def __init__(self, spdi: str = None, var_type: str = None, sample_ids: SampleSet = None):
        self.spdi = spdi
        self.var_type = var_type
        self.sample_ids = sample_ids

    @hybrid_property
    def spdi(self) -> str:
        return self.to_spdi(self.sequence, self.position, self.deletion, self.insertion)

    @spdi.setter
    def spdi(self, spdi_value: str):
        self.sequence, self.position, self.deletion, self.insertion = self.from_spdi(spdi_value)
        self._spdi = spdi_value

    @hybrid_property
    def sample_ids(self) -> SampleSet:
        if self._sample_ids is None:
            raise Exception('_sample_ids is not set')
        else:
            return SampleSet.from_bytes(self._sample_ids)

    @sample_ids.setter
    def sample_ids(self, sample_ids: SampleSet) -> None:
        if sample_ids is None:
            raise Exception('Cannot set sample_ids to None')
        else:
            self._sample_ids = sample_ids.get_bytes()

    @classmethod
    def to_spdi(cls, sequence_name: str, position: int, ref: str, alt: str) -> str:
        return f'{sequence_name}:{position}:{ref}:{alt}'

    @classmethod
    def from_spdi(cls, spdi: str) -> Tuple[str, int, str, str]:
        if spdi is None:
            raise Exception('Cannot parse value spdi=None')

        values = spdi.split(':')
        if len(values) != 4:
            raise Exception(f'Incorrect number of items for spdi=[{spdi}]')
        else:
            return str(values[0]), int(values[1]), str(values[2]), str(values[3])

    def __repr__(self):
        return (f'<NucleotideVariantsSamples(spdi={self.spdi}, var_type={self.var_type}, num_samples={len(self.sample_ids)})>')


class Reference(Base):
    __tablename__ = 'reference'
    id = Column(Integer, primary_key=True)
    name = Column(String(255))
    length = Column(Integer)
    _tree = Column('tree', UnicodeText(10 ** 6))
    tree_alignment_length = Column(Integer)

    sequences = relationship('ReferenceSequence')
    sample_nucleotide_variation = relationship('SampleNucleotideVariation', back_populates='reference')

    @hybrid_property
    def tree(self) -> Tree:
        if self._tree is None:
            raise Exception('Cannot convert an empty tree')
        else:
            return Tree(self._tree)

    @tree.setter
    def tree(self, in_tree: Tree):
        if in_tree is None:
            self._tree = None
        else:
            self._tree = in_tree.write()

    def __repr__(self):
        return f'<Reference(id={self.id}, name={self.name}, length={self.length})>'


class SampleNucleotideVariation(Base):
    __tablename__ = 'sample_nucleotide_variation'
    sample_id = Column(Integer, ForeignKey('sample.id'), primary_key=True)
    reference_id = Column(Integer, ForeignKey('reference.id'), primary_key=True)
    _masked_regions_file = Column('masked_regions_file', String(255))
    _nucleotide_variants_file = Column('nucleotide_variants_file', String(255))

    @hybrid_property
    def nucleotide_variants_file(self) -> Path:
        if self._nucleotide_variants_file is None:
            raise Exception('Empty _nucleotide_variants_file')
        else:
            return Path(self._nucleotide_variants_file)

    @nucleotide_variants_file.setter
    def nucleotide_variants_file(self, file: Path) -> None:
        if file is None:
            self._nucleotide_variants_file = None
        else:
            self._nucleotide_variants_file = str(file)

    @property
    def masked_regions(self) -> MaskedGenomicRegions:
        return MaskedGenomicRegions.from_file(self.masked_regions_file)

    @hybrid_property
    def masked_regions_file(self) -> Path:
        if self._masked_regions_file is None:
            raise Exception('Empty _masked_regions_file')
        else:
            return Path(self._masked_regions_file)

    @masked_regions_file.setter
    def masked_regions_file(self, file: Path) -> None:
        if file is None:
            self._masked_regions_file = None
        else:
            self._masked_regions_file = str(file)

    sample = relationship('Sample', back_populates='sample_nucleotide_variation')
    reference = relationship('Reference', back_populates='sample_nucleotide_variation')


class ReferenceSequence(Base):
    __tablename__ = 'reference_sequence'
    id = Column(Integer, primary_key=True)
    reference_id = Column(Integer, ForeignKey('reference.id'))
    sequence_name = Column(String(255))
    sequence_length = Column(Integer)

    def __repr__(self):
        return (f'<ReferenceSequence(id={self.id}, sequence_name={self.sequence_name},'
                f'sequence_length={self.sequence_length}, reference_id={self.reference_id})>')


class Sample(Base):
    __tablename__ = 'sample'
    id = Column(Integer, primary_key=True)
    name = Column(String(255))

    sample_nucleotide_variation = relationship('SampleNucleotideVariation', back_populates='sample')
    sample_kmer_index = relationship('SampleKmerIndex', uselist=False, back_populates='sample')

    def __repr__(self):
        return f'<Sample(id={self.id}, name={self.name})>'


class SampleKmerIndex(Base):
    __tablename__ = 'sample_kmer_index'
    sample_id = Column(Integer, ForeignKey('sample.id'), primary_key=True)
    _kmer_index_path = Column('kmer_index_path', String(255))

    sample = relationship('Sample', uselist=False, back_populates='sample_kmer_index')

    def __init__(self, sample: Sample, kmer_index_path: Union[str, Path]):
        self.sample = sample
        self.kmer_index_path = kmer_index_path

    @hybrid_property
    def kmer_index_path(self) -> Path:
        return Path(self._kmer_index_path)

    @kmer_index_path.setter
    def kmer_index_path(self, file: Union[str, Path]):
        if isinstance(file, Path):
            self._kmer_index_path = str(file)
        else:
            self._kmer_index_path = file
