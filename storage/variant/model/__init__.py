from pathlib import Path
from ete3 import Tree
from sqlalchemy import Column, Integer, String, ForeignKey, Table, LargeBinary, UnicodeText
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from storage.variant.CoreBitMask import CoreBitMask

Base = declarative_base()


class NucleotideVariantsSamples(Base):
    __tablename__ = 'nucleotide_variants_samples'
    spdi = Column(String(255), primary_key=True)
    var_type = Column(String(255))
    _samples_bitmap = Column(LargeBinary(length=100 * 10 ** 6))  # Max of 100 million bytes

    @classmethod
    def to_spdi(cls, sequence_name: str, position: int, ref: str, alt: str) -> str:
        return f'{sequence_name}:{position}:{ref}:{alt}'

    def __repr__(self):
        return (f'<NucleotideVariantsSamples(spdi={self.spdi}, var_type={self.var_type})>')


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
    _consensus_file = Column(String(255))
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
    kmer_index_path = Column(String(255))

    sample = relationship('Sample', uselist=False, back_populates='sample_kmer_index')
