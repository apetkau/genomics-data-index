from ete3 import Tree
from sqlalchemy import Column, Integer, String, ForeignKey, Table, LargeBinary, UnicodeText
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship

from storage.variant.CoreBitMask import CoreBitMask

Base = declarative_base()

sample_variation_association = Table('sample_variation_allele', Base.metadata,
                                     Column('sample_id', Integer, ForeignKey('sample.id')),
                                     Column('variantion_allele_id', String(255), ForeignKey('variation_allele.id')),
                                     )


class VariationAllele(Base):
    __tablename__ = 'variation_allele'
    id = Column(String(255), primary_key=True)
    sequence_id = Column(Integer, ForeignKey('reference_sequence.id'))
    position = Column(Integer)
    ref = Column(String(255))
    alt = Column(String(255))
    var_type = Column(String(255))

    samples = relationship('Sample', secondary=sample_variation_association, back_populates='variants')
    sequence = relationship('ReferenceSequence', back_populates='variants')

    def __init__(self, sequence=None, position: int = -1, ref: str = None, alt: str = None,
                 var_type: str = None):
        self.sequence = sequence
        self.position = position
        self.ref = ref
        self.alt = alt
        self.var_type = var_type

        self.id = self.to_spdi()

    def to_spdi(self):
        return VariationAllele.spdi(sequence_name=self.sequence.sequence_name,
                                    position=self.position,
                                    ref=self.ref,
                                    alt=self.alt)

    @classmethod
    def spdi(cls, sequence_name: str, position: int, ref: str, alt: str) -> str:
        return f'{sequence_name}:{position}:{ref}:{alt}'

    def __repr__(self):
        return (f'<VariationAllele(sequence_name={self.sequence.sequence_name}'
                f', position={self.position}, ref={self.ref}, alt={self.alt}, var_type={self.var_type})>')


class Reference(Base):
    __tablename__ = 'reference'
    id = Column(Integer, primary_key=True)
    name = Column(String(255))
    length = Column(Integer)
    _tree = Column('tree', UnicodeText(10 ** 6))
    tree_alignment_length = Column(Integer)
    sequences = relationship('ReferenceSequence')

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


class SampleSequence(Base):
    __tablename__ = 'sample_sequence'
    sample_id = Column(Integer, ForeignKey('sample.id'), primary_key=True)
    sequence_id = Column(Integer, ForeignKey('reference_sequence.id'), primary_key=True)
    _core_mask = Column(LargeBinary(length=100 * 10 ** 6))  # Max of 100 million
    flag = Column(String(255))

    sequence = relationship('ReferenceSequence', back_populates='sample_sequences')
    sample = relationship('Sample', back_populates='sample_sequences')

    @hybrid_property
    def core_mask(self) -> CoreBitMask:
        if self._core_mask is None:
            raise Exception('core_mask is not set')
        else:
            return CoreBitMask.from_bytes(self._core_mask, self.sequence.sequence_length)

    @core_mask.setter
    def core_mask(self, core_mask: CoreBitMask) -> None:
        if core_mask is None:
            raise Exception('Cannot set core_mask to None')
        elif self.sequence is None or self.sequence.sequence_length is None:
            raise Exception(f'Cannot set core_mask without the corresponding sequence (and length) set')
        elif self.sequence.sequence_length != len(core_mask):
            raise Exception(f'Cannot set core_mask, len(core_mask)=[{len(core_mask)}] '
                            f'is not the same as sequence_lenght=[{self.sequence.sequence_length}]')
        else:
            self._core_mask = core_mask.get_bytes()

    def __repr__(self):
        return f'<SampleSequence(sample_id={self.sample_id}, sequence_id={self.sequence_id}, flag={self.flag})>'


class ReferenceSequence(Base):
    __tablename__ = 'reference_sequence'
    id = Column(Integer, primary_key=True)
    reference_id = Column(Integer, ForeignKey('reference.id'))
    sequence_name = Column(String(255))
    sequence_length = Column(Integer)

    variants = relationship('VariationAllele', back_populates='sequence')
    sample_sequences = relationship('SampleSequence', back_populates='sequence')

    def __repr__(self):
        return (f'<ReferenceSequence(id={self.id}, sequence_name={self.sequence_name},'
                f'sequence_length={self.sequence_length}, reference_id={self.reference_id})>')


class Sample(Base):
    __tablename__ = 'sample'
    id = Column(Integer, primary_key=True)
    name = Column(String(255))

    variants = relationship('VariationAllele', secondary=sample_variation_association, back_populates='samples')
    sample_sequences = relationship('SampleSequence', back_populates='sample')
    sample_kmer_index = relationship('SampleKmerIndex', uselist=False, back_populates='sample')

    def __repr__(self):
        return f'<Sample(id={self.id}, name={self.name})>'


class SampleKmerIndex(Base):
    __tablename__ = 'sample_kmer_index'
    sample_id = Column(Integer, ForeignKey('sample.id'), primary_key=True)
    kmer_index_path = Column(String(255))

    sample = relationship('Sample', uselist=False, back_populates='sample_kmer_index')
