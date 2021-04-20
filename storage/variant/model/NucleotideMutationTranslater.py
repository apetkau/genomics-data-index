import re
from typing import Union, Tuple

from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation


class NucleotideMutationTranslater:

    @classmethod
    def convert_deletion(cls, deletion: Union[str, int]) -> int:
        if isinstance(deletion, str):
            if re.match(r'^\d+$', deletion):
                deletion = int(deletion)
            else:
                if not set(deletion).issubset({'A', 'T', 'C', 'G'}):
                    raise Exception('Deletion must either be an integer or a string with alphabet {A,T,C,G}'
                                    f': {deletion}')
                deletion = len(deletion)
        elif isinstance(deletion, int):
            if deletion < 0:
                raise Exception(f'ref=[{deletion}] must be a non-negative integer')
        else:
            raise Exception(f'ref=[{deletion}] must be either a string or a non-negative integer')

        return deletion

    @classmethod
    def from_spdi(cls, spdi: str, convert_deletion = True) -> Tuple[str, int, Union[int, str], str]:
        if spdi is None:
            raise Exception('Cannot parse value spdi=None')

        values = spdi.split(':')
        if len(values) != 4:
            raise Exception(f'Incorrect number of items for spdi=[{spdi}]')
        else:
            position = int(values[1])

            if convert_deletion:
                deletion = cls.convert_deletion(values[2])
            else:
                deletion = values[2]

            if position < 0:
                raise Exception(f'Position must be non-negative: {position}')

            return str(values[0]), position, deletion, str(values[3])

    @classmethod
    def to_spdi(cls, sequence_name: str, position: int, ref: Union[str, int], alt: str) -> str:
        if position < 0:
            raise Exception(f'Position must be non-negative: {position}')

        ref = cls.convert_deletion(ref)
        return f'{sequence_name}:{position}:{ref}:{alt}'

    @classmethod
    def to_db_feature(cls, feature: QueryFeatureMutation) -> QueryFeatureMutation:
        new_id = cls.to_spdi(feature.scope, feature.position, feature.ref, feature.alt)
        return QueryFeatureMutation(new_id)
