from typing import List, Dict, Set, Any, Union
from pathlib import Path
import logging
import json
from collections import Counter
import re

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.kind.isa.typing.SamplesTypingIsaKind import SamplesTypingIsaKind
from genomics_data_index.storage.SampleSet import SampleSet

logger = logging.getLogger(__name__)


class MutationParsingError(Exception):

    def __init__(self, msg):
        super().__init__(msg)


class ExperimentalSARSCov2ConstellationsTyper(SamplesTypingIsaKind):
    """
    This is an experimental/test class to try out typing SARS-CoV-2 based on constellations of variants/mutations.
    This parses constellation definitions from https://github.com/cov-lineages/constellations/blob/main/constellations/definitions
    """

    def __init__(self, constellation_files: Union[List[Path], List[str]],
                 sequence_name: str = 'MN908947.3', valid_gene_identifiers: Set[str] = None):
        super().__init__()

        self._sequence_name = sequence_name
        if valid_gene_identifiers is None:
            self._valid_gene_identifiers = {'orf1ab', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF8', 'N', 'ORF10'}
        self._typing_definitions = self._parse_definitions(constellation_files)

    def _mutation_to_identifier(self, mutation: str) -> str:
        values = mutation.split(':')
        if len(values) == 2:
            gene, mutation = values
        else:
            raise MutationParsingError(f'Cannot parse mutation [{mutation}], cannot be split by ":".')

        if gene == 'nuc':
            match = re.match(r'([ATCG]+)(\d+)([ATCG]+)', mutation)
            if not match:
                raise MutationParsingError(f'Could not parse mutation={mutation}, does not match pattern like [A150T].')
            sequence_id = f'n.{match.group(2)}{match.group(1)}>{match.group(3)}'
        else:
            # Curation of different possible gene identifiers
            gene_curation_map = {
                's': 'S',
                '1ab': 'orf1ab',
                'ORF1a': 'orf1ab',
                'ORF1ab': 'orf1ab',
                'ORF1b': 'orf1ab',
                '8': 'ORF8'
            }
            if gene in gene_curation_map:
                gene = gene_curation_map[gene]
            if gene not in self._valid_gene_identifiers:
                raise MutationParsingError(
                    f'gene=[{gene}] for mutation=[{mutation}] is not one of the valid '
                    f"gene identifiers={self._valid_gene_identifiers}.")
            sequence_id = f'{gene}:p.{mutation}'

        return f'hgvs_gn:{self._sequence_name}:{sequence_id}'

    def _parse_definitions(self, constellation_files: Union[List[Path], List[str]]) -> Dict[str, Any]:
        typing_definitions = dict()
        for file in constellation_files:
            if isinstance(file, str):
                file = Path(file)

            try:
                constellation_info = self._load_definitions(file)
                signature_mutations_raw = constellation_info['sites']
                signature_mutation_ids = {self._mutation_to_identifier(x) for x in signature_mutations_raw}
                logger.debug(f'constellation_file=[{file}], mutation_ids={signature_mutation_ids}')

                min_alt = constellation_info['rules'].get('min_alt', len(signature_mutation_ids))
                max_ref = constellation_info['rules'].get('max_ref', 0)
                labels = constellation_info['tags']
                labels = labels + [constellation_info['label']]

                for label in labels:
                    if label in typing_definitions:
                        previous_definition = typing_definitions[label]
                        previous_file = previous_definition['file']
                        logger.warning(f'Attempting to set typing definition for tag/label=[{label}], '
                                       f'file=[{file.name}], '
                                       f'but it has already been set in file [{previous_file.name}]. '
                                       f'Will ignore definition for tag/label=[{label}], file=[{file.name}].')
                    else:
                        typing_definitions[label] = {
                            'signature_mutations': signature_mutation_ids,
                            'min_alt': min_alt,
                            'max_ref': max_ref,
                            'file': file,
                            'tags': labels,
                        }
            except MutationParsingError as e:
                logger.error(f'Skipping file [{file}]: {e}')
        return typing_definitions

    def _load_definitions(self, constellation_file: Path) -> Dict[str, Any]:
        with open(constellation_file, 'r') as f:
            data = json.load(f)
        return data

    @property
    def name(self) -> str:
        return 'constellations_typer'

    @property
    def version(self) -> str:
        return 'experimental'

    def _validate_type_name(self, type_name: str) -> None:
        if type_name not in self._typing_definitions:
            raise Exception(f'type_name=[{type_name}] is not a valid type.')

    @property
    def type_names(self) -> List[str]:
        return list(self._typing_definitions.keys())

    def perfect_matches(self, data: str, query: SamplesQuery) -> SamplesQuery:
        self._validate_type_name(data)
        signature_mutations = self._typing_definitions[data]['signature_mutations']

        q = query
        for mutation in signature_mutations:
            q = q.hasa(mutation, kind='mutation')

        return q

    def _imperfect_pass_criteria(self, mutation_present_count: int, mutation_unknown_count: int,
                                 total_mutations: int,
                                 min_alt: int, max_ref: int) -> bool:
        number_ref = total_mutations - mutation_present_count - mutation_unknown_count
        return (mutation_present_count >= min_alt) and (number_ref <= max_ref)

    def imperfect_matches(self, data: str, query: SamplesQuery) -> SamplesQuery:
        self._validate_type_name(data)
        signature_mutations = self._typing_definitions[data]['signature_mutations']
        total_mutations = len(signature_mutations)
        min_alt = self._typing_definitions[data]['min_alt']
        max_ref = self._typing_definitions[data]['max_ref']

        # Use a Counter to keep track of number of mutations associated with each sample
        samples_present_counter = Counter(query.sample_set)
        sample_unknown_counter = Counter(query.unknown_set)

        for mutation in signature_mutations:
            mutation_query = query.hasa(mutation, kind='mutation')
            samples_present_counter.update(mutation_query.sample_set)
            sample_unknown_counter.update(mutation_query.unknown_set)

        sample_ids_pass = {s for s in samples_present_counter if
                           self._imperfect_pass_criteria(samples_present_counter[s],
                                                         sample_unknown_counter[s],
                                                         total_mutations=total_mutations,
                                                         min_alt=min_alt, max_ref=max_ref)}
        sample_set_present = SampleSet(sample_ids_pass)

        return query.intersect(sample_set_present)

    def isa_type(self, data: str, query: SamplesQuery) -> SamplesQuery:
        self._validate_type_name(data)
        signature_mutations = self._typing_definitions[data]['signature_mutations']
        min_alt = self._typing_definitions[data]['min_alt']

        # I first look for perfect matches to all mutations to avoid having to count mutations for these
        query_perfect_match = self.perfect_matches(data, query)
        if min_alt < len(signature_mutations):
            query_imperfect_matches = query.reset_universe() & (~(query_perfect_match.reset_universe()))
            # Now I have to look for imperfect matches to mutations and include them in my final result
            query_imperfect_matches = self.imperfect_matches(data, query_imperfect_matches)
            query_type_results = query_perfect_match | query_imperfect_matches
        else:
            query_type_results = query_perfect_match

        return query_type_results
