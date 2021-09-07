import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import List, Dict, Set, Any, Union

from Bio.SeqRecord import SeqRecord

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
    TODO: I need to review code in here https://github.com/cov-lineages/scorpio/blob/main/scorpio/scripts/type_constellations.py
          to make sure I match what's going on there.
    TODO: There are still some bugs with this. In particular, I've noticed typing B.1.1.7 doesn't work due to different standards in naming
          mutations between HGVS notation and names given in the constellation files (S:p.Y145del vs. s:Y144-).
    """

    def __init__(self, constellation_files: Union[List[Path], List[str]],
                 sequence: SeqRecord, perfect_match_only: bool = False,
                 valid_gene_identifiers: Set[str] = None):
        super().__init__()

        sequence_name = str(sequence.id)
        if sequence_name not in ['NC_045512.2', 'NC_045512']:
            logger.warning(
                f"sequence_name=[{sequence_name}] is not 'NC_045512.2' or 'NC_045512'. This is what the constellations are with respect"
                "to <https://github.com/cov-lineages/constellations/blob/main/constellations/data/SARS-CoV-2.json>")

        self._sequence_name = sequence_name
        self._sequence = sequence
        if valid_gene_identifiers is None:
            self._valid_gene_identifiers = {'ORF1ab', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N',
                                            'ORF10'}
        self._typing_definitions = self._parse_definitions(constellation_files)
        self._perfect_match_only = perfect_match_only

    def _handle_deletion_identifier(self, mutation_values: List[str], mutation: str) -> str:
        if mutation_values[0] == 'del':
            # Based on looking at some real-world data (for a single deletion identifier), it looks like the
            # position here is the position of the unchanged base, with the deletion occuring right after.
            position = int(mutation_values[1])
            position_0coord = position - 1
            length = int(mutation_values[2])

            if position_0coord < 0:
                raise MutationParsingError(f'Cannot parse mutation [{mutation}], position is negative')
            else:
                # I need to specify an insertion that's not of length 0 (for my use of these identifiers).
                # For example, say the sequence is "ATCG" and position=1 (position right before deletion),
                # and deletion length = 2. Then this means deleting "A[TC]G" and inserting nothing.
                # But to convert this to an identifier where I have at least 1 insertion, I need to delete
                # and re-insert one of the reference characters. I do this with the left-most base, so the actual
                # deletion becomes "[ATC]G" followed by an insertion of "A". That is, the identifier
                # becomes "sequence:ATC:A".
                deletion_start = position_0coord
                deletion_stop = deletion_start + (length + 1)
                position_start = position

            deletion_sequence = self._sequence.seq[deletion_start:deletion_stop]
            insertion_sequence = self._sequence.seq[deletion_start:deletion_start + 1]

            return f'{self._sequence_name}:{position_start}:{deletion_sequence}:{insertion_sequence}'
        else:
            raise MutationParsingError(f'Cannot parse mutation [{mutation}], has three parts but first '
                                       f'is not "del".')

    def _handle_gene_nucleotide_identifier(self, mutation_values: List[str], mutation: str) -> str:
        gene, mutation = mutation_values
        if gene == 'nuc':
            match = re.match(r'([ATCG]+)(\d+)([ATCG]+)', mutation)
            if not match:
                raise MutationParsingError(f'Could not parse mutation={mutation}, does not match pattern like [A150T].')
            mutation_id = f'{self._sequence_name}:{match.group(2)}:{match.group(1)}:{match.group(3)}'
        else:
            # Curation of different possible gene identifiers
            gene_curation_map = {
                's': 'S',
                '1ab': 'ORF1ab',
                'orf1ab': 'ORF1ab',
                'ORF1a': 'ORF1ab',
                'ORF1b': 'ORF1ab',
                'NSP2': 'ORF1ab',
                'nsp3': 'ORF1ab',
                'nsp4': 'ORF1ab',
                'nsp5': 'ORF1ab',
                'nsp6': 'ORF1ab',
                'nsp7': 'ORF1ab',
                'nsp12': 'ORF1ab',
                'nsp13': 'ORF1ab',
                'nsp15': 'ORF1ab',
                '8': 'ORF8'
            }
            if gene in gene_curation_map:
                gene = gene_curation_map[gene]
            if gene not in self._valid_gene_identifiers:
                raise MutationParsingError(
                    f'gene=[{gene}] for mutation=[{mutation}] is not one of the valid '
                    f"gene identifiers={self._valid_gene_identifiers}.")
            if mutation.endswith('-'):
                mutation = self._convert_aa_deletion(mutation)
            mutation_id = f'hgvs_gn:{self._sequence_name}:{gene}:p.{mutation}'
        return mutation_id

    def mutation_to_identifier(self, mutation: str) -> str:
        values = mutation.split(':')
        if len(values) == 2:
            return self._handle_gene_nucleotide_identifier(values, mutation=mutation)
        elif len(values) == 3:
            return self._handle_deletion_identifier(values, mutation=mutation)
        else:
            raise MutationParsingError(f'Cannot parse mutation [{mutation}], cannot be split by ":".')

    def _convert_aa_deletion(self, mutation: str) -> str:
        match = re.match(r'(\D+)(\d+)-$', mutation)
        if not match:
            raise MutationParsingError(f'Could not parse mutation={mutation}, does not match pattern like [HV69-].')
        aa_values = match.group(1)
        position = match.group(2)
        if len(aa_values) == 1:
            return f'{aa_values}{position}del'
        else:
            first_aa = aa_values[0]
            last_aa = aa_values[-1]
            last_position = int(position) + len(aa_values) - 1
            return f'{first_aa}{position}_{last_aa}{last_position}del'

    def _parse_definitions(self, constellation_files: Union[List[Path], List[str]]) -> Dict[str, Any]:
        typing_definitions = dict()
        for file in constellation_files:
            if isinstance(file, str):
                file = Path(file)

            try:
                constellation_info = self._load_definitions(file)
                signature_mutations_raw = constellation_info['sites']
                signature_mutation_ids = {self.mutation_to_identifier(x) for x in signature_mutations_raw}
                logger.debug(f'constellation_file=[{file}], mutation_ids={signature_mutation_ids}')

                rules = constellation_info['rules']
                min_alt = rules.get('min_alt', len(signature_mutation_ids))
                max_ref = rules.get('max_ref', 0)
                other_rules = set(rules.keys()) - {'min_alt', 'max_ref'}

                # Handle any other special-rules
                must_have_mutations = set()
                for rule_key in other_rules:
                    if rules[rule_key] == 'alt':
                        mutation_id = self.mutation_to_identifier(rule_key)
                        must_have_mutations.add(mutation_id)
                    else:
                        logger.warning(f'Skipping unknown rule {rule_key}={rules[rule_key]}, file=[{file.name}]')

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
                            'signature_mutations_original': signature_mutations_raw,
                            'must_have_mutations': must_have_mutations,
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

    def perfect_matches(self, type_name: str, query: SamplesQuery,
                        signature_mutations: Set[str] = None) -> SamplesQuery:
        self._validate_type_name(type_name)
        if signature_mutations is None:
            signature_mutations = self._typing_definitions[type_name]['signature_mutations']

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
        must_have_mutations = self._typing_definitions[data]['must_have_mutations']
        min_alt = self._typing_definitions[data]['min_alt']

        # If any must_have mutations, handle these first to try to eliminate as much as possible ahead of time
        if len(must_have_mutations) > 0:
            query = self.perfect_matches(type_name=data, query=query, signature_mutations=must_have_mutations)

        if self._perfect_match_only or min_alt == len(signature_mutations):
            # Separate case for all perfect matches since I can take advantage of tracking unknowns
            return self.perfect_matches(type_name=data, query=query)
        else:
            # SARS-CoV-2 typing with constellations does not have the concept of an "unknown" sample and
            # instead tries to classify these using e.g., min_alt and max_ref, and so on
            # So, for the purposes of typing, I want to ignore any "unknown" samples for typing and add
            # the unknown set back in at the end. The is because no matter what the result of the sequence-typing
            # code below, the statement 'query.isa("Type")' will still be unknown if the sample was already unknown
            # in the input query (for example, if Sample A is already "Unknown" in the input query and the result
            # below is "True" for Sample A, then the final result for Sample A is still "Unknown" since
            # "Unknown" AND "True" = "Unknown").
            query_no_unknowns = query.select_present().reset_universe()

            # Look for perfect matches
            query_perfect_match = self.perfect_matches(type_name=data,
                                                       query=query_no_unknowns,
                                                       signature_mutations=signature_mutations)

            # Look for imperfect matches in what's left over
            query_imperfect_matches = query_perfect_match.select_unknown() | query_perfect_match.select_absent()
            query_imperfect_matches = self.imperfect_matches(data, query_imperfect_matches)

            # and query to add back any existing unknown samples
            return query & (query_perfect_match.select_present() | query_imperfect_matches.select_present())
