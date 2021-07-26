import logging
from typing import Dict, Any, List

import pandas as pd
from pandas.api.types import CategoricalDtype

from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger(__name__)


class InvalidSnpEffVcfError(Exception):

    def __init__(self, message):
        super().__init__(message)


class VcfSnpEffAnnotationParser:
    ANNOTATION_COLUMNS = ['ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                          'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']
    IMPACT_TYPE = CategoricalDtype(categories=['HIGH', 'MODERATE', 'LOW', 'MODIFIER'], ordered=True)

    # Order of categories derived from list provided in <http://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf>
    # But has been modified a bit
    ANNOTATION_TYPE = CategoricalDtype(categories=[
        "chromosome_number_variation",
        "exon_loss_variant",
        "frameshift_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "rare_amino_acid_variant",
        "missense_variant",
        "disruptive_inframe_insertion",
        "conservative_inframe_insertion",
        "disruptive_inframe_deletion",
        "conservative_inframe_deletion",
        "5_prime_UTR_truncation+exon_loss_variant",
        "3_prime_UTR_truncation+exon_loss",
        "splice_branch_variant",
        "splice_region_variant",
        "stop_retained_variant",
        "initiator_codon_variant",
        "synonymous_variant",
        "initiator_codon_variant+non_canonical_start_codon",
        "coding_sequence_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "5_prime_UTR_premature_start_codon_gain_variant",
        "intragenic_variant",
        "intergenic_region",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TF_binding_site_variant",
        "regulatory_region_variant",
        "miRNA",
        "custom",
        "sequence_feature",
        "conserved_intron_variant",
        "intron_variant",
        "conserved_intergenic_variant",
        "non_coding_exon_variant",
        "nc_transcript_variant",
        "gene_variant",
        "chromosome",
    ], ordered=True)

    def __init__(self):
        pass

    def parse_annotation_headers(self, vcf_info: Dict[str, Any]) -> List[str]:
        expected_start_str = "Functional annotations: '"

        if 'ANN' in vcf_info:
            ann_info = vcf_info['ANN'].desc
            if not ann_info.startswith(expected_start_str):
                raise InvalidSnpEffVcfError(f"Found 'ANN' in VCF information but description does "
                                            f"not start with [{expected_start_str}]. desc=[{ann_info}]")
            else:
                ann_info = ann_info[len(expected_start_str):]

                if ann_info.endswith("'"):
                    ann_info = ann_info[0:len(ann_info) - 1]

                return [x.strip() for x in ann_info.split('|')]
        else:
            logger.log(TRACE_LEVEL, f"VCF does not contain 'ANN' in vcf_info, keys=[{vcf_info.keys()}]")
            return []

    def _setup_vcf_df_index(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        vcf_df_with_keys = vcf_df.copy()
        if len(vcf_df_with_keys) == 0:
            vcf_df_with_keys = vcf_df_with_keys.reindex(columns=vcf_df_with_keys.columns.tolist() + ['VARIANT_ID'])
        else:
            vcf_df_with_keys['VARIANT_ID'] = vcf_df_with_keys.apply(
                lambda x: f"{x['CHROM']}:{x['POS']}:{x['REF']}:{x['ALT']}",
                axis='columns')
            vcf_df_with_keys = vcf_df_with_keys.reset_index().rename({'index': 'original_index'}, axis='columns')
            vcf_df_with_keys = vcf_df_with_keys.set_index('VARIANT_ID')

        return vcf_df_with_keys

    def _extract_ann_from_info(self, vcf_df_with_keys: pd.DataFrame, vcf_ann_headers: List[str]) -> pd.DataFrame:
        ann_groups = vcf_df_with_keys['INFO'].map(lambda x: x['ANN']).explode()
        ann_split_fields = ann_groups.map(lambda x: x.split('|') if not pd.isna(x) else [pd.NA] * len(vcf_ann_headers))
        ann_split_fields = ann_split_fields.map(lambda x: dict(zip(vcf_ann_headers, x))).to_frame()
        ann_split_fields = ann_split_fields.rename({'INFO': 'ANN'}, axis='columns').reset_index()

        return ann_split_fields

    def vcf_has_annotations(self, vcf_df: pd.DataFrame) -> bool:
        return len(vcf_df) > 0 and 'INFO' in vcf_df and 'ANN' in vcf_df.iloc[0]['INFO']

    def parse_annotation_entries(self, vcf_ann_headers: List[str], vcf_df: pd.DataFrame) -> pd.DataFrame:
        """
        Given a list of snpeff VCF 'ANN' headers and the VCF datqframe, splits up the snpeff effects into separate rows.
        Returns a new dataframe with all the snpeff 'ANN' effects, one per row. The index is the same as the input 'vcf_df'
        dataframe which means you can merge the returned value to 'vcf_df' afterwards by the indexes.
        :param vcf_ann_headers: A list of headers for the VCF 'ANN' entries (order matters here).
        :param vcf_df: The dataframe containing the VCF information. The 'ANN' information is assumed to exist
                       in the 'INFO' column.
        :return: A dataframe with the snpeff annotation entries, one entry per line. The index is the same as the input
                 vcf_df, which means that you can merge the returned dataframe with 'vcf_df' by the index. If there are
                 no snpeff annotations then this will return an empty dataframe with the appropriate snpeff ANN columns
                 names.
        """
        vcf_df_with_keys = self._setup_vcf_df_index(vcf_df)

        if self.vcf_has_annotations(vcf_df):
            ann_split_fields = self._extract_ann_from_info(vcf_df_with_keys, vcf_ann_headers)

            def insert_key_to_dictionary(x: pd.Series):
                """
                I use this to insert the 'VARIANT_KEY' into the 'ANN' dictionary
                so that when I run pd.json_normalize I can keep everything linked by the same key
                (since json_normalize removes my index in the Series).
                :param x: A series of columns.
                :return: A series of columns with the 'VARIANT_KEY' inserted into the 'ANN' dictionary.
                """
                x['ANN']['VARIANT_ID'] = x['VARIANT_ID']
                return x

            ann_series = ann_split_fields.apply(insert_key_to_dictionary, axis='columns')['ANN']
            ann_df = pd.json_normalize(ann_series).add_prefix('ANN.').rename(
                {'ANN.VARIANT_ID': 'VARIANT_ID'}, axis='columns').set_index('VARIANT_ID')

            # Copy index to column 'VARIANT_ID' so it remains after merging
            ann_df['VARIANT_ID'] = ann_df.index

            # Now join ann_df back to vcf_df_with_keys to recover the original index
            vcf_df_with_keys = vcf_df_with_keys.merge(ann_df, left_index=True, right_index=True).set_index(
                'original_index')

            vcf_df_with_keys = vcf_df_with_keys[self.ANNOTATION_COLUMNS + ['VARIANT_ID']]

            # Fill NA with '' and convert '' back to NA because a straight conversion of '' to NA
            # fails in cases where a column has a mixture of '' and NA values.
            vcf_df_with_keys = vcf_df_with_keys.fillna('').replace('', pd.NA)
        else:
            logger.log(TRACE_LEVEL, 'vcf_df has no snpeff annotations, will set all annotations as NA')
            vcf_df_with_keys = vcf_df_with_keys.reset_index().reindex(
                columns=self.ANNOTATION_COLUMNS + ['original_index', 'VARIANT_ID']).set_index('original_index')

        return vcf_df_with_keys

    def select_variant_annotations(self, vcf_df: pd.DataFrame) -> pd.DataFrame:
        """
        Given a dataframe containing variants from a VCF file, select the appropriate annotations.
        :param vcf_df: The dataframe of variants from a VCF file.
        :return: A new dataframe containing only rows with the selected annotations.
        """
        non_annotation_columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID']
        expected_columns_set = set(non_annotation_columns + self.ANNOTATION_COLUMNS)
        actual_columns_set = set(vcf_df.columns.tolist())
        if expected_columns_set != actual_columns_set:
            raise Exception(f'Invalid columns for passed vcf_df.'
                            f' Unique expected={expected_columns_set - actual_columns_set}.'
                            f' Unique actual = {actual_columns_set - expected_columns_set}')

        # Use .copy() to get rid of warnings like "A value is trying to be set on a copy of a slice from a DataFrame"
        vcf_df_with_annotations = vcf_df[~vcf_df['ANN.Allele'].isna()].copy()

        # Set ordered types
        vcf_df_with_annotations['ANN.Annotation_Impact'] = vcf_df_with_annotations['ANN.Annotation_Impact'].astype(
            self.IMPACT_TYPE
        )
        vcf_df_with_annotations['ANN.Annotation'] = vcf_df_with_annotations['ANN.Annotation'].astype(
            self.ANNOTATION_TYPE
        )

        # Remove annotations from being considered
        vcf_df_with_annotations = vcf_df_with_annotations[
            vcf_df_with_annotations['ANN.Allele'] == vcf_df_with_annotations['ALT']]

        # Order and select first entry (use nth() instead of first() so that NA values are handled properly)
        vcf_df_with_annotations = vcf_df_with_annotations.sort_values(
            ['ANN.Annotation_Impact', 'ANN.Annotation', 'ANN.HGVS.c']).groupby('VARIANT_ID').nth(0).reset_index()

        # Merge back with original dataframe of variants to make sure I include those without annotations
        all_vcf_entries_grouped = vcf_df[non_annotation_columns].groupby('VARIANT_ID').first()
        return_vcf_df = all_vcf_entries_grouped.merge(vcf_df_with_annotations, how='left', left_on='VARIANT_ID',
                                                      right_on='VARIANT_ID',
                                                      suffixes=(None, '_with_annotations'))

        return_vcf_df = return_vcf_df[non_annotation_columns + self.ANNOTATION_COLUMNS].sort_values(
            'POS')

        return return_vcf_df
