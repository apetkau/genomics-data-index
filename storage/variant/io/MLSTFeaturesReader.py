from typing import Set

from storage.variant.io.FeaturesReader import FeaturesReader


class MLSTFeaturesReader(FeaturesReader):

    def __init__(self):
        super().__init__()

    def _minimal_expected_columns(self) -> Set[str]:
        return {'Sample', 'Scheme', 'Locus', 'Allele'}

    def get_or_create_feature_file(self, sample_name: str):
        raise Exception('Not implemented')

    def get_scheme_for_sample(self, sample_name: str):
        mlst_df = self.get_features_table()
        scheme_set = set(mlst_df.loc[mlst_df['Sample'] == sample_name, 'Scheme'].values)
        if len(scheme_set) == 0:
            raise Exception(f'No found schemes in the MLST table for sample [{sample_name}]')
        elif len(scheme_set) == 1:
            return scheme_set.pop()
        else:
            raise Exception(f'More than one schemes found in the MLST table for sample '
                            f'[{sample_name}]: [{scheme_set}]')
