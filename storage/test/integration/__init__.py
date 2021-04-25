from os import path, listdir
from pathlib import Path
from typing import List

root_data_dir = Path(path.dirname(__file__), 'data')
mlst_dir = root_data_dir / 'mlst'
data_dir = root_data_dir / 'snippy'
extra_snippy_dir = root_data_dir / 'extra-snippy-samples'
regular_vcf_dir = root_data_dir / 'regular-vcf'
variation_dir = root_data_dir / 'variation'
consensus_dir = root_data_dir / 'consensus'
sourmash_dir = root_data_dir / 'sourmash'
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
reference_file = data_dir / 'genome.fasta.gz'
tree_file = data_dir / 'tree.tre'
basic_mlst_file = mlst_dir / 'mlst-basic.tsv'
mlst_file_unknown = mlst_dir / 'mlst-unknown.tsv'
sistr_mlst_file = mlst_dir / 'mlst-sistr.csv'
chewbbaca_mlst_file = mlst_dir / 'mlst-chewbbaca-small.tsv'
mlst_file_single_scheme = mlst_dir / 'mlst-single-scheme.tsv'
mlst_snippy_file = mlst_dir / 'mlst-snippy-data.csv'

data_dir_empty = root_data_dir / 'empty_vcfs'

test_project_dir = root_data_dir / 'test_project_dir'

sourmash_signatures = {
    'SampleA': sourmash_dir / 'SampleA.sig.gz',
    'SampleB': sourmash_dir / 'SampleB.sig.gz',
    'SampleC': sourmash_dir / 'SampleC.sig.gz',
}

snippy_snps_dataframes = {
    'SampleA': data_dir / 'SampleA' / 'mutations-dataframe.snps.tsv',
    'SampleB': data_dir / 'SampleB' / 'mutations-dataframe.snps.tsv',
    'SampleC': data_dir / 'SampleC' / 'mutations-dataframe.snps.tsv',
}

snippy_all_dataframes = {
    'SampleA': data_dir / 'SampleA' / 'mutations-dataframe.all.tsv',
    'SampleB': data_dir / 'SampleB' / 'mutations-dataframe.all.tsv',
    'SampleC': data_dir / 'SampleC' / 'mutations-dataframe.all.tsv',
}


def expand_list_by(list_in: List[str], number: int) -> List[str]:
    new_list = []
    for value in list_in:
        new_list.extend([value] * number)

    return new_list
