from os import path, listdir
from pathlib import Path
from typing import List

root_data_dir = Path(path.dirname(__file__), 'data')
mlst_dir = root_data_dir / 'mlst'
data_dir = root_data_dir / 'snippy'
snpeff_data_dir = root_data_dir / 'snpeff-vcf'
snpeff_database_dir = root_data_dir / 'snpeff-database'
extra_snippy_dir = root_data_dir / 'extra-snippy-samples'
regular_vcf_dir = root_data_dir / 'regular-vcf'
variation_dir = root_data_dir / 'variation'
consensus_dir = root_data_dir / 'consensus'
sourmash_dir = root_data_dir / 'sourmash'
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
sample_dirs_AB = [data_dir / 'SampleA', data_dir / 'SampleB']
sample_dirs_C = [data_dir / 'SampleC']
sample_dirs_A = [data_dir / 'SampleA']
sample_dirs_BC = [data_dir / 'SampleB', data_dir / 'SampleC']
reference_file = data_dir / 'genome.fasta.gz'
reference_file_snpeff = snpeff_data_dir / 'NC_011083.fasta.gz'
reference_file_5000_snpeff = snpeff_database_dir / 'NC_011083-5000.gbk.gz'
reference_file_5000_snpeff_2 = snpeff_database_dir / 'NC_011083_CP001602-5000.gbk.gz'
snpeff_vcf_file = snpeff_database_dir / 'SampleA.vcf.gz'
tree_file = data_dir / 'tree.tre'
snpeff_tree_file = snpeff_data_dir / 'tree.txt'
basic_mlst_file = mlst_dir / 'mlst-basic.tsv'
mlst_file_unknown = mlst_dir / 'mlst-unknown.tsv'
sistr_mlst_file = mlst_dir / 'mlst-sistr.csv'
chewbbaca_mlst_file = mlst_dir / 'mlst-chewbbaca-small.tsv'
mlst_file_single_scheme = mlst_dir / 'mlst-single-scheme.tsv'
mlst_file_single_scheme2 = mlst_dir / 'mlst-single-scheme2.tsv'
mlst_file_single_scheme3 = mlst_dir / 'mlst-single-scheme3.tsv'
mlst_snippy_file = mlst_dir / 'mlst-snippy-data.csv'

data_dir_empty = root_data_dir / 'empty_vcfs'

test_project_dir = root_data_dir / 'test_project_dir'

snippy_sample_vcfs_dict = {
    'SampleA': data_dir / 'SampleA' / 'snps.vcf.gz',
    'SampleB': data_dir / 'SampleB' / 'snps.vcf.gz',
    'SampleC': data_dir / 'SampleC' / 'snps.vcf.gz',
}

snippy_sample_mask_sequences_dict = {
    'SampleA': data_dir / 'SampleA' / 'snps.aligned.fa',
    'SampleB': data_dir / 'SampleB' / 'snps.aligned.fa',
    'SampleC': data_dir / 'SampleC' / 'snps.aligned.fa',
}

snippy_sample2_vcfs_dict = {
    'SampleA_2': data_dir / 'SampleA' / 'snps.vcf.gz',
    'SampleB_2': data_dir / 'SampleB' / 'snps.vcf.gz',
    'SampleC_2': data_dir / 'SampleC' / 'snps.vcf.gz',
}

snippy_sample2_mask_sequences_dict = {
    'SampleA_2': data_dir / 'SampleA' / 'snps.aligned.fa',
    'SampleB_2': data_dir / 'SampleB' / 'snps.aligned.fa',
    'SampleC_2': data_dir / 'SampleC' / 'snps.aligned.fa',
}

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

snpeff_sample_vcfs = {
    'SH10-014': snpeff_data_dir / 'SH10-014.vcf.gz',
    'SH14-001': snpeff_data_dir / 'SH14-001.vcf.gz',
    'SH14-014': snpeff_data_dir / 'SH14-014.vcf.gz',
}

snpeff_sample_vcfs_fake_dup = {
    'SH10-014-dup-gene-variant': snpeff_data_dir / 'SH10-014-dup-gene-variant.vcf.gz',
    'SH10-014-dup-gene-variant-2': snpeff_data_dir / 'SH10-014-dup-gene-variant-2.vcf.gz',
}


def expand_list_by(list_in: List[str], number: int) -> List[str]:
    new_list = []
    for value in list_in:
        new_list.extend([value] * number)

    return new_list
