from os import path, listdir
from pathlib import Path

root_data_dir = Path(path.dirname(__file__), '..', 'data')
mlst_dir = root_data_dir / 'mlst'
data_dir = root_data_dir / 'snippy'
regular_vcf_dir = root_data_dir / 'regular-vcf'
variation_dir = root_data_dir / 'variation'
consensus_dir = root_data_dir / 'consensus'
sourmash_dir = root_data_dir / 'sourmash'
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
reference_file = data_dir / 'genome.fasta.gz'
tree_file = data_dir / 'tree.tre'
basic_mlst_file = mlst_dir / 'mlst-basic.tsv'

data_dir_empty = root_data_dir / 'empty_vcfs'

sourmash_signatures = {
    'SampleA': sourmash_dir / 'SampleA.sig.gz',
    'SampleB': sourmash_dir / 'SampleB.sig.gz',
    'SampleC': sourmash_dir / 'SampleC.sig.gz',
}
