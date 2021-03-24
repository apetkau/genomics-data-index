from os import path, listdir
from pathlib import Path

root_data_dir = Path(path.dirname(__file__), '..', 'data')
data_dir = root_data_dir / 'snippy'
variation_dir = root_data_dir / 'variation'
consensus_dir = root_data_dir / 'consensus'
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
reference_file = data_dir / 'genome.fasta.gz'
tree_file = data_dir / 'tree.tre'

data_dir_empty = root_data_dir / 'empty_vcfs'
