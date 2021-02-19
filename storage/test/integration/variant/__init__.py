from os import path, listdir
from pathlib import Path

data_dir = Path(path.dirname(__file__), '..', 'data', 'snippy')
sample_dirs = [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]
reference_file = data_dir / 'genome.fasta.gz'
tree_file = data_dir / 'tree.tre'

data_dir_empty = Path(path.dirname(__file__), '..', 'data', 'empty_vcfs')
