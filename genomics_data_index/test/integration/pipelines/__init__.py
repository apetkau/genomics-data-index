from pathlib import Path
import os


assemblies_dir = Path(__file__).parent / 'data' / 'assemblies'
samples_dir = assemblies_dir / 'samples'
expected_mutations_dir = assemblies_dir / 'expected'
assemblies_reference = assemblies_dir / 'genome.fasta'
assemblies_samples = {os.path.splitext(s)[0]: samples_dir / s for s in os.listdir(samples_dir) if s.startswith('Sample')}
expected_mutations = {os.path.splitext(s)[0]: expected_mutations_dir / f'{s}-mutations.txt' for s in os.listdir(samples_dir) if s.startswith('Sample')}
