from pathlib import Path
import os


assemblies_dir = Path(__file__).parent / 'data' / 'assemblies'
samples_dir = assemblies_dir / 'samples'
expected_mutations_dir = assemblies_dir / 'expected'
assemblies_reference = assemblies_dir / 'genome.fasta'
samples = [os.path.splitext(s)[0] for s in os.listdir(samples_dir)]
assemblies_samples = {s: samples_dir / f'{s}.fasta' for s in samples}
expected_mutations = {s: expected_mutations_dir / f'{s}-mutations.txt' for s in samples}
