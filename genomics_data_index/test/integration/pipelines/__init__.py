from pathlib import Path
import os


assemblies_dir = Path(__file__).parent / 'data' / 'assemblies'
assemblies_reference = assemblies_dir / 'reference.fasta'
assemblies_samples = [assemblies_dir / s for s in os.listdir(assemblies_dir) if s.startswith('Sample')]
