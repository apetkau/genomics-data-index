import os
from pathlib import Path

assemblies_dir = Path(__file__).parent / 'data' / 'assemblies'
reads_dir = Path(__file__).parent / 'data' / 'reads'
samples_dir = assemblies_dir / 'samples'
expected_mutations_dir = assemblies_dir / 'expected'
assemblies_reference = assemblies_dir / 'genome.fasta'
samples = [os.path.splitext(s)[0] for s in os.listdir(samples_dir)]
assemblies_samples = {s: samples_dir / f'{s}.fasta' for s in samples}
expected_mutations = {s: expected_mutations_dir / f'{s}-mutations.txt' for s in samples}

snpeff_input_sampleA = Path(__file__).parent / '..' / 'data' / 'snpeff-database' / 'input' / 'SampleA.fasta.gz'
snpeff_reference_genome = Path(__file__).parent / '..' / 'data' / 'snpeff-database' / 'NC_011083-5000.gbk.gz'

snpeff_reads_paired = [reads_dir / 'SampleA-snpeff_1.fastq.gz', reads_dir / 'SampleA-snpeff_2.fastq.gz']
snpeff_reads_single = [reads_dir / 'SampleA-single-snpeff.fastq.gz']
