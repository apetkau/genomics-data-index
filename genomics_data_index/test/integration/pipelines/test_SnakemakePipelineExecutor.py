from typing import List

import pandas as pd
import vcf
from pathlib import Path
from tempfile import TemporaryDirectory

from genomics_data_index.test.integration.pipelines import assemblies_samples, assemblies_reference, expected_mutations
from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor


def vcf_to_mutations_list(vcf_file: Path) -> List[str]:
    reader = vcf.Reader(filename=str(vcf_file))
    df = pd.DataFrame([vars(r) for r in reader])
    df['ALT'] = df['ALT'].apply(lambda x: x[0])
    mutations = df.apply(lambda x: f'{x["CHROM"]}:{x["POS"]}:{x["REF"]}:{x["ALT"]}', axis='columns')

    return mutations.tolist()


def read_expected_mutations(file: Path) -> List[str]:
    with open(file, 'r') as fh:
        mutations = [l.strip() for l in fh.readlines()]
        return mutations


def test_create_fofn_file_single_sample():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        actual_mutations_file = tmp_dir / 'variant' / 'SampleA.vcf.gz'
        actual_consensus_file = tmp_dir / 'consensus' / 'SampleA.fasta.gz'
        input_samples = [assemblies_samples['SampleA']]

        sampleA_expected_mutations_file = expected_mutations['SampleA']
        sampleA_expected_mutations = read_expected_mutations(sampleA_expected_mutations_file)

        pipeline_executor = SnakemakePipelineExecutor(working_directory=tmp_dir, use_conda=False)

        input_fofn = pipeline_executor.execute(input_files=input_samples,
                                               reference_file=assemblies_reference,
                                               ncores=1)

        assert input_fofn.exists()
        assert actual_mutations_file.exists()
        assert actual_consensus_file.exists()

        sampleA_actual_mutations = vcf_to_mutations_list(actual_mutations_file)

        assert len(sampleA_actual_mutations) == len(sampleA_expected_mutations)
        assert sampleA_actual_mutations == sampleA_expected_mutations
