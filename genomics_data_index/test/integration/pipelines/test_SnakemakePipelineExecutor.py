from typing import List

import pandas as pd
import vcf
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from genomics_data_index.test.integration.pipelines import assemblies_samples, assemblies_reference, expected_mutations
from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor
from genomics_data_index.storage.util import parse_sequence_file


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


def get_consensus_sequences(file: Path) -> List[SeqRecord]:
    ref_name, sequences = parse_sequence_file(file)
    return sequences


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

        # Verify input file of file names for rest of gdi software (used as input to the indexing component)
        fofn_df = pd.read_csv(input_fofn, sep='\t')
        print(fofn_df)
        assert ['Sample', 'VCF', 'Mask File'] == fofn_df.columns.tolist()

        assert 1 == len(fofn_df)
        assert ['SampleA'] == fofn_df['Sample'].tolist()
        actual_mutations_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['VCF'].tolist()[0])
        actual_consensus_file_from_df = Path(fofn_df[fofn_df['Sample'] == 'SampleA']['Mask File'].tolist()[0])

        assert actual_mutations_file == actual_mutations_file_from_df
        assert actual_consensus_file == actual_consensus_file_from_df

        # Verify mutations
        sampleA_actual_mutations = vcf_to_mutations_list(actual_mutations_file)
        assert len(sampleA_actual_mutations) == len(sampleA_expected_mutations)
        assert sampleA_actual_mutations == sampleA_expected_mutations

        # Verify consensus/mask sequence (position of Ns and -s are positions that should be masked out
        sampleA_consensus_records = get_consensus_sequences(actual_consensus_file)
        assert 1 == len(sampleA_consensus_records)

        sampleA_consensus_record = sampleA_consensus_records[0]
        assert 5180 == len(sampleA_consensus_record)
        assert 'reference' == sampleA_consensus_record.id
        assert 0 == sampleA_consensus_record.upper().seq.count('N')
        assert 0 == sampleA_consensus_record.upper().seq.count('-')
