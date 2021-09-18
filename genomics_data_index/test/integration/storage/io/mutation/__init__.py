from os import path
from pathlib import Path
from typing import Dict


def vcf_and_mask_files(sample_dirs) -> Dict[str, Dict[str, Path]]:
    sample_vcf_map = {}
    sample_genomic_files_mask = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        genomic_file_mask = Path(d, 'snps.aligned.fa')

        sample_vcf_map[sample_name] = vcf_file
        sample_genomic_files_mask[sample_name] = genomic_file_mask

    return {
        'vcfs': sample_vcf_map,
        'masks': sample_genomic_files_mask
    }


def vcf_and_bed_mask_files(sample_dirs) -> Dict[str, Dict[str, Path]]:
    sample_vcf_map = {}
    sample_genomic_files_mask = {}
    sample_genomic_files_mask_minus_vcf = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        genomic_file_mask = Path(d, 'snps.aligned.bed.gz')
        genomic_file_mask_minus_vcf = Path(d, 'snps.aligned.minus-vcf.bed.gz')

        sample_vcf_map[sample_name] = vcf_file
        sample_genomic_files_mask[sample_name] = genomic_file_mask
        sample_genomic_files_mask_minus_vcf[sample_name] = genomic_file_mask_minus_vcf

    return {
        'vcfs': sample_vcf_map,
        'masks': sample_genomic_files_mask,
        'masks_minus_vcf': sample_genomic_files_mask_minus_vcf,
    }
