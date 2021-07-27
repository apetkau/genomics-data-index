import logging
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Tuple

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import FastTreeCommandline
from ete3 import Tree

from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.service.ReferenceService import ReferenceService

logger = logging.getLogger(__name__)


class TreeService:
    TREE_BUILD_TYPES = ['fasttree', 'iqtree']
    ALIGN_TYPES = ['core', 'full']

    def __init__(self, database_connection: DatabaseConnection, reference_service: ReferenceService,
                 core_alignment_service: CoreAlignmentService):
        self._database = database_connection
        self._reference_service = reference_service
        self._core_alignment_service = core_alignment_service

    def build_tree(self, alignment: MultipleSeqAlignment,
                   tree_build_type: str = 'fasttree', num_cores: int = 1, align_type: str = 'core',
                   extra_params: str = None) -> Tuple[Tree, str]:
        if not tree_build_type in self.TREE_BUILD_TYPES:
            raise Exception(
                f'tree_type=[{tree_build_type}] is not one of the valid tree builder types {self.TREE_BUILD_TYPES}')

        if num_cores < 1:
            raise Exception(f'num_cores=[{num_cores}] is not supported')

        if extra_params is not None and tree_build_type == 'fasttree':
            raise Exception(f'extra_params is not supported for tree_build_type=[{tree_build_type}]')

        if align_type not in self.ALIGN_TYPES:
            raise Exception(f'align_type=[{align_type}] is not supported, must be one of {self.ALIGN_TYPES}')

        with TemporaryDirectory() as tmp_dir:
            input_file = Path(tmp_dir, 'input.fasta')
            with open(input_file, 'w') as f:
                AlignIO.write(alignment, f, 'fasta')

            if tree_build_type == 'fasttree':
                output_file = Path(tmp_dir, 'fasttree.tre')
                command = FastTreeCommandline(input=str(input_file), out=str(output_file))
                out, err = command()
                logger.debug('Output from FastTree')
                logger.debug(out)

                tree = Tree(str(output_file))
                return tree, out
            elif tree_build_type == 'iqtree':
                output_file = f'{str(input_file)}.treefile'
                command = ['iqtree', '--terrace', '--threads-max', str(num_cores), '-T', 'AUTO',
                           '-s', str(input_file)]

                extra_params_list = None
                if extra_params is not None:
                    extra_params_list = extra_params.split()
                    command.extend(extra_params_list)

                # Add ascertainment bias correction for core alignment
                # But only add if '-m' is not included in extra_params
                if extra_params_list is None or not ('-m' in extra_params_list):
                    if align_type == 'core':
                        command.extend(['-m', 'MFP+ASC'])
                    elif align_type == 'full':
                        command.extend(['-m', 'MFP'])

                try:
                    logger.debug(f'Running: {" ".join(command)}')
                    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                               check=True, text=True)
                    out = str(completed.stdout)
                    tree = Tree(str(output_file))
                    return tree, out
                except subprocess.CalledProcessError as e:
                    err_msg = str(e.stderr.strip())
                    raise Exception(f'Could not run iqtree on alignment: error {err_msg}')
            else:
                raise Exception(f'tree_type=[{tree_build_type}] is invalid')

    def rebuild_tree(self, reference_name: str, num_cores: int = 1, tree_build_type='iqtree',
                     align_type='core', extra_params=None):
        logger.debug('Building alignment')
        alignment = self._core_alignment_service.construct_alignment(reference_name=reference_name,
                                                                     include_reference=True,
                                                                     align_type=align_type)

        logger.debug('Building tree')
        tree, out = self.build_tree(alignment=alignment,
                                    tree_build_type=tree_build_type,
                                    num_cores=num_cores,
                                    align_type=align_type,
                                    extra_params=extra_params)

        logger.debug(f'Updating tree for reference genome [{reference_name}]')
        self._reference_service.update_tree(reference_name=reference_name,
                                            tree=tree,
                                            alignment_length=alignment.get_alignment_length())
