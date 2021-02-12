import logging
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import FastTreeCommandline
from ete3 import Tree

from storage.variant.service import DatabaseConnection

logger = logging.getLogger(__name__)


class TreeService:
    TREE_BUILD_TYPES = ['fasttree', 'iqtree']

    def __init__(self, database_connection: DatabaseConnection):
        self._database = database_connection

    def build_tree(self, alignment: MultipleSeqAlignment, tree_build_type: str = 'fasttree') -> Tree:
        if not tree_build_type in self.TREE_BUILD_TYPES:
            raise Exception(
                f'tree_type=[{tree_build_type}] is not one of the valid tree builder types {self.TREE_BUILD_TYPES}')

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
                return tree
            elif tree_build_type == 'iqtree':
                output_file = f'{str(input_file)}.treefile'
                command = ['iqtree', '-s', str(input_file)]
                logger.debug(' '.join(command))
                try:
                    subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                except subprocess.CalledProcessError as e:
                    err_msg = str(e.stderr.strip())
                    raise Exception(f'Could not run iqtree on alignment: error {err_msg}')

                tree = Tree(str(output_file))
                return tree
            else:
                raise Exception(f'tree_type=[{tree_build_type}] is invalid')
