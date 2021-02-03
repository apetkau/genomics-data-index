from tempfile import TemporaryDirectory
from pathlib import Path
import logging

from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio import Phylo, AlignIO

from storage.variant.service import DatabaseConnection
from Bio.Phylo.Applications import FastTreeCommandline

logger = logging.getLogger(__file__)


class TreeService:
    TREE_BUILD_TYPES = ['fasttree']

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
                print(out)

                with open(output_file, 'r') as f:
                    trees = list(Phylo.parse(f, 'newick'))
                    if len(trees) != 1:
                        raise Exception(f'Error, read more than one tree [{len(trees)}] from {output_file}')
                    return trees[0]
            else:
                raise Exception(f'tree_type=[{tree_build_type}] is invalid')