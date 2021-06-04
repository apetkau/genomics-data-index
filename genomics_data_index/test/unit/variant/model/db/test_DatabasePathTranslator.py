from os import mkdir
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from genomics_data_index.storage.model.db.DatabasePathTranslator import DatabasePathTranslator


def test_translate_from_database():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        expected1 = tmp_dir / 'file1'
        expected2 = tmp_dir / 'dir' / 'file2'
        expected3 = tmp_dir / 'dir' / 'subdir' / 'file3.txt'

        # Create files
        mkdir(tmp_dir / 'dir')
        mkdir(tmp_dir / 'dir' / 'subdir')
        open(expected1, 'w').close()
        open(expected2, 'w').close()
        open(expected3, 'w').close()

        dpt = DatabasePathTranslator(tmp_dir)

        actual_file = dpt.from_database('file1')
        assert actual_file == expected1

        actual_file = dpt.from_database('dir/file2')
        assert actual_file == expected2

        actual_file = dpt.from_database('dir/subdir/file3.txt')
        assert actual_file == expected3


def test_translate_from_database_fail():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        dpt = DatabasePathTranslator(tmp_dir)

        with pytest.raises(Exception) as execinfo:
            dpt.from_database('file1')
        assert 'does not exist for relative path' in str(execinfo.value)


def test_translate_to_database():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        file1 = tmp_dir / 'file1'
        file2 = tmp_dir / 'dir' / 'file2'
        file3 = tmp_dir / 'dir' / 'subdir' / 'file3.txt'

        dpt = DatabasePathTranslator(tmp_dir)

        assert 'file1' == dpt.to_database(file1)
        assert 'dir/file2' == dpt.to_database(file2)
        assert 'dir/subdir/file3.txt' == dpt.to_database(file3)
