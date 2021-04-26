from pathlib import Path


class DatabasePathTranslator:
    """
    Used to translate between absolute and relative paths when persisting Paths in the database.
    """

    def __init__(self, root_path):
        self._root_path = root_path

    def to_database(self, file: Path) -> str:
        """
        Translates the given file to the relative representation to store in the database.

        :param file: The file to translate.
        :return: A string with the relative representation to store in the database.
        """
        try:
            relative_path = file.relative_to(self._root_path)
            return str(relative_path)
        except ValueError as e:
            raise Exception(f'File to save [{file}] is not relative to root data path [{self._root_path}]', e)

    def from_database(self, file: str) -> Path:
        """
        Translates the given file representation in the database to an absolute Path on the filesystem.

        :param file: The string representation of the file in the database.
        :return: A Path representing this file on the filesystem.
        """
        file_path = Path(file)
        if file_path.is_absolute():
            raise Exception(f'Path from database [{file_path}] is absolute')

        absolute_file_path = self._root_path / file_path

        if not absolute_file_path.exists():
            raise Exception(f'absolute_file_path=[{absolute_file_path}] does not exist for relative path [{file_path}]')

        return absolute_file_path
