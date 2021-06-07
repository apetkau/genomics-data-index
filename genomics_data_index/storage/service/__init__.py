import logging
from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import genomics_data_index.storage.model.db
from genomics_data_index.storage.model.db import Base
from genomics_data_index.storage.model.db.DatabasePathTranslator import DatabasePathTranslator

logger = logging.getLogger(__name__)


class DatabaseConnection:

    def __init__(self, connection_string: str, database_path_translator: DatabasePathTranslator):
        self._setup_database_type(connection_string)

        engine = create_engine(connection_string, echo=False)

        Session = sessionmaker(bind=engine)
        self._session = Session()
        self._database_path_translator = database_path_translator

        # Sets global variable here for translating between relative/absolute paths in the database
        # I don't like using global variables and have to look into some other method to set this later
        if genomics_data_index.storage.model.db.database_path_translator is not None:
            logger.warning(f'Attempting to set global database_path_translator={database_path_translator}'
                           ' but it is already set')
        genomics_data_index.storage.model.db.database_path_translator = database_path_translator

        Base.metadata.create_all(engine)

    def _setup_database_type(self, connection_string: str) -> None:
        logger.debug(f'Setup database type: {connection_string}')
        self._is_sqlite = connection_string.startswith('sqlite')

        if self._is_sqlite:
            self._sqlite_path = Path(connection_string[len('sqlite:///'):])

    def get_session(self):
        return self._session

    def get_database_size(self) -> pd.DataFrame:
        if self._is_sqlite:
            db_size = self._sqlite_path.stat().st_size
            size_df = pd.DataFrame([['Database', str(self._sqlite_path.name), 'All', db_size, pd.NA, 1]],
                                   columns=['Type', 'Name', 'Division', 'Data Size', 'Index Size', 'Number of Items'])
            return size_df
        else:
            database_name = self._session.bind.url.database

            table_names = self._session.bind.table_names()

            # TODO: This could maybe be improved to use something more appropriate form SQLAlchemy
            statement = f'''
            SELECT table_name AS `Table`,
            data_length as `Data size`,
            index_length as `Index size`
            FROM information_schema.TABLES
            WHERE table_schema = "{database_name}";
            '''
            result = self._session.execute(statement)
            size_of_table = result.fetchall()

            size_df = pd.DataFrame(size_of_table, columns=['Division', 'Data Size', 'Index Size'])
            size_df.insert(0, 'Name', database_name)
            size_df.insert(0, 'Type', 'Database')

            # Get number of rows from tables in database
            size_df['Number of Items'] = pd.NA
            for name in table_names:
                count_statement = f'''
                SELECT count(*) from {name};
                '''
                result = self._session.execute(count_statement)
                table_counts = result.fetchall()
                assert 1 == len(table_counts)

                size_df.loc[size_df['Division'] == name, 'Number of Items'] = table_counts[0][0]

            return size_df


class EntityExistsError(Exception):

    def __init__(self, msg):
        super().__init__(msg)
