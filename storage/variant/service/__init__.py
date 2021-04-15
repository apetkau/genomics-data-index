import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from storage.variant.model.db import Base


class DatabaseConnection:

    def __init__(self, connection_string: str):
        engine = create_engine(connection_string, echo=False)

        Session = sessionmaker(bind=engine)
        self._session = Session()

        Base.metadata.create_all(engine)

    def get_session(self):
        return self._session

    def get_database_size(self) -> pd.DataFrame:
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
