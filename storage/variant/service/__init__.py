from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from storage.variant.model import Base


class DatabaseConnection:

    def __init__(self, connection_string: str):
        engine = create_engine(connection_string, echo=False)

        Session = sessionmaker(bind=engine)
        self._session = Session()

        Base.metadata.create_all(engine)

    def get_session(self):
        return self._session
