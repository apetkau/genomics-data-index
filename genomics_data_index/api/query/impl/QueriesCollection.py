from __future__ import annotations

from typing import List, Union

from genomics_data_index.storage.model.QueryFeature import QueryFeature


class QueriesCollection:

    def __init__(self, queries_list: List[Union[QueryFeature, str]]):
        self._queries = queries_list

    def append(self, query_obj: Union[QueryFeature, str]) -> QueriesCollection:
        new_queries_list = self._queries.copy()
        new_queries_list.append(query_obj)
        return QueriesCollection(new_queries_list)

    @property
    def last(self):
        return self._queries[-1]

    def query_expression(self) -> str:
        query_strs = [str(q) for q in self._queries]
        return ' AND '.join(query_strs)

    @classmethod
    def create_empty(cls):
        return QueriesCollection(queries_list=[])
