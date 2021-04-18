from __future__ import annotations

from storage.variant.service.MutationQueryService import MutationQueryService
from storage.variant.service.QueryService import QueryService
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation


class QueryConnection:

    def __init__(self):
        self._query_services = {}

    def register_feature_service(self, feature_type, feature_service) -> None:
        self._query_services[feature_type] = feature_service

    def get_feature_service(self, feature_type) -> QueryService:
        return self._query_services[feature_type]

    @classmethod
    def create(cls) -> QueryConnection:
        connection = cls()
        connection.register_feature_service(QueryFeatureMutation, MutationQueryService)
        return connection
