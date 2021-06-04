from __future__ import annotations

import abc
from typing import Union, List, Set, Tuple

import numpy as np
import pandas as pd
from ete3 import Tree

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class SamplesQuery(abc.ABC):
    """
    The base class for queries related to samples. The different query methods in a SamplesQuery
    operate on a set of samples (specifically a set of Sample identifiers,
    see :py:class:`genomics_data_index.storage.SampleSet`) and returns a new SamplesQuery consisting of the
    subset of samples matching the defined critera.

    That is to say, if A = SamplesQuery(), then B = A.method() will consist of a subset of A matching the criteria
    defined by A.method().

    Multiple subtypes of SamplesQuery exist which may provided access to data joined to a set of Samples for
    implementing different queries (e.g., a :py:class:`genomics_data_index.api.query.impl.TreeSamplesQuery` provides
    a tree joined to a SampleSet and methods to query based on the tree.
    """

    ALL_SAMPLES = SampleSet.create_all()

    def __init__(self):
        pass

    @property
    @abc.abstractmethod
    def universe_set(self) -> SampleSet:
        """
        The set of samples :py:class:`genomics_data_index.storage.SampleSet` defining the universe under consideration
        for this SamplesQuery. Useful for e.g., complement() functionality.

        :returns: A SampleSet.
        """
        pass

    @property
    @abc.abstractmethod
    def sample_set(self) -> SampleSet:
        """
        The set of selected samples :py:class:`genomics_data_index.storage.SampleSet` in this query.

        :returns: A SampleSet defining the selected samples of this query.
        """
        pass

    @abc.abstractmethod
    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None, default_isa_kind: str = 'names',
             default_isa_column: str = None) -> SamplesQuery:
        """
        Joins the passed dataframe onto the current query using the passed column name in the dataframe.
        The column can either contain sample IDs (in sample_ids_column) or sample names (sample_names_column).
        This will modify the universe set to the subset of samples found within the passed data frame.

        :param data_frame: The data frame to join on.
        :param sample_ids_column: The column name in the data frame containing internal sample IDs used for joining.
        :param sample_names_column: The column name in the data frame containing sample names to join on.
                                    Internally these will be mapped to Sample IDs.
        :param default_isa_kind: The default 'kind' parameter to isa(data, kind=kind, ...) queries. This is used
                                 to simplify queries after joining with a dataframe in cases where you want
                                 isa() to select by a particular column in the dataframe.
        :param default_isa_column: If default_isa_kind is set to 'dataframe' sets the default column to use for
                                   isa() queries.
        :return: A SamplesQuery representing the query joined with the data frame.
        """
        pass

    @abc.abstractmethod
    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        """
        Converts the selected set of samples to a DataFrame with one row per sample. By default only samples selected
        by this query will be returned (setting exclude_absent to False will return all samples in the defined universe
        as rows in the DataFrame with a column in the DataFrame used to define if the sample is present or absent).

        :param exclude_absent: Whether or not samples absent in this query (but in the universe of samples) should be
        included.

        :return: A DataFrame of the samples in this query, one row per sample.
        """
        pass

    @abc.abstractmethod
    def summary(self) -> pd.DataFrame:
        """
        Summarizes the selected samples in a DataFrame. This includes the number present and absent in the selection
        as well as percents.

        :return: A DataFrame summarizing the selected samples.
        """
        pass

    @abc.abstractmethod
    def summary_features(self, kind: str = 'mutations', **kwargs) -> pd.DataFrame:
        """
        Summarizes the selected features in a DataFrame. Please specify the kind of features with the kind parameter.

        :param kind: The kind of feature to summarize. By default this is *mutations*.
        :param **kwargs: Additional keyword arguments. Please see the documentation for the underlying implementation.

        :return: A DataFrame summarizing the features within the selected samples.
        """
        pass

    @abc.abstractmethod
    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all',
                      ncores: int = 1) -> Set[str]:
        """
        Returns all features as a set of strings for the selected samples.

        :param
        """
        pass

    @abc.abstractmethod
    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        """
        Intersects this SamplesQuery with the given :py:class:`genomics_data_index.storage.SampleSet`.
        Appends the given *query_message* to the list of query messages.

        :param sample_set: The :py:class:`genomics_data_index.storage.SampleSet` to intersect.
        :param query_message: A message to include as part of the query messages.

        :return: A new SamplesQuery which is the intersection of this query and the passed set.
        """
        pass

    @abc.abstractmethod
    def and_(self, other: SamplesQuery) -> SamplesQuery:
        """
        Performs an AND operation (intersection) on this query with another query.

        That is if A and B are sample queries then C = A.and_(B) means that C consists of samples in the intersection
        of A and B (C consists of samples that are in A AND B).

        :param other: The other SamplesQuery to use.

        :return: A new SamplesQuery which is the intersection of this query and the other query.
        """
        pass

    @abc.abstractmethod
    def or_(self, other: SamplesQuery) -> SamplesQuery:
        """
        Performs an OR operation (union) on this query with another query.

        That is if A and B are sample queries then C = A.or_(B) means that C consists of samples in the union
        of A and B (C consists of samples that are in A OR B).

        :param other: The other SamplesQuery to use.

        :return: A new SamplesQuery which is the union of this query and the other query.
        """
        pass

    @abc.abstractmethod
    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        """
        Builds a tree from the samples contained in this query using the passed *kind* of features.

        :param kind: The *kind* of features to use.
        :param **kwargs: See documentation for the underlying implementation.

        :return: A new SamplesQuery which has a tree attached to it to be used for additional querying.
        """
        pass

    @abc.abstractmethod
    def join_tree(self, tree: Tree, kind='mutation', **kwargs) -> SamplesQuery:
        """
        Joins a tree to this SamplesQuery representing a tree defined by the passed *kind* of features.

        :param tree: The tree to join to this query.
        :param kind: The *kind* of features that were used to build the tree
                      (used to define what type of tree we are joining/the distance units).
        :param **kwargs: See documentation for the underlying implementation.

        :return: A new SamplesQuery which has the passed tree attached to it to be used for additional querying.
        """
        pass

    @abc.abstractmethod
    def reset_universe(self) -> SamplesQuery:
        """
        Resets the *universe* set to be the set of currently selected samples.
        That is, if `A` is a SamplesQuery consisting of some selected samples, then
        `B = A.reset_universe()` implies that `B.universe_set == A`.

        :return: A SamplesQuery with the universe reset to whatever is selected.
        """
        pass

    @abc.abstractmethod
    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        """
        Queries for samples that have a (**hasa**) particular feature/property. That is if `A` is a SamplesQuery
        and `m` is a mutation, then A.hasa(m) selects all those samples of A that have the mutation m.

        :param property: The particular property/feature to query by. This can be either QueryFeature defining
        the particular feature, a string defining the feature, or a pandas Series consisting of boolean values
        which is the result of a DataFrame selection expression.
        :param kind: The kind of *property* that was passed.

        :return: A new SamplesQuery consisting of those samples that have the passed property.
        """
        pass

    def has(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        """
        Queries for samples that have a particular property. Synonym for hasa().
        """
        return self.hasa(property=property, kind=kind)

    @abc.abstractmethod
    def complement(self):
        """
        Returns the complement of the selected set of samples in the SamplesQuery.
        That is if A is a SamplesQuery with U = A.universe_set, then B = A.complement() will result in a
        SamplesQuery, B, consisting of all those samples in U not in A.

        :return: The complement of this SamplesQuery.
        """
        pass

    def within(self, data: Union[str, List[str], SamplesQuery, SampleSet], **kwargs) -> SamplesQuery:
        """
        Queries for samples within a particular distance. This is identical to calling
        isin(data, kind='distance', ...). Used to help make code easier to read.

        :param data: The data to use for selecting samples by (e.g., the samples being used to measure distance against).
        :param **kwargs: Other arguments for the the internal implementations.
        :return: A SamplesQuery with the matched samples.
        """
        return self.isin(data=data, kind='distance', **kwargs)

    @abc.abstractmethod
    def _within_distance(self, data: Union[str, List[str], SamplesQuery, SampleSet], distance: float, units: str,
                         **kwargs) -> SamplesQuery:
        pass

    @abc.abstractmethod
    def _distance_units(self) -> List[str]:
        pass

    @abc.abstractmethod
    def isa(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str = None,
            **kwargs) -> SamplesQuery:
        """
        Queries for samples which are a particular type/belong to a particular category.
        Read as "subset samples which are a (isa) particular type defined by 'data'".
        The default implementation will select samples by sample name but this type of query is also useful when used
        with an attached dataframe (at which point it selects based on matches to a column,
        see documentation for :py:class:`genomics_data_index.api.query.impl.DataFrameSamplesQuery`).

        :param data: The data to match.
        :param kind: The particular kind of data passed.
        :**kwargs: Arguments for the underlying implementation.
        :return: A SamplesQuery with the matched samples.
        """
        pass

    def isan(self, data: Union[str, List[str], SamplesQuery, SampleSet], kind: str = None,
             **kwargs) -> SamplesQuery:
        """
        Synonym for isa()
        """
        return self.isa(data=data, kind=kind, **kwargs)

    @abc.abstractmethod
    def isin(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str = None,
             **kwargs) -> SamplesQuery:
        """
        Queries for samples which are in (isin) the passed data. This can have a number of interpretations
        depending on the passed kind
        1. If kind == 'sample' (or 'samples') then isin(['Name']) returns a set of samples with one of the passed names.
        2. If kind == 'distance' (or 'distances') then isin(['Name']) returns a set of samples that are within
            a particular distance of the passed samples.

        :param data: The data to match. This could be a list of strings (representing sample names) or a SampleSet or
                    SamplesQuery representing a set of samples.
        :param kind: The particular kind of data passed.
        :**kwargs: Arguments for the underlying implementation.
        :return: A SamplesQuery with the matched samples.
        """
        pass

    @abc.abstractmethod
    def _get_sample_names_query_infix_from_data(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet]
                                                ) -> Tuple[Set[str], str]:
        pass

    @abc.abstractmethod
    def _can_handle_isin_kind(self, kind: str) -> bool:
        pass

    @abc.abstractmethod
    def _can_handle_isa_kind(self, kind: str) -> bool:
        pass

    @abc.abstractmethod
    def _isa_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def _isin_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def to_distances(self, kind: str, **kwargs) -> Tuple[np.ndarray, List[str]]:
        """
        Returns a distance-matrix representing the pairwise distances between all selected samples.

        :param kind: The features kind used to generate the distance matrix (e.g,. 'mutations', or 'kmer').
        :param **kwargs: Additional arguments depending on the implementation.

        :return: A Tuple consisting of a pairwise distance matrix (as a NumPy array) and a list of labels of the distance
                 matrix (sample names).
        """
        pass

    @abc.abstractmethod
    def _get_has_kinds(self) -> List[str]:
        pass

    @abc.abstractmethod
    def is_empty(self):
        """
        Whether or not the selected set of samples is empty.
        :return: True if the selected set of samples is empty, False otherwise.
        """
        pass

    @abc.abstractmethod
    def query_expression(self) -> str:
        """
        A string representing the series of queries performed to select the given samples in this SamplesQuery.
        For example if A is a SamplesQuery then `A.query_expression() == "hasa('mutation') AND hasa('mutation2')"
        would mean that the SamplesQuery A represents samples that have both 'mutation' and 'mutation2'.

        :return: A string representing the set of queries performed to generate this set.
        """
        pass

    @property
    @abc.abstractmethod
    def tree(self):
        """
        If this SamplesQuery has a joined tree, then returns the tree.

        :return: An ete3.Tree that has been joined to this SamplesQuery (or raises an exception if no tree).
        """
        pass

    @abc.abstractmethod
    def has_tree(self) -> bool:
        """
        Whether or not this SamplesQuery has a joined tree.

        :return: True if this SamplesQuery has a joined tree, False otherwise.
        """
        pass

    @abc.abstractmethod
    def tolist(self, names=True) -> List[Union[int, str]]:
        """
        Converts the set of selected samples into a list of either sample names or sample IDs.

        :param names: If True (default) return a list of sample names as strings, if False return a list of sample IDs.
        :return: A list of sample names or IDs.
        """
        pass

    @abc.abstractmethod
    def __and__(self, other):
        """
        Performs an **and** operation (intersection) between two different SamplesQuery objects.
        If A and B are two SamplesQuery objects then `A and B` is equivalent to `A.and_(B)`.

        :return: The and (intersection) of two SamplesQuery objects.
        """
        pass

    @abc.abstractmethod
    def __or__(self, other):
        """
        Performs an **or** operation (union) between two different SamplesQuery objects.
        If A and B are two SamplesQuery objects then `A or B` is equivalent to `A.or_(B)`.

        :return: The or (union) of two SamplesQuery objects.
        """
        pass

    @abc.abstractmethod
    def __len__(self):
        """
        The number of selected samples.

        :return: The number of selected samples.
        """
        pass
