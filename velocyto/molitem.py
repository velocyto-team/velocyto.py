from typing import *
import velocyto as vcy
import numpy as np
from collections import defaultdict


def dictionary_union(d1: DefaultDict[Any, List], d2: DefaultDict[Any, List]) -> DefaultDict[Any, List]:
    """Set union (|) operation on default dicitonary

    Arguments
    ---------
    d1: defaultdict
        First default dict
    d2: defaultdict
        Second default dict
    Returns
    -------
    A dictionary with the key the set union of the keys.
    If same key is present the entry will be combined using __add__
    """
    keys_set = set(d1) | set(d2)
    return defaultdict(list, {k: d1[k] + d2[k] for k in keys_set})


def dictionary_intersect(d1: DefaultDict[Any, List], d2: DefaultDict[Any, List]) -> DefaultDict[Any, List]:
    """Set intersection (&) operation on default dicitonary

    Arguments
    ---------
    d1: defaultdict
        First default dict
    d2: defaultdict
        Second default dict
    
    Returns
    -------
    A dictionary with the key the set intersection of the keys.
    If same key is present the entry will be combined using __add__
    """
    keys_set = set(d1) & set(d2)
    return defaultdict(list, ((k, d1[k] + d2[k]) for k in keys_set))


class Molitem:
    """Object that represents a molecule in the counting pipeline"""
    __slots__ = ["mappings_record"]  # , "final_report"]

    def __init__(self) -> None:
        self.mappings_record: DefaultDict[vcy.TranscriptModel, List[vcy.SegmentMatch]] = None
        # self.final_report: Tuple[int, int, int, int, int, int] = None

    def add_mappings_record(self, mappings_record: DefaultDict[vcy.TranscriptModel, List[vcy.SegmentMatch]]) -> None:
        if self.mappings_record is None:
            self.mappings_record = mappings_record
        else:
            self.mappings_record = dictionary_intersect(self.mappings_record, mappings_record)
