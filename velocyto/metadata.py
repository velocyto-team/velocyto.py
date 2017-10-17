import csv
import logging
from typing import *


class Metadata:
	def __init__(self, keys: List, values: List, types: Iterable) -> None:
		self.types = dict(zip(keys, types))
		self.dict = dict(zip(keys, values))
		for ix in range(len(keys)):
			setattr(self, keys[ix], values[ix])


class MetadataCollection:
	def __init__(self, filename: str) -> None:
		self.items = []  # type: List
		self.load(filename)

	def load(self, filename: str) -> None:
		keys = None
		types = None

		with open(filename, "r", encoding='MacRoman') as f:
			reader = csv.reader(f, delimiter=';')
			for row in reader:
				if len(row) == 0:
					continue
				if keys is None:
					keys = [r.split(':')[0] for r in row]
					types = [r.split(':')[1] for r in row]  # NOTE: I don't use type anymore
				else:
					self.items.append(Metadata(keys, row, types))

	def where(self, key: Any, value: Any) -> List:
		result = []  # type: List
		for item in self.items:
			if getattr(item, key) == value:
				result.append(item)
		return result
