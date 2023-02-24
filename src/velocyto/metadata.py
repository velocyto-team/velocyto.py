import csv
from typing import Any, Iterable


class Metadata:
    def __init__(self, keys: list, values: list, types: Iterable) -> None:
        self.types = dict(zip(keys, types, strict=True))
        self.dict = dict(zip(keys, values, strict=True))
        for ix in range(len(keys)):
            setattr(self, keys[ix], values[ix])


class MetadataCollection:
    def __init__(self, filename: str) -> None:
        self.items = []  # type: list
        self.load(filename)

    def load(self, filename: str) -> None:
        keys = None
        types = None

        with open(filename, newline="") as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read())
            csvfile.seek(0)
            reader = csv.reader(csvfile, dialect)
            for row in reader:
                if len(row) == 0:
                    continue
                if keys is None:
                    if len(row[0].split(":")) == 2:
                        keys = [r.split(":")[0] for r in row]
                        types = [r.split(":")[1] for r in row]  # NOTE: I don't use type anymore
                    else:
                        keys = row
                        types = ["None" for _ in row]
                else:
                    self.items.append(Metadata(keys, row, types))

    def where(self, key: Any, value: Any) -> list:
        return [item for item in self.items if getattr(item, key) == value]
