from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property
from typing import ClassVar
from urllib.parse import urljoin

import requests
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.generate import SchemaIDs

from constants import ABSTRACT_MAPPING, EXCLUDED_SCHEMAS


class NodeColor(Enum):
    Donor = "#cfe2f3ff"
    Biosample = "#d9d2e9ff"
    Library = "#fff2ccff"
    SequenceFileSet = "#efd5e7ff"
    SequenceFile = "#f4b3b6ff"
    MatrixFileSet = "#b3c4f4ff"
    RawMatrixFile = "#d9ead3ff"
    ProcessedMatrixFile = "#fce5cdff"
    TabularFile = "#cdf6fcff"
    ExperimentalCondition = "#fcf2cdff"
    GeneticModification = "#d1cdfcff"
    Treatment = "#cde0fcff"

    @classmethod
    def _missing_(cls, name):
        ids = SchemaIDs(name)

        name_check = (
            ids.api_name
            if "File" in ids.api_name
            else ABSTRACT_MAPPING.get(ids.api_name)
        )
        for member in cls:
            if member.name == name_check:
                return member
        return None


@dataclass
class LatticeNode:
    type_and_uuid: str
    excluded: set[str] = field(default_factory=lambda: EXCLUDED_SCHEMAS, repr=False)
    _cache: ClassVar[dict] = {}
    connection: ClassVar[Connection] = Connection("db2_prod")

    def __post_init__(self):
        self.type_and_uuid = self.type_and_uuid.strip("/")
        schema_ids, self.uuid = self.type_and_uuid.split("/")
        self.schema_ids = SchemaIDs(schema_ids)
        self.node_type = self.schema_ids.url_prefix

    def get_object_json(self, refetch: bool = False) -> dict:
        if refetch:
            del self._cache[self.type_and_uuid]

        if self.type_and_uuid in self._cache:
            return self._cache[self.type_and_uuid]

        if self.node_type in self.excluded:
            entry = {"@id": None}
            self._cache[self.type_and_uuid] = entry
            return entry

        print(f"API request for {self.type_and_uuid}")
        response = requests.get(
            urljoin(
                self.connection.server,
                self.type_and_uuid,
            )
        ).json()
        self._cache[self.type_and_uuid] = response
        return response

    @cached_property
    def object_json(self, refetch: bool = False) -> dict:
        return self.get_object_json(refetch)

    @property
    def alias(self) -> str | None:
        aliases = self.object_json.get("aliases")
        return aliases[0] if aliases is not None else None

    @property
    def object_path(self) -> str:
        return f"/{self.type_and_uuid}/"

    @staticmethod
    def get_ids(obj, excluded_schemas):
        """
        Yield all string values starting with '/' anywhere in a nested structure.

        Recurses through dicts (values only) and lists/tuples/sets.
        """
        excluded_values = {"/terms/", "audits"}

        if isinstance(obj, str):
            if obj.startswith("/") and obj.endswith("/") and obj not in excluded_values:
                node = LatticeNode(obj)
                if node.node_type not in excluded_schemas:
                    yield obj
        elif isinstance(obj, dict):
            for value in obj.values():
                yield from LatticeNode.get_ids(value, excluded_schemas)
        elif isinstance(obj, (list, tuple, set)):
            for item in obj:
                yield from LatticeNode.get_ids(item, excluded_schemas)

    @property
    def neighbors(self) -> set[str]:
        neighbor_set = set(self.get_ids(self.object_json, self.excluded))
        neighbor_set.discard(self.object_path)
        return neighbor_set
