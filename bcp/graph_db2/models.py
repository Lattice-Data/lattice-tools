from dataclasses import dataclass, field
from enum import Enum
from typing import ClassVar
from urllib.parse import urljoin

import requests
from connection import get_connection
from constants import ABSTRACT_MAPPING, DEFAULT_MODE, EXCLUDED_SCHEMAS
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.generate import SchemaIDs


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
    mode: str = DEFAULT_MODE
    excluded: set[str] = field(default_factory=lambda: EXCLUDED_SCHEMAS, repr=False)
    _cache: ClassVar[dict] = {}

    def __post_init__(self):
        self.type_and_uuid = self.type_and_uuid.strip("/")
        schema_ids, self.uuid = self.type_and_uuid.split("/")
        self.schema_ids = SchemaIDs(schema_ids)
        self.node_type = self.schema_ids.url_prefix
        self.uuid_path = f"/{self.type_and_uuid}/"

    @property
    def connection(self) -> Connection:
        return get_connection(self.mode)

    def fetch_object_json(self) -> dict:
        print(f"API request for {self.uuid_path}")
        return requests.get(
            urljoin(
                self.connection.server,
                self.uuid_path,
            )
        ).json()

    @property
    def object_json(self) -> dict:
        if self.uuid_path in self._cache:
            return self._cache[self.uuid_path]
        response = self.fetch_object_json()
        self.object_json = response
        return response

    @object_json.setter
    def object_json(self, value: dict) -> None:
        self._cache[self.uuid_path] = value

    @staticmethod
    def batch_cache_update(result_response: list[dict]) -> None:
        """
        Add batched report results to the cache to prevent individual API calls
        Expects list of dicts from DB2Gatherer.chunk_and_fetch()
        """
        for index, json_object in enumerate(result_response):
            uuid_path = json_object.get("@id")
            if uuid_path is None:
                print(f"Could not find '@id' at index: {index} of batch response")
                continue
            LatticeNode._cache[uuid_path] = json_object

    @property
    def alias(self) -> str | None:
        aliases = self.object_json.get("aliases")
        return aliases[0] if aliases is not None else None

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
        neighbor_set.discard(self.uuid_path)
        return neighbor_set
