from dataclasses import dataclass
from functools import cache, cached_property
from typing import Protocol

from db2_flattener.gather.gatherer import DB2Gatherer
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.constants import Configs
from db2_flattener.schema.generate import YAML_PATH, load_yaml_config


@cache
def load_config():
    return load_yaml_config(YAML_PATH)
