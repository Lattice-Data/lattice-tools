from functools import cache

from constants import DEFAULT_MODE
from db2_flattener.gather.lattice import Connection
from db2_flattener.schema.constants import Configs
from db2_flattener.schema.generate import (
    YAML_PATH,
    load_and_return_constant_dicts,
    load_yaml_config,
)


@cache
def create_configs(mode: str = DEFAULT_MODE, fetch_new: bool = False) -> Configs:
    if fetch_new:
        field_types, object_config = load_and_return_constant_dicts(mode)
        return Configs(
            FIELD_TYPES=field_types,
            OBJECT_CONFIG=object_config,
        )

    endpoint = Connection(mode).server
    config = load_yaml_config(YAML_PATH)
    if endpoint not in config:
        raise KeyError(
            f"No constants.yaml entry for {endpoint}. "
            f"Available: {sorted(config.keys())}"
            "Or --fetch_new to True"
        )

    return Configs(
        FIELD_TYPES=config[endpoint].field_types,
        OBJECT_CONFIG=config[endpoint].object_config,
    )
