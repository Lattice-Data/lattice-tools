from functools import cache

from db2_flattener.schema.generate import YAML_PATH, load_yaml_config


@cache
def load_config():
    return load_yaml_config(YAML_PATH)
