import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from os import PathLike
from pathlib import Path
from typing import Any, TypeAlias

import yaml
from constants import FIELD_TYPES, OBJECT_CONFIG
from DB2lattice import Connection
from extract_lattice_profiles import (
    LatticeProfileClient,
    ProfileSummary,
    summarize_profile,
)

Hierarchy: TypeAlias = dict[str, dict[str, dict]]
JSONProfile: TypeAlias = dict[str, Any]
FieldType: TypeAlias = dict[str, dict[str, str]]
ObjectConfig: TypeAlias = dict[str, dict[str, Any]]


YAML_PATH = Path("constants.yaml")
CONFIGS_TO_SAVE = {
    "https://api.data.lattice-data.org/",
    "https://lattice-api-dev.demo.lattice-data.org/",
    "https://api.sandbox.lattice-data.org",
}


def snake_to_camel(snake_case_str: str) -> str:
    return "".join(split.capitalize() for split in snake_case_str.split("_"))


def camel_to_snake(text: str) -> str:
    # Insert an underscore before any capital letter preceded by a lowercase letter or digit
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', text)
    # Insert an underscore between a lowercase letter/digit and a capital letter
    s2 = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1)
    return s2.lower()


def make_plural(input_str: str) -> str:
    vowels = "aeiou"
    last_char = input_str[-1]
    penultimate_char = input_str[-2]

    if last_char == "y" and penultimate_char not in vowels:
        return input_str[:-1] + "ies"
    return input_str + "s"


def make_singular(input_str: str) -> str:
    if not input_str.endswith("s"):
        return input_str
    if input_str.endswith("ies"):
        return input_str[:-3] + "y"
    return input_str[:-1]


@dataclass
class SchemaIDs:
    """
    Takes any input in the form of API name/title, URL prefix, or profile slug
    
    Parses correct ids into attributes, input_str is also available for other checks/asserts 
    """
    input_str: str
    
    def __post_init__(self):
        self.slug, self.url_prefix, self.api_name = self._generate_ids(self.input_str)
    
    def _generate_ids(self, input_str: str) -> tuple[str, str, str]:
        if input_str.endswith("s"):
            url_prefix = input_str
            slug = make_singular(input_str)
            api_name = snake_to_camel(slug)
        elif input_str[0].isupper():
            api_name = input_str
            slug = camel_to_snake(input_str)
            url_prefix = make_plural(slug)
        else:
            slug = input_str
            api_name = snake_to_camel(slug)
            url_prefix = make_plural(slug)

        return slug, url_prefix, api_name


def get_concrete_classes(input_name: str, hierarchy_dict: Hierarchy) -> list[str] | str:
    """Return str or list[str] of API names/Titles if concrete classes exist, else return initial input"""
    schema_id = SchemaIDs(input_name)
    if (
        schema_id.api_name in hierarchy_dict 
        and hierarchy_dict[schema_id.api_name]
    ):
        return [SchemaIDs(subclass).url_prefix for subclass in hierarchy_dict[schema_id.api_name]]

    return schema_id.url_prefix


@dataclass
class ConstantYAML:
    date_parsed: str
    field_types: FieldType
    object_config: ObjectConfig
    field_types_sha256: str
    object_config_sha256: str


def create_field_types(profiles: JSONProfile) -> FieldType:
    built_types: FieldType = {}

    for schema in profiles["profiles"]:
        for field in schema["properties"]:
            name = field["name"]
            data_type = field["type"]
            element_type = field["element_type"]
            values = {"type": data_type}

            if element_type is not None:
                values["elements"] = element_type 
            if name in built_types and built_types[name] != values:
                print(f"{schema}.{name} has different values than current parse")
                print(name, data_type, element_type)

            built_types[name] = values

    built_types = {key: value for key, value in sorted(built_types.items())}

    print("Created FIELD_TYPES")
    
    return built_types


def create_object_config(profiles: JSONProfile) -> ObjectConfig:
    object_config: ObjectConfig = {}
    hierarchy = profiles["hierarchy"]

    for schema in profiles["profiles"]:
        references_values: dict[str, str | list[str]] = {}
        slug = schema["slug"]
        schema_ids = SchemaIDs(slug)
        fields = [field["name"] for field in schema["properties"]]

        for field in schema["properties"]:
            field_name = field["name"]
            link_to = field["link_to"]
            if link_to is not None:
                references = get_concrete_classes(link_to, hierarchy)
                references_values[field_name] = references

        config_values = {
            "api_type": schema_ids.api_name,
            "fields": fields,
            "references": references_values
        }
        object_config[schema_ids.url_prefix] = config_values

    print("Created OBJECT_CONFIG")
            
    return object_config


def create_json_profiles(connection: Connection) -> JSONProfile:
    client = LatticeProfileClient(
        endpoint=connection.server,
        access_key=connection.authid,
        secret_key=connection.authpw,
    )

    slugs: list[str] = client.fetch_profile_slug_list()
    summaries: list[ProfileSummary] = []

    print(f"Fetching JSON profiles from: {client.endpoint}")

    for i, slug in enumerate(slugs, start=1):
        print(f"[{i}/{len(slugs)}] Fetching {slug}", file=sys.stderr)
        profile = client.fetch_profile(slug)
        summaries.append(summarize_profile(slug, profile))

    hierarchy = client.fetch_profile_hierarchy()

    return {
        "endpoint": client.endpoint,
        "profiles": [asdict(summary) for summary in summaries],
        "profile_count": len(summaries),
        "hierarchy": hierarchy
    }


def create_yaml_config(configs: dict[str, ConstantYAML], yaml_path: PathLike | str) -> None:
    output = {endpoint: asdict(constants) for endpoint, constants in configs.items()}
    with open(yaml_path, "w") as file:
        yaml.dump(output, file, default_flow_style=False)


def load_yaml_config(yaml_path: PathLike | str) -> dict[str, ConstantYAML]:
    with open(yaml_path, "rt") as f:
        full_yaml_config = yaml.safe_load(f.read())

    return {
        endpoint: ConstantYAML(**values)
        for endpoint, values in full_yaml_config.items()
    }


def hash_constant_dict(constant_dict: dict[str, Any]) -> str:
    """Dictionary cannot be directly hashed, but the JSON version can"""
    constant_in_json = json.dumps(constant_dict, sort_keys=True).encode("utf-8")
    return hashlib.sha256(constant_in_json).hexdigest()


def compare_field_types(old_types: FieldType, new_types: FieldType) -> None:
    old_fields_set = {field for field in old_types}
    new_fields_set = {field for field in new_types}
    new_fields = new_fields_set - old_fields_set
    removed_fields = old_fields_set - new_fields_set
    unchanged_fields = new_fields_set.intersection(old_fields_set)
    print("Comparing old and new fields...")

    print("New fields added:")
    for field in new_fields:
        print(f"\t{field}")

    print("Fields removed:")
    for field in removed_fields:
        print(f"\t{field}")

    for field in unchanged_fields:
        if old_types[field] != new_types[field]:
            print(f"Field type change for {field}:")
            print(f"\tOld type: '{old_types[field]}'")
            print(f"\tNew type: '{new_types[field]}'")


def compare_object_configs(old_configs: ObjectConfig, new_configs: ObjectConfig) -> None:
    old_objects_set = {schema for schema in old_configs}
    new_objects_set = {schema for schema in new_configs}
    new_objects = new_objects_set - old_objects_set
    removed_objects = old_objects_set - new_objects_set
    unchanged_objects = new_objects_set.intersection(old_objects_set)
    print("Comparing old and new schemas/objects...")

    print("New schemas/objects added:")
    for schema in new_objects:
        print(f"\t{schema}")

    print("Schemas/objects removed:")
    for schema in removed_objects:
        print(f"\t{schema}")

    # TODO: this needs a better way to display, or look to something else for the diff
    for schema in unchanged_objects:
        if old_configs[schema]["references"] != new_configs[schema]["references"]:
            print(f"Reference change for {schema}:")


def main() -> None:
    loaded_config: dict[str, ConstantYAML] = load_yaml_config(YAML_PATH)
    update_config = False

    connection = Connection("demo")
    endpoint = connection.server
    profiles: JSONProfile = create_json_profiles(connection)

    field_types = create_field_types(profiles)
    object_config = create_object_config(profiles)
    field_hash = hash_constant_dict(field_types)
    object_hash = hash_constant_dict(object_config)

    new_config = {
        endpoint: ConstantYAML(
            date_parsed=datetime.today().strftime("%Y-%m-%d"),
            field_types=field_types,
            object_config=object_config,
            field_types_sha256=field_hash,
            object_config_sha256=object_hash,
        )
    }

    if endpoint in loaded_config:
        print(f"Found prior constant config for '{endpoint}'")
        field_hash_mismatch = loaded_config[endpoint].field_types_sha256 != new_config[endpoint].field_types_sha256
        object_hash_mismatch = loaded_config[endpoint].object_config_sha256 != new_config[endpoint].object_config_sha256

        if not all([field_hash_mismatch, object_hash_mismatch]):
            print("Hashes are the same for old and new configs")
        else:
            print("Hash mismatch")
            update_config = True

        if field_hash_mismatch:
            print("New and old FIELD_TYPE hashes differ")
            compare_field_types(
                loaded_config[endpoint].field_types,
                new_config[endpoint].field_types,
            )
        if object_hash_mismatch:
            print("New and old OBJECT_CONFIG hashes differ")
            compare_object_configs(
                loaded_config[endpoint].object_config,
                new_config[endpoint].object_config,
            )

    # will only save major instances for now instead of any give demo instance
    if (
        update_config 
        or (endpoint not in loaded_config and endpoint in CONFIGS_TO_SAVE)
    ):
        print(f"Updating {YAML_PATH} with new config")
        loaded_config.update(new_config)
        create_yaml_config(
            configs=loaded_config,
            yaml_path="old_constants.yaml",
        )


    print(f"field start hash: {field_hash}")
    print(f"field import hash: {hash_constant_dict(loaded_config[connection.server].field_types)}")
    print(f"object config start hash: {object_hash}")
    print(f"object config import hash: {hash_constant_dict(loaded_config[connection.server].object_config)}")


def make_current_constant_yaml():
    endpoint = "https://lattice-api-dev.demo.lattice-data.org/"
    field_hash = hash_constant_dict(FIELD_TYPES)
    object_hash = hash_constant_dict(OBJECT_CONFIG)
    create_yaml_config(
        configs={
            endpoint: ConstantYAML(
                field_types_sha256=field_hash,
                object_config_sha256=object_hash,
                date_parsed=datetime.today().strftime("%Y-%m-%d"),
                field_types=FIELD_TYPES,
                object_config=OBJECT_CONFIG,
            )
        },
        yaml_path=Path("old_constants.yaml"),
    )


if __name__ == "__main__":
    main()
