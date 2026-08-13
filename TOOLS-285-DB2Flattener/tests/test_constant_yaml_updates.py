import pytest

import generate_constants
from generate_constants import (
    ConstantYAML,
    create_field_types,
    create_object_config,
    hash_constant_dict,
    create_yaml_config,
    load_and_return_constant_dicts,
    load_yaml_config,
)


ENDPOINT = "https://fake.test.instance/"


class FakeConnection:
    """Stands in for DB2lattice.Connection, which needs DB2_* env vars"""
    def __init__(self, mode):
        self.mode = mode
        self.server = ENDPOINT
        self.authid = "key"
        self.authpw = "secret"


def make_profiles(field_type="string", link_to=None, slug="human_donor"):
    """
    Minimal JSONProfile with a single schema and a single field.

    field_type only affects the FIELD_TYPES hash; slug/link_to only affect
    the OBJECT_CONFIG hash, so each can be changed independently.
    """
    return {
        "endpoint": ENDPOINT,
        "profiles": [
            {
                "slug": slug,
                "properties": [
                    {
                        "name": "sex",
                        "type": field_type,
                        "element_type": None,
                        "link_to": link_to,
                    },
                ],
            }
        ],
        "profile_count": 1,
        "hierarchy": {},
    }


def constants_from(profiles):
    field_types = create_field_types(profiles)
    object_config = create_object_config(profiles)
    return ConstantYAML(
        date_parsed="2020-01-01",
        field_types=field_types,
        object_config=object_config,
        field_types_sha256=hash_constant_dict(field_types),
        object_config_sha256=hash_constant_dict(object_config),
    )


@pytest.fixture
def yaml_path(tmp_path, monkeypatch):
    """Point the module at a throwaway yaml and stub out the network"""
    path = tmp_path / "constants.yaml"
    monkeypatch.setattr(generate_constants, "YAML_PATH", path)
    monkeypatch.setattr(generate_constants, "Connection", FakeConnection)
    return path


def seed(yaml_path, profiles):
    create_yaml_config(configs={ENDPOINT: constants_from(profiles)}, yaml_path=yaml_path)


def run_against(monkeypatch, profiles):
    """Run the generator as if the live API returned these profiles"""
    monkeypatch.setattr(
        generate_constants, "create_json_profiles", lambda connection: profiles
    )
    return load_and_return_constant_dicts("db2_fake")


def test_field_types_only_change_is_saved(yaml_path, monkeypatch):
    """A change to FIELD_TYPES alone must still be written to the yaml"""
    seed(yaml_path, make_profiles(field_type="string"))
    before = load_yaml_config(yaml_path)[ENDPOINT]

    run_against(monkeypatch, make_profiles(field_type="array"))
    after = load_yaml_config(yaml_path)[ENDPOINT]

    assert after.field_types["sex"] == {"type": "array"}
    assert after.field_types_sha256 != before.field_types_sha256
    # object config was untouched, so this half of the config should not move
    assert after.object_config_sha256 == before.object_config_sha256


def test_object_config_only_change_is_saved(yaml_path, monkeypatch):
    """A change to OBJECT_CONFIG alone must still be written to the yaml"""
    seed(yaml_path, make_profiles(link_to=None))
    before = load_yaml_config(yaml_path)[ENDPOINT]

    run_against(monkeypatch, make_profiles(link_to="ControlledTerm"))
    after = load_yaml_config(yaml_path)[ENDPOINT]

    assert after.object_config["human_donors"]["references"] == {
        "sex": "controlled_terms"
    }
    assert after.object_config_sha256 != before.object_config_sha256
    assert after.field_types_sha256 == before.field_types_sha256


def test_both_changed_is_saved(yaml_path, monkeypatch):
    seed(yaml_path, make_profiles(field_type="string", slug="human_donor"))
    before = load_yaml_config(yaml_path)[ENDPOINT]

    run_against(monkeypatch, make_profiles(field_type="array", slug="tissue"))
    after = load_yaml_config(yaml_path)[ENDPOINT]

    assert after.field_types_sha256 != before.field_types_sha256
    assert after.object_config_sha256 != before.object_config_sha256
    assert "tissues" in after.object_config


def test_unchanged_config_is_left_alone(yaml_path, monkeypatch):
    """Identical profiles should not rewrite the file, including date_parsed"""
    seed(yaml_path, make_profiles())
    before = yaml_path.read_text()

    run_against(monkeypatch, make_profiles())

    assert yaml_path.read_text() == before


def test_unsaved_endpoint_does_not_raise(yaml_path, monkeypatch, capsys):
    """
    A custom mode pointing at an instance that is neither already stored nor in
    CONFIGS_TO_SAVE has no entry to read back, and must not KeyError.
    """
    # an existing yaml with no entry for our endpoint, and ENDPOINT is not in
    # CONFIGS_TO_SAVE, so nothing gets written and there is nothing to read back
    create_yaml_config(configs={}, yaml_path=yaml_path)
    assert ENDPOINT not in load_yaml_config(yaml_path)

    field_types, object_config = run_against(monkeypatch, make_profiles())

    # the freshly built constants are still returned to the caller
    assert field_types["sex"] == {"type": "string"}
    assert "human_donors" in object_config
    assert "no import hashes to compare" in capsys.readouterr().out


def test_returns_new_constants_not_stored_ones(yaml_path, monkeypatch):
    """The return value reflects the live API, not whatever was on disk"""
    seed(yaml_path, make_profiles(field_type="string"))

    field_types, _ = run_against(monkeypatch, make_profiles(field_type="array"))

    assert field_types["sex"] == {"type": "array"}
