"""Tests for graph_db2.cli argument parsing and seed defaulting."""

from __future__ import annotations

import pytest

from graph_db2 import cli
from graph_db2.constants import DEFAULT_MODE, DEFAULT_PORT, FETCH_NEW
from graph_db2.explorer import SAMPLE_SEED

from tests.graph_db2_helpers import clean_graph_state  # noqa: F401  (autouse)

DEMO = "db2_demo"


def parse(*argv: str):
    return cli.build_parser().parse_args(list(argv))


def resolved_seed(*argv: str) -> str:
    """Mirror the seed defaulting in cli.main()."""
    args = parse(*argv)
    return args.seed or (SAMPLE_SEED if args.mode == DEFAULT_MODE else "")


# --------------------------------------------------------------------------
# defaults
# --------------------------------------------------------------------------


def test_defaults() -> None:
    args = parse()
    assert args.mode == DEFAULT_MODE
    assert args.port == DEFAULT_PORT
    assert args.fetch_new is FETCH_NEW
    assert args.debug is False
    assert args.seed is None


def test_mode_and_port_are_overridable() -> None:
    args = parse("--mode", DEMO, "--port", "9001")
    assert args.mode == DEMO
    assert args.port == 9001


def test_port_is_an_int() -> None:
    assert parse("--port", "8080").port == 8080


def test_bad_port_is_rejected() -> None:
    with pytest.raises(SystemExit):
        parse("--port", "not-a-number")


def test_fetch_new_and_debug_are_flags() -> None:
    args = parse("--fetch-new", "--debug")
    assert args.fetch_new is True
    assert args.debug is True


def test_unknown_flag_is_rejected() -> None:
    with pytest.raises(SystemExit):
        parse("--nope")


# --------------------------------------------------------------------------
# seed defaulting - a prod path 404s on any other deployment
# --------------------------------------------------------------------------


def test_default_mode_gets_the_sample_seed() -> None:
    assert resolved_seed() == SAMPLE_SEED


def test_other_mode_without_seed_starts_empty() -> None:
    """SAMPLE_SEED only exists on DEFAULT_MODE's server. Carrying it to demo
    used to 404 into an empty canvas with only a small grey status line."""
    assert resolved_seed("--mode", DEMO) == ""


def test_explicit_seed_wins_on_any_mode() -> None:
    seed = "/matrix_file_sets/f4aed931-fdcd-456a-a387-02cfad8473c7/"
    assert resolved_seed("--mode", DEMO, "--seed", seed) == seed
    assert resolved_seed("--seed", seed) == seed


@pytest.mark.parametrize(
    "raw",
    [
        "/matrix_file_sets/abc/",
        "matrix_file_sets/abc",
        "matrix_file_sets/abc/",
        "/matrix_file_sets/abc",
    ],
)
def test_seed_is_taken_verbatim_and_normalized_downstream(raw: str) -> None:
    """The parser does not normalize; build_app does, so every spelling has to
    survive argument parsing untouched."""
    from graph_db2.cyto_elements import normalize_path

    assert parse("--seed", raw).seed == raw
    assert normalize_path(raw) == "/matrix_file_sets/abc/"


# --------------------------------------------------------------------------
# main() wiring
# --------------------------------------------------------------------------


def test_main_passes_resolved_arguments_to_build_app(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict = {}

    class FakeApp:
        def run(self, **kwargs) -> None:
            captured["run"] = kwargs

    def fake_build_app(seed: str, mode: str, fetch_new: bool):
        captured["build"] = (seed, mode, fetch_new)
        return FakeApp()

    monkeypatch.setattr(cli, "build_app", fake_build_app)
    cli.main(["--mode", DEMO, "--port", "9999", "--fetch-new"])

    assert captured["build"] == ("", DEMO, True)
    assert captured["run"] == {"debug": False, "port": 9999}


def test_main_uses_sample_seed_on_default_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict = {}

    class FakeApp:
        def run(self, **kwargs) -> None:
            pass

    def fake_build_app(seed: str, mode: str, fetch_new: bool):
        captured["build"] = (seed, mode, fetch_new)
        return FakeApp()

    monkeypatch.setattr(cli, "build_app", fake_build_app)
    cli.main([])

    assert captured["build"] == (SAMPLE_SEED, DEFAULT_MODE, FETCH_NEW)


def test_help_names_the_right_module() -> None:
    """cli.py's docstring is what --help prints; it said `python -m db2_graph`
    while the package is graph_db2."""
    help_text = cli.build_parser().format_help()
    assert "python -m graph_db2" in help_text
    assert "db2_graph" not in help_text


def test_prog_name_matches_the_package() -> None:
    assert cli.build_parser().prog == "graph_db2"
