"""Tests for the guidesig --guidescan exact-match report.

Neither guidescan nor S3 is used: the pipeline run is replaced by a stub that
reads the input file guidesig wrote, and S3 by a fake client that serves bytes.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from guidesig import guidescancheck
from guidesig.cli import main
from guidesig.guidescancheck import (
    REFERENCE_FASTA,
    REFERENCE_SUFFIXES,
    GuidescanCheckError,
    GuidescanReport,
    TargetingGuide,
    check_species,
    ensure_references,
    format_report,
    parse_pam,
    pipeline_pam,
    read_targeting_guides,
    run_guidescan_check,
)
from tests.guidesig_helpers import write_delimited

HEADER = [
    "property",
    "guide_id",
    "guide_protospacer",
    "guide_role",
    "guide_target_gene_id",
    "guide_PAM",
]
MOUSE_GENE = "ENSMUSG00000003923"
HUMAN_GENE = "ENSG00000141510"
TARGETING = "Targeting a Gene"


def _sheet(tmp_path: Path, rows: list[list[str]], header=HEADER) -> Path:
    return write_delimited(tmp_path / "sheet.tsv", [header, *rows], file_format="tsv")


def _row(gid: str, seq: str, gene: str = MOUSE_GENE, role: str = TARGETING):
    return ["", gid, seq, role, gene, "3' NGG"]


# --- PAM parsing ---------------------------------------------------------------


@pytest.mark.parametrize("raw", ["3' NGG", "NGG", "3’ NGG", " 3'ngg ", "3'NGG"])
def test_parse_pam_accepts_3prime_ngg(raw: str) -> None:
    assert parse_pam(raw) == "NGG"


@pytest.mark.parametrize(
    "raw, message",
    [
        ("5' TTTV", "5' PAM"),
        ("", "unrecognized PAM"),
        ("3' N-G", "unrecognized PAM"),
    ],
)
def test_parse_pam_rejects(raw: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_pam(raw)


@pytest.mark.parametrize("pam", ["NGG", "AGG", "CGG", "GGG", "TGG"])
def test_pipeline_pam_maps_ngg_family_to_ngg(pam: str) -> None:
    assert pipeline_pam(pam) == "NGG"


@pytest.mark.parametrize("pam", ["TAG", "NAG", "GGA", "NNG", "NNGRRT"])
def test_pipeline_pam_passes_other_patterns_through(pam: str) -> None:
    assert pipeline_pam(pam) == pam


def test_parse_pam_accepts_any_length() -> None:
    assert parse_pam("3' NNGRRT") == "NNGRRT"


@pytest.mark.parametrize("pam", ["3' NAG", "TAG", "3' NNGRRT"])
def test_sheet_with_one_other_pattern_runs_unchanged(tmp_path: Path, pam: str) -> None:
    rows = [
        ["", "g1", "AAAA", TARGETING, MOUSE_GENE, pam],
        ["", "g2", "TTTA", TARGETING, MOUSE_GENE, pam],
    ]
    guides, _targets = read_targeting_guides(_sheet(tmp_path, rows), file_format="tsv")
    assert {guide.pam for guide in guides} == {parse_pam(pam)}


def test_sheet_with_actual_pam_bases_runs_as_ngg(tmp_path: Path) -> None:
    """Per-guide bases and the pattern can be mixed; all are sent as NGG."""
    rows = [
        ["", "g1", "AAAA", TARGETING, MOUSE_GENE, "AGG"],
        ["", "g2", "TTTA", TARGETING, MOUSE_GENE, "3' TGG"],
        ["", "g3", "ATAT", TARGETING, MOUSE_GENE, "cgg"],
        ["", "g4", "TATA", TARGETING, MOUSE_GENE, "3' NGG"],
    ]
    guides, _targets = read_targeting_guides(_sheet(tmp_path, rows), file_format="tsv")
    assert {guide.pam for guide in guides} == {"NGG"}
    seen: dict = {}
    run_guidescan_check(guides, Path("idx"), runner=_stub_runner(set(), seen))
    assert {row["pam"] for row in seen["rows"]} == {"NGG"}


def test_sheet_mixing_pam_patterns_raises(tmp_path: Path) -> None:
    rows = [
        ["", "g1", "AAAA", TARGETING, MOUSE_GENE, "AGG"],
        ["", "g2", "TTTA", TARGETING, MOUSE_GENE, "3' NAG"],
    ]
    with pytest.raises(GuidescanCheckError) as exc:
        read_targeting_guides(_sheet(tmp_path, rows), file_format="tsv")
    assert "targeting guides use more than one PAM pattern (NAG, NGG)" in str(exc.value)


# --- reading targeting guides --------------------------------------------------


def test_read_targeting_guides_skips_controls_comments_and_empty(
    tmp_path: Path,
) -> None:
    path = _sheet(
        tmp_path,
        [
            ["#example", "gx", "TTTT", TARGETING, MOUSE_GENE, "3' NGG"],
            _row("g1", "aaaa"),
            _row("ctl", "CCCA", gene="NA", role="control"),
            _row("g2", "", gene=MOUSE_GENE),
            _row("r1", "TTTA", gene="NA", role="Targeting a Region"),
        ],
    )
    guides, targets = read_targeting_guides(path, file_format="tsv")
    assert guides == [
        TargetingGuide("g1", "AAAA", "NGG"),
        TargetingGuide("r1", "TTTA", "NGG"),
    ]
    assert targets == {MOUSE_GENE, "NA"}


def test_read_targeting_guides_requires_pam_column(tmp_path: Path) -> None:
    path = _sheet(tmp_path, [["", "g1", "AAAA", TARGETING, MOUSE_GENE]], HEADER[:-1])
    with pytest.raises(GuidescanCheckError, match="'guide_PAM' absent from header"):
        read_targeting_guides(path, file_format="tsv")


def test_read_targeting_guides_bad_pam_names_row(tmp_path: Path) -> None:
    path = _sheet(tmp_path, [["", "g1", "AAAA", TARGETING, MOUSE_GENE, "5' TTTV"]])
    with pytest.raises(GuidescanCheckError, match="row 2: 5' PAM"):
        read_targeting_guides(path, file_format="tsv")


def test_read_targeting_guides_requires_guide_id(tmp_path: Path) -> None:
    path = _sheet(tmp_path, [_row("", "AAAA")])
    with pytest.raises(GuidescanCheckError, match="row 2: targeting guide has no"):
        read_targeting_guides(path, file_format="tsv")


def test_repeated_guide_id_with_same_sequence_collapses(tmp_path: Path) -> None:
    path = _sheet(tmp_path, [_row("g1", "AAAA"), _row("g1", "AAAA")])
    guides, _targets = read_targeting_guides(path, file_format="tsv")
    assert guides == [TargetingGuide("g1", "AAAA", "NGG")]


def test_repeated_guide_id_with_different_sequence_raises(tmp_path: Path) -> None:
    path = _sheet(tmp_path, [_row("g1", "AAAA"), _row("g1", "TTTA")])
    with pytest.raises(GuidescanCheckError, match="row 3: guide_id 'g1' appears"):
        read_targeting_guides(path, file_format="tsv")


# --- species cross-check ---------------------------------------------------------


@pytest.mark.parametrize(
    "ids, species",
    [
        ({MOUSE_GENE, "NA"}, "mouse"),
        ({HUMAN_GENE}, "human"),
        ({f"{HUMAN_GENE}.18"}, "human"),
    ],
)
def test_species_matches(ids: set[str], species: str) -> None:
    assert check_species(ids, species) is None


def test_species_mismatch_raises() -> None:
    with pytest.raises(GuidescanCheckError, match="--species human but .* mouse"):
        check_species({MOUSE_GENE}, "human")


def test_mixed_species_raises() -> None:
    with pytest.raises(GuidescanCheckError, match="mixes human"):
        check_species({MOUSE_GENE, HUMAN_GENE}, "mouse")


def test_placeholder_ids_skip_cross_check() -> None:
    """ENSGFAKE... is not a human ID; with nothing recognizable, a note results."""
    note = check_species({"ENSGFAKE0000001", "NA"}, "human")
    assert note is not None
    assert "cross-check skipped" in note


# --- reference cache -----------------------------------------------------------


class FakeS3:
    """Serves fixed-size objects; records every download."""

    def __init__(self, sizes: dict[str, int], fail_download: bool = False) -> None:
        self.sizes = sizes
        self.fail_download = fail_download
        self.downloads: list[str] = []

    def head_object(self, Bucket: str, Key: str) -> dict:
        if Key not in self.sizes:
            raise RuntimeError("Not Found")
        return {"ContentLength": self.sizes[Key]}

    def download_file(self, bucket: str, key: str, filename: str) -> None:
        self.downloads.append(key)
        Path(filename).write_bytes(b"x" * (self.sizes[key] // 2))
        if self.fail_download:
            raise RuntimeError("connection reset")
        Path(filename).write_bytes(b"x" * self.sizes[key])


def _reference_sizes(species: str) -> dict[str, int]:
    fasta = REFERENCE_FASTA[species]
    return {
        f"{guidescancheck.REFERENCE_PREFIX}{fasta}{suffix}": 10 + i
        for i, suffix in enumerate(REFERENCE_SUFFIXES)
    }


def test_ensure_references_downloads_then_reuses_cache(tmp_path: Path) -> None:
    client = FakeS3(_reference_sizes("mouse"))
    prefix = ensure_references("mouse", cache_dir=tmp_path, s3_client=client)
    fasta = REFERENCE_FASTA["mouse"]
    assert prefix == tmp_path / "mouse" / f"{fasta}.index"
    assert len(client.downloads) == len(REFERENCE_SUFFIXES)
    assert not list((tmp_path / "mouse").glob("*.partial"))

    client.downloads.clear()
    ensure_references("mouse", cache_dir=tmp_path, s3_client=client)
    assert client.downloads == []


def test_ensure_references_redownloads_on_size_mismatch(tmp_path: Path) -> None:
    sizes = _reference_sizes("human")
    client = FakeS3(sizes)
    ensure_references("human", cache_dir=tmp_path, s3_client=client)
    stale = tmp_path / "human" / (REFERENCE_FASTA["human"] + ".index.gs")
    stale.write_bytes(b"short")
    client.downloads.clear()
    ensure_references("human", cache_dir=tmp_path, s3_client=client)
    assert client.downloads == [
        f"{guidescancheck.REFERENCE_PREFIX}{REFERENCE_FASTA['human']}.index.gs"
    ]


def test_interrupted_download_leaves_no_cached_file(tmp_path: Path) -> None:
    client = FakeS3(_reference_sizes("mouse"), fail_download=True)
    with pytest.raises(GuidescanCheckError, match="could not download"):
        ensure_references("mouse", cache_dir=tmp_path, s3_client=client)
    assert list((tmp_path / "mouse").iterdir()) == []


def test_missing_reference_raises(tmp_path: Path) -> None:
    with pytest.raises(GuidescanCheckError, match="could not find reference"):
        ensure_references("mouse", cache_dir=tmp_path, s3_client=FakeS3({}))


# --- running and reporting -----------------------------------------------------


def _stub_runner(matching_sequences: set[str], seen: dict):
    """Stand in for GuidescanPipeline: NA for any sequence not in the set."""

    def runner(index_prefix: Path, guides_file: Path, output_dir: Path) -> dict:
        with guides_file.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        seen["rows"] = rows
        seen["index_prefix"] = index_prefix
        return {
            "id": [row["id"] for row in rows],
            "chromosome": [
                "1" if row["sequence"] in matching_sequences else "NA" for row in rows
            ],
        }

    return runner


def test_run_counts_exact_matches_and_maps_ids_back() -> None:
    guides = [
        TargetingGuide("0033", "AAAA", "NGG"),
        TargetingGuide("g2", "TTTA", "NGG"),
        TargetingGuide("g3", "ATAT", "NGG"),
    ]
    seen: dict = {}
    report = run_guidescan_check(
        guides, Path("idx"), runner=_stub_runner({"AAAA", "ATAT"}, seen)
    )
    assert (report.total, report.matched, report.unmatched) == (3, 2, ["g2"])
    assert [row["id"] for row in seen["rows"]] == ["g000001", "g000002", "g000003"]
    assert {row["pam"] for row in seen["rows"]} == {"NGG"}
    assert list(seen["rows"][0]) == [
        "id",
        "sequence",
        "pam",
        "chromosome",
        "position",
        "sense",
    ]


def test_run_with_no_guides_skips_guidescan() -> None:
    def runner(*args):
        raise AssertionError("guidescan must not run")

    report = run_guidescan_check([], None, runner=runner)
    assert report.total == 0
    assert report.notes == ["no targeting guides; guidescan not run"]
    assert format_report(report) == ["note: no targeting guides; guidescan not run"]


def test_report_line_and_no_warning_at_threshold() -> None:
    report = GuidescanReport(total=10, matched=9, unmatched=["g7"])
    assert format_report(report) == [
        "guidescan: 9 of 10 targeting guides (90.0%) have an exact match",
        "guidescan: no exact match: g7",
    ]


def test_warning_below_threshold_and_examples_capped() -> None:
    unmatched = [f"g{i:02d}" for i in range(12)]
    report = GuidescanReport(total=100, matched=88, unmatched=unmatched)
    lines = format_report(report)
    assert lines[0] == (
        "guidescan: 88 of 100 targeting guides (88.0%) have an exact match"
    )
    assert lines[1] == "guidescan: no exact match: " + ", ".join(unmatched[:10])
    assert lines[2] == "warning: guidescan exact-match rate 88.0% is below 90%"


# --- CLI -----------------------------------------------------------------------


@pytest.fixture
def fake_guidescan(monkeypatch):
    """Replace reference download and the pipeline run for CLI tests."""
    state: dict = {"matching": set(), "downloaded": False}

    def fake_refs(species: str, **kwargs) -> Path:
        state["downloaded"] = True
        state["species"] = species
        return Path("idx")

    def fake_pipeline(index_prefix, guides_file, output_dir):
        return _stub_runner(state["matching"], {})(
            index_prefix, guides_file, output_dir
        )

    monkeypatch.setattr(guidescancheck, "ensure_references", fake_refs)
    monkeypatch.setattr(guidescancheck, "_run_pipeline", fake_pipeline)
    return state


def _cli_sheet(tmp_path: Path, gene: str = MOUSE_GENE) -> Path:
    return _sheet(
        tmp_path,
        [
            _row("g1", "AAAA", gene=gene),
            _row("g2", "TTTA", gene=gene),
            _row("ctl", "ATAT", gene="NA", role="control"),
        ],
    )


def test_cli_reports_rate_and_prints_signature(
    tmp_path: Path, fake_guidescan, capsys
) -> None:
    fake_guidescan["matching"] = {"AAAA", "TTTA"}
    sheet = _cli_sheet(tmp_path)
    main(["--format", "tsv", "--guidescan", "--species", "mouse", str(sheet)])
    captured = capsys.readouterr()
    assert captured.out.startswith("gsig1:set:n=3:")
    assert "guidescan: 2 of 2 targeting guides (100.0%)" in captured.err
    assert "warning: guidescan" not in captured.err
    assert fake_guidescan["species"] == "mouse"


def test_cli_low_rate_warns_but_still_prints_signature(
    tmp_path: Path, fake_guidescan, capsys
) -> None:
    fake_guidescan["matching"] = {"AAAA"}
    sheet = _cli_sheet(tmp_path)
    main(["--format", "tsv", "--guidescan", "--species", "mouse", str(sheet)])
    captured = capsys.readouterr()
    assert captured.out.startswith("gsig1:set:n=3:")
    assert "guidescan: 1 of 2 targeting guides (50.0%)" in captured.err
    assert "guidescan: no exact match: g2" in captured.err
    assert "warning: guidescan exact-match rate 50.0% is below 90%" in captured.err


def test_cli_species_mismatch_fails_before_download(
    tmp_path: Path, fake_guidescan, capsys
) -> None:
    sheet = _cli_sheet(tmp_path, gene=HUMAN_GENE)
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", "--guidescan", "--species", "mouse", str(sheet)])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "--species mouse but guide_target_gene_id holds human" in captured.err
    assert fake_guidescan["downloaded"] is False


@pytest.mark.parametrize(
    "args, message",
    [
        (["--guidescan", "a.tsv"], "--guidescan requires --species"),
        (["--species", "mouse", "a.tsv"], "--species requires --guidescan"),
        (
            ["--guidescan", "--species", "mouse", "--compare", "a.tsv", "b.tsv"],
            "--guidescan is not allowed with --compare",
        ),
        (
            ["--guidescan", "--species", "mouse", "a.tsv", "b.tsv"],
            "--guidescan checks exactly one FILE",
        ),
    ],
)
def test_cli_guidescan_usage_errors(capsys, args: list[str], message: str) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--format", "tsv", *args])
    assert exc.value.code == 2
    assert message in capsys.readouterr().err
