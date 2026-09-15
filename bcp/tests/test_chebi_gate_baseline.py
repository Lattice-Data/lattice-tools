"""Regression anchor: the findings the gate produces over the batch already sent.

283 of the 290 records in ``chebi_run_2026_08b`` were deposited with ChEBI, so this
is not a pre-send check -- it is the fixed point that stops a later change to the
checks from quietly altering what the gate says. The handoff's terms: the original
two-file validator produced 821 findings and any divergence from them has to be a
recorded, reasoned choice rather than a drift. The current figures live in
``findings_baseline.json`` and are asserted from it -- deliberately not repeated
here, because a number written into prose goes stale and this docstring already
did, claiming 814 and "two deliberate choices" after a third was added.

The inputs cannot be committed: bcp/.gitignore excludes ``chebi_run_*/`` because a
run directory is regenerable network cache and data derived from one spreadsheet.
So these tests skip with a reason when the directory is absent, and the expected
numbers live in ``fixtures/chebi_gate/findings_baseline.json`` alongside the
SHA256 of each input, so a run that does happen is known to be judging the same
bytes. Every check itself is covered by fixture-based tests that always run; this
file adds scale and the historical anchor, not coverage.
"""

from __future__ import annotations

import collections
import hashlib
import json
from pathlib import Path

import pytest

from chebi_gate import checks, sdf, structure

BCP = Path(__file__).resolve().parent.parent
BASELINE = Path(__file__).parent / "fixtures" / "chebi_gate" / "findings_baseline.json"


@pytest.fixture(scope="module")
def baseline() -> dict:
    return json.loads(BASELINE.read_text())


@pytest.fixture(scope="module")
def run_dir(baseline) -> Path:
    path = BCP.parent / baseline["generated_from"]
    if not path.is_dir():
        pytest.skip(
            f"{baseline['generated_from']} is not present (it is gitignored); "
            "the per-check tests cover the logic, this file adds the 290-record anchor"
        )
    missing = [name for name in baseline["inputs"] if not (path / name).exists()]
    if missing:
        pytest.skip(f"run directory present but missing {missing}")
    return path


def _findings(path: Path, baseline: dict) -> collections.Counter:
    counts: collections.Counter = collections.Counter()
    for name, role in baseline["roles"].items():
        parsed = sdf.parse_file(path / name)
        ctxs = [
            checks.RecordContext(record=r, structure=structure.analyse(r), role=role)
            for r in parsed.records
        ]
        for ctx in ctxs:
            for finding in checks.run_record_checks(ctx):
                counts[f"{finding.check}/{finding.severity}"] += 1
        for finding in checks.run_file_checks(
            checks.FileContext(contexts=ctxs, raw=parsed.raw)
        ):
            counts[f"{finding.check}/{finding.severity}"] += 1
    return counts


def test_the_inputs_are_the_bytes_the_baseline_was_built_from(run_dir, baseline):
    """Without this, a divergence below could be a different input, not a regression."""
    for name, expected in baseline["inputs"].items():
        data = (run_dir / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == expected["sha256"], name
        assert len(data) == expected["bytes"], name


def test_record_counts_match(run_dir, baseline):
    for name, expected in baseline["inputs"].items():
        parsed = sdf.parse_file(run_dir / name)
        assert len(parsed.records) == expected["records"], name


def test_every_record_round_trips_byte_identically(run_dir, baseline):
    """Handoff invariant 19, at scale: 290 real records, not a fixture."""
    for name in baseline["inputs"]:
        data = (run_dir / name).read_bytes()
        parsed = sdf.parse_bytes(data)
        assert parsed.dumps() == data, name
        rejoined = b"".join(r.raw + sdf.RECORD_TERMINATOR for r in parsed.records)
        assert rejoined == data, name


def test_findings_match_the_recorded_baseline(run_dir, baseline):
    counts = _findings(run_dir, baseline)
    expected = collections.Counter(baseline["counts_by_check_severity"])
    assert dict(counts) == dict(expected)
    assert sum(counts.values()) == baseline["total_findings"]


def test_the_v1_divergence_is_only_the_recorded_choices(run_dir, baseline):
    """Every divergence from v1 is one the fixture records, with its reason."""
    counts = _findings(run_dir, baseline)
    v1 = baseline["v1_baseline"]
    original = collections.Counter(v1["counts_by_check_severity"])

    diverged = {
        key
        for key in set(original) | set(counts)
        if original.get(key, 0) != counts.get(key, 0)
    }
    assert diverged == {d["key"] for d in v1["deliberate_divergences"]}

    for entry in v1["deliberate_divergences"]:
        assert original.get(entry["key"], 0) == entry["baseline"]
        assert counts.get(entry["key"], 0) == entry["gate"]
        assert entry["why"]

    reproduced = sum(v for k, v in original.items() if counts.get(k) == v)
    assert reproduced == v1["rows_reproduced_identically"]


def test_the_baseline_records_why_each_divergence_is_intended(baseline):
    """Runs without the inputs: a gap with no stated reason is just an unfixed bug."""
    divergences = baseline["v1_baseline"]["deliberate_divergences"]
    # Not a hard-coded count. Handoff case 18 is about exactly this: a figure
    # written into prose goes stale, and an assertion carrying a magic number
    # fails for the right reason once and then gets "fixed" by bumping it. What
    # matters is that every divergence names a real check and states a reason.
    assert divergences
    for entry in divergences:
        assert entry["key"] in baseline["v1_baseline"]["counts_by_check_severity"]
        assert entry["baseline"] != entry["gate"], (
            f"{entry['key']} is listed as a divergence but the counts agree"
        )
        assert len(entry["why"]) > 40


def test_the_baseline_covers_every_check_that_fired(baseline):
    """A check ID in the baseline that the registry no longer has is a silent loss."""
    for key in baseline["counts_by_check_severity"]:
        check_id, severity = key.split("/")
        assert check_id in checks.CHECKS, check_id
        assert severity in checks.CHECKS[check_id].severities, key
