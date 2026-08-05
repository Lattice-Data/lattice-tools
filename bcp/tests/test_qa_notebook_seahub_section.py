"""
Executes the notebook's trailing SeaHub section against the REF3 fixture.

The layout test next door pins cell *order*. This one pins that the cells
actually run: it binds the namespace the shared pipeline leaves behind, executes
every SeaHub code cell in order, and asserts on the files they write. That is the
gap between "the notebook parses" and "the notebook works" -- and the only way to
catch a cell referring to a name the modules no longer export.
"""

from __future__ import annotations

import ast
import contextlib
import io
import json
import os

import pandas as pd
import pytest

from tests.qa_seahub_helpers import (
    BUCKET,
    ref3_sizes,
    ref3_trimmed_keys,
    ref3_vendor_keys,
    vendor_uri,
)
from tests.test_qa_gather import MockS3Client

NOTEBOOK_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "qa.ipynb")
VENDOR_READ_COUNT = 260527531


def _notebook_cells():
    with open(NOTEBOOK_PATH) as fh:
        return json.load(fh)["cells"]


def _import_namespace(cell_zero) -> dict:
    """Execute only cell 0's imports -- no boto3 client, no AWS configuration."""
    tree = ast.parse("".join(cell_zero["source"]))
    imports = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    namespace: dict = {}
    exec(compile(ast.Module(body=imports, type_ignores=[]), "cell0", "exec"), namespace)
    return namespace


def _run_section(tmp_path, *, raw_assay="seahub_sci", sources=None, keys=None):
    cells = _notebook_cells()
    namespace = _import_namespace(cells[0])

    trimmed_keys = ref3_trimmed_keys() if keys is None else keys
    vendor_keys = ref3_vendor_keys()
    sidecars = {
        key: json.dumps({"read_count": VENDOR_READ_COUNT})
        for key in vendor_keys
        if key.endswith(".cram-metadata.json")
    }
    if sources is None:
        sources = [vendor_uri("NVUS0000000000-11"), vendor_uri("NVUS0000000000-12")]

    namespace.update(
        {
            "pd": pd,
            "display": lambda *a, **k: None,
            "s3client": MockS3Client(keys=vendor_keys, file_contents=sidecars),
            "raw_assay": raw_assay,
            "bucket": BUCKET,
            "order": "REF3",
            "output_dir": str(tmp_path),
            "validate_raw": True,
            "all_raw_files": trimmed_keys,
            "read_metadata": {},
            "raw_file_sizes": ref3_sizes(),
            "seahub_fail_counts": {},
            "sop_violations": [],
            "discovered_wafers": set(),
            "untrimmed_s3_paths": list(sources),
        }
    )

    seahub = [c for c in cells if "seahub" in c.get("metadata", {}).get("tags", [])]
    assert seahub, "no seahub-tagged cells found"

    stdout = io.StringIO()
    for cell in seahub:
        if cell["cell_type"] != "code":
            continue
        with contextlib.redirect_stdout(stdout):
            exec(compile("".join(cell["source"]), "seahub-cell", "exec"), namespace)
    return stdout.getvalue(), namespace


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("qa_outs")
    output, namespace = _run_section(tmp_path)
    return output, namespace, tmp_path


class TestSectionRuns:
    def test_every_expected_output_is_written(self, run):
        _output, _ns, tmp_path = run

        assert sorted(os.listdir(tmp_path)) == [
            "REF3_errors.txt",
            "REF3_raw_sop_violations.csv",
            "REF3_seahub_rename_mapping.csv",
            "REF3_seahub_source_coverage.csv",
            "REF3_seahub_well_status.csv",
            "REF3_trimming_completeness.csv",
        ]

    def test_the_headline_is_the_per_well_status(self, run):
        output, _ns, _tmp = run

        assert (
            "PER-WELL STATUS: COMPLIANT=1, RENAMEABLE=3, DATA_GAP=1, UNKNOWN=1"
            in output
        )

    def test_the_well_status_csv_has_one_row_per_well(self, run):
        _output, _ns, tmp_path = run
        frame = pd.read_csv(tmp_path / "REF3_seahub_well_status.csv")

        assert len(frame) == 6
        assert set(frame["verdict"]) == {
            "COMPLIANT",
            "RENAMEABLE",
            "DATA_GAP",
            "UNKNOWN",
        }

    def test_the_rename_mapping_never_reuses_a_destination(self, run):
        _output, _ns, tmp_path = run
        frame = pd.read_csv(tmp_path / "REF3_seahub_rename_mapping.csv")
        moveable = frame[frame["status"] == "rename"]

        assert len(moveable) == len(set(moveable["expected_s3_uri"]))
        assert moveable["expected_s3_uri"].str.startswith(f"s3://{BUCKET}/").all()

    def test_junk_is_reported_without_a_destination(self, run):
        output, _ns, tmp_path = run
        frame = pd.read_csv(tmp_path / "REF3_seahub_rename_mapping.csv")
        junk = frame[frame["status"] == "not_data"]

        assert len(junk) == 5
        assert junk["expected_s3_uri"].isna().all()
        assert "likely files unrelated to the CRO sequencing" in output

    def test_the_unsourced_wafer_is_called_out(self, run):
        """Wafer 439000 has no vendor order, mirroring REF3_P05_1 in the real data."""
        output, _ns, tmp_path = run
        errors = (tmp_path / "REF3_errors.txt").read_text()

        assert "have no listed untrimmed source: 439000" in output
        assert "SEAHUB SOURCE COVERAGE: wafer(s) 439000" in errors

    def test_the_data_gap_reaches_the_errors_file_with_its_vendor_key(self, run):
        _output, _ns, tmp_path = run
        errors = (tmp_path / "REF3_errors.txt").read_text()

        assert "SEAHUB DATA GAP: wafer 438514 Z0305" in errors
        assert "the vendor delivered it as" in errors

    def test_orphans_are_not_written_to_the_errors_file(self, run):
        """With an order missing they are noise; the coverage line covers it."""
        _output, _ns, tmp_path = run
        errors = (tmp_path / "REF3_errors.txt").read_text()

        assert "ORPHAN_TRIMMED" not in errors

    def test_the_sop_type_is_filled_from_the_vendor_delivery(self, run):
        _output, _ns, tmp_path = run
        frame = pd.read_csv(tmp_path / "REF3_raw_sop_violations.csv")
        typed = frame[frame["type"] == "invalid_sublibrary_type"]

        assert typed["detail"].str.contains("untrimmed vendor delivery").any()

    def test_the_folder_rule_is_reported_per_sublibrary(self, run):
        _output, _ns, tmp_path = run
        frame = pd.read_csv(tmp_path / "REF3_raw_sop_violations.csv")
        folder = frame[frame["type"] == "sublibrary_folder_truncated"]

        # Four truncated folders in the fixture, whatever the well count beneath.
        assert len(folder) == 4
        assert set(folder["expected_folder"]) == {
            "REF3_P04_1",
            "REF3_P05_1",
            "REF3_P06_1",
            "REF3_P07_1",
        }


class TestSectionSkips:
    def test_a_non_seahub_assay_skips_every_cell_and_writes_nothing(self, tmp_path):
        output, _ns = _run_section(tmp_path, raw_assay="10x")

        assert os.listdir(tmp_path) == []
        assert "raw_assay='10x' — skipped" in output

    def test_no_untrimmed_source_names_what_will_not_run(self, tmp_path):
        output, _ns = _run_section(tmp_path, sources=[])

        assert "NO untrimmed source given" in output
        assert "trimming completeness" in output
        assert "untrimmed_s3_paths = [" in output
        errors = (tmp_path / "REF3_errors.txt").read_text()
        assert "SEAHUB SOURCE COVERAGE: no untrimmed_s3_paths set" in errors

    def test_without_a_source_the_data_gap_is_still_reported(self, tmp_path):
        """DATA_GAP outranks the un-nameable UNKNOWN, by design."""
        _output, _ns = _run_section(tmp_path, sources=[])
        frame = pd.read_csv(tmp_path / "REF3_seahub_well_status.csv")

        assert "DATA_GAP" in set(frame["verdict"])
