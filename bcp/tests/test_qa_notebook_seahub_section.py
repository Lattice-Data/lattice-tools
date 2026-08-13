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
from types import SimpleNamespace

import pandas as pd
import pytest

from tests.qa_seahub_helpers import (
    BUCKET,
    PROJECT,
    VENDOR_BUCKET,
    ref3_sizes,
    ref3_trimmed_keys,
    ref3_vendor_keys,
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


def _run_section(
    tmp_path,
    *,
    raw_assay="seahub_sci",
    sources=None,
    keys=None,
    experiment_id="REF3",
    output_label="REF3",
    bucket=BUCKET,
):
    """Run the SeaHub block over a synthetic namespace.

    ``experiment_id`` and ``output_label`` are separate because the notebook's
    ``order`` is ``ctx.output_label``, which ``run_label`` overrides. Defaulting
    them to the same value mirrors an unset ``run_label``.
    """
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
        # A project root, not the two orders: discovery finds both from the
        # wafers in the upload, which is the point of the change.
        sources = [f"s3://{VENDOR_BUCKET}/{PROJECT}"]

    namespace.update(
        {
            "pd": pd,
            "display": lambda *a, **k: None,
            "s3client": MockS3Client(
                buckets={VENDOR_BUCKET: vendor_keys}, file_contents=sidecars
            ),
            "raw_assay": raw_assay,
            "bucket": bucket,
            "ctx": SimpleNamespace(order=experiment_id),
            "order": output_label,
            "output_dir": str(tmp_path),
            "validate_raw": True,
            "all_raw_files": trimmed_keys,
            "read_metadata": {},
            "raw_file_sizes": ref3_sizes(),
            "seahub_fail_counts": {},
            "sop_violations": [],
            "discovered_wafers": set(),
            "untrimmed_search_roots": list(sources),
            "untrimmed_search_siblings": True,
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

    def test_a_run_with_no_wells_still_writes_the_status_csv(self, tmp_path):
        _output, _ns = _run_section(tmp_path, keys=[], sources=[])

        frame = pd.read_csv(tmp_path / "REF3_seahub_well_status.csv")
        assert len(frame) == 0
        assert list(frame.columns)

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


class TestTheVendorIndexIsIndependentOfTheLabel:
    """No name the operator chooses can change which vendor wells are indexed.

    This used to be a narrow guarantee -- the ExperimentID filter read
    ``ctx.order`` rather than the notebook's ``order``, which ``run_label``
    overrides, and filtering on the wrong one matched no vendor well and orphaned
    the whole upload. Now that deliveries are located by the wafers inside them,
    the property is stronger and worth keeping under a name that says so: neither
    label reaches the vendor index at all.
    """

    def test_a_run_label_changes_nothing_about_the_vendor_index(self, tmp_path):
        plain_dir, labelled_dir = tmp_path / "plain", tmp_path / "labelled"
        plain_dir.mkdir()
        labelled_dir.mkdir()
        _output, plain = _run_section(plain_dir)
        _output, labelled = _run_section(
            labelled_dir, experiment_id="REF3", output_label="REF3-rerun"
        )

        assert labelled["untrimmed_index"], "vendor index was filtered away"
        assert sorted(labelled["untrimmed_index"]) == sorted(plain["untrimmed_index"])

    def test_the_orphan_count_is_unchanged_by_a_run_label(self, tmp_path):
        plain_dir, labelled_dir = tmp_path / "plain", tmp_path / "labelled"
        plain_dir.mkdir()
        labelled_dir.mkdir()
        _run_section(plain_dir)
        _run_section(labelled_dir, experiment_id="REF3", output_label="REF3-rerun")

        def orphans(directory, label):
            frame = pd.read_csv(directory / f"{label}_trimming_completeness.csv")
            return (frame["category"] == "orphan_trimmed").sum()

        assert orphans(tmp_path / "labelled", "REF3-rerun") == orphans(
            tmp_path / "plain", "REF3"
        )


class TestSectionSkips:
    def test_a_non_seahub_assay_skips_every_cell_and_writes_nothing(self, tmp_path):
        output, _ns = _run_section(tmp_path, raw_assay="10x")

        assert os.listdir(tmp_path) == []
        assert "raw_assay='10x' — skipped" in output

    def test_no_untrimmed_source_names_what_will_not_run(self, tmp_path):
        output, _ns = _run_section(tmp_path, sources=[])

        assert "NO untrimmed search root given" in output
        assert "trimming completeness" in output
        assert "untrimmed_search_roots = [" in output
        errors = (tmp_path / "REF3_errors.txt").read_text()
        assert "SEAHUB SOURCE COVERAGE: no untrimmed_search_roots set" in errors

    def test_without_a_source_the_data_gap_is_still_reported(self, tmp_path):
        """DATA_GAP outranks the un-nameable UNKNOWN, by design."""
        _output, _ns = _run_section(tmp_path, sources=[])
        frame = pd.read_csv(tmp_path / "REF3_seahub_well_status.csv")

        assert "DATA_GAP" in set(frame["verdict"])


class TestAnUploadInTheWrongPlaceIsNotSilent:
    """Upload-scope defects have to be reported when there is nothing to move.

    They are filtered out of the per-object residual check, so an upload that is
    otherwise SOP-clean produces no rename rows -- and the banner used to be
    gated on there being some. The result was inverted severity: a clean upload
    in the wrong bucket reported every well COMPLIANT and wrote nothing to
    errors.txt, while adding a trivial naming defect made the location defect
    audible.
    """

    STEM = "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT"
    SUFFIXES = (
        ".trim.cram",
        ".trim.csv",
        ".trim.stderr",
        ".trim.stdout",
        ".trim_fail.csv",
    )

    def _clean_keys(self) -> list[str]:
        return [
            f"labalpha-seahub-bcp/REF3/raw/REF3_P05_2/436830/{self.STEM}{s}"
            for s in self.SUFFIXES
        ]

    def _run(self, tmp_path, bucket):
        return _run_section(
            tmp_path, keys=self._clean_keys(), sources=[], bucket=bucket
        )

    def test_a_clean_upload_in_the_right_place_says_nothing(self):
        """The control: no banner when the location is fine."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            out, _ = self._run(tmp, "czi-labalpha")

        assert "not in a location the SOP allows" not in out

    @pytest.mark.parametrize(
        "bucket,rule",
        [
            ("czi-labbeta", "lab_project_mismatch"),
            ("labalpha-data", "bad_bucket"),
        ],
    )
    def test_the_banner_prints_with_nothing_to_move(self, tmp_path, bucket, rule):
        out, namespace = self._run(tmp_path, bucket)

        assert namespace["rename_mapping"].moveable() == []
        assert "not in a location the SOP allows" in out
        assert rule in out

    @pytest.mark.parametrize("bucket", ["czi-labbeta", "labalpha-data"])
    def test_it_reaches_errors_txt(self, tmp_path, bucket):
        self._run(tmp_path, bucket)

        errors = (tmp_path / "REF3_errors.txt").read_text()
        assert "SEAHUB UPLOAD LOCATION" in errors

    def test_the_wording_says_the_well_status_reads_clean(self, tmp_path):
        """Because it does -- and that is the trap this banner exists for."""
        out, _ = self._run(tmp_path, "czi-labbeta")

        assert "PER-WELL STATUS: COMPLIANT=1" in out
        assert "the right files" in out and "in the wrong place" in out
