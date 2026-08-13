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
from collections import Counter
from types import SimpleNamespace

import pandas as pd
import pytest

from tests.qa_seahub_helpers import (
    BUCKET,
    PROJECT,
    RAW,
    VENDOR_BUCKET,
    VENDOR_ORDER,
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
    vendor_keys=None,
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
    vendor_keys = ref3_vendor_keys() if vendor_keys is None else vendor_keys
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


class TestTheErrorFileStaysReadable:
    """errors.txt has to stay bounded as well as complete."""

    def test_the_data_gap_lines_are_bounded(self, tmp_path):
        """An upload that misspells everything makes every well a gap.

        480 near-identical lines on one real listing, for a defect the SOP cell
        already reports once with a FAIL banner. The count and the CSV keep it
        audible without the file becoming unreadable.
        """
        stem = "437120-REF3_P04_1_A{}_GEX_hash_oligo-Z{:04d}-CAGCTCGAATGCGAT"
        keys = [
            f"{RAW}/REF3_P04_1/437120/{stem.format(i, i)}.trimmed.ucram"
            for i in range(1, 26)
        ]
        _output, _ns = _run_section(tmp_path, keys=keys, sources=[])
        errors = (tmp_path / "REF3_errors.txt").read_text()

        gap_lines = [
            ln for ln in errors.splitlines() if ln.startswith("SEAHUB DATA GAP")
        ]
        assert len(gap_lines) <= 11, len(gap_lines)
        assert any("more of" in ln for ln in gap_lines)


class TestTheReconCellBoundsItsErrorFileWrites:
    """The recon cell's own errors.txt loop must be bounded like the console above it.

    metadata_unavailable is in practice a whole-upload fact -- vendor sidecars are
    generated for the upload as a whole, so a run where CZI never made them
    reports it for every matched well. Unbounded, that write loop repeats near
    identically once per well, exactly what the console output five lines above
    it, and the DATA_GAP loop in the well-status cell, both already bound.
    """

    N_WELLS = 12  # > RECON_EXAMPLES_PER_CATEGORY (5), so the cap is exercised

    def _many_matched_wells_with_no_sidecar(self):
        """N trimmed wells, each with a vendor CRAM but no ``.cram-metadata.json``.

        No sidecar means no read_count, so every one of them degrades to
        ``metadata_unavailable`` rather than being checked against the trimmer's
        declared totals.
        """
        trimmed, vendor = [], []
        for i in range(self.N_WELLS):
            wafer = f"5000{i:02d}"
            stem = f"{wafer}-REF3_P05_1_A{i}_GEX_hash_oligo-Z0{i:03d}-CAGTCAGTTGCAGAT"
            trimmed += [
                f"{RAW}/REF3_P05_1/{wafer}/{stem}{suffix}"
                for suffix in (
                    ".trim.cram",
                    ".trim.csv",
                    ".trim.stderr",
                    ".trim.stdout",
                    ".trim_fail.csv",
                )
            ]
            vendor.append(f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/{wafer}/{stem}.cram")
        return trimmed, vendor

    def test_metadata_unavailable_is_capped_in_the_error_file(self, tmp_path):
        trimmed, vendor = self._many_matched_wells_with_no_sidecar()
        _output, ns = _run_section(tmp_path, keys=trimmed, vendor_keys=vendor)

        report = ns["trimming_report"]
        unavailable = [
            r for r in report.rows if r["category"] == "metadata_unavailable"
        ]
        assert len(unavailable) == self.N_WELLS, "fixture must exceed the cap"

        errors = (tmp_path / "REF3_errors.txt").read_text()
        lines = [
            ln
            for ln in errors.splitlines()
            if ln.startswith("TRIMMING METADATA_UNAVAILABLE")
        ]
        assert len(lines) <= 6, len(
            lines
        )  # RECON_EXAMPLES_PER_CATEGORY + 1 "more" line
        assert any("more of" in ln for ln in lines)

    def test_a_category_under_the_cap_is_not_truncated(self, tmp_path):
        """The bound must not drop rows that fit inside it."""
        _output, ns = _run_section(tmp_path)
        report = ns["trimming_report"]

        errors = (tmp_path / "REF3_errors.txt").read_text()
        categories = ns["ERRORS_TXT_CATEGORIES"]
        for category, count in Counter(r["category"] for r in report.rows).items():
            if category not in categories or count > 5:
                continue
            lines = [
                ln
                for ln in errors.splitlines()
                if ln.startswith(f"TRIMMING {category.upper()}")
            ]
            assert len(lines) == count, category


class TestASearchRootWithNoWafersDoesNotCrash:
    """The two cells must agree about whether a vendor index exists.

    The index cell needs a root *and* a wafer to search for; the recon cell
    checked only the root. With a root set and an upload yielding no seed,
    reconcile_trimming received None and raised TypeError mid-run -- aborting the
    rename cell and the headline per-well table with it.

    Zero seeds is reachable on exactly the input this mode exists to catch: every
    seed reading needs a six-segment well path, so an upload delivered flat under
    raw/ has none. That is the case the SOP block reports as bad_path_depth, and
    it must fail loudly rather than with a stack trace.
    """

    FLAT = f"{RAW}/438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT.trim.cram"

    def test_the_section_completes_and_writes_the_headline(self, tmp_path):
        output, _ns = _run_section(tmp_path, keys=[self.FLAT])

        assert "Skipping trimming reconciliation" in output
        assert "no wafer to search for" in output
        assert "REF3_seahub_well_status.csv" in os.listdir(tmp_path)

    def test_it_says_so_in_the_error_file(self, tmp_path):
        _run_section(tmp_path, keys=[self.FLAT])
        errors = (tmp_path / "REF3_errors.txt").read_text()

        assert "TRIMMING NOT RECONCILED" in errors

    def test_an_empty_root_list_still_says_that_instead(self, tmp_path):
        output, _ns = _run_section(tmp_path, sources=[])

        assert "untrimmed_search_roots is empty" in output
