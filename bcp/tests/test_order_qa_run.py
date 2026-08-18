"""Tests for order_qa.registry, ledger and runner.

The load-bearing test here is
``TestRegistryInvariant::test_every_declared_check_is_reported``. The tool exists
because the notebook can silently not run a check; if that invariant ever breaks,
this CLI has the same defect the notebook has and is worth less than the notebook.
"""

from datetime import datetime, timedelta, timezone

import pytest


from qa_gather import QAGatheredData

from order_qa import exits
from order_qa.ledger import Ledger, diff_against_ledger, load_ledger, save_ledger
from order_qa.registry import (
    REGISTRY,
    CheckContext,
    CheckResult,
    Scope,
    Status,
    run_registry,
)
from order_qa.report import RunReport, order_root, status_line
from order_qa.verify import S3Object
from order_qa.runner import RunOptions, _qa_verdict, run
from order_qa.spec import resolve_spec

NOW = datetime(2026, 8, 18, 12, 0, tzinfo=timezone.utc)
SPEC = resolve_spec("czi-psomagen/orglab-alpha/AN00000001")


def listing_page(keys, minutes_old=60, multipart=True):
    etag = '"abc-12"' if multipart else '"abc"'
    return [
        {
            "Contents": [
                {
                    "Key": key,
                    "Size": 1000,
                    "LastModified": NOW - timedelta(minutes=minutes_old),
                    "ETag": etag,
                }
                for key in keys
            ]
        }
    ]


class FakeS3:
    def __init__(self, pages, versions=None):
        self._pages = pages
        self._versions = versions or {}

    def get_paginator(self, operation):
        if operation == "list_object_versions":
            raise RuntimeError("ListObjectVersions denied")
        pages = self._pages
        outer = self

        class _P:
            def paginate(self, **kwargs):
                if operation == "list_objects":
                    return iter(outer._pages)
                return iter(pages)

        return _P()

    def head_object(self, Bucket, Key, PartNumber=None):
        if PartNumber is not None:
            return {"PartsCount": 12, "ContentLength": 8 << 20}
        return {"ContentLength": 1000, "VersionId": self._versions.get(Key, "v1")}

    def get_object_attributes(self, **kwargs):
        raise RuntimeError("GetObjectAttributes denied")

    def list_object_versions(self, **kwargs):
        raise RuntimeError("ListObjectVersions denied")

    def download_file(self, *args, **kwargs):
        raise RuntimeError("no processed files in this fixture")


def make_context(scope, data=None):
    return CheckContext(
        data=data or QAGatheredData(),
        scope=scope,
        bucket="czi-psomagen",
        order_label="AN00000001",
        s3_client=FakeS3([]),
        output_dir="",
    )


class TestRegistryInvariant:
    def test_every_declared_check_is_reported(self):
        """One result per declaration, always.

        This is the whole product: a caller can diff declared against reported and
        get nothing, so a check cannot go missing between the registry and the
        report the way a notebook cell goes missing between running and reading.
        """
        for assay in ("10x", "scale", "seahub_sci", "sci_plex", "10x_cram"):
            for raw, processed in ((True, True), (False, False), (True, False)):
                scope = Scope.resolve(
                    raw_assay=assay, validate_raw=raw, validate_processed=processed
                )
                results = run_registry(make_context(scope))
                assert len(results) == len(REGISTRY), f"{assay} raw={raw}"
                assert [r.name for r in results] == [c.name for c in REGISTRY]

    def test_no_check_name_is_duplicated(self):
        names = [c.name for c in REGISTRY]
        assert len(names) == len(set(names))

    def test_a_raising_check_does_not_end_the_run(self):
        """One broken check must not cost the other eleven their results."""
        scope = Scope.resolve(
            raw_assay="10x", validate_raw=True, validate_processed=False
        )
        results = run_registry(make_context(scope))
        assert len(results) == len(REGISTRY)
        assert all(isinstance(r, CheckResult) for r in results)


class TestApplicability:
    def test_scale_checks_skip_on_a_10x_order(self):
        scope = Scope.resolve(
            raw_assay="10x", validate_raw=True, validate_processed=True
        )
        by_name = {c.name: c for c in REGISTRY}
        applicable, reason = by_name["scale_cb_tag"].applicability(scope)
        assert applicable is False
        assert "scale" in reason

    def test_cellranger_skips_on_scale_because_scale_has_its_own_checks(self):
        scope = Scope.resolve(
            raw_assay="scale", validate_raw=True, validate_processed=True
        )
        by_name = {c.name: c for c in REGISTRY}
        assert by_name["processed_cellranger"].applicability(scope)[0] is False

    def test_seahub_is_raw_only_qa(self):
        """The notebook forces validate_processed=False for seahub_sci."""
        scope = Scope.resolve(
            raw_assay="seahub_sci", validate_raw=True, validate_processed=True
        )
        assert scope.validate_processed is False
        by_name = {c.name: c for c in REGISTRY}
        assert by_name["processed_cellranger"].applicability(scope)[0] is False

    def test_unimplemented_seahub_checks_report_not_run_not_skipped(self):
        """A gap must be declared, not omitted.

        Otherwise a seahub_sci order passes on the checks that are wired up while
        the cross-bucket reconciliation -- the point of QA-ing a trimmed upload --
        silently never happens.
        """
        scope = Scope.resolve(
            raw_assay="seahub_sci", validate_raw=True, validate_processed=False
        )
        results = {r.name: r for r in run_registry(make_context(scope))}
        recon = results["seahub_trimming_completeness"]
        assert recon.status is Status.NOT_RUN
        assert recon.status.denies_verdict is True
        assert "qa.ipynb" in recon.reason

    def test_skipped_does_not_deny_a_verdict_but_not_run_does(self):
        assert Status.SKIPPED.denies_verdict is False
        assert Status.PASS.denies_verdict is False
        assert Status.FINDINGS.denies_verdict is False
        assert Status.NOT_RUN.denies_verdict is True
        assert Status.ERROR.denies_verdict is True


class TestChecksAgainstGatheredData:
    def test_q30_findings_come_from_the_real_validator(self):
        data = QAGatheredData(pct_q30_values={"sub_a": 90.0, "sub_b": 40.0})
        scope = Scope.resolve(
            raw_assay="10x", validate_raw=True, validate_processed=False
        )
        results = {r.name: r for r in run_registry(make_context(scope, data))}
        assert results["raw_pct_q30"].status is Status.FINDINGS
        assert results["raw_pct_q30"].findings

    def test_q30_with_nothing_gathered_is_not_run_not_pass(self):
        """No Q30 values means Q30 is unknown, which is not the same as fine."""
        scope = Scope.resolve(
            raw_assay="10x", validate_raw=True, validate_processed=False
        )
        results = {r.name: r for r in run_registry(make_context(scope))}
        assert results["raw_pct_q30"].status is Status.NOT_RUN

    def test_extra_files_check_needs_the_expected_files_check_to_have_run(self):
        """Without the recognised-file set, an extra file cannot be identified."""
        scope = Scope.resolve(
            raw_assay="10x", validate_raw=True, validate_processed=False
        )
        results = {r.name: r for r in run_registry(make_context(scope))}
        assert results["raw_extra_files"].status is Status.NOT_RUN


class TestLedger:
    def test_silent_reupload_is_detected_by_version_identity(self, tmp_path):
        """A re-uploaded key leaves count and bytes identical.

        Versioning means the old version still exists, so only the current
        VersionId changes -- which is invisible to any count-and-sum check.
        """
        previous = Ledger(
            order_key=SPEC.key,
            run_at="2026-08-17T09:00:00+00:00",
            verdict="ok",
            file_count=2,
            version_by_key={"p/a": "v1", "p/b": "v1"},
            version_method="head_object",
        )
        save_ledger(tmp_path, previous)
        diff = diff_against_ledger(
            load_ledger(tmp_path),
            current_keys={"p/a", "p/b"},
            current_versions={"p/a": "v1", "p/b": "v2-RESEQUENCED"},
        )
        assert diff.reuploaded_keys == ["p/b"]
        assert diff.any_change is True
        assert "re-uploaded" in diff.summary

    def test_new_and_removed_keys(self, tmp_path):
        save_ledger(tmp_path, Ledger(run_at="t", version_by_key={"p/a": "v1"}))
        diff = diff_against_ledger(
            load_ledger(tmp_path),
            current_keys={"p/b"},
            current_versions={"p/b": "v1"},
        )
        assert diff.new_keys == ["p/b"]
        assert diff.removed_keys == ["p/a"]

    def test_first_run_is_not_a_change(self, tmp_path):
        diff = diff_against_ledger(
            load_ledger(tmp_path), current_keys={"p/a"}, current_versions={"p/a": "v"}
        )
        assert diff.compared is False
        assert diff.any_change is False

    def test_unversioned_previous_run_cannot_be_compared(self, tmp_path):
        save_ledger(tmp_path, Ledger(run_at="t", version_by_key={}))
        diff = diff_against_ledger(
            load_ledger(tmp_path), current_keys={"p/a"}, current_versions={"p/a": "v"}
        )
        assert diff.compared is False
        assert "cannot be detected" in diff.reason

    def test_keys_not_version_checked_are_counted_not_assumed_unchanged(self, tmp_path):
        """A bounded sweep must not report unchecked keys as identical."""
        save_ledger(
            tmp_path, Ledger(run_at="t", version_by_key={"p/a": "v1", "p/b": "v1"})
        )
        diff = diff_against_ledger(
            load_ledger(tmp_path),
            current_keys={"p/a", "p/b"},
            current_versions={"p/a": "v1"},
        )
        assert diff.unchanged == 1
        assert diff.uncomparable_keys == 1

    def test_corrupt_ledger_is_ignored_rather_than_fatal(self, tmp_path):
        (tmp_path / "ledger.json").write_text("{not json")
        assert load_ledger(tmp_path) is None

    def test_ledger_write_is_atomic_and_leaves_no_temp_files(self, tmp_path):
        save_ledger(tmp_path, Ledger(run_at="t", version_by_key={"p/a": "v"}))
        save_ledger(tmp_path, Ledger(run_at="t2", version_by_key={"p/a": "v2"}))
        assert sorted(p.name for p in tmp_path.iterdir()) == ["ledger.json"]
        assert load_ledger(tmp_path).run_at == "t2"


class TestExitCodePrecedence:
    def _report(self, checks=(), warnings=(), errors=()):
        report = RunReport(
            order_key=SPEC.key,
            bucket=SPEC.bucket,
            project=SPEC.project,
            order=SPEC.order,
            order_shape=SPEC.order_shape,
            assay="10x",
            started_at=NOW,
        )
        report.checks = list(checks)
        report.gathering_warnings = list(warnings)
        report.gathering_errors = list(errors)
        return report

    def test_all_pass_is_ok(self):
        report = self._report([CheckResult("a", "A", Status.PASS)])
        assert _qa_verdict(report) == exits.OK

    def test_findings_are_non_zero(self):
        """A finding here means the order is not releasable."""
        report = self._report([CheckResult("a", "A", Status.FINDINGS, findings=["x"])])
        assert _qa_verdict(report) == exits.QA_FINDINGS

    def test_unanswered_outranks_findings(self):
        """Findings are not worth triaging until the gaps are closed."""
        report = self._report(
            [
                CheckResult("a", "A", Status.FINDINGS, findings=["x"]),
                CheckResult("b", "B", Status.NOT_RUN, reason="unavailable"),
            ]
        )
        assert _qa_verdict(report) == exits.DEGRADED

    def test_skipped_checks_alone_are_still_ok(self):
        report = self._report(
            [
                CheckResult("a", "A", Status.PASS),
                CheckResult("b", "B", Status.SKIPPED, reason="scale only"),
            ]
        )
        assert _qa_verdict(report) == exits.OK

    def test_gathering_error_degrades(self):
        report = self._report([CheckResult("a", "A", Status.PASS)], errors=["boom"])
        assert _qa_verdict(report) == exits.DEGRADED

    def test_coverage_warning_degrades_even_when_checks_passed(self):
        """LISTING TRUNCATED means a check did not see everything it should."""
        report = self._report(
            [CheckResult("a", "A", Status.PASS)],
            warnings=["LISTING TRUNCATED for raw/"],
        )
        assert _qa_verdict(report) == exits.DEGRADED

    def test_advisory_warning_does_not_degrade(self):
        report = self._report(
            [CheckResult("a", "A", Status.PASS)],
            warnings=["some advisory note"],
        )
        assert _qa_verdict(report) == exits.OK


class TestRunnerEndToEnd:
    def _options(self, tmp_path, **overrides):
        defaults = dict(
            spec=SPEC,
            assay="10x",
            output_root=tmp_path,
            quiescence_minutes=15,
            now=NOW,
        )
        defaults.update(overrides)
        return RunOptions(**defaults)

    def test_still_uploading_stops_before_qa(self, tmp_path):
        """QA against a live upload reports absences as defects."""
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"], minutes_old=1))
        report = run(self._options(tmp_path), s3)
        assert report.exit_code == exits.STILL_UPLOADING
        assert report.checks == []
        assert report.report_dir.exists()

    def test_force_proceeds_past_quiescence_and_records_it(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"], minutes_old=1))
        report = run(self._options(tmp_path, force=True), s3)
        assert report.exit_code != exits.STILL_UPLOADING
        assert report.quiescence["forced"] is True

    def test_empty_prefix_is_verification_failure(self, tmp_path):
        report = run(self._options(tmp_path), FakeS3([{"Contents": []}]))
        assert report.exit_code == exits.VERIFICATION_FAILED

    def test_dry_run_verifies_without_gathering(self, tmp_path, monkeypatch):
        def explode(*args, **kwargs):
            raise AssertionError("--dry-run must not gather")

        monkeypatch.setattr("order_qa.runner.gather_qa_data", explode)
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"]))
        report = run(self._options(tmp_path, dry_run=True), s3)
        assert report.exit_code == exits.OK
        assert report.checks == []

    def test_report_and_ledger_are_written(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"]))
        report = run(self._options(tmp_path), s3)
        assert (report.report_dir / "report.md").exists()
        assert (report.report_dir / "summary.json").exists()
        root = order_root(tmp_path, SPEC.bucket, SPEC.project, SPEC.order)
        assert load_ledger(root) is not None
        assert load_ledger(root).version_by_key

    def test_rerunning_does_not_corrupt_the_prior_report(self, tmp_path, monkeypatch):
        """A resequence re-triggers QA on an order already reported on."""
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        keys = ["orglab-alpha/AN00000001/raw/a"]
        first = run(self._options(tmp_path), FakeS3(listing_page(keys)))
        first_body = (first.report_dir / "report.md").read_text()

        second = run(
            self._options(tmp_path),
            FakeS3(listing_page(keys), versions={keys[0]: "v2-RESEQUENCED"}),
        )
        assert second.report_dir != first.report_dir
        assert (first.report_dir / "report.md").read_text() == first_body
        assert second.rerun["reuploaded_keys"] == keys
        assert second.rerun["any_change"] is True

    def test_denied_capabilities_are_reported_not_claimed(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"]))
        report = run(self._options(tmp_path), s3)
        caps = report.capabilities["capabilities"]
        assert caps["get_object_attributes"]["status"] != "available"
        assert report.versions["method"] == "head_object"

    def test_multipart_objects_are_counted_so_no_md5_claim_is_made(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(
            listing_page(["orglab-alpha/AN00000001/raw/a.fastq.gz"], multipart=True)
        )
        report = run(self._options(tmp_path), s3)
        assert report.listing_summary["multipart_files"] == 1
        assert report.listing_summary["etag_is_md5_files"] == 0

    def test_status_line_is_one_line_of_valid_json(self, tmp_path, monkeypatch):
        import json

        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"]))
        report = run(self._options(tmp_path), s3)
        line = status_line(report)
        assert "\n" not in line
        parsed = json.loads(line)
        assert parsed["order"] == SPEC.key
        assert parsed["verdict"] == exits.name_for(report.exit_code)


class TestNothingValidatedIsNotAPass:
    """An order that no check applied to has not been QA'd.

    Found by review: with neither raw/ nor processed/ detected -- a vendor using
    an unexpected layout, or the wrong --assay -- every check skips legitimately
    and the run summed those skips into exit 0. That is precisely the false clean
    pass this tool exists to prevent, arrived at from the other direction.
    """

    def _report(self, checks):
        report = RunReport(
            order_key=SPEC.key,
            bucket=SPEC.bucket,
            project=SPEC.project,
            order=SPEC.order,
            order_shape=SPEC.order_shape,
            assay="10x",
            started_at=NOW,
        )
        report.checks = list(checks)
        return report

    def test_all_skipped_is_degraded_not_ok(self):
        report = self._report(
            [
                CheckResult(f"c{i}", "T", Status.SKIPPED, reason="not applicable")
                for i in range(12)
            ]
        )
        assert _qa_verdict(report) == exits.DEGRADED

    def test_one_real_pass_among_skips_is_ok(self):
        report = self._report(
            [
                CheckResult("a", "A", Status.PASS, summary="fine"),
                CheckResult("b", "B", Status.SKIPPED, reason="scale only"),
            ]
        )
        assert _qa_verdict(report) == exits.OK

    def test_one_real_finding_among_skips_is_findings(self):
        report = self._report(
            [
                CheckResult("a", "A", Status.FINDINGS, findings=["x"]),
                CheckResult("b", "B", Status.SKIPPED, reason="scale only"),
            ]
        )
        assert _qa_verdict(report) == exits.QA_FINDINGS

    def test_no_recognizable_layout_degrades_end_to_end(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )
        s3 = FakeS3(listing_page(["orglab-alpha/AN00000001/raw/a"]))
        report = run(
            RunOptions(
                spec=SPEC,
                assay="10x",
                output_root=tmp_path,
                quiescence_minutes=15,
                now=NOW,
            ),
            s3,
        )
        assert report.exit_code == exits.DEGRADED


class TestLedgerBaselineIsProtected:
    """A run that could not see the whole order must not become the baseline.

    Found by review: the ledger was replaced on every path. A failed listing then
    overwrote a good baseline with a partial one, so the next run reported every
    unseen object as new -- and the version identities that would have caught a
    silent re-upload were gone.
    """

    def _options(self, tmp_path, **overrides):
        defaults = dict(
            spec=SPEC,
            assay="10x",
            output_root=tmp_path,
            quiescence_minutes=15,
            now=NOW,
        )
        defaults.update(overrides)
        return RunOptions(**defaults)

    @pytest.fixture(autouse=True)
    def _no_gather(self, monkeypatch):
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )

    def test_good_baseline_survives_a_failed_listing(self, tmp_path):
        keys = ["orglab-alpha/AN00000001/raw/a", "orglab-alpha/AN00000001/raw/b"]
        run(self._options(tmp_path), FakeS3(listing_page(keys)))
        root = order_root(tmp_path, SPEC.bucket, SPEC.project, SPEC.order)
        good = load_ledger(root)
        assert set(good.version_by_key) == set(keys)

        class BrokenS3(FakeS3):
            def get_paginator(self, operation):
                raise RuntimeError("AccessDenied listing")

        report = run(self._options(tmp_path), BrokenS3([]))
        assert report.exit_code == exits.DEGRADED
        after = load_ledger(root)
        assert set(after.version_by_key) == set(keys), "baseline was overwritten"
        assert after.run_at == good.run_at

    def test_still_uploading_run_does_not_become_the_baseline(self, tmp_path):
        keys = ["orglab-alpha/AN00000001/raw/a"]
        run(self._options(tmp_path), FakeS3(listing_page(keys)))
        root = order_root(tmp_path, SPEC.bucket, SPEC.project, SPEC.order)
        first = load_ledger(root)

        in_flight = ["orglab-alpha/AN00000001/raw/a", "orglab-alpha/AN00000001/raw/b"]
        report = run(
            self._options(tmp_path), FakeS3(listing_page(in_flight, minutes_old=1))
        )
        assert report.exit_code == exits.STILL_UPLOADING
        assert load_ledger(root).run_at == first.run_at

    def test_the_failed_run_still_writes_its_own_report(self, tmp_path):
        """The report is always written; only the baseline is withheld."""

        class BrokenS3(FakeS3):
            def get_paginator(self, operation):
                raise RuntimeError("AccessDenied listing")

        report = run(self._options(tmp_path), BrokenS3([]))
        assert (report.report_dir / "report.md").exists()
        assert "Listing incomplete" in (report.report_dir / "report.md").read_text()

    def test_forced_run_may_become_the_baseline(self, tmp_path):
        """--force is an explicit decision that the observation counts."""
        keys = ["orglab-alpha/AN00000001/raw/a"]
        report = run(
            self._options(tmp_path, force=True),
            FakeS3(listing_page(keys, minutes_old=1)),
        )
        root = order_root(tmp_path, SPEC.bucket, SPEC.project, SPEC.order)
        assert report.quiescence["forced"] is True
        assert load_ledger(root) is not None


class TestManifestSizeAttribution:
    """A size pinned to the wrong key invents failures.

    Found by review: taking the first bare integer on the line read the lane out
    of ``sampleA,1,s3://...,1048576`` and reported every file as the wrong size.
    """

    def test_unambiguous_size_is_used(self, tmp_path):
        from order_qa.verify import read_manifest_keys

        manifest = tmp_path / "m.tsv"
        manifest.write_text("s3://b/p/a.fastq.gz\t1048576\n")
        assert read_manifest_keys(manifest) == {"p/a.fastq.gz": 1048576}

    def test_ambiguous_size_is_refused_not_guessed(self, tmp_path):
        from order_qa.verify import read_manifest_keys

        manifest = tmp_path / "m.csv"
        manifest.write_text("sampleA,1,s3://b/p/a.fastq.gz,1048576\n")
        assert read_manifest_keys(manifest) == {"p/a.fastq.gz": None}

    def test_refused_size_does_not_produce_a_mismatch(self, tmp_path):
        from order_qa.verify import Listing, compare_to_manifest

        manifest = tmp_path / "m.csv"
        manifest.write_text("sampleA,1,s3://b/p/a.fastq.gz,1048576\n")
        listing = Listing(
            "b",
            "p/",
            [
                S3Object(
                    key="p/a.fastq.gz",
                    size_bytes=999,
                    last_modified=NOW,
                    etag='"x"',
                )
            ],
        )
        comparison = compare_to_manifest(listing, manifest)
        assert comparison.size_mismatches == []
        assert comparison.ok is True


class TestNoInputIsNotAPass:
    """A check with nothing to work on has not passed.

    Found by running the tool: an order whose files matched no recognised naming
    reported PASS on four checks -- including ``raw_fastq_counts``, whose own
    summary read "nothing to compare", and ``raw_expected_files`` at
    "0/1 file groups complete". Four green checks on an unreadable order.
    """

    def _results(self, data, assay="10x"):
        scope = Scope.resolve(
            raw_assay=assay, validate_raw=True, validate_processed=False
        )
        return {r.name: r for r in run_registry(make_context(scope, data))}

    def test_unrecognised_naming_is_not_complete(self):
        data = QAGatheredData(
            all_raw_files=["p/o/raw/nothing_recognisable.txt"], has_raw=True
        )
        result = self._results(data)["raw_expected_files"]
        assert result.status is Status.NOT_RUN
        assert "recognised" in result.reason

    def test_extra_files_cascades_to_not_run(self):
        """Without a recognised set, every file looks equally unexpected."""
        data = QAGatheredData(
            all_raw_files=["p/o/raw/nothing_recognisable.txt"], has_raw=True
        )
        assert self._results(data)["raw_extra_files"].status is Status.NOT_RUN

    def test_empty_fastq_log_is_not_run(self):
        result = self._results(QAGatheredData(has_raw=True))["raw_fastq_counts"]
        assert result.status is Status.NOT_RUN
        assert "nothing to compare" in result.reason

    def test_empty_read_metadata_is_not_run(self):
        result = self._results(QAGatheredData(has_raw=True))["raw_read_metadata"]
        assert result.status is Status.NOT_RUN

    def test_empty_q30_is_not_run(self):
        result = self._results(QAGatheredData(has_raw=True))["raw_pct_q30"]
        assert result.status is Status.NOT_RUN

    def test_an_order_with_no_usable_input_degrades(self):
        """The whole point: four not_run beats four pass."""
        data = QAGatheredData(
            all_raw_files=["p/o/raw/nothing_recognisable.txt"], has_raw=True
        )
        results = self._results(data)
        assert all(
            results[name].status is Status.NOT_RUN
            for name in (
                "raw_expected_files",
                "raw_extra_files",
                "raw_fastq_counts",
                "raw_read_metadata",
                "raw_pct_q30",
            )
        )

    def test_real_input_still_passes(self):
        """The guards must not turn a genuinely clean order into a degraded one."""
        data = QAGatheredData(pct_q30_values={"s1": 91.0, "s2": 88.0}, has_raw=True)
        result = self._results(data)["raw_pct_q30"]
        assert result.status is Status.PASS
        assert "2/2" in result.summary
