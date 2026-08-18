"""Tests for order_qa.cli: argument handling, exit codes, and stream discipline.

stdout must carry only the status line. A tracker consuming
``order-qa ... | jq`` breaks the moment a log message lands on stdout, so that
separation is asserted rather than assumed.
"""

import json
import logging
import sys
from datetime import datetime, timedelta, timezone

import pytest

from qa_gather import QAGatheredData

from order_qa import exits
from order_qa.cli import build_parser, main

NOW = datetime.now(timezone.utc)


class FakeS3:
    def __init__(self, keys=(), minutes_old=60):
        self._keys = list(keys)
        self._minutes_old = minutes_old

    def get_paginator(self, operation):
        if operation == "list_object_versions":
            raise RuntimeError("denied")
        pages = [
            {
                "Contents": [
                    {
                        "Key": key,
                        "Size": 1000,
                        "LastModified": NOW - timedelta(minutes=self._minutes_old),
                        "ETag": '"abc-12"',
                    }
                    for key in self._keys
                ]
            }
        ]

        class _P:
            def paginate(self, **kwargs):
                return iter(pages)

        return _P()

    def head_object(self, Bucket, Key, PartNumber=None):
        if PartNumber is not None:
            return {"PartsCount": 12, "ContentLength": 8 << 20}
        return {"ContentLength": 1000, "VersionId": "v1"}

    def get_object_attributes(self, **kwargs):
        raise RuntimeError("denied")

    def list_object_versions(self, **kwargs):
        raise RuntimeError("denied")


@pytest.fixture
def stub_aws(monkeypatch):
    """Install a fake S3 client and a gather that returns nothing to check."""

    def install(keys=("orglab-alpha/AN00000001/raw/a.fastq.gz",), minutes_old=60):
        monkeypatch.setattr(
            "order_qa.cli._make_client",
            lambda: FakeS3(keys, minutes_old=minutes_old),
        )
        monkeypatch.setattr(
            "order_qa.runner.gather_qa_data",
            lambda ctx, client: QAGatheredData(has_raw=False, has_processed=False),
        )

    return install


class TestArgumentHandling:
    def test_assay_is_required_for_qa(self, caplog):
        assert main(["czi-psomagen/orglab-alpha/AN00000001"]) == exits.USAGE
        assert "--assay is required" in caplog.text

    def test_dry_run_needs_no_assay(self, stub_aws, tmp_path, caplog):
        """Verification is assay-independent, so --dry-run must not demand one."""
        stub_aws()
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--dry-run",
                "--output-dir",
                str(tmp_path),
            ]
        )
        assert code == exits.OK
        assert "--assay is required" not in caplog.text

    def test_unresolvable_spec_exits_resolution_failed(self):
        assert main(["czi-psomagen/orglab-alpha"]) == exits.RESOLUTION_FAILED

    def test_non_order_path_is_refused(self):
        assert main(["czi-psomagen/orglab-beta/jsonlds"]) == exits.RESOLUTION_FAILED

    def test_missing_manifest_is_a_usage_error(self, tmp_path):
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "10x",
                "--manifest",
                str(tmp_path / "nope.tsv"),
            ]
        )
        assert code == exits.USAGE

    def test_unknown_assay_is_a_usage_error(self, stub_aws, tmp_path):
        """resolve_qa_run_context rejects it; that is usage, not a fact about S3."""
        stub_aws()
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "not_an_assay",
                "--output-dir",
                str(tmp_path),
            ]
        )
        assert code == exits.USAGE

    def test_raw_and_no_raw_are_mutually_exclusive(self):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["spec", "--raw", "--no-raw"])

    def test_scope_flags_default_to_autodetect(self):
        args = build_parser().parse_args(["spec"])
        assert args.validate_raw is None
        assert args.validate_processed is None


class TestStreamDiscipline:
    def test_stdout_carries_only_the_status_line(
        self, stub_aws, tmp_path, capsys, caplog
    ):
        stub_aws()
        # basicConfig is a no-op once caplog has installed a root handler, so the
        # root logger stays at WARNING and INFO records would be dropped.
        caplog.set_level(logging.INFO)
        main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "10x",
                "--output-dir",
                str(tmp_path),
                "--verbose",
            ]
        )
        captured = capsys.readouterr()
        lines = [line for line in captured.out.splitlines() if line.strip()]
        assert len(lines) == 1
        payload = json.loads(lines[0])
        assert payload["order"] == "czi-psomagen/orglab-alpha/AN00000001"
        # Progress was emitted, but through logging rather than onto stdout, so
        # it cannot corrupt the line a tracker parses. (pytest's logging plugin
        # intercepts the records, hence caplog rather than the captured stderr.)
        assert caplog.text
        assert "Listing s3://" in caplog.text

    def test_logging_is_wired_to_stderr_not_stdout(self):
        """The stream discipline above depends on this, so assert it directly."""
        from order_qa.cli import _configure_logging

        root = logging.getLogger()
        saved = root.handlers[:]
        for handler in saved:
            root.removeHandler(handler)
        try:
            _configure_logging(False)
            streams = [
                getattr(h, "stream", None)
                for h in logging.getLogger().handlers
                if isinstance(h, logging.StreamHandler)
            ]
            assert streams, "no stream handler installed"
            assert sys.stdout not in streams
            assert sys.stderr in streams
        finally:
            for handler in logging.getLogger().handlers[:]:
                logging.getLogger().removeHandler(handler)
            for handler in saved:
                root.addHandler(handler)

    def test_still_uploading_reports_its_own_code(self, stub_aws, tmp_path, capsys):
        stub_aws(minutes_old=1)
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "10x",
                "--output-dir",
                str(tmp_path),
            ]
        )
        assert code == exits.STILL_UPLOADING
        payload = json.loads(capsys.readouterr().out.strip())
        assert payload["verdict"] == "still_uploading"

    def test_force_overrides_the_quiescence_gate(self, stub_aws, tmp_path):
        stub_aws(minutes_old=1)
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "10x",
                "--force",
                "--output-dir",
                str(tmp_path),
            ]
        )
        assert code != exits.STILL_UPLOADING

    def test_quiescence_minutes_is_honoured(self, stub_aws, tmp_path):
        stub_aws(minutes_old=30)
        code = main(
            [
                "czi-psomagen/orglab-alpha/AN00000001",
                "--assay",
                "10x",
                "--quiescence-minutes",
                "60",
                "--output-dir",
                str(tmp_path),
            ]
        )
        assert code == exits.STILL_UPLOADING


class TestProbe:
    def test_probe_prints_capability_json_and_exits_ok(
        self, stub_aws, tmp_path, capsys
    ):
        stub_aws()
        code = main(["czi-psomagen/orglab-alpha/AN00000001", "--probe"])
        assert code == exits.OK
        payload = json.loads(capsys.readouterr().out)
        capabilities = payload["capabilities"]
        assert capabilities["head_object"]["status"] == "available"
        # The two that matter: unproven on these buckets, so they must be
        # reported rather than assumed either way.
        assert capabilities["get_object_attributes"]["status"] != "available"
        assert capabilities["list_object_versions"]["status"] != "available"

    def test_probe_needs_no_assay(self, stub_aws):
        stub_aws()
        assert main(["czi-psomagen/orglab-alpha/AN00000001", "--probe"]) == exits.OK

    def test_probe_on_empty_prefix_fails_verification(self, stub_aws):
        stub_aws(keys=())
        code = main(["czi-psomagen/orglab-alpha/AN00000001", "--probe"])
        assert code == exits.VERIFICATION_FAILED
