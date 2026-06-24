"""
Tests for multi-wafer-per-sublibrary trimmer gathering (Psomagen / Novogene).

Existing tests in test_qa_gather.py are unchanged; scenarios here cover
multiple RunIDs under one library folder and provider-specific merged stats.
"""

from __future__ import annotations

import csv
import os

import pytest

from qa_checks import build_wafer_failure_stats
from qa_gather import gather_qa_data
from qa_mods import QARunContext, resolve_wafer_run_id, trimmer_failure_storage_key

from tests.test_qa_gather import MockS3Client, _make_ctx

QA_FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "qa")
PSOM_LISTING = os.path.join(QA_FIXTURES_DIR, "psomagen_multi_wafer_s3_listing.tsv")

TRIMMER_CSV = (
    "code,format,segment,reason,failed read count,total read count\n"
    "8,trim,preamble,rsq file,100,1000\n"
    "101,trim,insert,sequence was too short,50,1000\n"
)


def _psom_ctx() -> QARunContext:
    return _make_ctx(
        bucket="czi-synthetic",
        provider="synthetic",
        proj="synthetic-project-alpha",
        order="SYN00000001",
        listing_prefix="synthetic-project-alpha/SYN00000001/",
        output_label="SYN00000001",
    )


def _listing_prefixes_from_fixture() -> list[str]:
    groups: set[str] = set()
    with open(PSOM_LISTING, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            top = row["File_Path"].split("/")[0]
            if top.endswith("_L01") or top.endswith("_L02"):
                groups.add(top)
    return [f"synthetic-project-alpha/SYN00000001/{g}/" for g in sorted(groups)]


def _trimmer_keys_for_groups(groups: list[str]) -> list[str]:
    group_names = {g.rstrip("/").split("/")[-1] for g in groups}
    keys: list[str] = []
    with open(PSOM_LISTING, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            path = row["File_Path"]
            parts = path.split("/")
            if len(parts) < 3 or parts[1] != "raw":
                continue
            if parts[0] not in group_names:
                continue
            if "trimmer-failure_codes" not in row["File_Name"]:
                continue
            keys.append(f"synthetic-project-alpha/SYN00000001/{path}")
    return keys


def _minimal_raw_stub(group: str, run_id: str, assay: str) -> str:
    return (
        f"synthetic-project-alpha/SYN00000001/{group}/raw/"
        f"{run_id}-{group}_{assay}-Z5001-CTCAGCATA.csv"
    )


class TestTrimmerFailureStorageKey:
    def test_storage_key_uses_run_id_when_present(self):
        key = (
            "proj/order/LIB01/raw/441001-LIB01_GEX-Z5001-BC01_trimmer-failure_codes.csv"
        )
        storage_key, run_id = trimmer_failure_storage_key(key)
        assert storage_key == "441001"
        assert run_id == "441001"

    def test_storage_key_falls_back_to_experiment(self):
        key = "proj/order/LIB01/raw/badname_trimmer-failure_codes.csv"
        storage_key, run_id = trimmer_failure_storage_key(key)
        assert storage_key == "order/LIB01"
        assert run_id is None

    def test_resolve_wafer_run_id_accepts_run_id_native_key(self):
        assert resolve_wafer_run_id("441001", {}) == "441001"

    def test_resolve_wafer_run_id_uses_legacy_exp_map(self):
        assert (
            resolve_wafer_run_id("order/LIB01", {"order/LIB01": "441001"}) == "441001"
        )


class TestGatherPsomagenMultiWafer:
    def test_trimmer_stats_keyed_by_run_id_not_sublibrary(self):
        """Three wafers per sublibrary must not collapse to one stats bucket."""
        groups = _listing_prefixes_from_fixture()
        trimmer_keys = _trimmer_keys_for_groups(groups)
        assert len(trimmer_keys) == 6  # 2 libraries x 3 trimmer files each

        raw_stubs = [
            _minimal_raw_stub("SYNLIB_L01", "441010", "GEX"),
            _minimal_raw_stub("SYNLIB_L02", "441013", "GEX"),
        ]
        order_prefix = "synthetic-project-alpha/SYN00000001/"
        keys = [order_prefix] + [g for g in groups] + trimmer_keys + raw_stubs
        file_contents = {k: TRIMMER_CSV for k in trimmer_keys}

        ctx = _psom_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        run_ids_in_stats = {
            k for k in data.trimmer_failure_stats if k.isdigit() and len(k) == 6
        }
        assert len(run_ids_in_stats) == 5  # 441049 shared across L01 + L02
        assert len(data.trimmer_failure_stats) == 5

        assert set(data.group_failure_stats) == {
            "SYN00000001/SYNLIB_L01",
            "SYN00000001/SYNLIB_L02",
        }
        for group_key in data.group_failure_stats:
            assert not (group_key.isdigit() and len(group_key) == 6)
            stats = data.group_failure_stats[group_key]
            assert len(stats["rsq"]) == 3
            assert len(stats["trimmer_fail"]) == 3

    def test_wafer_failure_stats_covers_all_run_ids_without_merged_files(self):
        """Psomagen-style orders may lack merged trimmer CSVs; wafer stats still populate."""
        groups = _listing_prefixes_from_fixture()
        trimmer_keys = _trimmer_keys_for_groups(groups)
        order_prefix = "synthetic-project-alpha/SYN00000001/"
        keys = [order_prefix] + list(groups) + trimmer_keys
        file_contents = {k: TRIMMER_CSV for k in trimmer_keys}

        ctx = _psom_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        assert data.merged_wafer_stats
        assert all(
            stats.get("rsq_total_reads") for stats in data.merged_wafer_stats.values()
        )
        wafer_stats = build_wafer_failure_stats(
            data.trimmer_failure_stats, data.exp_to_run_map
        )
        assert len(wafer_stats) == 5
        assert "441049" in wafer_stats
        assert len(wafer_stats["441049"]["rsq"]) == 2  # both sublibraries
        for run_id, stats in wafer_stats.items():
            assert stats["rsq"]
            assert stats["trimmer_fail"]


class TestGatherNovogeneMultiWaferWithMerged:
    def test_multi_wafer_per_library_plus_order_level_merged(self):
        """Novogene: per-wafer trimmer stats plus order-level merged CSVs."""
        g1_trimmers = [
            "testproj/ORD01/G1/raw/439047-G1_GEX-Z0001-BC01_trimmer-failure_codes.csv",
            "testproj/ORD01/G1/raw/439048-G1_CRI-Z0002-BC02_trimmer-failure_codes.csv",
        ]
        g2_trimmers = [
            "testproj/ORD01/G2/raw/439049-G2_GEX-Z0003-BC03_trimmer-failure_codes.csv",
        ]
        merged_keys = [
            "testproj/ORD01/439047_merged_trimmer-failure_codes.csv",
            "testproj/ORD01/439049_merged_trimmer-failure_codes.csv",
        ]
        merged_csv = (
            "read group,code,format,segment,reason,"
            "failed read count,total read count\n"
            "TT,8,trim,preamble,rsq file,1000,10000\n"
            "TT,101,trim,insert,sequence was too short,100,10000\n"
        )
        raw_stub = "testproj/ORD01/G1/raw/439047-G1_GEX-Z0001-BC01.csv"
        keys = (
            ["testproj/ORD01/", "testproj/ORD01/G1/", "testproj/ORD01/G2/"]
            + g1_trimmers
            + g2_trimmers
            + merged_keys
            + [raw_stub]
        )
        file_contents = {k: TRIMMER_CSV for k in g1_trimmers + g2_trimmers}
        file_contents.update({k: merged_csv for k in merged_keys})

        ctx = _make_ctx(raw_assay="10x")
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        assert set(data.trimmer_failure_stats) == {"439047", "439048", "439049"}
        assert set(data.merged_wafer_stats) == {"439047", "439048", "439049"}
        assert data.merged_wafer_stats["439047"]["tt_total_reads"] == 10000
        assert data.merged_wafer_stats["439049"]["tt_total_reads"] == 10000
        assert data.merged_wafer_stats["439048"]["rsq_total_reads"] == 1000
        assert "tt_total_reads" not in data.merged_wafer_stats["439048"]
        wafer_stats = build_wafer_failure_stats(
            data.trimmer_failure_stats, data.exp_to_run_map
        )
        assert set(wafer_stats) == {"439047", "439048", "439049"}

        assert set(data.group_failure_stats) == {"ORD01/G1", "ORD01/G2"}
        assert len(data.group_failure_stats["ORD01/G1"]["rsq"]) == 2
        assert len(data.group_failure_stats["ORD01/G1"]["trimmer_fail"]) == 2
        assert len(data.group_failure_stats["ORD01/G2"]["rsq"]) == 1
        assert len(data.group_failure_stats["ORD01/G2"]["trimmer_fail"]) == 1


class TestGatherPsomagenPerSampleTrimmerStats:
    def test_per_sample_files_populate_rsq_but_no_tt(self):
        """Psomagen (no merged files): RSQ populates; TT and Q30 stay blank."""
        flex_path = os.path.join(
            QA_FIXTURES_DIR, "psomagen_trimmer_stats_flex_v2_sample.csv"
        )
        failure_path = os.path.join(
            QA_FIXTURES_DIR, "psomagen_per_sample_trimmer_failure_codes.csv"
        )
        with open(flex_path) as fh:
            flex_csv = fh.read()
        with open(failure_path) as fh:
            failure_csv = fh.read()

        stats_key = (
            "synthetic-project-alpha/SYN00000001/SYNLIB_L10/raw/"
            "441010-SYNLIB_L10_GEX-Z5019-CACATCACA_trimmer-stats.csv"
        )
        failure_key = (
            "synthetic-project-alpha/SYN00000001/SYNLIB_L10/raw/"
            "441010-SYNLIB_L10_GEX-Z5019-CACATCACA_trimmer-failure_codes.csv"
        )
        raw_stub = (
            "synthetic-project-alpha/SYN00000001/SYNLIB_L10/raw/"
            "441010-SYNLIB_L10_GEX-Z5019-CACATCACA.csv"
        )
        keys = [
            "synthetic-project-alpha/SYN00000001/",
            "synthetic-project-alpha/SYN00000001/SYNLIB_L10/",
            stats_key,
            failure_key,
            raw_stub,
        ]
        file_contents = {
            stats_key: flex_csv,
            failure_key: failure_csv,
        }
        ctx = _psom_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        assert "441010" in data.merged_wafer_stats
        stats = data.merged_wafer_stats["441010"]
        # No TT: whole-sample volume must not be mislabeled as template reads.
        assert "tt_total_reads" not in stats
        # RSQ comes from the per-sample failure_codes rsq-file row.
        assert stats["rsq_total_reads"] == 10695752970
        assert stats["rsq_failed_reads"] == 3261649289
        assert stats["rsq_fail_pct"] == pytest.approx(
            100.0 * 3261649289 / 10695752970, rel=1e-9
        )
        # Per-sample trimmer-stats.csv is no longer parsed for Q30.
        assert "sample_q30_pct" not in stats

    def test_shared_run_id_aggregates_rsq_across_sublibraries(self):
        """A wafer reused across libraries sums its per-sample RSQ totals."""
        failure_path = os.path.join(
            QA_FIXTURES_DIR, "psomagen_per_sample_trimmer_failure_codes.csv"
        )
        with open(failure_path) as fh:
            failure_csv = fh.read()

        base = "synthetic-project-alpha/SYN00000001"
        f_l01 = f"{base}/SYNLIB_L01/raw/441049-SYNLIB_L01_CRI-Z5019-CACATCACA_trimmer-failure_codes.csv"
        f_l02 = f"{base}/SYNLIB_L02/raw/441049-SYNLIB_L02_CRI-Z5019-CACATCACA_trimmer-failure_codes.csv"
        keys = [
            f"{base}/",
            f"{base}/SYNLIB_L01/",
            f"{base}/SYNLIB_L02/",
            f_l01,
            f_l02,
        ]
        file_contents = {f_l01: failure_csv, f_l02: failure_csv}
        ctx = _psom_ctx()
        s3 = MockS3Client(keys=keys, file_contents=file_contents)
        data = gather_qa_data(ctx, s3)

        assert "441049" in data.merged_wafer_stats
        # Same wafer in two sublibraries → RSQ totals summed.
        assert data.merged_wafer_stats["441049"]["rsq_total_reads"] == 2 * 10695752970
        assert data.merged_wafer_stats["441049"]["rsq_failed_reads"] == 2 * 3261649289


class TestBuildWaferFailureStatsRunIdKeys:
    def test_run_id_native_keys_aggregate_shared_wafer_across_groups(self):
        """Same RunID in multiple stats buckets rolls up to one wafer row."""
        trimmer_failure_stats = {
            "441099": {"rsq": [1.0], "trimmer_fail": [5.0]},
            "441100": {"rsq": [2.0], "trimmer_fail": [6.0]},
        }
        wafer_stats = build_wafer_failure_stats(trimmer_failure_stats, {})
        assert set(wafer_stats) == {"441099", "441100"}
