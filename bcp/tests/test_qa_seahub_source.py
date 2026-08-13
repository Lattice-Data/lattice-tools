"""
Tests for cross-bucket SeaHub trimming reconciliation (qa_seahub_source).

The vendor layout mirrored here is the real one measured from
``czi-novogene/labalpha-seahub-bcp/NVUS0000000000-02``: the group folder sits
*before* ``raw``, giving a different path depth from the lab's trimmed upload.
"""

from __future__ import annotations

import json

import pytest

from qa_seahub_recon import reconcile_trimming
from qa_seahub_source import (
    WaferSeeds,
    _descend_to_raw_prefixes,
    _discover_untrimmed_sources,
    _normalize_search_roots,
    _wafer_seeds,
    derive_source_experiment,
    derive_source_order,
    index_trimmed_upload,
    index_untrimmed_sources,
    load_source_read_counts,
    normalize_source_uris,
    parse_source_uri,
    source_experiment_matches,
    source_order_by_wafer,
)

from tests.qa_seahub_helpers import (
    PROJECT,
    UNTRIMMED_WAFER,
    VENDOR_ORDER,
    VENDOR_ORDER_2,
    ref3_trimmed_keys,
    ref3_vendor_keys,
    vendor_uri,
)
from tests.test_qa_gather import MockS3Client

SOURCE_URI = "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-02"
SOURCE_DIR = "labalpha-seahub-bcp/NVUS0000000000-02/REF5_P01/raw/448642"
TRIMMED_DIR = "labalpha-seahub-bcp/REF5/raw/REF5_P01/448642"

VENDOR_STEM_A = "448642-REF5_P01_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
VENDOR_STEM_B = "448642-REF5_P01_GEX_hash_oligo-Z0002-AGCTCGAATGCGATC"
TRIM_STEM_A = "448642-REF5_P01_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
TRIM_STEM_B = "448642-REF5_P01_A2_GEX_hash_oligo-Z0002-AGCTCGAATGCGATC"


def _vendor_keys(*stems: str) -> list[str]:
    keys: list[str] = []
    for stem in stems:
        keys.extend(
            [
                f"{SOURCE_DIR}/{stem}.cram",
                f"{SOURCE_DIR}/{stem}.cram-metadata.json",
                f"{SOURCE_DIR}/{stem}.csv",
                f"{SOURCE_DIR}/{stem}.json",
                f"{SOURCE_DIR}/{stem}_FlowQ.metric",
                f"{SOURCE_DIR}/{stem}_trimmer-stats.csv",
            ]
        )
    keys.extend(
        [
            f"{SOURCE_DIR}/448642_LibraryInfo.xml",
            f"{SOURCE_DIR}/448642_UploadCompleted.json",
            f"{SOURCE_DIR}/448642_merged_trimmer-stats.csv",
            f"{SOURCE_DIR}/448642_run_SecondaryAnalysis.txt",
            f"{SOURCE_DIR}/448642-unmatched.cram",
            "labalpha-seahub-bcp/NVUS0000000000-02/file-manifest.json",
        ]
    )
    return keys


def _trimmed_keys(*stems: str, omit_cram_for: str = "") -> list[str]:
    keys: list[str] = []
    for stem in stems:
        suffixes = [".trim.csv", ".trim.stderr", ".trim.stdout", ".trim_fail.csv"]
        if stem != omit_cram_for:
            suffixes += [".trim.cram", ".trim.cram-metadata.json"]
        keys.extend(f"{TRIMMED_DIR}/{stem}{suffix}" for suffix in suffixes)
    return keys


def _source_index(keys: list[str]) -> dict:
    return index_untrimmed_sources(MockS3Client(keys=keys), [SOURCE_URI]).index


class TestParseSourceUri:
    def test_splits_bucket_and_prefix(self):
        assert parse_source_uri(SOURCE_URI) == (
            "czi-novogene",
            "labalpha-seahub-bcp/NVUS0000000000-02/",
        )

    def test_rejects_non_s3_uri(self):
        with pytest.raises(ValueError):
            parse_source_uri("/local/path")


class TestIndexUntrimmedSource:
    def test_indexes_one_entry_per_well(self):
        index = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        assert sorted(index) == [("448642", "Z0001"), ("448642", "Z0002")]

    def test_captures_identity_fields_and_sidecar(self):
        entry = _source_index(_vendor_keys(VENDOR_STEM_A))[("448642", "Z0001")]
        assert entry.barcode == "CAGCTCGAATGCGAT"
        assert entry.group == "REF5_P01"
        assert entry.assay == "GEX_hash_oligo"
        assert entry.metadata_key.endswith(".cram-metadata.json")

    def test_skips_wafer_level_and_unmatched_files(self):
        index = _source_index(_vendor_keys(VENDOR_STEM_A))
        assert len(index) == 1
        assert "unmatched" not in index[("448642", "Z0001")].cram_key

    def test_reads_counts_from_sidecars(self):
        keys = _vendor_keys(VENDOR_STEM_A)
        contents = {
            f"{SOURCE_DIR}/{VENDOR_STEM_A}.cram-metadata.json": json.dumps(
                {"read_count": 260527531}
            )
        }
        s3 = MockS3Client(keys=keys, file_contents=contents)
        index = index_untrimmed_sources(s3, [SOURCE_URI]).index
        load_source_read_counts(s3, index, "czi-novogene")
        assert index[("448642", "Z0001")].read_count == 260527531


class TestOneBadSidecarDoesNotKillTheRun:
    """A sidecar that cannot be read leaves one count empty, not the whole run.

    ``future.result()`` re-raises in the caller's thread, so a single missing,
    truncated or non-numeric sidecar out of several hundred took the notebook
    cell down with it. The reconciliation already handles a well with no vendor
    count -- it reports ``metadata_unavailable`` -- so an unreadable one degrades
    into a category the report understands.
    """

    def _run(self, bad_payload: str | None) -> dict:
        """Index two wells, give A a good sidecar and B ``bad_payload``."""
        keys = _vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B)
        contents = {
            f"{SOURCE_DIR}/{VENDOR_STEM_A}.cram-metadata.json": json.dumps(
                {"read_count": 260527531}
            )
        }
        if bad_payload is not None:
            contents[f"{SOURCE_DIR}/{VENDOR_STEM_B}.cram-metadata.json"] = bad_payload
        s3 = MockS3Client(keys=keys, file_contents=contents)
        index = _source_index(keys)
        load_source_read_counts(s3, index, "czi-novogene")
        return index

    @pytest.mark.parametrize(
        "payload,label",
        [
            (None, "object absent"),
            ("{not json at all", "malformed json"),
            ('{"read_count": "not a number"}', "non-numeric count"),
            ("{}", "no read_count key"),
        ],
    )
    def test_the_good_sidecar_is_still_read(self, payload, label):
        index = self._run(payload)

        assert index[("448642", "Z0001")].read_count == 260527531, label

    @pytest.mark.parametrize(
        "payload",
        [None, "{not json at all", '{"read_count": "not a number"}', "{}"],
    )
    def test_the_bad_one_is_left_as_no_count(self, payload):
        index = self._run(payload)

        assert index[("448642", "Z0002")].read_count is None

    def test_the_reconciliation_degrades_to_metadata_unavailable(self):
        """The category that already exists for a well with no vendor count."""
        index = self._run("{not json at all")
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A, TRIM_STEM_B))
        fail_counts = {
            (TRIMMED_DIR, stem): {"JumboSciGEX": {"failed": 1, "total": 260527531}}
            for stem in (TRIM_STEM_A, TRIM_STEM_B)
        }

        report = reconcile_trimming(index, trimmed, fail_counts)

        unavailable = [
            r for r in report.rows if r["category"] == "metadata_unavailable"
        ]
        assert [r["ug"] for r in unavailable] == ["Z0002"]


class TestIndexTrimmedUpload:
    def test_indexes_by_wafer_and_ug(self):
        index = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A, TRIM_STEM_B))
        assert sorted(index) == [("448642", "Z0001"), ("448642", "Z0002")]
        assert index[("448642", "Z0001")].sublibrary == "REF5_P01"
        assert index[("448642", "Z0001")].has_cram is True

    def test_flags_a_well_whose_cram_is_absent(self):
        index = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A, omit_cram_for=TRIM_STEM_A)
        )
        assert index[("448642", "Z0001")].has_cram is False

    def test_indexes_the_bare_family_too(self):
        bare = "448642-448642-REF5_P01_A1-Z0001-CAGCTCGAATGCGAT"
        keys = [f"{TRIMMED_DIR}/{bare}{s}" for s in (".cram", ".csv", "_fail.csv")]
        index = index_trimmed_upload(keys)
        entry = index[("448642", "Z0001")]
        assert entry.family == "bare"
        assert entry.assay is None
        assert entry.has_cram is True


class TestReconcileTrimming:
    def _fail_counts(self, stem: str, total: int) -> dict:
        # Keyed by (raw_dir, stem), as the gatherer writes it.
        return {(TRIMMED_DIR, stem): {"JumboSciGEX": {"failed": 100, "total": total}}}

    def test_consistent_pair_is_clean(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(TRIM_STEM_A, 260527531)
        )
        assert report.rows == []
        assert report.matched == 1

    def test_vendor_cram_with_no_trim_output(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed)
        not_trimmed = [r for r in report.rows if r["category"] == "not_trimmed"]
        assert [r["ug"] for r in not_trimmed] == ["Z0002"]

    def test_trimmed_well_whose_cram_is_absent_counts_as_not_trimmed(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A, omit_cram_for=TRIM_STEM_A)
        )
        report = reconcile_trimming(source, trimmed)
        assert [r["category"] for r in report.rows] == ["not_trimmed"]

    def test_orphan_trimmed_well(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A, TRIM_STEM_B))
        report = reconcile_trimming(source, trimmed)
        orphans = [r for r in report.rows if r["category"] == "orphan_trimmed"]
        assert [r["ug"] for r in orphans] == ["Z0002"]

    def test_identity_mismatch_on_barcode(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 100
        wrong_barcode = "448642-REF5_P01_A1_GEX_hash_oligo-Z0001-TTTTTTTTTTTTTTT"
        trimmed = index_trimmed_upload(_trimmed_keys(wrong_barcode))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(wrong_barcode, 100)
        )
        assert [r["category"] for r in report.rows] == ["identity_mismatch"]
        assert "barcode differs" in report.rows[0]["detail"]

    def test_read_count_mismatch_when_trimmer_saw_fewer_reads(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(
            source, trimmed, self._fail_counts(TRIM_STEM_A, 158233602)
        )
        assert [r["category"] for r in report.rows] == ["read_count_mismatch"]
        assert "260527531" in report.rows[0]["detail"]

    def test_any_matching_format_total_is_accepted(self):
        # Format blocks declare different totals for the same well, so matching
        # any one of them means the trimmer consumed the delivered file.
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        source[("448642", "Z0001")].read_count = 260527531
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        fail_counts = {
            (TRIMMED_DIR, TRIM_STEM_A): {
                "JumboSciGEX": {"failed": 1, "total": 158233602},
                "JumboSciHash": {"failed": 2, "total": 260527531},
            }
        }
        report = reconcile_trimming(source, trimmed, fail_counts)
        assert report.rows == []

    def test_missing_vendor_read_count_is_reported_as_unavailable(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed, self._fail_counts(TRIM_STEM_A, 1))
        assert [r["category"] for r in report.rows] == ["metadata_unavailable"]

    def test_verdict_summarizes_both_sides(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed)
        assert "2 vendor CRAMs" in report.verdict()
        assert "1 trimmed wells" in report.verdict()
        assert report.counts["not_trimmed"] == 1

    def test_verdict_appends_identity_mismatch_when_present(self):
        wrong_barcode = TRIM_STEM_A.replace("CAGCTCGAATGCGAT", "AAAAAAAAAAAAAAAA")
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(wrong_barcode))
        report = reconcile_trimming(source, trimmed)

        assert "1 identity mismatches" in report.verdict()

    def test_verdict_appends_metadata_unavailable_when_present(self):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))
        report = reconcile_trimming(source, trimmed, self._fail_counts(TRIM_STEM_A, 1))

        assert "1 metadata unavailable" in report.verdict()


class TestReconcileDoesNotMutateItsInput:
    """Reconciling twice must give the same answer twice.

    ``list(untrimmed.coverage)`` duplicates the list and not the SourceCoverage
    objects in it, so tallying ``entry.matched += 1`` wrote back into the
    caller's UntrimmedSources. Re-running the recon cell without re-running the
    source-index cell -- how a notebook is ordinarily used -- doubled ``matched``,
    and ``unmatched`` is ``max(0, indexed - matched)``, so a source with genuine
    unmatched wells reported 0 of them.
    """

    def _sources(self):
        return index_untrimmed_sources(
            MockS3Client(keys=_vendor_keys(VENDOR_STEM_A, VENDOR_STEM_B)),
            [SOURCE_URI],
        )

    def test_a_second_run_reports_the_same_coverage(self):
        sources = self._sources()
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))

        first = reconcile_trimming(sources, trimmed)
        # Snapshot before the second run: under the bug both reports alias the
        # same SourceCoverage objects, so comparing them afterwards compares two
        # views of one mutated value and passes either way.
        before = [(c.matched, c.unmatched) for c in first.source_coverage]

        second = reconcile_trimming(sources, trimmed)

        assert [(c.matched, c.unmatched) for c in second.source_coverage] == before

    def test_a_second_run_still_sees_the_unmatched_well(self):
        """The consequence: unmatched floors at 0 once matched is double-counted."""
        sources = self._sources()
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))

        reconcile_trimming(sources, trimmed)
        second = reconcile_trimming(sources, trimmed)

        assert [c.matched for c in second.source_coverage] == [1]
        assert [c.unmatched for c in second.source_coverage] == [1]

    def test_the_caller_object_is_left_alone(self):
        sources = self._sources()
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))

        reconcile_trimming(sources, trimmed)

        assert [c.matched for c in sources.coverage] == [0]

    def test_the_report_owns_its_coverage_rows(self):
        sources = self._sources()
        trimmed = index_trimmed_upload(_trimmed_keys(TRIM_STEM_A))

        report = reconcile_trimming(sources, trimmed)

        assert all(
            row is not original
            for row, original in zip(report.source_coverage, sources.coverage)
        )


# The vendor layout as measured on a real delivery: the segment before ``raw`` is
# the ExperimentID alone and carries no sublibrary.  SOURCE_DIR above predates
# that observation and uses ``{ExperimentID}_{sublibrary}``; both shapes must
# work, since order derivation is positional.  Order ids are sanitized.
ORDER_11 = "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-11"
ORDER_12 = "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-12"


def _vendor_layout_key(order: str, wafer: str, stem: str, suffix: str = ".cram") -> str:
    return f"labalpha-seahub-bcp/{order}/REF3/raw/{wafer}/{stem}{suffix}"


P04_STEM = "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
P07_STEM = "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT"


class TestNormalizeSourceUris:
    def test_a_bare_string_is_wrapped(self):
        assert normalize_source_uris(ORDER_11) == [ORDER_11]

    def test_blanks_are_dropped_and_repeats_collapse(self):
        assert normalize_source_uris([ORDER_11, "  ", ORDER_11, ORDER_12]) == [
            ORDER_11,
            ORDER_12,
        ]

    def test_none_and_empty_list_are_empty(self):
        assert normalize_source_uris(None) == []
        assert normalize_source_uris([]) == []

    def test_a_bucket_only_uri_is_rejected(self):
        """It would paginate an entire bucket."""
        with pytest.raises(ValueError, match="too broad"):
            normalize_source_uris(["s3://czi-novogene"])

    def test_a_project_only_uri_is_rejected(self):
        with pytest.raises(ValueError, match="too broad"):
            normalize_source_uris(["s3://czi-novogene/labalpha-seahub-bcp"])

    def test_a_non_s3_entry_names_itself(self):
        with pytest.raises(ValueError, match="/local/path"):
            normalize_source_uris([ORDER_11, "/local/path"])

    def test_a_comma_is_never_split(self):
        """A prefix may legitimately contain a comma."""
        uri = "s3://czi-novogene/proj,x/NVUS0000000000-02"
        assert normalize_source_uris(uri) == [uri]


class TestDeriveSourceOrder:
    def test_reads_the_order_from_the_measured_layout(self):
        key = _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM)
        assert derive_source_order(key) == "NVUS0000000000-11"

    def test_reads_the_order_from_the_fixture_layout(self):
        assert derive_source_order(f"{SOURCE_DIR}/{VENDOR_STEM_A}.cram") == (
            "NVUS0000000000-02"
        )

    def test_an_unexpected_depth_degrades_rather_than_guessing(self):
        assert derive_source_order("a/b/c/raw/d/e/f.cram") == ""
        assert derive_source_order("f.cram") == ""


class TestIndexUntrimmedSources:
    def _s3(self, *specs):
        keys = [_vendor_layout_key(o, w, s) for o, w, s in specs]
        keys += [
            _vendor_layout_key(o, w, s, ".cram-metadata.json") for o, w, s in specs
        ]
        return MockS3Client(keys=keys)

    def test_several_prefixes_merge_into_one_index(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        assert sorted(result.index) == [("437120", "Z0001"), ("438514", "Z0305")]
        assert len(result) == 2

    def test_a_bare_string_is_accepted(self):
        s3 = self._s3(("NVUS0000000000-11", "437120", P04_STEM))
        assert len(index_untrimmed_sources(s3, ORDER_11)) == 1

    def test_each_entry_carries_its_provenance(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        first = result.index[("437120", "Z0001")]
        assert first.source_uri == ORDER_11
        assert first.source_order == "NVUS0000000000-11"
        assert first.bucket == "czi-novogene"
        assert result.index[("438514", "Z0305")].source_order == "NVUS0000000000-12"
        assert result.buckets == ("czi-novogene",)

    def test_sidecars_pair_with_their_own_cram(self):
        s3 = self._s3(("NVUS0000000000-11", "437120", P04_STEM))
        entry = index_untrimmed_sources(s3, [ORDER_11]).index[("437120", "Z0001")]

        assert (
            entry.metadata_key
            == f"{entry.cram_key[: -len('.cram')]}.cram-metadata.json"
        )

    def test_a_well_delivered_twice_keeps_the_newest_and_reports_the_other(self):
        """A re-delivery must not silently overwrite, and must not lose either.

        The previous single-prefix indexer assigned unconditionally, so one of
        the two orders vanished with no record. The order-taking rule that
        replaced it then took the lexicographic *minimum*, which keeps the
        original and discards the re-delivery -- backwards, since a well is
        mostly re-delivered because the first copy was wrong.
        """
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "437120", P04_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        assert len(result.index) == 1
        assert result.index[("437120", "Z0001")].source_order == "NVUS0000000000-12"
        duplicates = [
            f for f in result.findings if f["category"] == "duplicate_source_well"
        ]
        assert len(duplicates) == 1
        assert "NVUS0000000000-11" in duplicates[0]["detail"]

    def test_an_unplaceable_key_loses_to_a_well_formed_one(self):
        """An empty order sorts before every real one, so it used to win outright.

        The key here has no ``raw`` segment at the expected depth, so
        ``derive_source_order`` reads no order from it -- and taking the minimum
        then preferred the one copy whose provenance QA cannot state.
        """
        good = _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM)
        odd = f"labalpha-seahub-bcp/NVUS0000000000-11/REF3/{P04_STEM}.cram"
        result = index_untrimmed_sources(MockS3Client(keys=[good, odd]), [ORDER_11])

        entry = result.index[("437120", "Z0001")]
        assert entry.cram_key == good
        assert entry.source_order == "NVUS0000000000-11"
        duplicates = [
            f for f in result.findings if f["category"] == "duplicate_source_well"
        ]
        assert [f["source_key"] for f in duplicates] == [odd]

    def test_the_winner_does_not_depend_on_the_caller_list_order(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "437120", P04_STEM),
        )
        forward = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])
        reverse = index_untrimmed_sources(s3, [ORDER_12, ORDER_11])

        assert (
            forward.index[("437120", "Z0001")].cram_key
            == reverse.index[("437120", "Z0001")].cram_key
        )

    def test_a_prefix_inside_another_is_listed_once(self):
        s3 = self._s3(("NVUS0000000000-11", "437120", P04_STEM))
        result = index_untrimmed_sources(s3, [ORDER_11, f"{ORDER_11}/REF3/raw"])

        assert len(result.index) == 1
        skipped = [c for c in result.coverage if c.skipped_reason]
        assert len(skipped) == 1
        assert any(
            f["category"] == "overlapping_source_prefix" for f in result.findings
        )

    def test_coverage_reports_what_each_prefix_contributed(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        by_uri = {c.source_uri: c for c in result.coverage}
        assert by_uri[ORDER_11].indexed == 1
        assert by_uri[ORDER_11].cram_keys == 1
        assert by_uri[ORDER_11].source_order == "NVUS0000000000-11"
        assert by_uri[ORDER_12].indexed == 1

    def test_an_empty_uri_list_is_not_an_error(self):
        """The notebook default must not blow up."""
        result = index_untrimmed_sources(MockS3Client(keys=[]), [])

        assert len(result) == 0
        assert result.coverage == []
        assert result.findings == []

    def test_vendor_housekeeping_files_are_still_skipped(self):
        index = _source_index(_vendor_keys(VENDOR_STEM_A))
        assert len(index) == 1


class TestDeriveSourceExperiment:
    def test_reads_the_experiment_from_the_measured_layout(self):
        key = _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM)
        assert derive_source_experiment(key) == "REF3"

    def test_reads_the_older_experiment_sublibrary_shape(self):
        assert derive_source_experiment(f"{SOURCE_DIR}/{VENDOR_STEM_A}.cram") == (
            "REF5_P01"
        )

    def test_an_unexpected_depth_degrades_rather_than_guessing(self):
        assert derive_source_experiment("a/b/c/raw/d/e/f.cram") == ""
        assert derive_source_experiment("f.cram") == ""


class TestSourceExperimentMatches:
    def test_the_experiment_id_alone_matches(self):
        assert source_experiment_matches("REF3", "REF3")

    def test_the_older_experiment_sublibrary_shape_matches(self):
        """A bare equality test here would orphan every well of an old delivery."""
        assert source_experiment_matches("REF5_P01", "REF5")

    def test_a_suffixed_reupload_folder_matches(self):
        """Measured: one order held REF3 beside a GENE7_reupload folder.

        The suffix is not always a sublibrary, so the rule cannot assume that
        shape -- it is any ``{ExperimentID}_...`` folder of the same experiment.
        """
        assert source_experiment_matches("GENE7_reupload", "GENE7")
        assert not source_experiment_matches("GENE7_reupload", "REF3")

    def test_a_longer_experiment_id_is_not_a_match(self):
        assert not source_experiment_matches("REF50", "REF5")
        assert not source_experiment_matches("REF50_P01", "REF5")

    def test_another_experiment_is_not_a_match(self):
        assert not source_experiment_matches("REF7", "REF3")


class TestExperimentScoping:
    """One order holds several experiments, so an order prefix spans them."""

    def _s3(self, *specs):
        keys = [
            f"labalpha-seahub-bcp/{order}/{experiment}/raw/{wafer}/{stem}{suffix}"
            for order, experiment, wafer, stem in specs
            for suffix in (".cram", ".cram-metadata.json")
        ]
        return MockS3Client(keys=keys)

    def test_another_experiments_wells_are_excluded(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "REF3", "437120", P04_STEM),
            ("NVUS0000000000-11", "REF7", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11], experiment_id="REF3")

        assert sorted(result.index) == [("437120", "Z0001")]

    def test_the_exclusion_is_reported_once_per_foreign_experiment(self):
        s3 = self._s3(
            ("NVUS0000000000-11", "REF3", "437120", P04_STEM),
            ("NVUS0000000000-11", "REF7", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11], experiment_id="REF3")

        spans = [
            r
            for r in result.findings
            if r["category"] == "source_prefix_spans_experiments"
        ]
        assert len(spans) == 1
        assert "REF7" in spans[0]["detail"]

    def test_the_older_experiment_sublibrary_shape_is_kept(self):
        s3 = self._s3(("NVUS0000000000-11", "REF3_P04_1", "437120", P04_STEM))
        result = index_untrimmed_sources(s3, [ORDER_11], experiment_id="REF3")

        assert sorted(result.index) == [("437120", "Z0001")]
        assert result.findings == []

    def test_a_key_whose_shape_hides_the_experiment_is_kept(self):
        """Dropping it would report a delivered well as never delivered."""
        s3 = MockS3Client(
            keys=[f"labalpha-seahub-bcp/NVUS0000000000-11/{P04_STEM}.cram"]
        )
        result = index_untrimmed_sources(s3, [ORDER_11], experiment_id="REF3")

        assert sorted(result.index) == [("437120", "Z0001")]
        assert [r["category"] for r in result.findings] == [
            "source_experiment_unreadable"
        ]

    def test_no_experiment_id_indexes_everything(self):
        """The default has to stay backward compatible."""
        s3 = self._s3(
            ("NVUS0000000000-11", "REF3", "437120", P04_STEM),
            ("NVUS0000000000-11", "REF7", "438514", P07_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11])

        assert sorted(result.index) == [("437120", "Z0001"), ("438514", "Z0305")]
        assert result.findings == []

    def test_an_excluded_well_does_not_become_a_not_trimmed_row(self):
        """The whole point: a foreign well must not read as a completeness gap."""
        s3 = self._s3(
            ("NVUS0000000000-11", "REF3", "437120", P04_STEM),
            ("NVUS0000000000-11", "REF7", "438514", P07_STEM),
        )
        sources = index_untrimmed_sources(s3, [ORDER_11], experiment_id="REF3")
        trimmed = index_trimmed_upload(
            [
                f"labalpha-seahub-bcp/REF3/raw/REF3_P04_1/437120/{P04_STEM}{suffix}"
                for suffix in (".trim.cram", ".trim.csv")
            ]
        )
        report = reconcile_trimming(sources, trimmed)

        assert "not_trimmed" not in {r["category"] for r in report.rows}


class TestSourceOrderByWafer:
    def test_one_order_per_wafer(self):
        s3 = MockS3Client(
            keys=[
                _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM),
                _vendor_layout_key("NVUS0000000000-12", "438514", P07_STEM),
            ]
        )
        index = index_untrimmed_sources(s3, [ORDER_11, ORDER_12]).index

        assert source_order_by_wafer(index) == {
            "437120": "NVUS0000000000-11",
            "438514": "NVUS0000000000-12",
        }

    def test_a_wafer_split_across_orders_shows_both(self):
        s3 = MockS3Client(
            keys=[
                _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM),
                _vendor_layout_key(
                    "NVUS0000000000-12",
                    "437120",
                    "437120-REF3_P04_1_A2_GEX_hash_oligo-Z0002-AGCTCGAATGCGATC",
                ),
            ]
        )
        index = index_untrimmed_sources(s3, [ORDER_11, ORDER_12]).index

        assert source_order_by_wafer(index) == {
            "437120": "NVUS0000000000-11|NVUS0000000000-12"
        }


class TestSizeAndCoverageReconciliation:
    def _pair(self, source_size: int, trimmed_size: int):
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        entry = source[("448642", "Z0001")]
        entry.read_count = 100
        entry.size_bytes = source_size
        trimmed = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A),
            sizes={f"{TRIMMED_DIR}/{TRIM_STEM_A}.trim.cram": trimmed_size},
        )
        fail_counts = {
            (TRIMMED_DIR, TRIM_STEM_A): {"JumboSciGEX": {"failed": 1, "total": 100}}
        }
        return reconcile_trimming(source, trimmed, fail_counts)

    def test_a_smaller_trimmed_cram_is_clean(self):
        """The real ratio: 14.4 GB of trimmed output from a 35-45 GB source."""
        assert self._pair(35_760_085_711, 14_400_000_000).rows == []

    def test_a_trimmed_cram_at_or_above_its_source_is_reported(self):
        report = self._pair(14_400_000_000, 35_760_085_711)
        assert [r["category"] for r in report.rows] == ["size_not_reduced"]
        assert "untrimmed data" in report.rows[0]["detail"]

    def test_an_equal_size_is_reported(self):
        report = self._pair(1_000, 1_000)
        assert [r["category"] for r in report.rows] == ["size_not_reduced"]

    def test_an_unknown_size_stays_silent(self):
        """Sizes are only collected in S3 mode; 0 means unknown, not empty."""
        assert self._pair(35_760_085_711, 0).rows == []
        assert self._pair(0, 14_400_000_000).rows == []

    def test_a_wafer_with_no_listed_source_is_named_in_the_verdict(self):
        """The signature of a forgotten vendor order.

        REF3 makes this real: order NVUS0000000000-11 holds six of its seven
        sublibraries, so listing it alone leaves every REF3_P05_1 well
        unsourced. Those wells surface as orphan_trimmed, which alone reads as a
        completeness failure rather than as incomplete input.
        """
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        unlisted = "999999-REF5_P09_A1_GEX_hash_oligo-Z0009-CAGCTCGAATGCGAT"
        trimmed = index_trimmed_upload(
            _trimmed_keys(TRIM_STEM_A)
            + [f"labalpha-seahub-bcp/REF5/raw/REF5_P09/999999/{unlisted}.trim.cram"]
        )

        report = reconcile_trimming(source, trimmed)

        assert report.unsourced_wafers == ["999999"]
        assert "an untrimmed order is probably missing" in report.verdict()

    def test_coverage_carries_the_matched_count(self):
        s3 = MockS3Client(
            keys=[
                _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM),
                _vendor_layout_key(
                    "NVUS0000000000-11", "437120", P04_STEM, ".cram-metadata.json"
                ),
            ]
        )
        sources = index_untrimmed_sources(s3, [ORDER_11])
        trimmed = index_trimmed_upload(
            [f"labalpha-seahub-bcp/REF3/raw/REF3_P04_1/437120/{P04_STEM}.trim.cram"]
        )

        report = reconcile_trimming(sources, trimmed)

        assert report.source_coverage[0].matched == 1
        assert report.source_coverage[0].unmatched == 0

    def test_indexer_findings_are_seeded_into_the_report(self):
        s3 = MockS3Client(
            keys=[
                _vendor_layout_key("NVUS0000000000-11", "437120", P04_STEM),
                _vendor_layout_key("NVUS0000000000-12", "437120", P04_STEM),
            ]
        )
        sources = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        report = reconcile_trimming(sources, {})

        assert any(r["category"] == "duplicate_source_well" for r in report.rows)
        assert "duplicate source wells" in report.verdict()

    def test_a_bare_index_still_works(self):
        """Backwards compatibility with the single-prefix call shape."""
        source = _source_index(_vendor_keys(VENDOR_STEM_A))
        report = reconcile_trimming(source, {})

        assert report.source_total == 1
        assert report.source_coverage == []


class TestDuplicateTrimmedWell:
    def test_one_well_under_two_names_is_reported(self):
        keys = [
            f"labalpha-seahub-bcp/REF5/raw/REF5_P01/448642/{TRIM_STEM_A}.trim.cram",
            "labalpha-seahub-bcp/REF5/raw/P01/448642/"
            "448642-448642-REF5_P01_A1-Z0001-CAGCTCGAATGCGAT.cram",
        ]
        findings: list[dict] = []

        index = index_trimmed_upload(keys, findings=findings)

        assert len(index) == 1
        assert [f["category"] for f in findings] == ["duplicate_trimmed_well"]

    def test_the_artifacts_of_one_well_are_not_duplicates(self):
        """Identity reuse is normal -- five artifacts share one (wafer, UG)."""
        findings: list[dict] = []

        index_trimmed_upload(_trimmed_keys(TRIM_STEM_A), findings=findings)

        assert findings == []

    def test_duplicate_is_reported_once_per_alternate_name_not_per_artifact(self):
        """Each artifact of a duplicated well reused the same finding -- six identical rows."""
        alt_stem = "448642-448642-REF5_P01_A1-Z0001-CAGCTCGAATGCGAT"
        alt_dir = "labalpha-seahub-bcp/REF5/raw/P01/448642"
        suffixes = [
            ".trim.csv",
            ".trim.stderr",
            ".trim.stdout",
            ".trim_fail.csv",
            ".trim.cram",
            ".trim.cram-metadata.json",
        ]
        keys = _trimmed_keys(TRIM_STEM_A) + [
            f"{alt_dir}/{alt_stem}{suffix}" for suffix in suffixes
        ]
        findings: list[dict] = []

        index_trimmed_upload(keys, findings=findings)

        assert len(findings) == 1
        assert findings[0]["category"] == "duplicate_trimmed_well"


class TestAWholeUntrimmedWaferIsLoud:
    """A delivered wafer with nothing trimmed must reach both reporting paths.

    Pinned on the current order-prefix behaviour, ahead of the change that
    locates vendor prefixes from the trimmed upload's wafer list instead. That
    change makes ``wafer in trimmed_wafers`` true of every indexed identity by
    construction, so a wafer the lab never touched would not be listed, not be
    indexed, and produce no rows at all -- and the two paths it drives are the
    only ones that would say so. ``not_trimmed`` here, and the vendor-only
    ``DATA_GAP`` wells in test_qa_seahub_rename, which is what errors.txt
    writes; between them, a forgotten plate reading as a clean upload.

    Within-wafer completeness is a different question and is covered above:
    wafer 438514 is trimmed with its CRAM absent, and still reports.
    """

    def _index(self):
        s3 = MockS3Client(keys=ref3_vendor_keys(extra_wafer=True))
        return index_untrimmed_sources(
            s3, [vendor_uri(VENDOR_ORDER), vendor_uri(VENDOR_ORDER_2)]
        ).index

    def test_the_wafer_is_indexed_from_an_order_prefix(self):
        """The premise of the rest: nothing about the upload gates the listing."""
        index = self._index()

        assert [ug for wafer, ug in sorted(index) if wafer == UNTRIMMED_WAFER] == [
            "Z0500",
            "Z0501",
        ]

    def test_every_well_of_it_reports_not_trimmed(self):
        report = reconcile_trimming(
            self._index(), index_trimmed_upload(ref3_trimmed_keys())
        )

        not_trimmed = [
            r
            for r in report.rows
            if r["category"] == "not_trimmed" and r["wafer"] == UNTRIMMED_WAFER
        ]
        assert [r["ug"] for r in not_trimmed] == ["Z0500", "Z0501"]

    def test_the_verdict_counts_it(self):
        """The printed line is what an operator reads before opening any CSV."""
        report = reconcile_trimming(
            self._index(), index_trimmed_upload(ref3_trimmed_keys())
        )

        # Wafer 438514's absent CRAM is the third.
        assert "3 not trimmed" in report.verdict()

    def test_it_is_not_reported_as_an_unsourced_wafer(self):
        """unsourced_wafers is the trimmed-side signal, so it cannot cover this.

        It names wafers in the *upload* that no source delivered -- the opposite
        direction -- which is why this wafer needs a pin of its own.
        """
        report = reconcile_trimming(
            self._index(), index_trimmed_upload(ref3_trimmed_keys())
        )

        assert UNTRIMMED_WAFER not in report.unsourced_wafers


class TestWaferSeeds:
    """The vendor search terms, derived from the upload and nothing else.

    Each of the three readings is here because it rescues a case the others
    miss; the tests are named for the case rather than the reading, since the
    union is the contract and which member supplies a wafer is not.
    """

    WAFER_DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/436830"
    STEM = "436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT"

    def _seeds(self, *keys: str, discovered: set[str] | None = None):
        keys_list = list(keys)
        return _wafer_seeds(
            trimmed_keys=keys_list,
            trimmed_index=index_trimmed_upload(keys_list),
            discovered_wafers=discovered,
        )

    def test_the_ordinary_case_yields_one_wafer(self):
        seeds = self._seeds(f"{self.WAFER_DIR}/{self.STEM}.trim.cram")

        assert seeds.wafers == ("436830",)
        assert seeds.rejected == ()

    def test_a_wafer_whose_filenames_do_not_parse_is_still_found(self):
        """The measured ScaleBio delivery: 192 CRAMs, no Z#### UG anywhere.

        The identity index is empty for these, so the folder reading is the only
        one that sees the wafer -- and without it the vendor side would never be
        searched for a wafer that is plainly present.
        """
        key = "labalpha-seahub-bcp/REF3/raw/RNA3_098/426971/426971-RNA3-098C_GEX_QSR-7_10A.cram"

        assert index_trimmed_upload([key]) == {}
        assert self._seeds(key).wafers == ("426971",)

    def test_a_well_whose_folder_is_not_a_wafer_is_still_found(self):
        """Only the filename reading sees this one, and the folder is refused."""
        key = f"labalpha-seahub-bcp/REF3/raw/REF3_P05_1/oops/{self.STEM}.trim.cram"
        seeds = self._seeds(key)

        assert seeds.wafers == ("436830",)
        assert seeds.rejected == ("oops",)

    def test_a_folder_disagreeing_with_its_filename_seeds_both(self):
        """wafer_mismatch is an SOP rule because this happens.

        Which token is right is not decidable here, so both are searched: the
        wrong one costs one lookup and reports as not found, whereas picking
        wrong silently compares against the wrong delivery.
        """
        key = f"labalpha-seahub-bcp/REF3/raw/REF3_P05_1/999999/{self.STEM}.trim.cram"

        assert self._seeds(key).wafers == ("436830", "999999")

    def test_junk_in_a_wafer_folder_still_seeds_the_wafer(self):
        assert self._seeds(f"{self.WAFER_DIR}/login.html").wafers == ("436830",)

    def test_an_object_off_wafer_depth_seeds_nothing(self):
        """Neither reading has a segment to read, and that is reported elsewhere.

        Documented rather than worked around: bad_path_depth already names these
        objects, and inventing a wafer for one would be a guess.
        """
        shallow = f"labalpha-seahub-bcp/REF3/raw/{self.STEM}.trim.cram"
        deep = f"{self.WAFER_DIR}/extra/{self.STEM}.trim.cram"

        assert self._seeds(shallow, deep) == WaferSeeds()

    def test_a_sublibrary_name_is_never_searched_for(self):
        """discovered_wafers is unvalidated at source, so it is filtered here."""
        seeds = self._seeds(discovered={"REF3_P04_1", "437120_old", "437120"})

        assert seeds.wafers == ("437120",)
        assert seeds.rejected == ("437120_old", "REF3_P04_1")

    def test_manifest_mode_still_yields_the_uploaded_wafers(self):
        """discovered_wafers is populated by the folder walk, so s3 mode only."""
        key = f"{self.WAFER_DIR}/{self.STEM}.trim.cram"

        assert self._seeds(key, discovered=set()).wafers == ("436830",)

    def test_the_whole_fixture_yields_its_five_wafers_once_each(self):
        keys = ref3_trimmed_keys()
        seeds = _wafer_seeds(
            trimmed_keys=keys, trimmed_index=index_trimmed_upload(keys)
        )

        assert seeds.wafers == ("436830", "436831", "437120", "438514", "439000")
        assert UNTRIMMED_WAFER not in seeds.wafers

    def test_no_input_at_all_is_empty_rather_than_an_error(self):
        assert _wafer_seeds() == WaferSeeds()
        assert len(WaferSeeds()) == 0


class TestASidecarStaysInItsOwnBucket:
    """A CRAM must not be credited with a sidecar from a different bucket.

    ``sidecar_by_stem_key`` outlives the per-prefix loop, so keying it on the
    stem alone let a CRAM in one bucket resolve a sidecar that exists at the
    same path in another. The entry then claims metadata it does not have, and
    the fetch that follows misses and degrades to ``metadata_unavailable`` --
    indistinguishable from a vendor sidecar with no ``read_count`` in it.

    Driven through ``index_untrimmed_sources`` -- the observable effect is
    ``metadata_key``, which the indexer assigns and which decides whether
    ``load_source_read_counts`` has anything to fetch.
    """

    ORDER = "labalpha-seahub-bcp/NVUS0000000000-11"
    PATH = f"{ORDER}/REF3/raw/437120/{P04_STEM}"

    def _index(self, per_bucket: dict[str, list[str]]):
        """List the same order prefix in every named bucket, sidecar bucket first."""
        s3 = MockS3Client(buckets=per_bucket)
        return index_untrimmed_sources(
            s3, [f"s3://{bucket}/{self.ORDER}" for bucket in per_bucket]
        )

    def test_a_sidecar_in_another_bucket_is_not_claimed(self):
        result = self._index(
            {
                "czi-other": [f"{self.PATH}.cram-metadata.json"],
                "czi-novogene": [f"{self.PATH}.cram"],
            }
        )
        entry = result.index[("437120", "Z0001")]

        assert entry.bucket == "czi-novogene"
        assert entry.metadata_key == ""

    def test_a_sidecar_beside_its_own_cram_still_resolves(self):
        """The fix must not stop the ordinary same-bucket case working."""
        result = self._index(
            {"czi-novogene": [f"{self.PATH}.cram", f"{self.PATH}.cram-metadata.json"]}
        )
        entry = result.index[("437120", "Z0001")]

        assert entry.metadata_key == f"{self.PATH}.cram-metadata.json"

    def test_the_well_then_reconciles_as_metadata_unavailable_not_as_a_count(self):
        """What the operator actually sees, and why claiming the sidecar is wrong.

        With no metadata_key there is nothing to fetch, so the well has no vendor
        read_count and says so. Claiming a sidecar from another bucket reached the
        same category by a worse route -- a fetch that could not succeed.
        """
        result = self._index(
            {
                "czi-other": [f"{self.PATH}.cram-metadata.json"],
                "czi-novogene": [f"{self.PATH}.cram"],
            }
        )

        assert result.index[("437120", "Z0001")].read_count is None


class TestNormalizeSearchRoots:
    """A search root is descended with a delimiter, so a bucket root is legal.

    The paired assertion lives in TestNormalizeSourceUris and must stay red on
    the same input: the listing validator still refuses what this one accepts,
    because the two feed different kinds of S3 call.
    """

    BUCKET_ROOT = "s3://czi-novogene"
    PROJECT_ROOT = "s3://czi-novogene/labalpha-seahub-bcp"

    def test_a_bucket_root_is_accepted(self):
        assert _normalize_search_roots(self.BUCKET_ROOT) == ["s3://czi-novogene/"]

    def test_the_listing_validator_still_rejects_the_same_root(self):
        """The guard is repurposed, not relaxed -- pinned here as well as there."""
        with pytest.raises(ValueError, match="too broad"):
            normalize_source_uris(f"{self.BUCKET_ROOT}/")

    def test_a_trailing_slash_is_never_the_difference_between_two_roots(self):
        assert _normalize_search_roots([self.BUCKET_ROOT, f"{self.BUCKET_ROOT}/"]) == [
            "s3://czi-novogene/"
        ]

    def test_a_project_root_is_accepted_and_is_the_recommended_shape(self):
        assert _normalize_search_roots(self.PROJECT_ROOT) == [f"{self.PROJECT_ROOT}/"]

    def test_a_bare_string_is_wrapped(self):
        assert len(_normalize_search_roots(self.PROJECT_ROOT)) == 1

    def test_none_and_empty_entries_are_dropped(self):
        assert _normalize_search_roots(None) == []
        assert _normalize_search_roots(["", "  "]) == []

    def test_a_non_s3_root_raises(self):
        with pytest.raises(ValueError):
            _normalize_search_roots("/local/path")

    def test_a_comma_is_never_split(self):
        """An S3 prefix may contain a comma; splitting one lists the wrong thing."""
        odd = "s3://czi-novogene/project,with-comma"
        assert _normalize_search_roots(odd) == [f"{odd}/"]

    def test_a_root_inside_another_root_is_dropped(self):
        assert _normalize_search_roots([self.PROJECT_ROOT, self.BUCKET_ROOT]) == [
            "s3://czi-novogene/"
        ]

    def test_the_same_prefix_in_another_bucket_is_kept(self):
        roots = _normalize_search_roots(
            [self.PROJECT_ROOT, "s3://czi-psomagen/labalpha-seahub-bcp"]
        )

        assert len(roots) == 2

    def test_a_sibling_project_is_not_swallowed_by_a_prefix_match(self):
        """``.../REF3/`` must not read as covering ``.../REF3_P05_1/``."""
        roots = _normalize_search_roots(
            ["s3://czi-novogene/REF3", "s3://czi-novogene/REF3_P05_1"]
        )

        assert len(roots) == 2

    @pytest.mark.parametrize(
        "root",
        [
            "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw",
            "s3://czi-novogene/labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/437120",
        ],
    )
    def test_a_root_inside_a_raw_folder_is_rejected(self, root):
        """That is an untrimmed_s3_paths entry, not a root to search from.

        The descent stops at a child folder named raw, so starting inside one
        locates nothing -- which would read as "this bucket holds no vendor data"
        rather than "you gave me the wrong kind of path".
        """
        with pytest.raises(ValueError, match="inside a 'raw' folder"):
            _normalize_search_roots(root)

    def test_a_prefix_merely_containing_the_word_raw_is_fine(self):
        """Segment-wise, not substring: a project may legitimately be named this."""
        root = "s3://czi-novogene/rawlings-seahub-bcp"
        assert _normalize_search_roots(root) == [f"{root}/"]


class TestDescendToRawPrefixes:
    """Locating vendor deliveries by shape, without being told where they are.

    A root already inside a ``raw`` folder is not tested here: step 4's
    ``_normalize_search_roots`` refuses one at the door, so the descent cannot
    receive it.
    """

    FOLDERED = f"labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/437120/{P04_STEM}.cram"
    # Real inside czi-novogene: no wafer directory, wafer only in the filename.
    FLAT = (
        "project-killifish/NVUS0000000000-20/CH13/raw/"
        "438586-CH13_GEX_hash_oligo-Z0005-ACGTACGTACGTACG.cram"
    )

    def _scan(self, keys, root="s3://czi-novogene", **kwargs):
        s3 = MockS3Client(buckets={"czi-novogene": list(keys)})
        return s3, _descend_to_raw_prefixes(s3, root, **kwargs)

    def test_it_finds_the_foldered_layout_from_a_bucket_root(self):
        _s3, scan = self._scan([self.FOLDERED])

        assert len(scan.prefixes) == 1
        found = scan.prefixes[0]
        assert found.prefix == "labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/"
        assert found.experiment_segment == "REF3"
        assert found.wafer_folders == ("437120",)
        assert found.loose_objects == ()
        assert scan.complete is True

    def test_it_finds_the_flat_layout_by_its_loose_objects(self):
        """Four of six real vendor prefixes have no wafer directory at all."""
        _s3, scan = self._scan([self.FLAT])

        found = scan.prefixes[0]
        assert found.prefix == "project-killifish/NVUS0000000000-20/CH13/raw/"
        assert found.experiment_segment == "CH13"
        assert found.wafer_folders == ()
        assert found.loose_objects == (self.FLAT,)

    def test_a_prefix_can_carry_both_shapes_at_once(self):
        """A wafer folder and a wafer-level file are not mutually exclusive."""
        loose = "labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/437120_LibraryInfo.xml"
        _s3, scan = self._scan([self.FOLDERED, loose])

        found = scan.prefixes[0]
        assert found.wafer_folders == ("437120",)
        assert found.loose_objects == (loose,)

    def test_several_deliveries_under_one_root_are_all_found(self):
        _s3, scan = self._scan([self.FOLDERED, self.FLAT])

        assert sorted(p.experiment_segment for p in scan.prefixes) == ["CH13", "REF3"]

    def test_it_never_lists_without_a_delimiter(self):
        """The mechanical proof that this is not the flat bucket scan.

        index_untrimmed_sources paginates with no Delimiter, which is why a
        bucket root is refused there; if this walk ever did the same it would
        still pass every other test here while listing every object in the
        bucket. Nothing else would notice.
        """
        s3, _scan = self._scan([self.FOLDERED, self.FLAT])

        assert s3.listing_calls
        assert [call for call in s3.listing_calls if call[2] != "/"] == []

    def test_a_wide_level_is_paginated_rather_than_truncated(self):
        keys = [
            f"proj/order-{i:04d}/REF3/raw/437120/{P04_STEM}.cram" for i in range(1200)
        ]
        s3 = MockS3Client(buckets={"czi-novogene": keys}, list_limit=1000)

        scan = _descend_to_raw_prefixes(s3, "s3://czi-novogene", max_listings=10_000)

        assert len(scan.prefixes) == 1200

    def test_the_depth_cap_stops_a_decoy_tree_and_says_so(self):
        deep = "a/b/c/d/e/f/g/REF3/raw/437120/x.cram"
        _s3, scan = self._scan([deep])

        assert scan.prefixes == []
        assert scan.depth_capped
        assert scan.complete is False

    def test_the_depth_cap_still_reaches_a_bucket_root_delivery(self):
        """Three levels to raw/ from a bucket root, so the default must allow it."""
        _s3, scan = self._scan([self.FOLDERED])

        assert len(scan.prefixes) == 1
        assert scan.depth_capped == ()

    def test_a_denied_branch_is_recorded_and_the_rest_still_walked(self):
        s3 = MockS3Client(
            buckets={"czi-novogene": [self.FOLDERED, self.FLAT]},
            deny_prefixes={"project-killifish/": "AccessDenied"},
        )

        scan = _descend_to_raw_prefixes(s3, "s3://czi-novogene")

        assert [p.experiment_segment for p in scan.prefixes] == ["REF3"]
        assert scan.unreadable == (
            ("s3://czi-novogene/project-killifish/", "AccessDenied"),
        )
        assert scan.complete is False

    def test_a_denied_root_yields_nothing_rather_than_raising(self):
        s3 = MockS3Client(
            buckets={"czi-novogene": [self.FOLDERED]},
            deny_prefixes={"": "AccessDenied"},
        )

        scan = _descend_to_raw_prefixes(s3, "s3://czi-novogene")

        assert scan.prefixes == []
        assert scan.unreadable and scan.unreadable[0][1] == "AccessDenied"

    def test_the_listing_budget_stops_the_walk_and_is_reported(self):
        keys = [f"proj/order-{i:03d}/REF3/raw/437120/x.cram" for i in range(50)]
        s3 = MockS3Client(buckets={"czi-novogene": keys})

        scan = _descend_to_raw_prefixes(s3, "s3://czi-novogene", max_listings=5)

        assert scan.budget_exhausted is True
        assert scan.listings <= 5
        assert scan.complete is False

    def test_a_project_root_is_cheaper_than_a_bucket_root(self):
        """The documented recommendation, made checkable rather than asserted."""
        keys = [self.FOLDERED, self.FLAT]
        bucket_s3 = MockS3Client(buckets={"czi-novogene": keys})
        project_s3 = MockS3Client(buckets={"czi-novogene": keys})

        _descend_to_raw_prefixes(bucket_s3, "s3://czi-novogene")
        _descend_to_raw_prefixes(project_s3, "s3://czi-novogene/labalpha-seahub-bcp")

        assert len(project_s3.listing_calls) < len(bucket_s3.listing_calls)

    def test_no_root_is_no_work(self):
        s3 = MockS3Client(buckets={"czi-novogene": [self.FOLDERED]})

        scan = _descend_to_raw_prefixes(s3, [])

        assert scan.prefixes == [] and s3.listing_calls == []
        assert scan.complete is True


class TestDiscoverUntrimmedSources:
    """Vendor deliveries found from the upload's own wafers, not from a typed order."""

    ROOT = f"s3://czi-novogene/{PROJECT}"

    def _seeds(self, keys=None):
        keys = list(keys or ref3_trimmed_keys())
        return _wafer_seeds(trimmed_keys=keys, trimmed_index=index_trimmed_upload(keys))

    def _s3(self, vendor_keys=None, **kwargs):
        keys = list(
            vendor_keys
            if vendor_keys is not None
            else ref3_vendor_keys(extra_wafer=True)
        )
        return MockS3Client(buckets={"czi-novogene": keys}, **kwargs)

    def _discover(self, s3=None, root=None, **kwargs):
        return _discover_untrimmed_sources(
            s3 or self._s3(), root or self.ROOT, self._seeds(), **kwargs
        )

    def test_it_reproduces_the_order_prefix_index_exactly(self):
        """The headline: the same answer, without naming either order.

        The operator has to know and type both orders to get this index today,
        and REF3_P05_1 is in neither the obvious one nor findable by name.
        """
        s3 = self._s3()
        by_order = index_untrimmed_sources(
            s3, [vendor_uri(VENDOR_ORDER), vendor_uri(VENDOR_ORDER_2)]
        )
        discovered = self._discover(s3)

        assert sorted(discovered.index) == sorted(by_order.index)
        assert {k: e.cram_key for k, e in discovered.index.items()} == {
            k: e.cram_key for k, e in by_order.index.items()
        }

    @pytest.mark.parametrize(
        "segment",
        ["REF3-reupload", "reupload_REF3", "ref3", "REF2", "delivery_7"],
    )
    def test_a_renamed_experiment_folder_is_still_found(self, segment):
        """Every one of these orphans the whole upload under name matching.

        source_experiment_matches accepts REF3_reupload and refuses all of these,
        which is what made the ExperimentID a bad key: the wafer inside the folder
        is unchanged by whatever the folder is called.
        """
        renamed = [k.replace("/REF3/", f"/{segment}/") for k in ref3_vendor_keys()]
        discovered = self._discover(self._s3(renamed))

        assert len(discovered.index) == 4
        assert {c.experiment_segment for c in discovered.coverage} == {segment}

    def test_a_never_trimmed_sibling_wafer_is_still_indexed(self):
        """The step-1 pin, now green through discovery rather than a typed order.

        Wafer 440000 cannot be a seed -- nothing was uploaded for it -- so this
        works only because a delivery a seed *did* point at is indexed whole.
        """
        discovered = self._discover()

        assert UNTRIMMED_WAFER in {e.wafer for e in discovered.index.values()}

    def test_the_sibling_wafer_reports_not_trimmed_and_a_data_gap(self):
        """Both reporting paths, which is what makes a forgotten plate audible."""
        keys = ref3_trimmed_keys()
        discovered = self._discover()

        report = reconcile_trimming(discovered.index, index_trimmed_upload(keys))
        not_trimmed = [
            r
            for r in report.rows
            if r["category"] == "not_trimmed" and r["wafer"] == UNTRIMMED_WAFER
        ]

        assert [r["ug"] for r in not_trimmed] == ["Z0500", "Z0501"]

    def test_turning_off_expansion_loses_the_sibling_again(self):
        """States the cost of the option rather than leaving it implied."""
        discovered = self._discover(expand_siblings=False)

        assert UNTRIMMED_WAFER not in {e.wafer for e in discovered.index.values()}

    def test_another_experiment_under_the_same_root_is_not_indexed(self):
        """What replaces the ExperimentID filter: the upload's own contents.

        The foreign delivery is enumerated -- the descent walked past it -- and
        named, but none of its wells reach the comparison.
        """
        foreign = [
            f"{PROJECT}/NVUS0000000000-99/GENE7/raw/999999/"
            "999999-GENE7_P01_A1_GEX_hash_oligo-Z0900-ACGTACGTACGTACG.cram"
        ]
        discovered = self._discover(self._s3(ref3_vendor_keys() + foreign))

        assert "999999" not in {e.wafer for e in discovered.index.values()}
        unseeded = [
            f
            for f in discovered.findings
            if f["category"] == "unseeded_vendor_delivery"
        ]
        assert len(unseeded) == 1
        assert "GENE7" in unseeded[0]["detail"]
        assert ("s3://czi-novogene/" + foreign[0].rsplit("/", 2)[0] + "/") in [
            uri for uri, _seg in discovered.search.unseeded
        ]

    def test_a_seed_found_nowhere_is_reported(self):
        """439000 is uploaded but was never delivered -- the orphan direction."""
        discovered = self._discover()

        assert discovered.search.not_located == ("439000",)
        assert [
            f["wafer"]
            for f in discovered.findings
            if f["category"] == "wafer_not_found"
        ] == ["439000"]

    def test_a_wafer_in_two_deliveries_is_one_row_not_one_per_well(self):
        copy = [
            k.replace("NVUS0000000000-11", "NVUS0000000000-77")
            for k in ref3_vendor_keys()
        ]
        discovered = self._discover(self._s3(ref3_vendor_keys() + copy))

        multiple = [
            f
            for f in discovered.findings
            if f["category"] == "wafer_multiple_deliveries"
        ]
        assert {f["wafer"] for f in multiple} == set(discovered.search.collided)
        assert len(multiple) == len(discovered.search.collided)

    def test_a_located_delivery_that_yields_no_wells_is_reported_not_silent(self):
        """The measured ScaleBio delivery: 192 CRAMs, no Z#### UG anywhere.

        The wafer is seeded, so its delivery is wanted and indexed -- and every
        filename is refused by both stem patterns, so it indexes nothing. Without
        this the coverage row reads cram_keys=0, which is indistinguishable from
        the prefix being wrong. Present-and-unreadable is a different problem from
        absent, and only one of them is the operator's to fix.
        """
        wafer = "426971"
        trimmed = [
            f"labalpha-seahub-bcp/REF3/raw/RNA3_098/{wafer}/"
            f"{wafer}-RNA3_098_A1_GEX_hash_oligo-Z0700-ACGTACGTACGTACG.trim.cram"
        ]
        scale = [
            f"{PROJECT}/NVUS0000000000-04/RNA3_098/raw/{wafer}/"
            f"{wafer}-RNA3-098C_GEX_QSR-7_10A.cram"
        ]
        discovered = _discover_untrimmed_sources(
            self._s3(scale), self.ROOT, self._seeds(trimmed)
        )

        no_wells = [
            f for f in discovered.findings if f["category"] == "wafer_folder_no_wells"
        ]
        assert len(no_wells) == 1
        assert no_wells[0]["wafer"] == wafer
        assert discovered.index == {}
        assert [c.skipped_reason for c in discovered.coverage] == ["no parseable wells"]

    def test_a_flat_vendor_layout_is_discovered_by_filename(self):
        """No wafer directory: four of six real vendor prefixes look like this."""
        stem = "436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT"
        flat = [f"{PROJECT}/NVUS0000000000-33/REF3/raw/{stem}.cram"]
        discovered = self._discover(self._s3(flat))

        assert ("436830", "Z0169") in discovered.index

    def test_a_non_wafer_token_is_reported_rather_than_searched_for(self):
        keys = ref3_trimmed_keys() + [
            "labalpha-seahub-bcp/REF3/raw/REF3_P09_1/not_a_wafer/x.trim.cram"
        ]
        seeds = self._seeds(keys)
        discovered = _discover_untrimmed_sources(self._s3(), self.ROOT, seeds)

        assert "not_a_wafer" in seeds.rejected
        assert [
            f["category"]
            for f in discovered.findings
            if f["category"] == "wafer_seed_rejected"
        ] == ["wafer_seed_rejected"]

    def test_a_refused_listing_becomes_a_finding_not_an_exception(self):
        s3 = self._s3(deny_prefixes={f"{PROJECT}/{VENDOR_ORDER_2}/": "AccessDenied"})
        discovered = self._discover(s3)

        assert any(
            f["category"] == "search_prefix_unreadable" for f in discovered.findings
        )
        assert discovered.search.complete is False

    def test_an_exhausted_budget_is_a_finding_not_a_silent_short_answer(self):
        discovered = self._discover(max_listings=2)

        assert discovered.search.budget_exhausted is True
        assert any(
            f["category"] == "search_budget_exhausted" for f in discovered.findings
        )

    def test_coverage_is_one_row_per_delivery_with_its_provenance(self):
        """One row labelled UNKNOWN_ORDER for the whole search would be useless."""
        discovered = self._discover()

        assert len(discovered.coverage) == 2
        assert all(c.experiment_segment == "REF3" for c in discovered.coverage)
        assert {c.found_by for c in discovered.coverage} <= {"seed", "sibling"}
        assert sum(c.indexed for c in discovered.coverage) == len(discovered.index)

    def test_running_it_twice_gives_the_same_answer_twice(self):
        """A fresh container per run, so re-executing the cell cannot double rows."""
        s3 = self._s3()
        first = self._discover(s3)
        second = self._discover(s3)

        assert len(first.findings) == len(second.findings)
        assert sorted(first.index) == sorted(second.index)

    def test_no_root_or_no_seed_is_no_work(self):
        s3 = self._s3()

        assert _discover_untrimmed_sources(s3, [], self._seeds()).index == {}
        assert _discover_untrimmed_sources(s3, self.ROOT, WaferSeeds()).index == {}
        assert s3.listing_calls == []

    def test_the_summary_names_what_an_operator_has_to_act_on(self):
        summary = self._discover().search.summary()

        assert "located 4 of 5 seed wafer(s)" in summary
        assert "not found: 439000" in summary
