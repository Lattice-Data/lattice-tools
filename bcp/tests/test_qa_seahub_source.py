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
    derive_source_experiment,
    derive_source_order,
    index_trimmed_upload,
    index_untrimmed_source,
    index_untrimmed_sources,
    load_source_read_counts,
    normalize_source_uris,
    parse_source_uri,
    source_experiment_matches,
    source_order_by_wafer,
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
        index = index_untrimmed_source(s3, SOURCE_URI)
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
            stem: {"JumboSciGEX": {"failed": 1, "total": 260527531}}
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
        return {stem: {"JumboSciGEX": {"failed": 100, "total": total}}}

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
            TRIM_STEM_A: {
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

    def test_a_well_delivered_twice_keeps_one_and_reports_the_other(self):
        """A re-delivery must not silently overwrite.

        The previous single-prefix indexer assigned unconditionally, so one of
        the two orders vanished with no record.
        """
        s3 = self._s3(
            ("NVUS0000000000-11", "437120", P04_STEM),
            ("NVUS0000000000-12", "437120", P04_STEM),
        )
        result = index_untrimmed_sources(s3, [ORDER_11, ORDER_12])

        assert len(result.index) == 1
        assert result.index[("437120", "Z0001")].source_order == "NVUS0000000000-11"
        duplicates = [
            f for f in result.findings if f["category"] == "duplicate_source_well"
        ]
        assert len(duplicates) == 1
        assert "NVUS0000000000-12" in duplicates[0]["detail"]

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
        fail_counts = {TRIM_STEM_A: {"JumboSciGEX": {"failed": 1, "total": 100}}}
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
