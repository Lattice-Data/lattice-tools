"""
Tests for SeaHub SOP filename/path validation and family-aware completeness.

Cases are drawn from two real S3 listings reviewed alongside these changes: a
trimmed upload (2035 objects) in which one sublibrary follows the SOP and six
carry trim outputs with the ``.trim`` infix dropped, and a vendor delivery
(5389 objects).  Lab and order identifiers are sanitized here.
"""

from __future__ import annotations

import os

from qa_checks import check_expected_raw_files, check_extra_raw_files
from qa_constants import SEAHUB_STEM_RE
from qa_gather import gather_qa_data
from qa_mods import grab_seahub_trim_fail_csv, seahub_stem_and_family
from qa_seahub_sop import (
    sop_violation_summary,
    validate_seahub_key,
    validate_seahub_stems,
)

from tests.test_qa_gather import MockS3Client, _make_ctx

QA_FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "qa")
SEAHUB_TRIM_FAIL = os.path.join(QA_FIXTURES_DIR, "seahub_trim_fail_sample.csv")

# Reviewer-supplied known-good objects.  The labbeta name has no well token and
# a hyphenated ExperimentID; the labalpha name has a well token and an
# ExperimentID that legitimately appears inside the sublibrary name.
GOOD_LABBETA = (
    "labbeta-seahub-bcp/CHEM3-R100/raw/R100E/441389/"
    "441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram"
)
GOOD_LABALPHA = (
    "labalpha-seahub-bcp/REF3/raw/REF3_P05_2/436830/"
    "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
)
# A real object from the upload: four independent defects at once.
REAL_MISNAMED = (
    "labalpha-seahub-bcp/REF3/raw/P05_2/436830/"
    "436830-436830-REF3_P05_2_A10-Z0169-CTCGCAATAGATGAT.cram"
)

BARE_STEM = "438514-438514-REF3_P07_1_A3-Z0305-CACACACAACATGAT"
BARE_DIR = "labalpha-seahub-bcp/REF3/raw/P07_1/438514"
TRIM_STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
TRIM_DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479"


def _types(violations) -> set[str]:
    return {v.type for v in violations}


class TestStemFamily:
    def test_trim_family(self):
        assert seahub_stem_and_family(f"{TRIM_STEM}.trim.cram") == (
            TRIM_STEM,
            ".trim.cram",
            "trim",
        )

    def test_bare_family(self):
        assert seahub_stem_and_family(f"{BARE_STEM}.cram") == (
            BARE_STEM,
            ".cram",
            "bare",
        )

    def test_metadata_suffix_wins_over_cram(self):
        # Longest-first matching: ".cram" must not truncate the sidecar name.
        assert seahub_stem_and_family(f"{BARE_STEM}.cram-metadata.json") == (
            BARE_STEM,
            ".cram-metadata.json",
            "bare",
        )

    def test_unknown_suffix(self):
        assert seahub_stem_and_family("index.html") is None


class TestKnownGoodNamesAreClean:
    def test_labbeta_example_has_no_violations(self):
        assert validate_seahub_key("czi-labbeta", GOOD_LABBETA) == []

    def test_labalpha_example_has_no_violations(self):
        assert validate_seahub_key("czi-labalpha", GOOD_LABALPHA) == []

    def test_no_well_token_parses_longest_type_first(self):
        # Guards the greedy-parse trap: group R100E + type GEX_hash_oligo, not
        # group R100E_GEX + type hash_oligo.
        match = SEAHUB_STEM_RE.match(
            "441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
        )
        assert match is not None
        assert match.group("group") == "R100E"
        assert match.group("assay") == "GEX_hash_oligo"


class TestRealMisnamedObject:
    def test_reports_exactly_the_four_defects(self):
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        assert _types(violations) == {
            "missing_trim_infix",
            "repeated_token",
            "invalid_sublibrary_type",
            "sublibrary_mismatch",
        }

    def test_expected_name_restores_the_trim_infix(self):
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        infix = next(v for v in violations if v.type == "missing_trim_infix")
        assert infix.expected_name.endswith(".trim.cram")

    def test_sublibrary_mismatch_names_both_sides(self):
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        detail = next(v for v in violations if v.type == "sublibrary_mismatch").detail
        assert "P05_2" in detail


class TestPathRules:
    def test_bad_bucket(self):
        violations = validate_seahub_key("labalpha-data", GOOD_LABALPHA)
        assert "bad_bucket" in _types(violations)

    def test_lab_project_mismatch(self):
        violations = validate_seahub_key("czi-labbeta", GOOD_LABALPHA)
        assert "lab_project_mismatch" in _types(violations)

    def test_bad_path_depth(self):
        violations = validate_seahub_key(
            "czi-labalpha",
            "labalpha-seahub-bcp/REF3/raw/436830/"
            "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram",
        )
        assert "bad_path_depth" in _types(violations)

    def test_hyphenated_experiment_id_is_not_split(self):
        # CHEM3-R100 must survive as one path segment.
        assert "bad_path_depth" not in _types(
            validate_seahub_key("czi-labbeta", GOOD_LABBETA)
        )


class TestFilenameRules:
    def test_wafer_mismatch(self):
        key = (
            "labalpha-seahub-bcp/REF3/raw/REF3_P05_2/436830/"
            "999999-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
        )
        assert "wafer_mismatch" in _types(validate_seahub_key("czi-labalpha", key))

    def test_bad_well(self):
        key = (
            "labalpha-seahub-bcp/REF3/raw/REF3_P05_2/436830/"
            "436830-REF3_P05_2_Z99_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
        )
        assert "bad_well" in _types(validate_seahub_key("czi-labalpha", key))

    def test_sublibrary_mismatch_when_folder_lacks_prefix(self):
        key = (
            "labalpha-seahub-bcp/REF3/raw/P05_2/436830/"
            "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
        )
        assert "sublibrary_mismatch" in _types(validate_seahub_key("czi-labalpha", key))

    def test_unexpected_suffix_for_browser_junk(self):
        key = "labalpha-seahub-bcp/REF3/raw/P05_2/436830/index.html"
        assert "unexpected_suffix" in _types(validate_seahub_key("czi-labalpha", key))

    def test_experiment_id_in_sublibrary_is_not_a_repeated_token(self):
        # REF3 appears in both the ExperimentID folder and the sublibrary name;
        # that is the reviewer's own known-good form.  repeated_token is scoped
        # to the basename precisely so this does not false-positive.
        assert "repeated_token" not in _types(
            validate_seahub_key("czi-labalpha", GOOD_LABALPHA)
        )


class TestStemLevelReporting:
    def test_one_stem_reported_once_across_its_artifacts(self):
        keys = [
            f"{BARE_DIR}/{BARE_STEM}{suffix}"
            for suffix in (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
        ]
        violations = validate_seahub_stems("czi-labalpha", keys)
        summary = sop_violation_summary(violations)
        assert summary["repeated_token"] == 1
        assert summary["invalid_sublibrary_type"] == 1
        assert summary["missing_trim_infix"] == 1

    def test_junk_files_are_reported_individually(self):
        keys = [
            f"{BARE_DIR}/index.html",
            f"{BARE_DIR}/ug-icon.png",
        ]
        violations = validate_seahub_stems("czi-labalpha", keys)
        assert sop_violation_summary(violations)["unexpected_suffix"] == 2


class TestFamilyAwareCompleteness:
    def _bare_well(self, stem: str, suffixes: tuple[str, ...]) -> list[str]:
        return [f"{BARE_DIR}/{stem}{suffix}" for suffix in suffixes]

    def test_complete_bare_well_is_not_reported_missing(self):
        keys = self._bare_well(
            BARE_STEM,
            (".cram", ".cram-metadata.json", ".csv", ".stderr", ".stdout", "_fail.csv"),
        )
        all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
        assert all_good == 1
        assert lost == []

    def test_absent_cram_in_bare_well_is_reported(self):
        # The real P07_1/438514 A3 defect: every artifact but the CRAM.
        keys = self._bare_well(BARE_STEM, (".csv", ".stderr", ".stdout", "_fail.csv"))
        all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
        assert all_good == 0
        assert len(lost) == 1
        assert lost[0]["path"] == BARE_STEM
        assert lost[0][".cram"].endswith(f"{BARE_STEM}.cram")

    def test_bare_well_is_not_asked_for_trim_names(self):
        keys = self._bare_well(
            BARE_STEM, (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
        )
        _all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
        assert lost == []

    def test_trim_well_still_requires_trim_names(self):
        keys = [
            f"{TRIM_DIR}/{TRIM_STEM}{s}"
            for s in (".trim.csv", ".trim.stderr", ".trim.stdout", ".trim_fail.csv")
        ]
        _all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
        assert len(lost) == 1
        assert ".trim.cram" in lost[0]

    def test_mixed_family_well_requires_both_sets(self):
        keys = [f"{BARE_DIR}/{BARE_STEM}.cram", f"{BARE_DIR}/{BARE_STEM}.trim.csv"]
        _all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
        assert len(lost) == 1
        assert ".trim.cram" in lost[0]
        assert ".csv" in lost[0]

    def test_bare_metadata_sidecar_is_optional_not_extra(self):
        keys = self._bare_well(
            BARE_STEM,
            (".cram", ".cram-metadata.json", ".csv", ".stderr", ".stdout", "_fail.csv"),
        )
        _all_good, _lost, found = check_expected_raw_files(keys, "seahub_sci")
        assert check_extra_raw_files(keys, found, "seahub_sci") == []

    def test_browser_junk_stays_extra(self):
        keys = self._bare_well(
            BARE_STEM,
            (".cram", ".cram-metadata.json", ".csv", ".stderr", ".stdout", "_fail.csv"),
        ) + [f"{BARE_DIR}/index.html", f"{BARE_DIR}/objects_list_CRO-929.txt"]
        _all_good, _lost, found = check_expected_raw_files(keys, "seahub_sci")
        extra = check_extra_raw_files(keys, found, "seahub_sci")
        assert sorted(e.split("/")[-1] for e in extra) == [
            "index.html",
            "objects_list_CRO-929.txt",
        ]


class TestFailCsvSuffixes:
    def test_absolute_counts_recorded_per_format(self):
        stats: dict = {}
        counts: dict = {}
        grab_seahub_trim_fail_csv(
            stats, "438514", SEAHUB_TRIM_FAIL, fail_counts=counts, stem_key=BARE_STEM
        )
        assert set(counts[BARE_STEM]) == {"JumboSciHash", "JumboSciGEX"}
        gex = counts[BARE_STEM]["JumboSciGEX"]
        assert gex["total"] == 158233602
        assert 0 < gex["failed"] < gex["total"]

    def test_bare_fail_csv_feeds_the_wafer_table(self):
        # A wafer whose only failure CSV is "_fail.csv" (no ".trim" infix) must
        # still produce trimmer stats and appear as a discovered wafer; before
        # this, such wafers rendered as blank rows.
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P07_1/438514"
        keys = [
            f"{raw_dir}/{BARE_STEM}.cram",
            f"{raw_dir}/{BARE_STEM}.csv",
            f"{raw_dir}/{BARE_STEM}.stderr",
            f"{raw_dir}/{BARE_STEM}.stdout",
            f"{raw_dir}/{BARE_STEM}_fail.csv",
        ]
        with open(SEAHUB_TRIM_FAIL) as fh:
            fail_content = fh.read()
        pages = {
            ("labalpha-seahub-bcp/REF3/", "/"): [
                {"CommonPrefixes": [{"Prefix": "labalpha-seahub-bcp/REF3/raw/"}]}
            ],
            ("labalpha-seahub-bcp/REF3/raw/", "/"): [
                {"CommonPrefixes": [{"Prefix": "labalpha-seahub-bcp/REF3/raw/P07_1/"}]}
            ],
            ("labalpha-seahub-bcp/REF3/raw/P07_1/", "/"): [
                {"CommonPrefixes": [{"Prefix": f"{raw_dir}/"}]}
            ],
            (f"{raw_dir}/", ""): [{"Contents": [{"Key": k} for k in keys]}],
        }
        s3 = MockS3Client(
            paginated_pages=pages,
            file_contents={f"{raw_dir}/{BARE_STEM}_fail.csv": fail_content},
        )
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, s3)
        assert data.discovered_wafers == {"438514"}
        assert data.trimmer_failure_stats["438514"]["trimmer_fail"]
        assert BARE_STEM in data.seahub_fail_counts

    def test_gather_records_sop_violations_and_one_summary_warning(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P07_1/438514"
        keys = [f"{raw_dir}/{BARE_STEM}{s}" for s in (".cram", ".csv", ".stdout")]
        pages = {
            ("labalpha-seahub-bcp/REF3/", "/"): [
                {"CommonPrefixes": [{"Prefix": "labalpha-seahub-bcp/REF3/raw/"}]}
            ],
            ("labalpha-seahub-bcp/REF3/raw/", "/"): [
                {"CommonPrefixes": [{"Prefix": "labalpha-seahub-bcp/REF3/raw/P07_1/"}]}
            ],
            ("labalpha-seahub-bcp/REF3/raw/P07_1/", "/"): [
                {"CommonPrefixes": [{"Prefix": f"{raw_dir}/"}]}
            ],
            (f"{raw_dir}/", ""): [{"Contents": [{"Key": k} for k in keys]}],
        }
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, MockS3Client(paginated_pages=pages))
        types = {v["type"] for v in data.sop_violations}
        assert "missing_trim_infix" in types
        assert "repeated_token" in types
        sop_warnings = [
            w for w in data.gathering_warnings if w.startswith("SOP VIOLATIONS:")
        ]
        assert len(sop_warnings) == 1
