"""
Tests for corrected-name composition and the per-well roll-up (qa_seahub_rename).
"""

from __future__ import annotations

import pytest

from qa_constants import (
    SEAHUB_BARE_SUFFIXES,
    SEAHUB_BARE_TO_TRIM_SUFFIX,
    SEAHUB_RENAMEABLE_SOP_TYPES,
    SEAHUB_SOP_RULES,
    SEAHUB_WELL_VERDICTS,
)
from qa_seahub_rename import (
    RENAME_COLUMNS,
    build_rename_mapping,
    expected_trimmed_key,
    roll_up_wells,
    rollup_summary,
    source_sublibrary_segment,
)
from qa_seahub_source import SourceEntry, index_untrimmed_sources

from tests.qa_seahub_helpers import (
    BUCKET,
    JUNK_NAMES,
    RAW,
    VENDOR_ORDER,
    ref3_trimmed_keys,
    ref3_vendor_keys,
    vendor_uri,
)
from tests.test_qa_gather import MockS3Client

# Reviewer-supplied known-good objects, in both labs' shapes.
GOOD_LABALPHA = (
    "labalpha-seahub-bcp/REF3/raw/REF3_P05_2/436830/"
    "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
)
GOOD_LABBETA = (
    "labbeta-seahub-bcp/CHEM3-R100/raw/R100E/441389/"
    "441389-R100E_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram"
)

WELL_A = f"{RAW}/P04_1/437120/437120-437120-REF3_P04_1_A1-Z0001-CAGCTCGAATGCGAT"
WELL_B = (
    f"{RAW}/P05_1/436830/436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT"
)
WELL_G2 = f"{RAW}/P06_1/439000/439000-439000-REF3_P06_1_B3-Z0401-ACGTACGTACGTACA"

VENDOR_A = SourceEntry(
    wafer="437120",
    ug="Z0001",
    barcode="CAGCTCGAATGCGAT",
    group="REF3_P04_1_A1",
    assay="GEX_hash_oligo",
    cram_key="labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/437120/x.cram",
    bucket="czi-novogene",
)


def _vendor_index():
    s3 = MockS3Client(keys=ref3_vendor_keys())
    return index_untrimmed_sources(
        s3, [vendor_uri(VENDOR_ORDER), vendor_uri("NVUS0000000000-12")]
    ).index


class TestKnownGoodIsIdempotent:
    """The primary regression guard: a compliant upload proposes nothing."""

    @pytest.mark.parametrize(
        "bucket,key",
        [("czi-labalpha", GOOD_LABALPHA), ("czi-labbeta", GOOD_LABBETA)],
    )
    def test_a_compliant_object_proposes_itself(self, bucket, key):
        proposal = expected_trimmed_key(bucket, key)

        assert proposal.compliant is True
        assert proposal.defects == ()
        assert proposal.unresolved == ()
        assert proposal.expected_s3_uri == proposal.current_s3_uri

    @pytest.mark.parametrize(
        "bucket,key",
        [("czi-labalpha", GOOD_LABALPHA), ("czi-labbeta", GOOD_LABBETA)],
    )
    def test_a_compliant_object_is_still_compliant_with_a_vendor_match(
        self, bucket, key
    ):
        source = SourceEntry(
            wafer="436830",
            ug="Z0169",
            barcode="CTCGCAATAGATGAT",
            group="REF3_P05_2_A10",
            assay="GEX_hash_oligo",
            cram_key="v/REF3/raw/436830/x.cram",
        )
        assert expected_trimmed_key(bucket, key, source).compliant is True

    def test_a_compliant_listing_yields_an_empty_mapping(self):
        keys = [
            GOOD_LABALPHA.replace(".trim.cram", suffix)
            for suffix in (
                ".trim.cram",
                ".trim.csv",
                ".trim.stderr",
                ".trim.stdout",
                ".trim_fail.csv",
            )
        ]

        mapping = build_rename_mapping("czi-labalpha", keys)

        assert mapping.rows == []
        assert mapping.compliant_objects == len(keys)


class TestExpectedTrimmedKey:
    def test_all_four_defects_compose_into_one_key(self):
        proposal = expected_trimmed_key(BUCKET, f"{WELL_A}.cram", VENDOR_A)

        assert proposal.expected_s3_uri == (
            f"s3://{BUCKET}/{RAW}/REF3_P04_1/437120/"
            "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram"
        )
        assert proposal.defects == (
            "duplicated_wafer_token",
            "invalid_sublibrary_type",
            "missing_trim_infix",
            "sublibrary_folder_truncated",
        )
        assert proposal.name_source == "vendor"

    def test_a_truncated_folder_changes_only_the_folder(self):
        proposal = expected_trimmed_key(BUCKET, f"{WELL_B}.trim.cram")

        assert proposal.defects == ("sublibrary_folder_truncated",)
        assert proposal.expected_s3_uri.endswith(
            "/REF3_P05_1/436830/"
            "436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
        )
        assert proposal.name_source == "inferred"

    @pytest.mark.parametrize("bare,trim", sorted(SEAHUB_BARE_TO_TRIM_SUFFIX.items()))
    def test_every_bare_suffix_maps_to_its_trim_form(self, bare, trim):
        proposal = expected_trimmed_key(
            BUCKET,
            f"{RAW}/REF3_P05_2/436831/"
            f"436831-REF3_P05_2_A11_GEX_hash_oligo-Z0170-CTCGCAATAGATGAC{bare}",
        )

        assert proposal.defects == ("missing_trim_infix",)
        assert proposal.expected_s3_uri.endswith(trim)

    def test_the_suffix_map_is_total_over_the_bare_family(self):
        """Guards against a KeyError on a suffix nobody mapped."""
        assert set(SEAHUB_BARE_TO_TRIM_SUFFIX) == set(SEAHUB_BARE_SUFFIXES)

    def test_a_missing_type_tag_without_a_vendor_match_is_unresolved(self):
        proposal = expected_trimmed_key(BUCKET, f"{WELL_G2}.cram")

        assert proposal.expected_s3_uri == ""
        assert proposal.unresolved == ("invalid_sublibrary_type",)
        assert proposal.name_source == "inferred"

    def test_two_different_wafer_tokens_are_refused(self):
        """Not a repair QA may guess at."""
        key = (
            f"{RAW}/P04_1/437120/437120-438514-REF3_P04_1_A1-Z0001-CAGCTCGAATGCGAT.cram"
        )
        proposal = expected_trimmed_key(BUCKET, key, VENDOR_A)

        assert proposal.unresolved == ("conflicting_wafer_tokens",)
        assert proposal.expected_s3_uri == ""

    def test_a_vendor_type_conflicting_with_the_filename_is_refused(self):
        conflicting = SourceEntry(
            wafer="436830",
            ug="Z0169",
            barcode="CTCGCAATAGATGAT",
            group="REF3_P05_1_A10",
            assay="CRI",
            cram_key="v/REF3/raw/436830/x.cram",
        )
        proposal = expected_trimmed_key(BUCKET, f"{WELL_B}.trim.cram", conflicting)

        assert proposal.unresolved == ("conflicting_sublibrary_type",)

    def test_junk_gets_no_destination(self):
        proposal = expected_trimmed_key(BUCKET, f"{RAW}/P04_1/437120/login.html")

        assert proposal.defects == ("non_sequencing_artifact",)
        assert proposal.expected_s3_uri == ""

    def test_a_wrong_depth_key_is_unresolved(self):
        proposal = expected_trimmed_key(BUCKET, "labalpha-seahub-bcp/REF3/raw/x.cram")

        assert proposal.unresolved == ("bad_path_depth",)

    def test_every_proposal_is_itself_sop_clean(self):
        """A plan that produces another bad name is worse than no plan."""
        index = _vendor_index()
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), index)

        for row in mapping.moveable():
            assert row["expected_s3_uri"].startswith(f"s3://{BUCKET}/")
            assert row["unresolved"] == ""


class TestSourceSublibrarySegment:
    def test_the_measured_vendor_layout_yields_the_experiment_id(self):
        """Which is exactly why it is not the primary source of the sublibrary."""
        key = (
            "labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/438514/"
            "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT.cram"
        )
        assert source_sublibrary_segment(key) == "REF3"

    def test_an_unexpected_shape_yields_nothing(self):
        assert source_sublibrary_segment("a/b/c.cram") == ""


class TestBuildRenameMapping:
    def test_columns_and_ordering(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert tuple(mapping.rows[0]) == RENAME_COLUMNS
        assert [r["current_s3_uri"] for r in mapping.rows] == sorted(
            r["current_s3_uri"] for r in mapping.rows
        )

    def test_every_object_is_accounted_for(self):
        keys = ref3_trimmed_keys()
        mapping = build_rename_mapping(BUCKET, keys, _vendor_index())

        assert mapping.total_objects == len(keys)
        assert len(mapping.rows) + mapping.compliant_objects == len(keys)

    def test_the_compliant_well_contributes_no_rows(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert not any(
            "436831-REF3_P05_2_A11" in r["current_s3_uri"] for r in mapping.rows
        )

    def test_junk_is_reported_but_never_given_a_destination(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())
        junk = [r for r in mapping.rows if r["status"] == "not_data"]

        assert len(junk) == len(JUNK_NAMES)
        assert all(r["expected_s3_uri"] == "" for r in junk)

    def test_the_data_gap_well_maps_only_the_artifacts_that_exist(self):
        """A rename mapping moves objects; it cannot move an absent CRAM."""
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())
        rows = [r for r in mapping.rows if "438514-438514" in r["current_s3_uri"]]

        assert len(rows) == 4
        assert not any(r["current_s3_uri"].endswith(".cram") for r in rows)

    def test_name_source_is_vendor_only_where_a_vendor_key_matched(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())

        by_source = {r["name_source"] for r in mapping.moveable()}
        assert by_source <= {"vendor", "inferred"}
        assert any(
            r["name_source"] == "vendor"
            for r in mapping.rows
            if "437120-437120" in r["current_s3_uri"]
        )

    def test_destinations_are_unique_and_do_not_overwrite(self):
        keys = ref3_trimmed_keys()
        mapping = build_rename_mapping(BUCKET, keys, _vendor_index())
        destinations = [r["expected_s3_uri"] for r in mapping.moveable()]

        assert len(destinations) == len(set(destinations))
        assert not set(destinations) & {f"s3://{BUCKET}/{k}" for k in keys}

    def test_two_objects_claiming_one_destination_are_blocked(self):
        """The bare and trim forms of one artifact collapse onto one key."""
        stem = "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
        keys = [
            f"{RAW}/REF3_P04_1/437120/{stem}.cram",
            f"{RAW}/REF3_P04_1/437120/{stem}.trim.cram",
        ]

        mapping = build_rename_mapping(BUCKET, keys)

        assert mapping.moveable() == []
        assert [c["kind"] for c in mapping.collisions] == ["destination_exists"]

    def test_a_many_to_one_collision_blocks_both(self):
        base = f"{RAW}/P04_1/437120"
        stem = "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
        keys = [f"{base}/{stem}.cram", f"{base}/437120-{stem}.cram"]

        mapping = build_rename_mapping(BUCKET, keys)

        assert mapping.moveable() == []
        assert any(c["kind"] == "many_to_one" for c in mapping.collisions)

    def test_the_mapping_is_stable_across_input_order(self):
        keys = ref3_trimmed_keys()
        index = _vendor_index()

        forward = build_rename_mapping(BUCKET, keys, index).rows
        reverse = build_rename_mapping(BUCKET, list(reversed(keys)), index).rows

        assert forward == reverse


class TestRollUpWells:
    def test_every_verdict_is_exercised(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert rollup_summary(rollup.rows) == {
            "COMPLIANT": 1,
            "RENAMEABLE": 3,
            "DATA_GAP": 1,
            "UNKNOWN": 1,
        }

    def test_the_data_gap_well_names_the_recoverable_vendor_key(self):
        """So the action is 're-run trim', not 'hunt for missing data'."""
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())
        row = next(r for r in rollup.rows if r["verdict"] == "DATA_GAP")

        assert row["wafer"] == "438514"
        assert row["missing_artifacts"] == ".cram"
        assert "the vendor delivered it as" in row["detail"]

    def test_a_well_with_no_derivable_name_is_unknown(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())
        row = next(r for r in rollup.rows if r["verdict"] == "UNKNOWN")

        assert row["ug"] == "Z0401"
        assert "invalid_sublibrary_type" in row["detail"]

    def test_data_gap_outranks_an_underivable_name(self):
        """Without a vendor index the CRAM is still absent, and that wins."""
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), None)
        row = next(r for r in rollup.rows if r["wafer"] == "438514")

        assert row["verdict"] == "DATA_GAP"

    def test_a_renameable_well_carries_only_repairable_defects(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        for row in rollup.rows:
            if row["verdict"] != "RENAMEABLE":
                continue
            defects = set(filter(None, row["defects"].split("|")))
            assert defects <= SEAHUB_RENAMEABLE_SOP_TYPES

    def test_one_row_per_well(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        identities = [(r["wafer"], r["ug"]) for r in rollup.rows]
        assert len(identities) == len(set(identities))

    def test_a_vendor_well_with_nothing_uploaded_is_unknown(self):
        index = _vendor_index()
        rollup = roll_up_wells(BUCKET, [], index)

        assert {r["verdict"] for r in rollup.rows} == {"UNKNOWN"}
        assert all("nothing was uploaded" in r["detail"] for r in rollup.rows)

    def test_a_wrong_depth_object_is_counted_as_unaccounted(self):
        rollup = roll_up_wells(BUCKET, ["labalpha-seahub-bcp/REF3/raw/stray.cram"])

        assert rollup.rows == []
        assert rollup.unaccounted == 1

    def test_summary_keys_are_the_documented_vocabulary(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert tuple(rollup_summary(rollup.rows)) == SEAHUB_WELL_VERDICTS


class TestRuleVocabulary:
    def test_renameable_types_are_a_subset_of_the_rules(self):
        assert SEAHUB_RENAMEABLE_SOP_TYPES <= SEAHUB_SOP_RULES
