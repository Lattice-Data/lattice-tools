"""
Tests for corrected-name composition and the per-well roll-up (qa_seahub_rename).
"""

from __future__ import annotations

import ast
from dataclasses import replace
import json
import pathlib
import re

import pytest

import qa_seahub_rename

from qa_checks import check_expected_raw_files
from qa_constants import (
    SEAHUB_BARE_SUFFIXES,
    SEAHUB_BARE_TO_TRIM_SUFFIX,
    SEAHUB_RENAMEABLE_SOP_TYPES,
    SEAHUB_RENAME_NAME_SOURCES,
    SEAHUB_RENAME_STATUSES,
    SEAHUB_SOP_RULES,
    SEAHUB_VIOLATION_SCOPES,
    SEAHUB_WELL_VERDICTS,
    raw_expected,
)
from qa_seahub_rename import (
    RENAME_COLUMNS,
    WELL_COLUMNS,
    build_rename_mapping,
    expected_trimmed_key,
    roll_up_wells,
    rollup_summary,
)
from qa_seahub_sop import validate_seahub_key, validate_seahub_stems
from qa_seahub_source import SourceEntry, index_untrimmed_sources

from tests.qa_seahub_helpers import (
    BUCKET,
    JUNK_NAMES,
    PROJECT,
    RAW,
    TRIM_SUFFIXES,
    UNTRIMMED_WAFER,
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


# The fixture's folders as delivered: five truncated, one already SOP-correct.
# Pinned by value so a mutation that applies the vendor's SOP name to an uploaded
# row -- rewriting "P05_1" to "REF3_P05_1" -- cannot pass on truthiness.
UPLOADED_SUBLIBRARIES = {
    ("436830", "Z0169"): "P05_1",
    ("436831", "Z0170"): "REF3_P05_2",
    ("437120", "Z0001"): "P04_1",
    ("438514", "Z0305"): "P07_1",
    ("439000", "Z0400"): "P06_1",
    ("439000", "Z0401"): "P06_1",
}


def _vendor_index(extra_wafer: bool = False):
    s3 = MockS3Client(keys=ref3_vendor_keys(extra_wafer=extra_wafer))
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
    def test_columns_and_ordering(self):
        """The well-status CSV is the headline output, so its columns are a
        contract -- pinned the way RENAME_COLUMNS is, which WELL_COLUMNS was not.
        """
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert tuple(rollup.rows[0]) == WELL_COLUMNS
        identities = [(r["wafer"], r["ug"]) for r in rollup.rows]
        assert identities == sorted(identities)

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
        # The SOP name, though this well was delivered bare: what is absent is
        # the artifact kind, and the name it should arrive under is the SOP one.
        assert row["missing_artifacts"] == ".trim.cram"
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

    def test_a_vendor_well_with_nothing_uploaded_is_a_data_gap(self):
        """The largest gap there is, so it is labelled as one.

        As UNKNOWN it was the only kind of gap the notebook left out of
        errors.txt, which writes DATA_GAP rows alone -- and nothing about the
        well is unidentifiable: the vendor names it exactly.
        """
        index = _vendor_index()
        rollup = roll_up_wells(BUCKET, [], index)

        assert {r["verdict"] for r in rollup.rows} == {"DATA_GAP"}
        assert all("nothing was uploaded" in r["detail"] for r in rollup.rows)

    def test_a_well_with_nothing_uploaded_lists_every_missing_artifact(self):
        """Zeroing `required` when there were no keys printed this blank."""
        rollup = roll_up_wells(BUCKET, [], _vendor_index())
        row = rollup.rows[0]

        missing = set(filter(None, row["missing_artifacts"].split("|")))
        assert missing == set(raw_expected["seahub_sci"])

    def test_a_well_with_nothing_uploaded_names_the_vendor_key(self):
        """Otherwise the row says something is absent without saying from where."""
        rollup = roll_up_wells(BUCKET, [], _vendor_index())

        assert all("the vendor delivered it as" in r["detail"] for r in rollup.rows)

    def test_a_well_with_nothing_uploaded_is_locatable(self):
        """sublibrary and well come off the uploaded objects, and there are none.

        Leaving them blank made the highest-severity verdict the least locatable
        row in the headline CSV -- the notebook printed ``wafer 438514 Z0305 ():``
        -- so they are read off the vendor stem, which names both.
        """
        rollup = roll_up_wells(BUCKET, [], _vendor_index())

        assert rollup.rows
        for row in rollup.rows:
            assert row["sublibrary"], row

    @pytest.mark.parametrize(
        "group,sublibrary,well",
        [
            ("REF3_P07_1_A3", "REF3_P07_1", "A3"),
            ("REF3_P07_1_A10", "REF3_P07_1", "A10"),
            # No trailing well token: the whole group is the sublibrary, and
            # there is no well to invent. The labbeta shape does this.
            ("R100E", "R100E", ""),
            ("REF3_P07_1", "REF3_P07_1", ""),
        ],
    )
    def test_the_vendor_stem_is_split_the_way_the_rename_splits_it(
        self, group, sublibrary, well
    ):
        source = SourceEntry(
            wafer="438514",
            ug="Z0305",
            barcode="CACACACAACATGAT",
            group=group,
            assay="GEX_hash_oligo",
            cram_key=f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/438514/x.cram",
        )

        row = roll_up_wells(BUCKET, [], {("438514", "Z0305"): source}).rows[0]

        assert (row["sublibrary"], row["well"]) == (sublibrary, well)

    def test_an_uploaded_well_still_takes_its_names_from_the_upload(self):
        """The vendor split applies only to the row that has no objects.

        Pinned by value, not truthiness: the fixture's folders are truncated, so
        a mutation applying the vendor split to every well would quietly rewrite
        ``P05_1`` to ``REF3_P05_1`` here and still leave every cell non-empty.
        """
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        uploaded = {
            (r["wafer"], r["ug"]): r["sublibrary"] for r in rollup.rows if r["objects"]
        }
        assert uploaded == UPLOADED_SUBLIBRARIES

    def test_a_gap_row_names_the_sublibrary_the_vendor_used(self):
        """The two sources of `sublibrary` can disagree, and that is deliberate.

        An uploaded row carries the folder the upload actually used, which for
        the commonest real defect is truncated; a gap row has no folder, so it
        carries the vendor's SOP name for the same sublibrary. Filtering the CSV
        on one spelling therefore misses the other -- pinned here so the
        divergence stays a documented property rather than a surprise.
        """
        gap = SourceEntry(
            wafer="438515",
            ug="Z0306",
            barcode="CACACACAACATGAA",
            group="REF3_P07_1_A4",
            assay="GEX_hash_oligo",
            cram_key=f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/438515/x.cram",
        )
        index = {**_vendor_index(), ("438515", "Z0306"): gap}

        rows = {
            (r["wafer"], r["ug"]): r["sublibrary"]
            for r in roll_up_wells(BUCKET, ref3_trimmed_keys(), index).rows
        }

        assert rows[("438514", "Z0305")] == "P07_1"  # uploaded, truncated folder
        assert rows[("438515", "Z0306")] == "REF3_P07_1"  # gap, vendor SOP name

    def test_a_doubled_wafer_in_the_vendor_stem_is_not_a_sublibrary(self):
        """index_untrimmed_sources does not normalize a doubled wafer.

        The second wafer stays inside ``group``, and unstripped it reached the
        headline CSV as sublibrary ``438514-REF3_P07_1``.
        """
        source = SourceEntry(
            wafer="438514",
            ug="Z0305",
            barcode="CACACACAACATGAT",
            group="438514-REF3_P07_1_A3",
            assay="GEX_hash_oligo",
            cram_key=f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/438514/x.cram",
        )

        row = roll_up_wells(BUCKET, [], {("438514", "Z0305"): source}).rows[0]

        assert (row["sublibrary"], row["well"]) == ("REF3_P07_1", "A3")

    def test_no_vendor_index_means_no_such_row(self):
        """The control: these rows exist only because a vendor delivery does."""
        assert roll_up_wells(BUCKET, [], None).rows == []

    def test_a_wrong_depth_object_is_counted_as_unaccounted(self):
        rollup = roll_up_wells(BUCKET, ["labalpha-seahub-bcp/REF3/raw/stray.cram"])

        assert rollup.rows == []
        assert rollup.unaccounted == 1

    def test_summary_keys_are_the_documented_vocabulary(self):
        rollup = roll_up_wells(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert tuple(rollup_summary(rollup.rows)) == SEAHUB_WELL_VERDICTS


class TestAllMisspelledArtifactsStillAppearInRollup:
    """A well whose every artifact is misspelled used to vanish from the headline CSV."""

    DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479"
    STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
    MISSPELLED = (
        ".trimmed.ucram",
        ".trimmer_stats.csv",
        ".trim.stderr.txt",
        ".trim.stdout.txt",
        ".trim_failures.csv",
    )

    def _keys(self):
        return [f"{self.DIR}/{self.STEM}{suffix}" for suffix in self.MISSPELLED]

    def test_one_well_with_only_misspelled_artifacts_gets_a_row(self):
        rollup = roll_up_wells(BUCKET, self._keys())

        assert len(rollup.rows) == 1
        assert rollup.rows[0]["verdict"] == "DATA_GAP"
        assert rollup.rows[0]["objects"] == len(self.MISSPELLED)
        assert rollup.unaccounted == 0

    def test_it_does_not_double_report_when_the_vendor_index_also_has_the_well(self):
        source = SourceEntry(
            wafer="430479",
            ug="Z0001",
            barcode="CAGCTCGAATGCGAT",
            group="REF3_P05_1_A1",
            assay="GEX_hash_oligo",
            cram_key="labalpha-seahub-bcp/NVUS0000000000-11/REF3/raw/430479/x.cram",
            bucket="czi-novogene",
        )
        rollup = roll_up_wells(BUCKET, self._keys(), {("430479", "Z0001"): source})

        assert len(rollup.rows) == 1
        assert rollup.rows[0]["verdict"] == "DATA_GAP"
        assert "nothing was uploaded" not in rollup.rows[0]["detail"]


_COMPLETE_DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479"
_COMPLETE_STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
_BARE_FIVE = (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
_TRIM_FIVE = (
    ".trim.cram",
    ".trim.csv",
    ".trim.stderr",
    ".trim.stdout",
    ".trim_fail.csv",
)


def _well_keys(*suffixes):
    return [f"{_COMPLETE_DIR}/{_COMPLETE_STEM}{s}" for s in suffixes]


def _both_paths(keys):
    """(roll-up verdict, inventory-missing suffixes) for one well."""
    verdict = roll_up_wells(BUCKET, keys).rows[0]["verdict"]
    _good, lost, _found = check_expected_raw_files(keys, "seahub_sci")
    missing = sorted(
        path.rsplit("/", 1)[-1].replace(_COMPLETE_STEM, "")
        for row in lost
        for key, path in row.items()
        if key != "path"
    )
    return verdict, missing


class TestCompletenessIsNamingAgnostic:
    """A well is complete when every artifact *kind* arrived, whatever its name.

    Judging by suffix family let one optional sidecar carrying the other
    family's name decide the whole requirement set, and the two completeness
    paths reacted to that differently -- the roll-up picking ``trim``, the
    inventory requiring both sets. A well with five correct ``.trim.*`` files
    plus a stray bare sidecar was reported complete by one and missing five
    files by the other.
    """

    @pytest.mark.parametrize(
        "label,suffixes",
        [
            ("five bare", _BARE_FIVE),
            ("five trim", _TRIM_FIVE),
            ("bare plus its own sidecar", _BARE_FIVE + (".cram-metadata.json",)),
            ("bare plus a trim sidecar", _BARE_FIVE + (".trim.cram-metadata.json",)),
            ("trim plus a bare sidecar", _TRIM_FIVE + (".cram-metadata.json",)),
            ("trim plus a stray bare csv", _TRIM_FIVE + (".csv",)),
            ("both families in full", _TRIM_FIVE + _BARE_FIVE),
            (
                "required artifacts split across families",
                (".trim.cram", ".trim.csv", ".trim.stderr", ".stdout", "_fail.csv"),
            ),
        ],
    )
    def test_a_well_with_every_kind_is_complete_however_named(self, label, suffixes):
        verdict, missing = _both_paths(_well_keys(*suffixes))

        assert verdict != "DATA_GAP", label
        assert missing == [], label

    @pytest.mark.parametrize(
        "label,suffixes,absent",
        [
            ("bare well missing its csv", _BARE_FIVE, ".csv"),
            ("trim well missing its cram", _TRIM_FIVE, ".trim.cram"),
        ],
    )
    def test_a_genuinely_absent_artifact_is_still_reported(
        self, label, suffixes, absent
    ):
        verdict, missing = _both_paths(
            _well_keys(*[s for s in suffixes if s != absent])
        )

        assert verdict == "DATA_GAP", label
        # Reported under its SOP name whichever spelling the upload used.
        assert missing == [SEAHUB_BARE_TO_TRIM_SUFFIX.get(absent, absent)], label

    def test_the_two_paths_never_disagree(self):
        """The invariant the split verdicts violated."""
        for suffixes in (
            _BARE_FIVE,
            _TRIM_FIVE,
            _TRIM_FIVE + (".cram-metadata.json",),
            _BARE_FIVE + (".trim.cram-metadata.json",),
            _TRIM_FIVE[:-1],
            _BARE_FIVE[:-1],
        ):
            verdict, missing = _both_paths(_well_keys(*suffixes))
            assert (verdict == "DATA_GAP") == bool(missing), suffixes


class TestAnEmptyCramIsNotCompliant:
    """A CRAM that arrived and is 0 bytes carries no data.

    Completeness asks only whether the artifact arrived, and the trimmed-vs-vendor
    size check skips a zero because it cannot tell "empty" from "size was never
    collected" -- so an interrupted `aws s3 cp`, which is exactly the shape that
    leaves a zero-byte object behind, read as COMPLIANT with the verdict line
    saying "1 matched, 0 not trimmed".

    Only the CRAM. An empty `.stdout` or `.stderr` is ordinary -- a step that
    printed nothing -- and the shared fixture ships `.stdout` at 0 bytes for
    exactly that reason.
    """

    DIR = f"{RAW}/REF3_P05_1/430479"
    STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
    SUFFIXES = (
        ".trim.cram",
        ".trim.csv",
        ".trim.stderr",
        ".trim.stdout",
        ".trim_fail.csv",
    )

    def _keys(self) -> list[str]:
        return [f"{self.DIR}/{self.STEM}{s}" for s in self.SUFFIXES]

    def _sizes(self, **override: int) -> dict[str, int]:
        sizes = {k: 100 for k in self._keys()}
        for suffix, value in override.items():
            sizes[f"{self.DIR}/{self.STEM}{suffix.replace('_', '.', 1)}"] = value
        return sizes

    def test_a_zero_byte_cram_is_a_data_gap(self):
        sizes = {
            **{k: 100 for k in self._keys()},
            f"{self.DIR}/{self.STEM}.trim.cram": 0,
        }

        rows = roll_up_wells(BUCKET, self._keys(), None, sizes=sizes).rows

        assert [r["verdict"] for r in rows] == ["DATA_GAP"]
        assert "0 bytes" in rows[0]["detail"]

    def test_a_normal_cram_is_still_compliant(self):
        sizes = {k: 100 for k in self._keys()}

        rows = roll_up_wells(BUCKET, self._keys(), None, sizes=sizes).rows

        assert [r["verdict"] for r in rows] == ["COMPLIANT"]

    def test_an_empty_stdout_is_not_a_gap(self):
        """A step that printed nothing is normal; the fixture ships this shape."""
        sizes = {
            **{k: 100 for k in self._keys()},
            f"{self.DIR}/{self.STEM}.trim.stdout": 0,
            f"{self.DIR}/{self.STEM}.trim.stderr": 0,
        }

        rows = roll_up_wells(BUCKET, self._keys(), None, sizes=sizes).rows

        assert [r["verdict"] for r in rows] == ["COMPLIANT"]

    def test_unknown_sizes_are_not_read_as_empty(self):
        """Manifest mode collects no sizes; absent must not mean zero."""
        assert [r["verdict"] for r in roll_up_wells(BUCKET, self._keys()).rows] == [
            "COMPLIANT"
        ]
        rows = roll_up_wells(BUCKET, self._keys(), None, sizes={}).rows
        assert [r["verdict"] for r in rows] == ["COMPLIANT"]

    def test_the_bare_family_cram_counts_too(self):
        keys = [
            f"{self.DIR}/{self.STEM}{s}"
            for s in (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
        ]
        sizes = {**{k: 100 for k in keys}, f"{self.DIR}/{self.STEM}.cram": 0}

        rows = roll_up_wells(BUCKET, keys, None, sizes=sizes).rows

        assert [r["verdict"] for r in rows] == ["DATA_GAP"]

    def test_it_fires_without_a_vendor_delivery(self):
        """The size check needs a vendor to compare against; this does not."""
        sizes = {
            **{k: 100 for k in self._keys()},
            f"{self.DIR}/{self.STEM}.trim.cram": 0,
        }

        rows = roll_up_wells(BUCKET, self._keys(), None, sizes=sizes).rows

        assert rows[0]["verdict"] == "DATA_GAP"
        assert "the vendor delivered it as" not in rows[0]["detail"]


class TestTheTwoIndexesKeyWellsTheSameWay:
    """``index_trimmed_upload`` and ``roll_up_wells`` must agree on identity.

    One keyed on the raw stem and the other on the doubled-wafer-normalized stem.
    That is the same answer for every real shape, and differs for exactly one:
    a doubled wafer with no group at all, which parses *before* normalization
    (wafer-group-UG-barcode) and not after (wafer-UG-barcode is three fields, not
    four). The roll-up then failed to identify a well the vendor index had
    matched, and emitted a "nothing was uploaded" row for a well that was
    uploaded -- beside the real row for the same five objects.

    Synthetic, and it failed safe by over-reporting, but Q1 promoted that row to
    DATA_GAP and DATA_GAP rows go to errors.txt, so a phantom one now costs an
    operator a lookup.
    """

    STEM = "123456-123456-Z0001-ACGTACGT"
    KEYS = [f"{RAW}/P01/123456/123456-123456-Z0001-ACGTACGT{s}" for s in TRIM_SUFFIXES]
    VENDOR = {
        ("123456", "Z0001"): SourceEntry(
            wafer="123456",
            ug="Z0001",
            barcode="ACGTACGT",
            group="REF3_P01_A1",
            assay="GEX_hash_oligo",
            cram_key=f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/123456/x.cram",
        )
    }

    def test_both_indexes_produce_the_same_identity(self):
        from qa_seahub_source import index_trimmed_upload

        indexed = set(index_trimmed_upload(self.KEYS, sizes={}))
        rolled = {(r["wafer"], r["ug"]) for r in roll_up_wells(BUCKET, self.KEYS).rows}

        assert indexed == rolled == {("123456", "Z0001")}

    def test_no_phantom_nothing_was_uploaded_row(self):
        rollup = roll_up_wells(BUCKET, self.KEYS, self.VENDOR)

        assert len(rollup.rows) == 1
        assert "nothing was uploaded" not in rollup.rows[0]["detail"]
        assert rollup.rows[0]["objects"] == len(self.KEYS)

    def test_the_well_still_reports_its_real_problem(self):
        """Identified, but with no group there is no corrected name to derive."""
        rollup = roll_up_wells(BUCKET, self.KEYS, self.VENDOR)

        assert rollup.rows[0]["verdict"] == "UNKNOWN"


class TestUploadScopeDoesNotBlackOutTheRenamePlan:
    """A wrong bucket is one finding, not a veto on every rename in the upload.

    The rename self-check fed each proposal back through ``validate_seahub_key``,
    which also runs ``bad_bucket`` and ``lab_project_mismatch``. Both are
    invariant under the proposed move -- same bucket, same project -- so a
    proposal that strictly reduces the violation set was withheld for a reason it
    could not possibly fix. Measured on this fixture before the filter: 20
    moveable objects to 0, and the wells from
    {COMPLIANT:1, RENAMEABLE:3, DATA_GAP:1, UNKNOWN:1} to {DATA_GAP:1, UNKNOWN:5}.
    """

    def _plan(self, bucket: str):
        keys = ref3_trimmed_keys()
        index = _vendor_index()
        mapping = build_rename_mapping(bucket, keys, index)
        rollup = roll_up_wells(bucket, keys, index)
        return mapping, rollup

    # "czi-wrong" is a well-formed bucket under the wrong lab, so it trips
    # lab_project_mismatch; "notczi-x" is not czi-<lab> at all, so it trips
    # bad_bucket. Both are upload-scope, and both used to veto the plan.
    @pytest.mark.parametrize("bad_bucket", ["czi-wrong", "notczi-x"])
    def test_a_wrong_bucket_changes_nothing_about_the_plan(self, bad_bucket):
        good, good_wells = self._plan(BUCKET)
        bad, bad_wells = self._plan(bad_bucket)

        assert dict(bad.counts) == dict(good.counts)
        assert len(bad.moveable()) == len(good.moveable())
        assert rollup_summary(bad_wells.rows) == rollup_summary(good_wells.rows)

    @pytest.mark.parametrize("bad_bucket", ["czi-wrong", "notczi-x"])
    def test_a_wrong_bucket_still_leaves_wells_moveable(self, bad_bucket):
        bad, bad_wells = self._plan(bad_bucket)

        assert len(bad.moveable()) > 0
        assert rollup_summary(bad_wells.rows)["RENAMEABLE"] > 0

    @pytest.mark.parametrize("bad_bucket", ["czi-wrong", "notczi-x"])
    def test_no_proposal_is_withheld_for_an_upload_scope_rule(self, bad_bucket):
        bad, _ = self._plan(bad_bucket)

        withheld = {
            u
            for row in bad.rows
            for u in row["unresolved"].split("|")
            if u.startswith("proposal_not_sop:")
        }
        assert not any("bad_bucket" in u for u in withheld)
        assert not any("lab_project_mismatch" in u for u in withheld)

    @pytest.mark.parametrize(
        "bad_bucket,rule",
        [("czi-wrong", "lab_project_mismatch"), ("notczi-x", "bad_bucket")],
    )
    def test_the_sop_table_still_reports_it(self, bad_bucket, rule):
        """Filtered only inside the rename gate; validate_seahub_key is unchanged."""
        violations = validate_seahub_stems(bad_bucket, ref3_trimmed_keys())

        rows = [v for v in violations if v.type == rule]
        assert rows, f"{rule} should still be reported for {bad_bucket}"
        assert all(v.scope == "upload" for v in rows)

    def test_a_non_upload_defect_still_withholds_the_proposal(self):
        """The gate keeps its job for facts a move can actually change."""
        residual = validate_seahub_key(
            BUCKET, f"{RAW}/P05_1/436830/not-a-sop-name.trim.cram"
        )

        assert residual
        assert any(v.scope != "upload" for v in residual)


class TestTheTwoHeadlineCsvsAgree:
    """The rename CSV and the well-status CSV must not contradict each other.

    They answer one question -- can this be fixed by moving objects -- and used
    to answer it two ways: ``expected_trimmed_key`` repairs a
    ``sublibrary_mismatch`` from the vendor group and returns a complete
    destination, while the roll-up tested the defect *names* against
    ``SEAHUB_RENAMEABLE_SOP_TYPES``, which did not list it. The operator was
    handed moves for a well the status CSV called impossible. The roll-up now
    asks ``RenameProposal.renameable``, so both read the same predicate.
    """

    MISMATCH_KEY = (
        f"{RAW}/P07_1/438514/"
        "438514-REF3_P09_9_A3_GEX_hash_oligo-Z0301-CAGCTCGAATGCGAA.trim.cram"
    )
    MISMATCH_VENDOR = SourceEntry(
        wafer="438514",
        ug="Z0301",
        barcode="CAGCTCGAATGCGAA",
        group="REF3_P07_1_A3",
        assay="GEX_hash_oligo",
        cram_key=f"{PROJECT}/{VENDOR_ORDER}/REF3/raw/438514/x.cram",
        bucket="czi-novogene",
    )

    def test_mismatch_with_a_vendor_group_is_renameable(self):
        proposal = expected_trimmed_key(BUCKET, self.MISMATCH_KEY, self.MISMATCH_VENDOR)

        assert proposal.unresolved == ()
        assert proposal.defects == ("sublibrary_mismatch",)
        assert proposal.renameable is True
        assert proposal.expected_s3_uri.endswith(
            "/REF3/raw/REF3_P07_1/438514/"
            "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0301-CAGCTCGAATGCGAA.trim.cram"
        )

    def test_mismatch_without_a_vendor_group_stays_unresolved(self):
        """The control: with no vendor there is no authority, so no proposal."""
        proposal = expected_trimmed_key(BUCKET, self.MISMATCH_KEY, None)

        assert proposal.unresolved == ("sublibrary_mismatch",)
        assert proposal.expected_s3_uri == ""
        assert proposal.renameable is False

    def test_a_mismatched_well_rolls_up_as_renameable(self):
        keys = [self.MISMATCH_KEY.replace(".trim.cram", s) for s in TRIM_SUFFIXES]
        index = {("438514", "Z0301"): self.MISMATCH_VENDOR}

        rollup = roll_up_wells(BUCKET, keys, index)

        assert [r["verdict"] for r in rollup.rows] == ["RENAMEABLE"]

    def test_a_mismatched_well_without_a_vendor_stays_unknown(self):
        keys = [self.MISMATCH_KEY.replace(".trim.cram", s) for s in TRIM_SUFFIXES]

        rollup = roll_up_wells(BUCKET, keys, None)

        assert [r["verdict"] for r in rollup.rows] == ["UNKNOWN"]

    @pytest.mark.parametrize("with_vendor", [True, False])
    def test_no_moveable_object_belongs_to_an_unknown_well(self, with_vendor):
        """The general invariant, in the direction the old test could not see.

        ``test_a_renameable_well_carries_only_repairable_defects`` reads
        well -> defects. This reads object -> well, which is the direction the
        contradiction appeared in.
        """
        index = _vendor_index() if with_vendor else None
        # All five artifacts, or the well is DATA_GAP and never reaches the
        # renameable test this invariant is about.
        keys = ref3_trimmed_keys() + [
            self.MISMATCH_KEY.replace(".trim.cram", s) for s in TRIM_SUFFIXES
        ]
        full_index = dict(index or {})
        full_index[("438514", "Z0301")] = self.MISMATCH_VENDOR

        mapping = build_rename_mapping(BUCKET, keys, full_index)
        verdicts = {
            (r["wafer"], r["ug"]): r["verdict"]
            for r in roll_up_wells(BUCKET, keys, full_index).rows
        }

        for row in mapping.moveable():
            identity = (row["wafer"], row["ug"])
            if identity not in verdicts:
                continue  # junk and odd-depth objects belong to no well
            assert verdicts[identity] != "UNKNOWN", (
                f"{row['current_s3_uri']} is offered a destination but its well "
                f"is UNKNOWN"
            )


def _module_level_constants(tree: ast.Module) -> set[str]:
    """Public constant names bound at a module's top level.

    Descends into top-level ``if``/``try``, since a name bound under a version
    or import guard is still module level, and flattens tuple targets. Anything
    lowercase is a variable by convention, not a constant, and anything
    underscore-prefixed is private.
    """

    def statements(body):
        for node in body:
            yield node
            if isinstance(node, ast.If):
                yield from statements(node.body)
                yield from statements(node.orelse)
            elif isinstance(node, ast.Try):
                yield from statements(node.body + node.orelse + node.finalbody)
                for handler in node.handlers:
                    yield from statements(handler.body)

    def bound(target):
        if isinstance(target, ast.Name):
            yield target.id
        elif isinstance(target, (ast.Tuple, ast.List)):
            for element in target.elts:
                yield from bound(element)

    names: set[str] = set()
    for node in statements(tree.body):
        targets = (
            [node.target]
            if isinstance(node, ast.AnnAssign)
            else node.targets
            if isinstance(node, ast.Assign)
            else []
        )
        names.update(n for t in targets for n in bound(t))
    return {n for n in names if not n.startswith("_") and not n.islower()}


def _identifiers_referenced_outside(bcp: pathlib.Path, exclude: str) -> set[str]:
    """Every identifier used by code under ``bcp``, ignoring one file.

    By AST, so a name that appears only in a comment or a docstring does not
    count as used -- including the docstrings of the guards below, which name
    the very constants they exist to catch and so exempted them permanently when
    this was a text search.

    ``qa.ipynb`` is scanned too: the notebook is the primary consumer of these
    modules, so a constant used only there is live. A cell that will not parse
    (a magic, a fragment) falls back to a word scan, which over-counts -- the
    safe direction, since the cost of a false positive here is a failing build
    on correct code.
    """
    names: set[str] = set()

    def collect(source: str) -> bool:
        try:
            tree = ast.parse(source)
        except SyntaxError:
            return False
        for node in ast.walk(tree):
            if isinstance(node, ast.Name):
                names.add(node.id)
            elif isinstance(node, ast.Attribute):
                names.add(node.attr)
            elif isinstance(node, ast.alias):
                names.add(node.name.split(".")[-1])
                if node.asname:
                    names.add(node.asname)
        return True

    for path in bcp.rglob("*.py"):
        if path.name != exclude:
            collect(path.read_text())

    notebook = bcp / "qa.ipynb"
    if notebook.exists():
        for cell in json.loads(notebook.read_text()).get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            if not collect(source):
                names.update(re.findall(r"\b\w+\b", source))
    return names


def _identifiers_referenced_by_production(bcp: pathlib.Path, exclude: str) -> set[str]:
    """Like the above, but blind to ``tests/``.

    A function nothing but a test calls is dead API however thoroughly it is
    tested. Constants are different -- ``WELL_COLUMNS`` is a CSV contract whose
    only legitimate reader is the test that pins it -- which is why the two
    guards below apply different rules rather than one.
    """
    names: set[str] = set()
    for path in bcp.rglob("*.py"):
        if path.name == exclude or "tests" in path.parts:
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Name):
                names.add(node.id)
            elif isinstance(node, ast.Attribute):
                names.add(node.attr)
            elif isinstance(node, ast.alias):
                names.add(node.name.split(".")[-1])
    notebook = bcp / "qa.ipynb"
    if notebook.exists():
        for cell in json.loads(notebook.read_text()).get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            try:
                tree = ast.parse(source)
            except SyntaxError:
                names.update(re.findall(r"\b\w+\b", source))
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Name):
                    names.add(node.id)
                elif isinstance(node, ast.Attribute):
                    names.add(node.attr)
                elif isinstance(node, ast.alias):
                    names.add(node.name.split(".")[-1])
    return names


class TestRuleVocabulary:
    """Makes the closed vocabularies do what their comments claim.

    ``SEAHUB_SOP_RULES`` says it is "kept explicit so a typo in a new rule shows
    up as a test failure"; nothing asserted that, so a misspelled rule type would
    have shipped silently. These walk the real fixture listings and check every
    value actually emitted against its declared vocabulary.
    """

    def test_no_seahub_constant_is_referenced_only_by_itself(self):
        """A constant nobody reads is a claim about the code that has stopped
        being true.

        ``SEAHUB_BARE_EXPECTED`` went dead when completeness moved to artifact
        kinds and sat there describing a required set nothing consulted, right
        beside ``SEAHUB_BARE_OPTIONAL``, which is live -- so its presence read as
        deliberate. Names derived inside ``qa_constants`` itself do not count as
        used: that is how the dead one looked alive.
        """
        bcp = pathlib.Path(qa_seahub_rename.__file__).parent
        constants = bcp / "qa_constants.py"
        declared = {
            name
            for name in _module_level_constants(ast.parse(constants.read_text()))
            if name.startswith("SEAHUB_")
        }
        assert declared, "no SEAHUB_* constants found"

        referenced = _identifiers_referenced_outside(bcp, exclude="qa_constants.py")
        unread = sorted(declared - referenced)
        assert unread == [], f"declared but read nowhere: {unread}"

    def test_no_seahub_module_function_is_production_dead(self):
        """A public function only the tests call is dead API.

        ``grab_seahub_trim_fail_csv`` and ``index_untrimmed_source`` were both
        thin wrappers with no caller outside the suite, and the constant guards
        could not see them: one is scoped to constants in ``qa_constants``, and
        the other counts a test reference as a use, which for a constant it
        should.

        Scoped to ``qa_seahub_*`` deliberately. ``qa_mods`` carries at least one
        such function that predates this work, and failing the build on
        unrelated code would make the guard something to switch off.
        """
        bcp = pathlib.Path(qa_seahub_rename.__file__).parent
        modules = sorted(bcp.glob("qa_seahub_*.py"))
        assert modules, "no qa_seahub_* modules found"

        unread: list[str] = []
        for module in modules:
            tree = ast.parse(module.read_text())
            public = {
                node.name
                for node in tree.body
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and not node.name.startswith("_")
            }
            called = {
                node.id
                for node in ast.walk(tree)
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
            }
            reachable = _identifiers_referenced_by_production(bcp, module.name)
            unread += [
                f"{module.name}:{name}" for name in sorted(public - called - reachable)
            ]

        assert unread == [], f"public but reachable only from tests: {unread}"

    def test_no_seahub_module_constant_is_unread(self):
        """The same guard, for the modules rather than the constants file.

        The test above is scoped to ``SEAHUB_*`` names in ``qa_constants``, so it
        could not see ``WELL_COLUMNS`` sitting dead a few lines from
        ``RENAME_COLUMNS``, which is pinned. The rule has to be looser here: a
        module constant used only inside its own module is doing its job, which
        is not true of a vocabulary in ``qa_constants``. So a name counts as read
        if it is *loaded* anywhere in its own module -- by AST, so a mention in a
        comment or docstring does not rescue it -- or referred to in any other
        file. Leading-underscore names are private and exempt.
        """
        bcp = pathlib.Path(qa_seahub_rename.__file__).parent
        modules = sorted(bcp.glob("qa_seahub_*.py"))
        assert modules, "no qa_seahub_* modules found"

        unread: list[str] = []
        for module in modules:
            tree = ast.parse(module.read_text())
            declared = _module_level_constants(tree)
            loaded = {
                node.id
                for node in ast.walk(tree)
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
            }
            referenced = _identifiers_referenced_outside(bcp, exclude=module.name)
            unread += [
                f"{module.name}:{name}"
                for name in sorted(declared - loaded - referenced)
            ]

        assert unread == [], f"declared but read nowhere: {unread}"

    def test_renameable_types_are_a_subset_of_the_rules(self):
        assert SEAHUB_RENAMEABLE_SOP_TYPES <= SEAHUB_SOP_RULES

    def test_renameable_set_matches_what_proposals_carry(self):
        """The set is documentation now, so something must keep it honest.

        The roll-up asks ``RenameProposal.renameable`` rather than this set, so
        nothing in production would notice the set going stale -- which is how it
        came to omit ``sublibrary_mismatch``. Read the defect literals straight
        out of the module, as ``TestSopRuleVocabulary`` does for rule types.
        """
        source = pathlib.Path(qa_seahub_rename.__file__).read_text()
        appended = {
            node.args[0].value
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Call)
            and getattr(node.func, "attr", None) == "append"
            and getattr(node.func.value, "id", None) == "defects"
            and node.args
            and isinstance(node.args[0], ast.Constant)
        }

        assert appended, "no defects.append('...') literals found"
        # non_sequencing_artifact is set on a junk object, which belongs to no
        # well and so never reaches the roll-up's renameable test.
        assert appended - {"non_sequencing_artifact"} == SEAHUB_RENAMEABLE_SOP_TYPES

    def test_every_emitted_violation_type_is_a_declared_rule(self):
        violations = validate_seahub_stems(BUCKET, ref3_trimmed_keys())

        assert {v.type for v in violations} <= SEAHUB_SOP_RULES

    def test_every_emitted_scope_is_a_declared_scope(self):
        violations = validate_seahub_stems(BUCKET, ref3_trimmed_keys())

        assert {v.scope for v in violations} <= set(SEAHUB_VIOLATION_SCOPES)

    def test_a_wholly_misnamed_upload_stays_inside_both_vocabularies(self):
        """The newest rules are the likeliest to be misspelled."""
        keys = [
            f"{RAW}/P04_1/437120/437120-REF3_P04_1_A1_GEX_hash_oligo"
            f"-Z000{i}-CAGCTCGAATGCGAT.trimmed.ucram"
            for i in range(1, 4)
        ]

        violations = validate_seahub_stems(BUCKET, keys)

        assert {v.type for v in violations} <= SEAHUB_SOP_RULES
        assert {v.scope for v in violations} <= set(SEAHUB_VIOLATION_SCOPES)

    def test_every_emitted_rename_status_is_declared(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())

        assert {r["status"] for r in mapping.rows} <= set(SEAHUB_RENAME_STATUSES)

    def test_every_emitted_name_source_is_declared(self):
        mapping = build_rename_mapping(BUCKET, ref3_trimmed_keys(), _vendor_index())
        sources = {r["name_source"] for r in mapping.rows if r["name_source"]}

        assert sources <= set(SEAHUB_RENAME_NAME_SOURCES)


class TestAWholeUntrimmedWaferIsADataGap:
    """The other half of the pin in test_qa_seahub_source.

    A delivered wafer with nothing trimmed is a gap per well, and ``DATA_GAP``
    is the only verdict the notebook writes to errors.txt -- so this is the path
    on which a forgotten plate is either audible or silent. Pinned on the
    current order-prefix behaviour, before vendor prefixes start being located
    from the trimmed upload's wafer list, which would drop these wells out of
    the index that seeds them.
    """

    def test_its_wells_roll_up_as_data_gaps(self):
        rollup = roll_up_wells(
            BUCKET, ref3_trimmed_keys(), _vendor_index(extra_wafer=True)
        )

        rows = [r for r in rollup.rows if r["wafer"] == UNTRIMMED_WAFER]
        assert [(r["ug"], r["verdict"]) for r in rows] == [
            ("Z0500", "DATA_GAP"),
            ("Z0501", "DATA_GAP"),
        ]

    def test_it_adds_to_the_gap_count_rather_than_replacing_it(self):
        """Wafer 438514's absent CRAM is the pre-existing gap; these are extra."""
        rollup = roll_up_wells(
            BUCKET, ref3_trimmed_keys(), _vendor_index(extra_wafer=True)
        )

        assert rollup_summary(rollup.rows)["DATA_GAP"] == 3

    def test_each_row_names_the_vendor_key_it_can_be_re_trimmed_from(self):
        """The action is 're-run trim on this', not 'find out what is missing'."""
        rollup = roll_up_wells(
            BUCKET, ref3_trimmed_keys(), _vendor_index(extra_wafer=True)
        )

        for row in (r for r in rollup.rows if r["wafer"] == UNTRIMMED_WAFER):
            assert "nothing was uploaded" in row["detail"]
            assert "the vendor delivered it as" in row["detail"]
            assert row["missing_artifacts"].split("|")
            assert row["sublibrary"] == "REF3_P08_1"

    def test_the_wafer_is_absent_from_the_trimmed_upload_entirely(self):
        """Guards the fixture: this must be a whole-wafer gap, not a partial one."""
        assert not [k for k in ref3_trimmed_keys() if UNTRIMMED_WAFER in k]


class TestADoubledVendorWaferDoesNotCostARepair:
    """Both readers of a vendor group must strip a leading wafer token.

    index_untrimmed_sources does not normalize a doubled vendor wafer, so the
    group can arrive as ``438514-REF3_P07_1_A3``. The gap-row path stripped it;
    the sublibrary_mismatch repair did not, so it proposed a key that re-validated
    as duplicated_wafer_token and was withheld -- a well that could have been
    renamed rolled up UNKNOWN instead. The self-check meant no bad key was ever
    proposed, so this cost a repair rather than producing a wrong one.
    """

    KEY = (
        f"{RAW}/OTHER_P09/438514/"
        "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT.trim.cram"
    )
    DOUBLED = SourceEntry(
        wafer="438514",
        ug="Z0305",
        barcode="CACACACAACATGAT",
        group="438514-REF3_P07_1_A3",
        assay="GEX_hash_oligo",
        cram_key="v/REF3/raw/438514/x.cram",
        bucket="czi-novogene",
    )

    def test_the_mismatch_repair_is_proposed_rather_than_withheld(self):
        proposal = expected_trimmed_key(BUCKET, self.KEY, self.DOUBLED)

        assert proposal.defects == ("sublibrary_mismatch",)
        assert proposal.unresolved == ()
        assert proposal.renameable is True
        assert proposal.expected_s3_uri.endswith(
            "/REF3_P07_1/438514/"
            "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT.trim.cram"
        )

    def test_a_vendor_group_with_no_doubled_wafer_is_unchanged(self):
        plain = replace(self.DOUBLED, group="REF3_P07_1_A3")
        proposal = expected_trimmed_key(BUCKET, self.KEY, plain)

        assert proposal.renameable is True
        assert proposal.expected_s3_uri == (
            expected_trimmed_key(BUCKET, self.KEY, self.DOUBLED).expected_s3_uri
        )
