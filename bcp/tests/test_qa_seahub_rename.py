"""
Tests for corrected-name composition and the per-well roll-up (qa_seahub_rename).
"""

from __future__ import annotations

import ast
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
)
from qa_seahub_rename import (
    RENAME_COLUMNS,
    build_rename_mapping,
    expected_trimmed_key,
    roll_up_wells,
    rollup_summary,
)
from qa_seahub_sop import validate_seahub_stems
from qa_seahub_source import SourceEntry, index_untrimmed_sources

from tests.qa_seahub_helpers import (
    BUCKET,
    JUNK_NAMES,
    PROJECT,
    RAW,
    TRIM_SUFFIXES,
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
            target.id
            for node in ast.parse(constants.read_text()).body
            for target in (
                [node.target]
                if isinstance(node, ast.AnnAssign)
                else node.targets
                if isinstance(node, ast.Assign)
                else []
            )
            if isinstance(target, ast.Name) and target.id.startswith("SEAHUB_")
        }
        assert declared, "no SEAHUB_* constants found"

        elsewhere = "".join(
            p.read_text()
            for p in [*bcp.glob("*.py"), *(bcp / "tests").rglob("*.py")]
            if p.name != "qa_constants.py"
        )
        unread = sorted(n for n in declared if not re.search(rf"\b{n}\b", elsewhere))
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
