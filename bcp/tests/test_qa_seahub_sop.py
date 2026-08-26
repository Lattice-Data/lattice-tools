"""
Tests for SeaHub SOP filename/path validation and family-aware completeness.

Cases are drawn from two real S3 listings reviewed alongside these changes: a
trimmed upload (2035 objects) in which one sublibrary follows the SOP and six
carry trim outputs with the ``.trim`` infix dropped, and a vendor delivery
(5389 objects).  Lab and order identifiers are sanitized here.
"""

from __future__ import annotations

import ast
import os
import pathlib

import pytest

from qa_checks import check_expected_raw_files, check_extra_raw_files
from qa_constants import (
    SEAHUB_SOP_RULES,
    SEAHUB_STEM_RE,
    SEAHUB_VIOLATION_SCOPES,
)
import qa_seahub_sop
from qa_gather import gather_qa_data
from qa_mods import (
    apply_seahub_trim_fail_blocks,
    parse_seahub_trim_fail_csv,
    seahub_stem_and_family,
)
from qa_seahub_sop import (
    _folder_candidates,
    _match_folder_candidates,
    group_seahub_keys,
    seahub_group_parts,
    sop_violation_summary,
    unrecognized_suffix,
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
# The same sublibrary with the redundant ExperimentID prefix elided from the
# folder, which is the shape every real trimmed upload uses.
GOOD_ELIDED_PREFIX = (
    "labalpha-seahub-bcp/REF3/raw/P05_2/436830/"
    "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
)
# A real object from the upload: three independent defects at once.
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

    def test_a_bare_suffix_needs_a_stem_that_parses(self):
        """Measured on a real 480-well upload that reported as 960 wells.

        ``.csv`` is generic, so ``<well>.trimmer_stats.csv`` matched the bare
        family and left ``<well>.trimmer_stats`` as a stem -- a second, phantom
        well alongside ``<well>.failure_codes``.
        """
        assert seahub_stem_and_family(f"{TRIM_STEM}.trimmer_stats.csv") is None
        assert seahub_stem_and_family(f"{TRIM_STEM}.failure_codes.csv") is None

    def test_a_doubled_wafer_bare_stem_still_parses(self):
        """The guard must not reject the real bare family it exists to support."""
        assert seahub_stem_and_family(f"{BARE_STEM}.cram") is not None

    def test_the_trim_family_is_not_gated(self):
        """`.trim.*` is distinctive, so a malformed stem stays a stem defect."""
        assert seahub_stem_and_family("not-a-seahub-name.trim.cram") == (
            "not-a-seahub-name",
            ".trim.cram",
            "trim",
        )


class TestUnrecognizedSuffix:
    def test_reads_the_extension_after_the_barcode(self):
        assert unrecognized_suffix(f"{TRIM_STEM}.trimmed.ucram") == ".trimmed.ucram"

    def test_a_name_with_no_hyphen_falls_back_to_the_first_dot(self):
        assert unrecognized_suffix("index.html") == ".html"

    def test_a_name_with_no_dot_reports_empty(self):
        assert unrecognized_suffix("436516") == ""


class TestUnknownSuffixReporting:
    """An upload that misspells its artifacts misspells every one of them."""

    def _keys(self, suffix, count):
        return [
            f"{TRIM_DIR}/430479-REF3_P05_1_A{i}_GEX_hash_oligo"
            f"-Z{i:04d}-CAGTCAGTTGCAGAT{suffix}"
            for i in range(1, count + 1)
        ]

    def test_one_row_per_distinct_suffix_not_per_object(self):
        keys = self._keys(".trimmed.ucram", 20)

        summary = sop_violation_summary(validate_seahub_stems("czi-labalpha", keys))

        assert summary["unexpected_suffix"] == 1

    def test_the_row_names_the_spelling_the_count_and_the_sop_artifacts(self):
        keys = self._keys(".trimmed.ucram", 20)

        row = next(
            v
            for v in validate_seahub_stems("czi-labalpha", keys)
            if v.type == "unexpected_suffix"
        )

        assert row.scope == "suffix"
        assert "20 object(s)" in row.detail
        assert "'.trimmed.ucram'" in row.detail
        assert ".trim.cram" in row.detail

    def test_distinct_spellings_report_separately(self):
        keys = self._keys(".trimmed.ucram", 5) + self._keys(".trimmer_stats.csv", 5)

        rows = [
            v
            for v in validate_seahub_stems("czi-labalpha", keys)
            if v.type == "unexpected_suffix"
        ]

        assert len(rows) == 2

    def test_an_upload_with_no_recognizable_artifact_fails_loudly(self):
        keys = self._keys(".trimmed.ucram", 5)

        rows = [
            v
            for v in validate_seahub_stems("czi-labalpha", keys)
            if v.type == "no_recognized_artifacts"
        ]

        assert len(rows) == 1
        assert rows[0].scope == "upload"
        assert "no well could be identified" in rows[0].detail

    def test_it_does_not_fire_when_some_artifact_is_recognized(self):
        keys = self._keys(".trimmed.ucram", 5) + [f"{TRIM_DIR}/{TRIM_STEM}.trim.cram"]

        summary = sop_violation_summary(validate_seahub_stems("czi-labalpha", keys))

        assert "no_recognized_artifacts" not in summary

    def test_an_empty_listing_is_not_a_failure(self):
        assert validate_seahub_stems("czi-labalpha", []) == []


class TestRuleSetIsClosed:
    """`SEAHUB_SOP_RULES` claims a typo in a new rule shows up as a test failure.

    Asserting only over the *emitted* types cannot deliver that: the fixtures do
    not exercise every rule, so a misspelled ``wafer_mismatch`` still passes.
    Read the literals straight out of the module instead, which covers rules no
    fixture reaches.
    """

    def _declared_literals(self, keyword: str) -> set[str]:
        source = pathlib.Path(qa_seahub_sop.__file__).read_text()
        found: set[str] = set()
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.Call):
                continue
            if getattr(node.func, "id", None) != "SopViolation":
                continue
            for kw in node.keywords:
                if kw.arg == keyword and isinstance(kw.value, ast.Constant):
                    found.add(kw.value.value)
        return found

    def test_every_rule_the_module_can_emit_is_declared(self):
        declared = self._declared_literals("type")

        assert declared, "no SopViolation(type=...) literals found"
        assert declared <= SEAHUB_SOP_RULES

    def test_no_declared_rule_is_unreachable(self):
        """The other direction: a rule nothing emits is dead vocabulary."""
        assert SEAHUB_SOP_RULES <= self._declared_literals("type")

    def test_every_scope_the_module_can_emit_is_declared(self):
        # "stem" is the dataclass default, so it need not appear as a literal.
        assert self._declared_literals("scope") <= set(SEAHUB_VIOLATION_SCOPES)


class TestGroupPartsCandidateOrder:
    """Candidate order must not affect the result, in any shape.

    It used to decide a verdict -- try the folder as-is first, or a correct
    folder reported as truncated. Then, once both spellings were accepted, it
    still decided which prefix got stripped, and in one shape it decided wrongly:
    a compliant object reported bad_well because the shorter candidate was tried
    first. The matcher now picks by the leftover rather than by position, so the
    order genuinely does not matter -- which is what these tests pin.
    """

    @staticmethod
    def _reversed(group, sublibrary, experiment_id):
        """The production matcher, fed the production candidates in reverse.

        Reversing the real list rather than re-implementing it: a hand-written
        mirror silently stops testing the ordering property the moment
        ``_folder_candidates`` grows a third entry or ``_match_folder_candidates``
        changes how it strips a prefix.
        """
        candidates = list(reversed(_folder_candidates(sublibrary, experiment_id)))
        return _match_folder_candidates(group, candidates, sublibrary)

    @pytest.mark.parametrize(
        "experiment_id,sublibrary,group",
        [
            # The three real shapes, in both folder spellings.
            ("REF3", "REF3_P05_2", "REF3_P05_2_A10"),
            ("REF3", "P05_2", "REF3_P05_2_A10"),
            ("CHEM16", "P03", "CHEM16_P03_A1"),
            ("CHEM3-R100", "R100E", "R100E"),
            # A genuine mismatch, and a trailing token that is not a well.
            ("REF3", "P05_2", "SOMETHING_ELSE_A10"),
            ("REF3", "P05", "REF3_P05_1_A10"),
        ],
    )
    def test_the_order_is_immaterial_on_every_realistic_input(
        self, experiment_id, sublibrary, group
    ):
        assert seahub_group_parts(group, sublibrary, experiment_id) == self._reversed(
            group, sublibrary, experiment_id
        )

    def test_candidate_order_does_not_affect_the_result(self):
        """The one shape where both candidates can explain the group.

        Folder ``P10`` under ExperimentID ``P10_2`` makes the full form
        ``P10_2_P10``, which itself starts with ``P10_``. This object is
        compliant: sublibrary ``P10_2_P10``, well ``A1``. Picking the first
        prefix match instead returned trailing ``2_P10_A1`` and reported
        ``bad_well`` on it, so the order was deciding correctness rather than
        breaking a tie. No real upload has this shape.
        """
        expected = ("P10_2_P10", "A1", True)

        assert seahub_group_parts("P10_2_P10_A1", "P10", "P10_2") == expected
        assert self._reversed("P10_2_P10_A1", "P10", "P10_2") == expected

    def test_the_compliant_degenerate_object_reports_no_bad_well(self):
        """The defect the order used to produce, at the rule level."""
        key = (
            "labalpha-seahub-bcp/P10_2/raw/P10/432640/"
            "432640-P10_2_P10_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram"
        )

        assert "bad_well" not in _types(validate_seahub_key("czi-labalpha", key))


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
    def test_reports_exactly_the_three_defects(self):
        """Its folder elides the ExperimentID prefix, which is not a defect."""
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        assert _types(violations) == {
            "missing_trim_infix",
            "duplicated_wafer_token",
            "invalid_sublibrary_type",
        }

    def test_expected_name_restores_the_trim_infix(self):
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        infix = next(v for v in violations if v.type == "missing_trim_infix")
        assert infix.expected_name.endswith(".trim.cram")

    def test_duplicated_wafer_expected_name_drops_the_repeat(self):
        violations = validate_seahub_key("czi-labalpha", REAL_MISNAMED)
        doubled = next(v for v in violations if v.type == "duplicated_wafer_token")
        assert doubled.expected_name == (
            "436830-REF3_P05_2_A10-Z0169-CTCGCAATAGATGAT.cram"
        )

    def test_vendor_index_fills_the_missing_sublibrary_type(self):
        violations = validate_seahub_key(
            "czi-labalpha",
            REAL_MISNAMED,
            assay_by_identity={("436830", "Z0169"): "GEX_hash_oligo"},
        )
        typed = next(v for v in violations if v.type == "invalid_sublibrary_type")
        assert typed.expected_name == (
            "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.cram"
        )
        assert "vendor" in typed.detail


class TestPathRules:
    def test_bad_bucket(self):
        violations = validate_seahub_key("labalpha-data", GOOD_LABALPHA)
        assert "bad_bucket" in _types(violations)

    def test_lab_project_mismatch(self):
        violations = validate_seahub_key("czi-labbeta", GOOD_LABALPHA)
        assert "lab_project_mismatch" in _types(violations)


class TestAnUnrecognizedTypeTokenIsNotAMissingOne:
    """The relaxed pattern matches either way; the report should not.

    A stem whose type token is outside SEAHUB_SUBLIBRARY_TYPES folds it into the
    group, so the rule said the stem "carries no sublibrary type" about a stem
    that plainly carries one -- and with a vendor match it proposed a corrected
    name with the unrecognised token still buried in the sublibrary.
    """

    DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P01/438514"
    UNRECOGNIZED = "438514-REF3_P01_ATAC-Z0305-CACACACAACATGAT"
    NO_TYPE = "438514-REF3_P01_A3-Z0305-CACACACAACATGAT"
    VENDOR = {("438514", "Z0305"): "GEX_hash_oligo"}

    def _row(self, stem: str, vendor=None):
        rows = [
            v
            for v in validate_seahub_key(
                "czi-labalpha",
                f"{self.DIR}/{stem}.trim.cram",
                assay_by_identity=vendor,
            )
            if v.type == "invalid_sublibrary_type"
        ]
        assert len(rows) == 1, rows
        return rows[0]

    def test_the_detail_names_the_token_that_is_there(self):
        row = self._row(self.UNRECOGNIZED)

        assert "carries 'ATAC' where the sublibrary type belongs" in row.detail
        assert "carries no sublibrary type" not in row.detail

    def test_no_corrected_name_is_proposed_around_it(self):
        """Appending the vendor type would leave ATAC inside the sublibrary."""
        row = self._row(self.UNRECOGNIZED, vendor=self.VENDOR)

        assert row.expected_name == ""
        assert "ATAC" not in row.expected_name

    def test_the_vendor_value_is_still_reported(self):
        row = self._row(self.UNRECOGNIZED, vendor=self.VENDOR)

        assert "GEX_hash_oligo" in row.detail

    def test_a_genuinely_absent_type_reads_as_before(self):
        row = self._row(self.NO_TYPE)

        assert "carries no sublibrary type" in row.detail

    def test_a_genuinely_absent_type_is_still_corrected_from_the_vendor(self):
        """The 288 REF3 rows go through here; withholding these would break them."""
        row = self._row(self.NO_TYPE, vendor=self.VENDOR)

        assert row.expected_name == (
            "438514-REF3_P01_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT.trim.cram"
        )


class TestLabProjectIsMatchedOnAnExactPrefix:
    """The lab name is a surname, and surnames contain hyphens.

    Comparing ``project.split("-")[0]`` to the lab only works while the lab name
    has none: ``czi-van-der-berg`` gives lab ``van-der-berg``, its own project
    ``van-der-berg-seahub-bcp`` splits to ``van``, and the rule fired on a
    correct upload. The rule is ``upload`` scope, so the cost is one row for the
    whole listing plus the rename cell's banner telling an operator to fix a
    location that is fine.
    """

    STEM = "436830-REF3_P05_2_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"

    def _fires(self, bucket: str, project: str) -> bool:
        key = f"{project}/REF3/raw/REF3_P05_2/436830/{self.STEM}"
        return "lab_project_mismatch" in _types(validate_seahub_key(bucket, key))

    @pytest.mark.parametrize(
        "bucket,project",
        [
            ("czi-van-der-berg", "van-der-berg-seahub-bcp"),
            ("czi-labalpha", "labalpha-seahub-bcp"),
            # The project may be the bare lab name, with nothing appended.
            ("czi-labalpha", "labalpha"),
        ],
    )
    def test_a_project_under_its_own_lab_is_clean(self, bucket, project):
        assert not self._fires(bucket, project)

    @pytest.mark.parametrize(
        "bucket,project",
        [
            ("czi-labalpha", "otherlab-seahub-bcp"),
            ("czi-van-der-berg", "van-seahub-bcp"),
            # The trailing hyphen in the prefix test is what stops a lab name
            # that is a prefix of another lab's from matching it.
            ("czi-lab", "labalpha-seahub-bcp"),
            # It has to be a *prefix*: the lab name appearing anywhere in the
            # project is not the same claim, and a containment test would let
            # this through.
            ("czi-labalpha", "notlabalpha-seahub-bcp"),
            ("czi-labalpha", "seahub-labalpha-bcp"),
        ],
    )
    def test_a_project_under_another_lab_still_fires(self, bucket, project):
        assert self._fires(bucket, project)

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

    def test_a_folder_that_elides_the_experiment_prefix_is_clean(self):
        """The accepted spelling every real trimmed upload uses.

        The ExperimentID is already an ancestor segment, so the prefix carries
        nothing the path does not. Demanding the full form reported every GENE7
        sublibrary as broken and proposed a move for all 5184 of its objects.
        """
        assert validate_seahub_key("czi-labalpha", GOOD_ELIDED_PREFIX) == []

    def test_the_real_reported_path_is_clean(self):
        """The CHEM16 object that reported sublibrary_folder_truncated in review.

        Folder ``P03`` under ExperimentID ``CHEM16``, filename ``CHEM16_P03``.
        Sanitized to the fixture lab, keeping the shape exactly.
        """
        key = (
            "labalpha-seahub-bcp/CHEM16/raw/P03/432640/"
            "432640-CHEM16_P03_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram"
        )
        assert validate_seahub_key("czi-labalpha", key) == []

    def test_sublibrary_mismatch_when_the_two_cannot_be_reconciled(self):
        """Neither the folder nor {ExperimentID}_{folder} explains the name."""
        key = (
            "labalpha-seahub-bcp/REF3/raw/P05_2/436830/"
            "436830-SOMETHING_ELSE_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT.trim.cram"
        )
        assert "sublibrary_mismatch" in _types(validate_seahub_key("czi-labalpha", key))

    def test_browser_junk_is_reported_as_non_sequencing(self):
        key = "labalpha-seahub-bcp/REF3/raw/P05_2/436830/index.html"
        types = _types(validate_seahub_key("czi-labalpha", key))
        assert types == {"non_sequencing_artifact"}

    def test_an_unrecognized_name_is_still_an_unexpected_suffix(self):
        """The junk rule must not become a catch-all for anything unparsed."""
        for name in ("README.txt", f"{TRIM_STEM}.trim.cram.1"):
            key = f"labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479/{name}"
            assert _types(validate_seahub_key("czi-labalpha", key)) == {
                "unexpected_suffix"
            }

    def test_experiment_id_in_sublibrary_is_not_a_repeated_token(self):
        # REF3 appears in both the ExperimentID folder and the sublibrary name;
        # that is the reviewer's own known-good form.  repeated_token is scoped
        # to the basename precisely so this does not false-positive.
        assert "repeated_token" not in _types(
            validate_seahub_key("czi-labalpha", GOOD_LABALPHA)
        )


class TestRepeatedTokenIsScopedToTheSublibrary:
    """The rule looks at the group slot, after the stem has matched.

    Splitting the whole stem on both ``-`` and ``_`` into one flat bag put the
    sublibrary type in the same bag as the sublibrary name, so a sublibrary whose
    name legitimately contains a type token collided with the type itself.
    """

    def _key(self, sublibrary: str, stem: str) -> str:
        return f"labalpha-seahub-bcp/REF3/raw/{sublibrary}/438514/{stem}.trim.cram"

    def test_a_type_token_inside_the_sublibrary_name_is_not_a_repeat(self):
        """The false positive: SOP-clean but for a rule with no real catch.

        Via the rename self-check this withheld the proposal and flipped the well
        to UNKNOWN -- for a whole sublibrary at once, with no filename the tool
        would have accepted.
        """
        key = self._key(
            "REF3_GEX_P01",
            "438514-REF3_GEX_P01_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT",
        )

        assert _types(validate_seahub_key("czi-labalpha", key)) == set()

    def test_a_genuine_repeat_inside_the_sublibrary_still_fires(self):
        """Narrowed, not deleted: a repeat within the group is a real error."""
        key = self._key(
            "REF3_P01_P01",
            "438514-REF3_P01_P01_GEX_hash_oligo-Z0305-CACACACAACATGAT",
        )

        violations = validate_seahub_key("czi-labalpha", key)

        assert "repeated_token" in _types(violations)
        row = next(v for v in violations if v.type == "repeated_token")
        assert "P01" in row.detail

    def test_a_repeated_wafer_belongs_to_its_own_rule(self):
        """Not double-reported: the wafer head is duplicated_wafer_token's."""
        key = self._key(
            "REF3_P01",
            "438514-438514-REF3_P01_GEX_hash_oligo-Z0305-CACACACAACATGAT",
        )

        types = _types(validate_seahub_key("czi-labalpha", key))

        assert "duplicated_wafer_token" in types
        assert "repeated_token" not in types

    def test_an_unparseable_stem_reports_only_that(self):
        """No slots to compare, so the rule cannot run and does not guess."""
        key = self._key("REF3_P01", "not-a-sop-stem-at-all")

        assert "repeated_token" not in _types(validate_seahub_key("czi-labalpha", key))


class TestMalformedUgAndBarcode:
    """No bad_ug / bad_barcode rule: both stem patterns already pin the tokens.

    A dedicated rule could never fire, since reaching it means the stem matched a
    pattern requiring ``Z\\d{4}`` and ``[ACGT]+``. A malformed token fails both
    patterns instead, so the fact still surfaces -- as ``unparseable_stem``.
    """

    def test_a_malformed_ug_is_an_unparseable_stem(self):
        key = (
            f"{TRIM_DIR}/430479-REF3_P05_1_A1_GEX_hash_oligo"
            "-ZZ169-CAGTCAGTTGCAGAT.trim.cram"
        )
        assert _types(validate_seahub_key("czi-labalpha", key)) == {"unparseable_stem"}

    def test_a_barcode_outside_acgt_is_an_unparseable_stem(self):
        key = (
            f"{TRIM_DIR}/430479-REF3_P05_1_A1_GEX_hash_oligo"
            "-Z0097-CAGTCAGTTGCAGAN.trim.cram"
        )
        assert _types(validate_seahub_key("czi-labalpha", key)) == {"unparseable_stem"}

    def test_a_well_formed_stem_is_still_clean(self):
        key = f"{TRIM_DIR}/{TRIM_STEM}.trim.cram"
        assert validate_seahub_key("czi-labalpha", key) == []

    def test_the_rule_set_does_not_advertise_them(self):
        assert "bad_ug" not in SEAHUB_SOP_RULES
        assert "bad_barcode" not in SEAHUB_SOP_RULES


class TestStemLevelReporting:
    def test_one_stem_reported_once_across_its_artifacts(self):
        keys = [
            f"{BARE_DIR}/{BARE_STEM}{suffix}"
            for suffix in (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
        ]
        violations = validate_seahub_stems("czi-labalpha", keys)
        summary = sop_violation_summary(violations)
        assert summary["duplicated_wafer_token"] == 1
        assert summary["invalid_sublibrary_type"] == 1
        assert summary["missing_trim_infix"] == 1

    def test_junk_files_are_reported_individually(self):
        keys = [
            f"{BARE_DIR}/index.html",
            f"{BARE_DIR}/ug-icon.png",
        ]
        violations = validate_seahub_stems("czi-labalpha", keys)
        assert sop_violation_summary(violations)["non_sequencing_artifact"] == 2

    def test_a_mixed_family_well_still_reports_its_bare_artifact(self):
        """Regression: every artifact is validated, not just the first sorted one.

        A well delivered as a compliant ``.trim.cram`` plus a bare ``_fail.csv``
        used to report nothing at all -- ``.trim.cram`` sorts first ('.' < '_'),
        it validated clean, and the rest of the well was skipped. The same
        ``_fail.csv`` on its own reported correctly, which is what made the gap
        invisible.
        """
        keys = [f"{TRIM_DIR}/{TRIM_STEM}.trim.cram", f"{TRIM_DIR}/{TRIM_STEM}_fail.csv"]

        summary = sop_violation_summary(validate_seahub_stems("czi-labalpha", keys))

        assert summary == {"missing_trim_infix": 1}

    def test_a_wrong_bucket_is_one_row_for_the_whole_upload(self):
        """A 288-well upload used to write 288 identical bad_bucket rows."""
        keys = [
            f"{TRIM_DIR}/{stem}{suffix}"
            for stem in (TRIM_STEM, TRIM_STEM.replace("Z0097", "Z0098"))
            for suffix in (".trim.cram", ".trim.csv", ".trim.stderr")
        ]

        summary = sop_violation_summary(validate_seahub_stems("labalpha-data", keys))

        assert summary["bad_bucket"] == 1

    def test_a_project_lab_mismatch_is_one_row_for_the_whole_upload(self):
        keys = [
            f"{TRIM_DIR}/{stem}{suffix}"
            for stem in (TRIM_STEM, TRIM_STEM.replace("Z0097", "Z0098"))
            for suffix in (".trim.cram", ".trim.csv", ".trim.stderr")
        ]

        summary = sop_violation_summary(validate_seahub_stems("czi-labbeta", keys))

        assert summary["lab_project_mismatch"] == 1

    def test_two_projects_under_one_bucket_each_report(self):
        """The fact is per distinct project, not one row for the listing."""
        other = TRIM_DIR.replace("labalpha-seahub-bcp", "labgamma-seahub-bcp")
        keys = [
            f"{TRIM_DIR}/{TRIM_STEM}.trim.cram",
            f"{other}/{TRIM_STEM}.trim.cram",
        ]

        summary = sop_violation_summary(validate_seahub_stems("czi-labbeta", keys))

        assert summary["lab_project_mismatch"] == 2

    def test_an_unparseable_object_does_not_re_report_the_bucket(self):
        """The dedup spans both the grouped and the ungrouped key loops."""
        keys = [
            f"{TRIM_DIR}/{TRIM_STEM}.trim.cram",
            f"{TRIM_DIR}/index.html",
        ]

        summary = sop_violation_summary(validate_seahub_stems("labalpha-data", keys))

        assert summary["bad_bucket"] == 1

    def test_the_upload_scope_is_a_known_scope(self):
        violations = validate_seahub_stems(
            "labalpha-data", [f"{TRIM_DIR}/{TRIM_STEM}.trim.cram"]
        )
        assert {v.scope for v in violations} <= set(SEAHUB_VIOLATION_SCOPES)
        bucket_rule = next(v for v in violations if v.type == "bad_bucket")
        assert bucket_rule.scope == "upload"

    def test_a_fully_compliant_well_reports_nothing(self):
        keys = [
            f"{TRIM_DIR}/{TRIM_STEM}{suffix}"
            for suffix in (
                ".trim.cram",
                ".trim.csv",
                ".trim.stderr",
                ".trim.stdout",
                ".trim_fail.csv",
                ".trim.cram-metadata.json",
            )
        ]
        assert validate_seahub_stems("czi-labalpha", keys) == []

    def test_a_whole_sublibrary_under_an_elided_folder_reports_nothing(self):
        """The shape of a real upload: many wells and wafers, one folder spelling.

        This used to be the folder-scope dedup test -- these same keys produced
        one ``sublibrary_folder_truncated`` row for the sublibrary, and reported
        per object would have produced one per well. Now the spelling is accepted
        and there is nothing to collapse, so the whole listing is clean.
        """
        keys = []
        for wafer, ugs in (
            ("437120", ("Z0001", "Z0002", "Z0003")),
            ("437121", ("Z0004",)),
        ):
            for ug in ugs:
                stem = f"{wafer}-REF3_P04_1_A1_GEX_hash_oligo-{ug}-CAGCTCGAATGCGAT"
                keys.append(
                    f"labalpha-seahub-bcp/REF3/raw/P04_1/{wafer}/{stem}.trim.cram"
                )

        summary = sop_violation_summary(validate_seahub_stems("czi-labalpha", keys))

        assert summary == {}


class TestStemGrouping:
    def test_a_well_groups_on_wafer_and_ug(self):
        """The grouping key matches the cross-bucket identity key.

        qa_seahub_source indexes both the vendor delivery and the trimmed upload
        on (wafer, UG), so grouping the same way lets a group join the vendor
        index without re-parsing.
        """
        keys = [
            f"{BARE_DIR}/{BARE_STEM}{suffix}"
            for suffix in (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
        ]

        groups, unparsed = group_seahub_keys("czi-labalpha", keys)

        assert unparsed == []
        assert len(groups) == 1
        group = groups[0]
        assert group.identity == ("438514", "Z0305")
        assert group.families == frozenset({"bare"})
        assert group.well_id == "A3"
        assert group.has_cram is True
        assert group.normalized_stem == "438514-REF3_P07_1_A3-Z0305-CACACACAACATGAT"

    def test_a_well_missing_its_cram_is_visible_on_the_group(self):
        keys = [
            f"{BARE_DIR}/{BARE_STEM}{suffix}"
            for suffix in (".csv", ".stderr", ".stdout", "_fail.csv")
        ]

        groups, _unparsed = group_seahub_keys("czi-labalpha", keys)

        assert groups[0].has_cram is False

    def test_junk_lands_in_unparsed_not_in_a_group(self):
        groups, unparsed = group_seahub_keys(
            "czi-labalpha", [f"{BARE_DIR}/login.html", f"{BARE_DIR}/{BARE_STEM}.cram"]
        )

        assert unparsed == [f"{BARE_DIR}/login.html"]
        assert len(groups) == 1

    def test_a_mixed_family_well_is_one_group_carrying_both_families(self):
        keys = [f"{TRIM_DIR}/{TRIM_STEM}.trim.cram", f"{TRIM_DIR}/{TRIM_STEM}_fail.csv"]

        groups, _unparsed = group_seahub_keys("czi-labalpha", keys)

        assert len(groups) == 1
        assert groups[0].families == frozenset({"trim", "bare"})


class TestCompletenessIsPerFolderNotPerStem:
    """A stem is unique within a sublibrary folder, not across them.

    ``beginnings`` was keyed on the stem alone while its value carried the
    folder every expected path is built from, so a well uploaded under two
    folders kept only the first folder's. The second copy's artifacts matched
    nothing, never reached ``raw_found``, and check_extra_raw_files takes
    everything not there -- so two complete wells reported as one complete well
    plus five unexpected files.

    That upload is already reported: index_trimmed_upload calls it
    duplicate_trimmed_well and the roll-up calls the well UNKNOWN. The inventory
    was the one output describing it wrongly rather than differently.
    """

    STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
    SUFFIXES = (
        ".trim.cram",
        ".trim.csv",
        ".trim.stderr",
        ".trim.stdout",
        ".trim_fail.csv",
    )
    # Two *spellings of one sublibrary*, not two sublibraries: the second folder
    # elides the ExperimentID prefix the first carries. That is the shape the
    # accepted-both-spellings rule makes possible, so it is the one the inventory
    # has to key per folder -- and the one whose duplicate is pinned by
    # test_one_well_under_both_folder_spellings_is_still_caught.
    DIRS = (
        "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479",
        "labalpha-seahub-bcp/REF3/raw/P05_1/430479",
    )

    def _keys(self, omit: tuple[str, str] | None = None) -> list[str]:
        return [
            f"{d}/{self.STEM}{s}"
            for d in self.DIRS
            for s in self.SUFFIXES
            if (d, s) != omit
        ]

    def test_both_copies_count_as_complete(self):
        keys = self._keys()

        all_good, raw_lost, raw_found = check_expected_raw_files(keys, "seahub_sci")

        assert (all_good, raw_lost) == (2, [])
        assert len(raw_found) == len(keys)

    def test_neither_copy_is_reported_as_extra(self):
        keys = self._keys()

        _good, _lost, raw_found = check_expected_raw_files(keys, "seahub_sci")

        assert check_extra_raw_files(keys, raw_found, "seahub_sci") == []

    def test_an_incomplete_second_copy_is_reported_missing(self):
        keys = self._keys(omit=(self.DIRS[1], ".trim.cram"))

        all_good, raw_lost, _found = check_expected_raw_files(keys, "seahub_sci")

        assert all_good == 1
        assert len(raw_lost) == 1
        assert raw_lost[0][".trim.cram"] == f"{self.DIRS[1]}/{self.STEM}.trim.cram"

    def test_an_object_outside_a_wafer_folder_is_not_a_well(self):
        """It belongs to no well, so it stays an extra file and invents nothing.

        The recursive raw/ walk collects wrong-depth objects on purpose so
        bad_path_depth can see them, and roll_up_wells counts them as
        unaccounted. Keyed per directory without this gate, such an object
        became a well of its own: four artifacts reported missing from a folder
        that must not exist, and the stray itself dropped from the extra list --
        the only place it should have appeared.
        """
        stray = f"{self.DIRS[0]}/reupload/{self.STEM}.trim.cram"
        keys = [f"{self.DIRS[0]}/{self.STEM}{s}" for s in self.SUFFIXES] + [stray]

        all_good, raw_lost, raw_found = check_expected_raw_files(keys, "seahub_sci")

        assert (all_good, raw_lost) == (1, [])
        assert check_extra_raw_files(keys, raw_found, "seahub_sci") == [stray]

    def test_a_stray_alone_invents_no_well_either(self):
        """True before this change only when some other object shared its stem."""
        stray = f"{self.DIRS[0]}/reupload/{self.STEM}.trim.cram"

        all_good, raw_lost, raw_found = check_expected_raw_files([stray], "seahub_sci")

        assert (all_good, raw_lost) == (0, [])
        assert check_extra_raw_files([stray], raw_found, "seahub_sci") == [stray]

    def test_the_row_shape_is_unchanged(self):
        """The folder is in the per-ending paths, so no new column is needed."""
        keys = self._keys(omit=(self.DIRS[1], ".trim.cram"))

        _good, raw_lost, _found = check_expected_raw_files(keys, "seahub_sci")

        assert set(raw_lost[0]) == {"path", ".trim.cram"}
        assert raw_lost[0]["path"] == self.STEM


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
        # Keyed on the SOP artifact even though this well was delivered bare:
        # what is absent is the kind, and it should arrive under the SOP name.
        assert lost[0][".trim.cram"].endswith(f"{BARE_STEM}.trim.cram")

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

    def test_a_mixed_family_well_is_judged_by_kind_not_by_spelling(self):
        """It used to require *both* spellings of every kind.

        A well delivering ``.cram`` and ``.trim.csv`` has two of the five kinds,
        so three are missing. Demanding both spellings instead reported the two
        it did deliver as missing too -- files that were never meant to exist.
        """
        keys = [f"{BARE_DIR}/{BARE_STEM}.cram", f"{BARE_DIR}/{BARE_STEM}.trim.csv"]

        _all_good, lost, _found = check_expected_raw_files(keys, "seahub_sci")

        assert len(lost) == 1
        missing = {k for k in lost[0] if k != "path"}
        assert missing == {".trim.stderr", ".trim.stdout", ".trim_fail.csv"}

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
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL),
            stats,
            "438514",
            fail_counts=counts,
            stem_key=BARE_STEM,
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
        s3 = MockS3Client(
            keys=keys,
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
        assert (BARE_DIR, BARE_STEM) in data.seahub_fail_counts

    def test_per_well_csvs_are_not_read_as_sublibrary_q30_stats(self):
        # Per-well ".csv" / "_fail.csv" are single-well trimmer output and carry
        # no PCT_PF_Q30_bases.  Treating them as sublibrary stats fetched one
        # object per well and emitted a "PCT_PF_Q30_bases missing" warning for
        # each.  Wells without a sublibrary type are the exposed case, since the
        # generic filter only skips names containing "hash_oligo".
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P07_1/438514"
        keys = [
            f"{raw_dir}/{BARE_STEM}.cram",
            f"{raw_dir}/{BARE_STEM}.csv",
            f"{raw_dir}/{BARE_STEM}_fail.csv",
            f"{raw_dir}/{TRIM_STEM}.trim.csv",
            f"{raw_dir}/{TRIM_STEM}.trim_fail.csv",
        ]
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        with open(SEAHUB_TRIM_FAIL) as fh:
            fail_content = fh.read()
        # Only the failure CSVs are registered.  A Q30 read of the per-well
        # ".csv" files would go through s3fs and find nothing; the warning
        # asserted below is what the reviewer saw hundreds of times.
        s3 = MockS3Client(
            keys=keys,
            file_contents={
                f"{raw_dir}/{BARE_STEM}_fail.csv": fail_content,
                f"{raw_dir}/{TRIM_STEM}.trim_fail.csv": fail_content,
            },
        )
        data = gather_qa_data(ctx, s3)
        assert data.pct_q30_values == {}
        assert not [
            w for w in data.gathering_warnings if "PCT_PF_Q30_bases missing" in w
        ]

    def test_gather_records_sop_violations_and_one_summary_warning(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P07_1/438514"
        keys = [f"{raw_dir}/{BARE_STEM}{s}" for s in (".cram", ".csv", ".stdout")]
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, MockS3Client(keys=keys))
        types = {v["type"] for v in data.sop_violations}
        assert "missing_trim_infix" in types
        assert "duplicated_wafer_token" in types
        sop_warnings = [
            w for w in data.gathering_warnings if w.startswith("SOP VIOLATIONS:")
        ]
        assert len(sop_warnings) == 1


class TestOneTokenIsClaimedByOneRuleOnly:
    """An unrecognised type token must not also be reported as a bad well.

    The relaxed stem pattern matches whether the type token is absent or merely
    unrecognised, and the trailing-token check ran unconditionally afterwards. So
    ``438514-REF3_P01_ATAC-...`` produced two rows about ``ATAC`` giving opposite
    advice -- one saying it should be a sublibrary type, the next saying it should
    be a well -- and an upload that misspells its type token on every well doubled
    its SOP table. This module is built on one row per distinct fact.
    """

    def _types(self, stem: str, folder: str = "REF3_P01") -> list[str]:
        key = f"labalpha-seahub-bcp/REF3/raw/{folder}/438514/{stem}.trim.cram"
        return [v.type for v in validate_seahub_key("czi-labalpha", key)]

    def test_an_unrecognised_type_token_is_reported_once(self):
        types = self._types("438514-REF3_P01_ATAC-Z0305-CACACACAACATGAT")

        assert types == ["invalid_sublibrary_type"]

    def test_a_genuinely_bad_well_on_a_typed_stem_still_fires(self):
        """The suppression must not swallow the rule it shares a token with."""
        types = self._types("438514-REF3_P01_Z9_GEX-Z0305-CACACACAACATGAT")

        assert "bad_well" in types
        assert "invalid_sublibrary_type" not in types

    def test_a_missing_type_token_is_still_reported(self):
        """A trailing *well* with no type is the other relaxed-branch case."""
        types = self._types("438514-REF3_P01_A3-Z0305-CACACACAACATGAT")

        assert "invalid_sublibrary_type" in types
        assert "bad_well" not in types
