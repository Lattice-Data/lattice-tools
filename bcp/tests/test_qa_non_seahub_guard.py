"""
Permanent guard: SeaHub work must not move the non-SeaHub assays.

Why this module exists
---------------------
A reviewer reported that this branch changed 10x behavior. The report traced to
notebook layout, not code: the SeaHub cells sat between the "extra raw files"
cell and the CellRanger cell, moving CellRanger validation from cell index 13 to
17, so a run that stopped early produced a different ``*_errors.txt``. The one
genuine non-SeaHub code change was in :func:`qa_checks.check_extra_raw_files`,
where the optional-suffix computation moved from ``f.replace(f"{raw_dir}/{b}",
"")`` to an offset slice. Those two expressions agree only when the
reconstructed beginning is a literal prefix of the key, and where they diverge
the slice *silently accepted* a file that belonged in the extra list.

Both halves are pinned here: the narrow offset-slice guard, and a frozen
baseline for every assay that shares :mod:`qa_checks` with SeaHub.

The baselines are inline constants on purpose. The failure mode to prevent is
someone regenerating expectations from current behavior, which turns a
regression into a passing test. If a value here needs to change, that is a
deliberate decision to state in the commit message, not a refresh.

Do not delete this module.
"""

from __future__ import annotations

import pytest

from qa_checks import (
    MIN_METADATA_READ_COUNT,
    _fastq_count_mode,
    check_expected_raw_files,
    check_extra_raw_files,
    validate_read_metadata,
)
from qa_constants import raw_expected, raw_optional
from qa_mods import parse_raw_filename

# A well-formed 10x beginning: parse_raw_filename's regex matches it, so the
# reconstructed beginning is a literal prefix of the key.
GOOD_10X_DIR = "proj/NVUS01/G1/raw"
GOOD_10X_BASE = f"{GOOD_10X_DIR}/439047-G1_GEX-Z0273-CAGTCAGTTGCAGAT"

# The same shape with the assay token missing. parse_raw_filename falls through
# to its hyphen-splitting fallback, which cannot tell group from assay and
# returns both as "G1" -- so it reconstructs "439047-G1_G1-Z0169-ACGT", three
# characters longer than the real "439047-G1-Z0169-ACGT" prefix.
DIVERGENT_10X_KEY = f"{GOOD_10X_DIR}/439047-G1-Z0169-ACGT_S1_Log.out"


class TestExtraRawSuffixOffsetGuard:
    """Pins the optional-suffix computation in check_extra_raw_files."""

    def test_beginning_that_does_not_prefix_the_key_is_reported_extra(self):
        """A misaligned tail must never be matched against optional endings.

        ``_Log.out`` is in ``raw_optional["10x"]``, and slicing this key by the
        length of the reconstructed beginning yields exactly that string --
        which would silently accept the file. The suffix is meaningless
        whenever the beginning is not a literal prefix, so the file is extra.
        """
        assert parse_raw_filename(DIVERGENT_10X_KEY, "10x") == (
            "439047",
            "G1",
            "G1",
            "Z0169",
            "ACGT",
        )
        assert check_extra_raw_files([DIVERGENT_10X_KEY], [], "10x") == [
            DIVERGENT_10X_KEY
        ]

    @pytest.mark.parametrize(
        "ending", raw_expected["10x"] + raw_optional["10x"], ids=lambda e: e
    )
    def test_offset_and_replace_agree_for_every_well_formed_ending(self, ending):
        """Documents when the two expressions are interchangeable.

        For a name whose reconstructed beginning *is* a literal prefix, the
        offset slice and the original ``.replace()`` both yield the ending. The
        guard added to check_extra_raw_files only affects the other case.
        """
        key = f"{GOOD_10X_BASE}{ending}"
        run, group, assay, ug, barcode = parse_raw_filename(key, "10x")
        beginning = f"{run}-{group}_{assay}-{ug}-{barcode}"
        prefix = f"{GOOD_10X_DIR}/{beginning}"

        assert key.startswith(prefix)
        assert key[len(prefix) :] == key.replace(prefix, "") == ending

    def test_optional_suffix_is_still_accepted_for_a_well_formed_name(self):
        assert check_extra_raw_files([f"{GOOD_10X_BASE}_Log.out"], [], "10x") == []

    def test_unknown_suffix_on_a_well_formed_name_is_still_extra(self):
        key = f"{GOOD_10X_BASE}.bam"
        assert check_extra_raw_files([key], [], "10x") == [key]

    def test_viral_orf_gex_well_borrows_the_10x_optional_set(self):
        """Pins the 10x_viral_ORF/GEX special case sharing this computation."""
        key = "proj/NVUS01/G1/raw/440002-G1_GEX-Z0101-CAGTCAGTTGCAGAT_Log.out"
        assert check_extra_raw_files([key], [], "10x_viral_ORF") == []

    @pytest.mark.parametrize("raw_assay", ["sci_plex", "sci_jumbo", "scale"])
    def test_assays_that_never_reach_the_offset_computation(self, raw_assay):
        """sci_plex and sci_jumbo have no raw_optional entry; scale's is empty.

        All three short-circuit before the suffix computation, so the guard
        cannot affect them. Asserted rather than asserted-in-a-comment.
        """
        assert not raw_optional.get(raw_assay)
        key = f"{GOOD_10X_BASE}.bam"
        assert check_extra_raw_files([key], [], raw_assay) == [key]


# base key, expected parse tuple, expected count of found paths for a complete well
BASELINES = {
    "10x": (
        GOOD_10X_BASE,
        ("439047", "G1", "GEX", "Z0273", "CAGTCAGTTGCAGAT"),
        20,
    ),
    "10x_cram": (
        "proj/NVUS01/LeS1867W11/raw/442488-LeS1867W11_ATAC-Z0027-CACTGTCAGCCAGAT",
        ("442488", "LeS1867W11", "ATAC", "Z0027", "CACTGTCAGCCAGAT"),
        8,
    ),
    "10x_viral_ORF": (
        "proj/NVUS01/G1/raw/440001-G1_viral_ORF-Z0100-CAGTCAGTTGCAGAT",
        ("440001", "G1", "viral_ORF", "Z0100", "CAGTCAGTTGCAGAT"),
        8,
    ),
    "sci_plex": (
        "proj/NVUS01/G1/raw/441000-G1_GEX_hash_oligo-Z0200-CAGTCAGTTGCAGAT",
        ("441000", "G1", "GEX_hash_oligo", "Z0200", "CAGTCAGTTGCAGAT"),
        6,
    ),
    "sci_jumbo": (
        "proj/NVUS01/G1/raw/441001-G1_GEX_hash_oligo-Z0201-CAGTCAGTTGCAGAT",
        ("441001", "G1", "GEX_hash_oligo", "Z0201", "CAGTCAGTTGCAGAT"),
        8,
    ),
}

EXPECTED_FASTQ_COUNT_MODES = {
    "10x": "10x",
    "10x_cram": "skip",
    "10x_viral_ORF": "skip",
    "sci_plex": "gex_hash",
    "sci_jumbo": "skip",
    "scale": "gex_hash",
    "seahub_sci": "gex_hash",
}


class TestNonSeahubAssayBaseline:
    """Frozen behavior for every assay that is not seahub_sci."""

    @pytest.mark.parametrize("raw_assay", list(BASELINES))
    def test_parse_raw_filename_is_unchanged(self, raw_assay):
        base, expected, _found = BASELINES[raw_assay]
        assert parse_raw_filename(f"{base}.csv", raw_assay) == expected

    @pytest.mark.parametrize("raw_assay", list(BASELINES))
    def test_a_complete_well_is_complete(self, raw_assay):
        base, _expected, found_count = BASELINES[raw_assay]
        keys = [f"{base}{ending}" for ending in raw_expected[raw_assay]]

        all_good, raw_lost, raw_found = check_expected_raw_files(keys, raw_assay)

        assert (all_good, raw_lost) == (1, [])
        assert len(raw_found) == found_count

    @pytest.mark.parametrize("raw_assay", list(BASELINES))
    def test_one_missing_ending_is_reported_once_with_its_full_path(self, raw_assay):
        base, _expected, _found = BASELINES[raw_assay]
        dropped = raw_expected[raw_assay][0]
        keys = [
            f"{base}{ending}" for ending in raw_expected[raw_assay] if ending != dropped
        ]

        _all_good, raw_lost, _raw_found = check_expected_raw_files(keys, raw_assay)

        assert raw_lost == [{"path": base.split("/")[-1], dropped: f"{base}{dropped}"}]

    @pytest.mark.parametrize("raw_assay", list(BASELINES))
    def test_an_unknown_artifact_is_extra_and_nothing_else_is(self, raw_assay):
        base, _expected, _found = BASELINES[raw_assay]
        keys = [f"{base}{ending}" for ending in raw_expected[raw_assay]]
        junk = f"{base}.bam"

        _all_good, _raw_lost, raw_found = check_expected_raw_files(
            keys + [junk], raw_assay
        )

        assert check_extra_raw_files(keys + [junk], raw_found, raw_assay) == [junk]

    @pytest.mark.parametrize(
        "raw_assay,expected_mode", sorted(EXPECTED_FASTQ_COUNT_MODES.items())
    )
    def test_fastq_count_mode_is_frozen(self, raw_assay, expected_mode):
        assert _fastq_count_mode(raw_assay) == expected_mode

    def test_expected_and_optional_ending_sets_are_frozen(self):
        """Catches an ending being added or removed as a side effect.

        Only the assays whose completeness is expressed through the suffix
        mechanism are listed; scale bypasses it entirely.
        """
        assert raw_expected["10x"] == [
            ".csv",
            ".json",
            "_trimmer-failure_codes.csv",
            "_trimmer-stats.csv",
            "_unmatched.cram",
            "_unmatched.cram-metadata.json",
            "_unmatched.csv",
            "_unmatched.json",
            "_S1_L001_R1_001.csv",
            "_S1_L001_R1_001.fastq.gz",
            "_S1_L001_R1_001.fastq.gz-metadata.json",
            "_S1_L001_R1_001.json",
            "_S1_L001_R1_001_sample.fastq.gz",
            "_S1_L001_R1_001_sample.fastq.gz-metadata.json",
            "_S1_L001_R2_001.csv",
            "_S1_L001_R2_001.fastq.gz",
            "_S1_L001_R2_001.fastq.gz-metadata.json",
            "_S1_L001_R2_001.json",
            "_S1_L001_R2_001_sample.fastq.gz",
            "_S1_L001_R2_001_sample.fastq.gz-metadata.json",
        ]
        assert raw_optional["10x"] == [
            ".scRNA.applicationQC.h5",
            ".scRNA.applicationQC.html",
            "_Log.final.out",
            "_Log.out",
            "_Log.progress.out",
            "_ReadsPerGene.out.tab",
            "_SJ.out.tab",
        ]
        assert raw_optional["10x_cram"] == ["_extract_stats.h5"]
        assert raw_expected["scale"] == []
        assert raw_optional["scale"] == []


class TestScaleBaseline:
    """Scale bypasses the suffix mechanism, so it is pinned separately."""

    WAFER_DIR = "proj/NVUS123/GENE9-R115/raw/440115"
    LISTING = [
        f"{WAFER_DIR}/440115-R115H_GEX_QSR-8-5B.cram",
        f"{WAFER_DIR}/440115-R115H_GEX_QSR-8_trimmer-stats.csv",
        f"{WAFER_DIR}/440115_LibraryInfo.xml",
    ]
    JUNK = f"{WAFER_DIR}/junk.bam"

    def test_wafer_listing_is_unchanged(self):
        listing = self.LISTING + [self.JUNK]

        all_good, raw_lost, raw_found = check_expected_raw_files(listing, "scale")

        # raw_expected["scale"] is empty, so all_good counts recognized
        # beginnings and nothing lands in raw_found.
        assert (all_good, raw_lost, raw_found) == (2, [], [])
        assert check_extra_raw_files(listing, raw_found, "scale") == [self.JUNK]


class TestReadMetadataBaseline:
    """R1/R2 pairing must stay on for every assay but 10x_cram and seahub_sci."""

    R1 = "439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
    R2 = "439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz"

    def _metadata(self):
        return {
            self.R1: {"read_count": MIN_METADATA_READ_COUNT, "errors": []},
            self.R2: {"read_count": MIN_METADATA_READ_COUNT, "errors": []},
        }

    def test_a_matched_pair_is_validated_for_10x(self):
        counts, errors, pairing = validate_read_metadata(self._metadata(), "10x")

        assert counts == {"G1": {"GEX": MIN_METADATA_READ_COUNT}}
        assert errors == []
        assert pairing == {"r1_without_r2_metadata": [], "r2_without_r1_metadata": []}

    def test_a_mismatched_pair_is_reported_for_10x(self):
        metadata = self._metadata()
        metadata[self.R2]["read_count"] = MIN_METADATA_READ_COUNT + 1

        _counts, errors, _pairing = validate_read_metadata(metadata, "10x")

        assert any(e.startswith("READ COUNT ERROR:") for e in errors)

    @pytest.mark.parametrize(
        "raw_assay,pairing_skipped",
        [
            ("10x", False),
            ("10x_viral_ORF", False),
            ("sci_plex", False),
            ("sci_jumbo", False),
            ("10x_cram", True),
        ],
    )
    def test_pairing_is_skipped_only_for_cram_only_layouts(
        self, raw_assay, pairing_skipped
    ):
        """Only R1 is summed into group counts when pairing is active.

        With pairing skipped, both reads are summed, so the group total
        doubles. That difference is the observable proxy for the flag.

        scale and seahub_sci are absent deliberately: both derive the group from
        path segments rather than the basename, so a fastq-style basename cannot
        be parsed for them at all, and neither delivery contains R1/R2 fastqs.
        Their own suites cover them.
        """
        counts, _errors, _pairing = validate_read_metadata(self._metadata(), raw_assay)

        total = sum(sum(assays.values()) for assays in counts.values())
        expected = MIN_METADATA_READ_COUNT * (2 if pairing_skipped else 1)
        assert total == expected
