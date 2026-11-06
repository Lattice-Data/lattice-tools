"""
Tests for Illumina delivery support (``sequencing_type="Illumina"``).

The QA pipeline grew up around Ultima (UG100): wafers, the Ultima trimmer, and
fastqs the trimmer emits under the single literal name ``_S1_L001_R1_001``.
An Illumina delivery is a different shape, and this module pins it:

    {proj}/{Order#}/{GroupID}/raw/        -- flat, one directory
        {FC}_CopyComplete.txt             -- eight per-flowcell run files
        ...
        {FC}_RunParameters.xml            -- says which reads to expect
        {FC}_{GroupID}_{Assay}_S{n}_L00#_{READ}_001.fastq.gz
        ...each with a -metadata.json sidecar

Two properties are load-bearing and are asserted here rather than assumed:

* Which reads are required comes from each flowcell's own RunParameters.xml.
  A run with no index reads plans two; a dual-indexed run plans four. QA must
  not hard-code either.
* Lanes are accepted as they arrive and only reported, so nothing may fail a
  delivery for having four lanes on one library and two on another.
"""

import pytest

from qa_checks import (
    check_expected_raw_files,
    format_read_count,
    check_extra_raw_files,
    summarize_illumina_flowcells,
    total_reads_by_flowcell,
    validate_fastq_counts,
    validate_read_metadata,
)
from qa_constants import (
    ALLOWED_RAW_ASSAYS,
    ILLUMINA_FLOWCELL_FILE_SUFFIXES,
    ILLUMINA_MODALITIES,
)
from qa_mods import (
    normalize_sequencing_type,
    parse_planned_reads,
    parse_raw_filename,
    resolve_qa_run_context,
)

RAW = "testproj/ORD01/G1/raw"
FC = "22YKTVLT4"

DUAL_INDEX_XML = """<?xml version="1.0"?>
<RunParameters>
  <PlannedReads>
    <Read ReadName="Read1" Cycles="151" />
    <Read ReadName="Index1" Cycles="10" />
    <Read ReadName="Index2" Cycles="10" />
    <Read ReadName="Read2" Cycles="151" />
  </PlannedReads>
</RunParameters>
"""

NO_INDEX_XML = """<?xml version="1.0"?>
<RunParameters>
  <PlannedReads>
    <Read ReadName="Read1" Cycles="151" />
    <Read ReadName="Read2" Cycles="151" />
  </PlannedReads>
</RunParameters>
"""


def flowcell_files(flowcell=FC, raw_dir=RAW, omit=()):
    return [
        f"{raw_dir}/{flowcell}{s}"
        for s in ILLUMINA_FLOWCELL_FILE_SUFFIXES
        if s not in omit
    ]


def fastqs(
    reads=("R1", "R2"),
    lanes=("001",),
    group="G1",
    assay="GEX",
    sample="1",
    flowcell=FC,
    raw_dir=RAW,
    sidecars=True,
    omit=(),
):
    out = []
    for lane in lanes:
        for read in reads:
            f = f"{raw_dir}/{flowcell}_{group}_{assay}_S{sample}_L{lane}_{read}_001.fastq.gz"
            if (read, lane) in omit:
                continue
            out.append(f)
            if sidecars:
                out.append(f"{f}-metadata.json")
    return out


class TestNormalizeSequencingType:
    """Mirrors TestNormalizeRawAssay: a hand-typed parameter must fail loudly."""

    def test_ultima_case_insensitive(self):
        assert normalize_sequencing_type("ultima") == "Ultima"
        assert normalize_sequencing_type("ULTIMA") == "Ultima"

    def test_illumina_case_insensitive(self):
        assert normalize_sequencing_type("illumina") == "Illumina"
        assert normalize_sequencing_type(" Illumina ") == "Illumina"

    def test_empty_raises(self):
        for value in ("", "   ", None):
            with pytest.raises(ValueError, match="sequencing_type"):
                normalize_sequencing_type(value)

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="not recognized"):
            normalize_sequencing_type("NovaSeq")


class TestSequencingTypeOnTheContext:
    """Mirrors the allow_truncated_stats_name trio, plus the 10x-only rule."""

    def _s3(self, **kw):
        # Components rather than an s3:// literal: a "czi-{lab}/{project}"
        # string in a tracked file is a project slot, and
        # test_sanitized_identifiers guards those against real names.
        return resolve_qa_run_context(
            data_source="s3",
            provider="novogene",
            proj="testproj",
            order="ORD01",
            **kw,
        )

    def test_default_is_ultima(self):
        assert self._s3(raw_assay="10x").sequencing_type == "Ultima"

    def test_propagated_s3(self):
        assert self._s3(
            raw_assay="10x", sequencing_type="Illumina"
        ).sequencing_type == ("Illumina")

    def test_propagated_manifest(self, tmp_path):
        manifest = tmp_path / "m.tsv"
        manifest.write_text("s3://czi-novogene/test-project/ORD01/G1/raw/a.fastq.gz\n")
        ctx = resolve_qa_run_context(
            data_source="manifest",
            raw_assay="10x",
            manifest_path=str(manifest),
            run_label="ORD01",
            sequencing_type="Illumina",
        )
        assert ctx.sequencing_type == "Illumina"

    @pytest.mark.parametrize("assay", sorted(ALLOWED_RAW_ASSAYS - {"10x"}))
    def test_illumina_rejects_every_other_assay(self, assay):
        # Parametrized over the live set so a newly added assay is covered
        # without anyone remembering to extend this list.
        with pytest.raises(ValueError, match="Illumina"):
            self._s3(raw_assay=assay, sequencing_type="Illumina")


class TestParsePlannedReads:
    def test_dual_index_run(self):
        assert parse_planned_reads(DUAL_INDEX_XML) == ["R1", "I1", "I2", "R2"]

    def test_run_without_index_reads(self):
        assert parse_planned_reads(NO_INDEX_XML) == ["R1", "R2"]

    def test_unknown_read_name_is_skipped_not_fatal(self):
        xml = DUAL_INDEX_XML.replace('ReadName="Index2"', 'ReadName="Read7"')
        assert parse_planned_reads(xml) == ["R1", "I1", "R2"]

    def test_malformed_xml_returns_none(self):
        assert parse_planned_reads("not xml at all") is None

    def test_missing_planned_reads_returns_none(self):
        assert parse_planned_reads("<RunParameters/>") is None


class TestParseRawFilenameIllumina:
    NAME = f"{RAW}/{FC}_G1_GEX_S3_L002_R1_001.fastq.gz"

    def test_parses_into_the_five_tuple(self):
        assert parse_raw_filename(self.NAME, "10x", sequencing_type="Illumina") == (
            FC,
            "G1",
            "GEX",
            "S3",
            "L002",
        )

    def test_group_id_containing_a_hyphen(self):
        name = f"{RAW}/{FC}_LAT-123_GEX_S1_L001_R2_001.fastq.gz"
        run, group, assay, _, _ = parse_raw_filename(
            name, "10x", sequencing_type="Illumina"
        )
        # The group must survive intact: _process_raw_file compares it against
        # the folder name and raises WRONG GROUP on a mismatch.
        assert (run, group, assay) == (FC, "LAT-123", "GEX")

    def test_the_delivery_spec_example(self):
        # A GroupID made of underscore-separated parts, as real deliveries are.
        name = f"{RAW}/B23YTVYLT4_BCP_1_U3_GEX_S3_L002_R1_001.fastq.gz"
        assert parse_raw_filename(name, "10x", sequencing_type="Illumina") == (
            "B23YTVYLT4",
            "BCP_1_U3",
            "GEX",
            "S3",
            "L002",
        )

    def test_an_assay_shaped_word_inside_the_group_id(self):
        # "CRI" is an assay token, but it is not followed by _S{n}_L{lane}, so
        # the group runs through it and GEX is the assay.
        name = f"{RAW}/{FC}_CRI_POOL_GEX_S2_L001_R1_001.fastq.gz"
        _, group, assay, _, _ = parse_raw_filename(
            name, "10x", sequencing_type="Illumina"
        )
        assert (group, assay) == ("CRI_POOL", "GEX")

    def test_the_hyphenated_form_is_not_illumina(self):
        name = f"{RAW}/{FC}-G1_GEX_S1_L001_R1_001.fastq.gz"
        assert parse_raw_filename(name, "10x", sequencing_type="Illumina") is None

    def test_combined_modality_token_is_not_split(self):
        name = f"{RAW}/{FC}_G1_ADT_HTO_S1_L001_R1_001.fastq.gz"
        _, group, assay, _, _ = parse_raw_filename(
            name, "10x", sequencing_type="Illumina"
        )
        # Not group "G1_ADT" + modality "HTO".
        assert (group, assay) == ("G1", "ADT_HTO")

    @pytest.mark.parametrize("modality", ILLUMINA_MODALITIES)
    def test_every_illumina_modality_parses(self, modality):
        name = f"{RAW}/{FC}_BCP_1_U3_{modality}_S3_L002_R1_001.fastq.gz"
        _, group, assay, _, _ = parse_raw_filename(
            name, "10x", sequencing_type="Illumina"
        )
        assert (group, assay) == ("BCP_1_U3", modality)

    def test_an_unlisted_modality_does_not_parse(self):
        # GEX_hash_oligo is an Ultima sublibrary type, not an Illumina modality.
        name = f"{RAW}/{FC}_G1_GEX_hash_oligo_S1_L001_R1_001.fastq.gz"
        assert parse_raw_filename(name, "10x", sequencing_type="Illumina") is None

    def test_flowcell_run_file_is_not_a_library(self):
        assert (
            parse_raw_filename(
                f"{RAW}/{FC}_RunInfo.xml", "10x", sequencing_type="Illumina"
            )
            is None
        )

    @pytest.mark.parametrize("suffix", ILLUMINA_FLOWCELL_FILE_SUFFIXES)
    def test_no_run_file_parses_as_a_library(self, suffix):
        """None of the eight may reach the Ultima fallback splitter.

        That splitter cannot tell group from assay, so it invents an assay
        token for whatever it is handed -- which would fire WRONG ASSAY on
        every flowcell run file in the delivery.
        """
        assert (
            parse_raw_filename(f"{RAW}/{FC}{suffix}", "10x", sequencing_type="Illumina")
            is None
        )

    def test_ultima_name_is_unchanged_by_the_new_keyword(self):
        name = "p/O/G1/raw/439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz"
        assert parse_raw_filename(name, "10x") == parse_raw_filename(
            name, "10x", sequencing_type="Ultima"
        )


class TestIlluminaExpectedFiles:
    def _check(self, keys, planned=None):
        return check_expected_raw_files(
            keys,
            "10x",
            sequencing_type="Illumina",
            planned_reads=planned
            if planned is not None
            else {FC: ["R1", "I1", "I2", "R2"]},
        )

    def test_complete_delivery_is_all_good(self):
        keys = flowcell_files() + fastqs(
            reads=("R1", "I1", "I2", "R2"), lanes=("001", "002")
        )
        all_good, raw_lost, raw_found = self._check(keys)
        assert raw_lost == []
        # one flowcell unit + one library-lane unit per lane
        assert all_good == 3
        assert set(raw_found) == set(keys)

    def test_missing_run_file_is_reported(self):
        keys = flowcell_files(omit=("_RTAExited.txt",)) + fastqs()
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        rows = [r for r in raw_lost if "_RTAExited.txt" in r]
        assert len(rows) == 1
        assert rows[0]["_RTAExited.txt"] == f"{RAW}/{FC}_RTAExited.txt"

    def test_missing_logs_gzip_is_reported(self):
        keys = flowcell_files(omit=("_Logs.gzip",)) + fastqs()
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        assert any(r.get("_Logs.gzip") == f"{RAW}/{FC}_Logs.gzip" for r in raw_lost)

    def test_read_missing_on_one_lane_only(self):
        keys = flowcell_files() + fastqs(lanes=("001", "002"), omit=(("R2", "002"),))
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        rows = [r for r in raw_lost if r["path"].endswith("_L002")]
        assert len(rows) == 1
        assert rows[0]["R2"].endswith("_L002_R2_001.fastq.gz")
        assert not [r for r in raw_lost if r["path"].endswith("_L001")]

    def test_index_reads_required_only_when_planned(self):
        keys = flowcell_files() + fastqs(reads=("R1", "R2"))
        _, lost_two, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        assert lost_two == []
        _, lost_four, _ = self._check(keys, planned={FC: ["R1", "I1", "I2", "R2"]})
        missing = {k for row in lost_four for k in row if k != "path"}
        assert missing == {"I1", "I2"}

    def test_missing_sidecar_is_reported(self):
        keys = flowcell_files() + fastqs(reads=("R1", "R2"))
        keys.remove(f"{RAW}/{FC}_G1_GEX_S1_L001_R2_001.fastq.gz-metadata.json")
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        rows = [r for r in raw_lost if "R2-metadata.json" in r]
        assert len(rows) == 1

    def test_an_absent_read_is_not_also_reported_as_an_absent_sidecar(self):
        keys = flowcell_files() + fastqs(reads=("R1",))
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        keys_reported = {k for row in raw_lost for k in row if k != "path"}
        assert keys_reported == {"R2"}

    def test_unreadable_xml_falls_back_to_paired_reads(self):
        keys = flowcell_files() + fastqs(reads=("R1", "R2"))
        all_good, raw_lost, _ = self._check(keys, planned={})
        assert raw_lost == []
        assert all_good == 2

    def test_two_flowcells_are_checked_against_their_own_planned_reads(self):
        other = "22ABCDEF9"
        keys = (
            flowcell_files()
            + fastqs(reads=("R1", "I1", "I2", "R2"))
            + flowcell_files(flowcell=other)
            + fastqs(reads=("R1", "R2"), flowcell=other, sample="2")
        )
        _, raw_lost, _ = self._check(
            keys, planned={FC: ["R1", "I1", "I2", "R2"], other: ["R1", "R2"]}
        )
        assert raw_lost == []

    def test_same_flowcell_under_two_groups_is_two_units(self):
        other_dir = "testproj/ORD01/G2/raw"
        keys = (
            flowcell_files()
            + fastqs()
            + flowcell_files(raw_dir=other_dir)
            + fastqs(raw_dir=other_dir, group="G2")
        )
        all_good, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        assert raw_lost == []
        assert all_good == 4  # two flowcell units + two library-lane units


class TestIlluminaExtraFiles:
    def test_expected_delivery_has_no_extras(self):
        keys = flowcell_files() + fastqs(
            reads=("R1", "I1", "I2", "R2"), lanes=("001", "002")
        )
        _, _, raw_found = check_expected_raw_files(
            keys,
            "10x",
            sequencing_type="Illumina",
            planned_reads={FC: ["R1", "I1", "I2", "R2"]},
        )
        assert (
            check_extra_raw_files(keys, raw_found, "10x", sequencing_type="Illumina")
            == []
        )

    def test_a_stray_file_is_flagged(self):
        stray = f"{RAW}/leftover.bam"
        ultima_stray = f"{RAW}/{FC}_G1_GEX_S1_L001_trimmer-stats.csv"
        keys = flowcell_files() + fastqs() + [stray, ultima_stray]
        _, _, raw_found = check_expected_raw_files(
            keys,
            "10x",
            sequencing_type="Illumina",
            planned_reads={FC: ["R1", "R2"]},
        )
        extra = check_extra_raw_files(
            keys, raw_found, "10x", sequencing_type="Illumina"
        )
        assert sorted(extra) == sorted([stray, ultima_stray])


class TestIlluminaFlowcellSummary:
    def test_reports_lanes_and_libraries_without_enforcing_them(self):
        keys = (
            flowcell_files()
            + fastqs(lanes=("001", "002", "003", "004"))
            + fastqs(lanes=("001", "002"), assay="CRI", sample="2")
        )
        rows = summarize_illumina_flowcells(keys, {FC: ["R1", "R2"]})
        assert len(rows) == 1
        assert rows[0]["flowcell"] == FC
        assert rows[0]["planned_reads"] == "R1, R2"
        assert rows[0]["lanes"] == "L001, L002, L003, L004"
        assert rows[0]["library_count"] == 2

    def test_columns_lead_with_group_id(self):
        """The table is read per GroupID, so that column comes first."""
        rows = summarize_illumina_flowcells(
            flowcell_files() + fastqs(), {FC: ["R1", "R2"]}
        )
        assert list(rows[0]) == [
            "group_id",
            "flowcell",
            "modality",
            "planned_reads",
            "read_types_found",
            "missing_read_type",
            "lanes",
            "library_count",
            "read_count",
        ]
        assert rows[0]["group_id"] == "G1"

    @staticmethod
    def _meta(reads_by_name, raw_dir=RAW):
        """read_metadata as the gather step builds it, keyed by basename."""
        return {
            name: {
                "read_count": n,
                "errors": [],
                "__source_key": f"{raw_dir}/{name}-metadata.json",
            }
            for name, n in reads_by_name.items()
        }

    def test_read_count_sums_r1_only_across_libraries_and_lanes(self):
        keys = (
            flowcell_files()
            + fastqs(reads=("R1", "I1", "R2"), lanes=("001", "002"))
            + fastqs(reads=("R1", "I1", "R2"), lanes=("001",), assay="CRI", sample="2")
        )
        meta = self._meta(
            {
                f"{FC}_G1_GEX_S1_L001_R1_001.fastq.gz": 10_000_000,
                f"{FC}_G1_GEX_S1_L002_R1_001.fastq.gz": 12_000_000,
                f"{FC}_G1_CRI_S2_L001_R1_001.fastq.gz": 3_000_000,
                # Same reads again -- must not be added.
                f"{FC}_G1_GEX_S1_L001_R2_001.fastq.gz": 10_000_000,
                f"{FC}_G1_GEX_S1_L001_I1_001.fastq.gz": 10_000_000,
            }
        )
        rows = summarize_illumina_flowcells(keys, {FC: ["R1", "I1", "R2"]}, meta)
        assert rows[0]["read_count"] == 25_000_000

    def test_read_count_is_per_group_and_flowcell(self):
        other_dir = "testproj/ORD01/G2/raw"
        keys = (
            flowcell_files()
            + fastqs()
            + flowcell_files(raw_dir=other_dir)
            + fastqs(raw_dir=other_dir, group="G2")
        )
        meta = {
            **self._meta({f"{FC}_G1_GEX_S1_L001_R1_001.fastq.gz": 4_000_000}),
            # Same basename shape in a different folder: must land in G2's row.
            "s3://czi-novogene/x/G2-r1": {
                "read_count": 6_000_000,
                "errors": [],
                "__source_key": f"{other_dir}/{FC}_G2_GEX_S1_L001_R1_001.fastq.gz-metadata.json",
            },
        }
        rows = {r["group_id"]: r for r in summarize_illumina_flowcells(keys, {}, meta)}
        assert rows["G1"]["read_count"] == 4_000_000
        assert rows["G2"]["read_count"] == 6_000_000

    def test_errored_or_countless_sidecars_are_skipped(self):
        keys = flowcell_files() + fastqs(lanes=("001", "002", "003"))
        meta = self._meta(
            {
                f"{FC}_G1_GEX_S1_L001_R1_001.fastq.gz": 5_000_000,
                f"{FC}_G1_GEX_S1_L002_R1_001.fastq.gz": 7_000_000,
            }
        )
        meta[f"{FC}_G1_GEX_S1_L002_R1_001.fastq.gz"]["errors"] = ["truncated"]
        meta[f"{FC}_G1_GEX_S1_L003_R1_001.fastq.gz"] = {
            "errors": [],
            "__source_key": f"{RAW}/{FC}_G1_GEX_S1_L003_R1_001.fastq.gz-metadata.json",
        }
        rows = summarize_illumina_flowcells(keys, {}, meta)
        assert rows[0]["read_count"] == 5_000_000

    def test_read_count_is_zero_without_metadata(self):
        rows = summarize_illumina_flowcells(flowcell_files() + fastqs(), {})
        assert rows[0]["read_count"] == 0

    FOUR = ["R1", "I1", "I2", "R2"]

    def test_complete_delivery_lists_every_type_and_nothing_missing(self):
        keys = (
            flowcell_files()
            + fastqs(reads=self.FOUR, lanes=("001", "002"))
            + fastqs(reads=self.FOUR, lanes=("002",), assay="CRI", sample="3")
        )
        row = summarize_illumina_flowcells(keys, {FC: self.FOUR})[0]
        # Library reads first, then index reads, whatever order the XML used.
        assert row["read_types_found"] == "R1, R2, I1, I2"
        assert row["missing_read_type"] == "none"

    def test_one_missing_read_in_one_lane_set_is_named(self):
        # The user's example: an ..._S3_L002 set without its I2 fastq. I2 is
        # still found on L001, so it is both found and missing.
        keys = flowcell_files() + fastqs(
            reads=self.FOUR, lanes=("001", "002"), sample="3", omit=(("I2", "002"),)
        )
        row = summarize_illumina_flowcells(keys, {FC: self.FOUR})[0]
        assert row["read_types_found"] == "R1, R2, I1, I2"
        assert row["missing_read_type"] == "I2"

    def test_several_missing_types_are_all_listed(self):
        keys = flowcell_files() + fastqs(
            reads=self.FOUR,
            lanes=("001", "002"),
            omit=(("R2", "001"), ("I1", "002"), ("I2", "002")),
        )
        row = summarize_illumina_flowcells(keys, {FC: self.FOUR})[0]
        assert row["missing_read_type"] == "R2, I1, I2"

    def test_reads_not_planned_are_not_missing(self):
        keys = flowcell_files() + fastqs(reads=("R1", "R2"))
        row = summarize_illumina_flowcells(keys, {FC: ["R1", "R2"]})[0]
        assert row["read_types_found"] == "R1, R2"
        assert row["missing_read_type"] == "none"
        row = summarize_illumina_flowcells(keys, {FC: self.FOUR})[0]
        assert row["missing_read_type"] == "I1, I2"

    def test_an_unplanned_read_is_found_but_not_missing(self):
        keys = flowcell_files() + fastqs(reads=("R1", "I1", "R2"))
        row = summarize_illumina_flowcells(keys, {FC: ["R1", "R2"]})[0]
        assert row["read_types_found"] == "R1, R2, I1"
        assert row["missing_read_type"] == "none"

    def test_a_sidecar_does_not_stand_in_for_its_fastq(self):
        keys = flowcell_files() + fastqs(reads=("R1", "R2"))
        keys.remove(f"{RAW}/{FC}_G1_GEX_S1_L001_R2_001.fastq.gz")
        row = summarize_illumina_flowcells(keys, {FC: ["R1", "R2"]})[0]
        assert row["read_types_found"] == "R1"
        assert row["missing_read_type"] == "R2"

    def test_a_flowcell_with_no_fastqs_misses_every_planned_read(self):
        row = summarize_illumina_flowcells(flowcell_files(), {FC: self.FOUR})[0]
        assert row["read_types_found"] == ""
        assert row["missing_read_type"] == "R1, R2, I1, I2"

    def test_unknown_planned_reads_fall_back_to_r1_r2(self):
        keys = flowcell_files() + fastqs(reads=("R1",))
        assert summarize_illumina_flowcells(keys, {})[0]["missing_read_type"] == "R2"

    def test_modality_lists_each_distinct_modality_once(self):
        keys = (
            flowcell_files()
            + fastqs(lanes=("001", "002"), assay="ADT_HTO", sample="2")
            + fastqs(lanes=("001", "002"), assay="GEX", sample="1")
            + fastqs(lanes=("001",), assay="TCR", sample="3")
            + fastqs(lanes=("001",), assay="GEX", sample="4")
        )
        row = summarize_illumina_flowcells(keys, {FC: ["R1", "R2"]})[0]
        # Unique, in the documented modality order, not file order.
        assert row["modality"] == "GEX, ADT_HTO, TCR"

    def test_modality_is_empty_without_fastqs(self):
        assert summarize_illumina_flowcells(flowcell_files(), {})[0]["modality"] == ""

    def test_unknown_planned_reads_are_named_as_such(self):
        rows = summarize_illumina_flowcells(flowcell_files() + fastqs(), {})
        assert rows[0]["planned_reads"].startswith("UNKNOWN")


class TestIlluminaFastqCounts:
    """Lane counts differ legitimately between modalities; that is not a defect."""

    def _log(self, gex_lanes, cri_lanes):
        return {
            "G1": {
                "GEX": [
                    f"{FC}_G1_GEX_S1_L{lane}_{r}_001.fastq.gz"
                    for lane in gex_lanes
                    for r in ("R1", "I1", "I2", "R2")
                ],
                "CRI": [
                    f"{FC}_G1_CRI_S2_L{lane}_{r}_001.fastq.gz"
                    for lane in cri_lanes
                    for r in ("R1", "I1", "I2", "R2")
                ],
            }
        }

    def test_unequal_lane_counts_are_not_a_mismatch(self):
        log = self._log(("001", "002", "003", "004"), ("001", "002"))
        assert validate_fastq_counts(log, "10x", sequencing_type="Illumina") == []

    def test_raw_file_counts_would_have_reported_a_mismatch(self):
        # Guards the fix: under the Ultima comparison the same delivery fails.
        log = self._log(("001", "002", "003", "004"), ("001", "002"))
        assert validate_fastq_counts(log, "10x") != []

    def test_a_genuinely_absent_modality_library_is_still_reported(self):
        log = self._log(("001",), ("001",))
        log["G1"]["GEX"] += [
            f"{FC}_G1_GEX_S9_L001_{r}_001.fastq.gz" for r in ("R1", "R2")
        ]
        errors = validate_fastq_counts(log, "10x", sequencing_type="Illumina")
        assert len(errors) == 1
        assert "MISMATCH FQ COUNTS" in errors[0]


class TestIndexReadsAreNotLibraryReads:
    """An index read is a sequencing artifact, not a library read."""

    def _metadata(self, reads):
        return {
            f"{RAW}/{FC}_G1_GEX_S1_L001_{r}_001.fastq.gz": {"read_count": 5_000_000}
            for r in reads
        }

    def test_index_reads_are_excluded_from_the_group_total(self):
        counts, errors, _ = validate_read_metadata(
            self._metadata(("R1", "I1", "I2", "R2")), "10x", sequencing_type="Illumina"
        )
        # R1 only: counting I1 and I2 too would report 15M and make every
        # group fail the raw-vs-CellRanger comparison.
        assert counts == {"G1": {"GEX": 5_000_000}}
        assert errors == []

    def test_the_r1_r2_pair_is_still_compared(self):
        meta = self._metadata(("R1", "I1", "I2", "R2"))
        meta[f"{RAW}/{FC}_G1_GEX_S1_L001_R2_001.fastq.gz"]["read_count"] = 4_999_999
        _, errors, _ = validate_read_metadata(meta, "10x", sequencing_type="Illumina")
        assert any("READ COUNT ERROR" in e for e in errors)

    def test_index_files_do_not_create_pairing_complaints(self):
        _, errors, pairing = validate_read_metadata(
            self._metadata(("R1", "I1", "I2", "R2")), "10x", sequencing_type="Illumina"
        )
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []

    def test_ultima_totals_are_unchanged(self):
        # The regression anchor for the shared code path: an Ultima delivery
        # has no index sidecars, so its totals must be exactly as before.
        meta = {
            "p/O/G1/raw/439047-G1_GEX-Z0273-BC01_S1_L001_R1_001.fastq.gz": {
                "read_count": 2_000_000
            },
            "p/O/G1/raw/439047-G1_GEX-Z0273-BC01_S1_L001_R2_001.fastq.gz": {
                "read_count": 2_000_000
            },
        }
        counts, errors, _ = validate_read_metadata(meta, "10x")
        assert counts == {"G1": {"GEX": 2_000_000}}
        assert errors == []


# ---------------------------------------------------------------------------
# Gather: the full S3 walk against a mock Illumina delivery
# ---------------------------------------------------------------------------

from qa_gather import gather_qa_data  # noqa: E402
from tests.test_qa_gather import MockS3Client, _make_ctx  # noqa: E402


class TestIlluminaGather:
    """End-to-end through QADataGatherer, not just the pure checkers."""

    def _delivery(self, flowcells=(FC,), xml=DUAL_INDEX_XML, with_xml=True):
        keys, contents = [], {}
        for fc in flowcells:
            keys += flowcell_files(flowcell=fc)
            if not with_xml:
                keys.remove(f"{RAW}/{fc}_RunParameters.xml")
            else:
                contents[f"{RAW}/{fc}_RunParameters.xml"] = xml
            keys += fastqs(
                reads=("R1", "I1", "I2", "R2"),
                lanes=("001", "002"),
                flowcell=fc,
                sidecars=False,
            )
        return keys, contents

    def _gather(self, keys, contents):
        return gather_qa_data(
            _make_ctx(raw_assay="10x", sequencing_type="Illumina"),
            MockS3Client(keys=keys, file_contents=contents),
        )

    def test_planned_reads_are_read_from_the_xml(self):
        data = self._gather(*self._delivery())
        assert data.illumina_planned_reads == {FC: ["R1", "I1", "I2", "R2"]}

    def test_lanes_are_collected_per_flowcell(self):
        data = self._gather(*self._delivery())
        assert data.illumina_lanes == {FC: {"L001", "L002"}}

    def test_a_well_formed_delivery_raises_no_assay_or_group_errors(self):
        # The regression guard for the suppressed Ultima fallback: without it
        # every flowcell run file parses as a library with an invented assay.
        data = self._gather(*self._delivery())
        assert data.gathering_errors == []

    def test_flowcell_run_files_are_not_counted_as_fastqs(self):
        data = self._gather(*self._delivery())
        assert len(data.fastq_log["G1"]["GEX"]) == 8  # 4 reads x 2 lanes

    def test_missing_xml_warns_and_leaves_planned_reads_empty(self):
        data = self._gather(*self._delivery(with_xml=False))
        assert data.illumina_planned_reads == {}
        assert any("RUN PARAMETERS MISSING" in w for w in data.gathering_warnings)

    def test_unreadable_xml_warns(self):
        keys, contents = self._delivery()
        contents[f"{RAW}/{FC}_RunParameters.xml"] = "<broken"
        data = self._gather(keys, contents)
        assert data.illumina_planned_reads == {}
        assert any("RUN PARAMETERS UNREADABLE" in w for w in data.gathering_warnings)

    def test_two_flowcells_each_get_an_entry(self):
        other = "22ABCDEF9"
        keys, contents = self._delivery(flowcells=(FC, other))
        contents[f"{RAW}/{other}_RunParameters.xml"] = NO_INDEX_XML
        data = self._gather(keys, contents)
        assert data.illumina_planned_reads == {
            FC: ["R1", "I1", "I2", "R2"],
            other: ["R1", "R2"],
        }

    def test_a_multi_part_group_id_raises_no_wrong_group(self):
        """The folder name and the filename's GroupID must agree end to end."""
        raw = "testproj/ORD01/BCP_1_U3/raw"
        fc = "B23YTVYLT4"
        keys = flowcell_files(flowcell=fc, raw_dir=raw) + fastqs(
            reads=("R1", "R2"),
            lanes=("001", "002"),
            group="BCP_1_U3",
            sample="3",
            flowcell=fc,
            raw_dir=raw,
            sidecars=False,
        )
        data = self._gather(keys, {f"{raw}/{fc}_RunParameters.xml": NO_INDEX_XML})
        assert data.gathering_errors == []
        assert len(data.fastq_log["BCP_1_U3"]["GEX"]) == 4

    def test_non_ultima_modalities_raise_no_wrong_assay(self):
        keys = flowcell_files() + [
            k
            for mod, n in (("ADT", "1"), ("HTO", "2"), ("BCR", "3"), ("TCR", "4"))
            for k in fastqs(assay=mod, sample=n, sidecars=False)
        ]
        data = self._gather(keys, {f"{RAW}/{FC}_RunParameters.xml": NO_INDEX_XML})
        assert data.gathering_errors == []
        assert set(data.fastq_log["G1"]) == {"ADT", "HTO", "BCR", "TCR"}

    def test_q30_is_skipped_but_not_silently(self):
        data = self._gather(*self._delivery())
        assert data.pct_q30_values == {}
        # Cell 2 writes warnings into {order}_errors.txt; a skipped check that
        # leaves a clean log is exactly what that rule exists to prevent.
        assert any("Q30 SKIPPED" in w for w in data.gathering_warnings)


class TestFormatReadCount:
    @pytest.mark.parametrize(
        "value, expected",
        [
            (0, "0"),
            (99_999, "99999"),
            (100_000, "0.100M"),
            (250_000, "0.250M"),
            (1_234_567, "1.23M"),
            (12_345_678, "12.3M"),
            (123_456_789, "123M"),
            (999_999_999, "1.00B"),  # rounds up across the unit boundary
            (999_500, "1.00M"),
            (1_234_567_890, "1.23B"),
            (45_678_901_234, "45.7B"),
        ],
    )
    def test_formats(self, value, expected):
        assert format_read_count(value) == expected

    def test_none_and_nan_stay_none(self):
        assert format_read_count(None) is None
        assert format_read_count(float("nan")) is None


class TestTotalReadsByFlowcell:
    def test_sums_every_group_on_a_flowcell(self):
        rows = [
            {"group_id": "G1", "flowcell": "FC_B", "read_count": 4_000_000},
            {"group_id": "G2", "flowcell": "FC_B", "read_count": 6_000_000},
            {"group_id": "G1", "flowcell": "FC_A", "read_count": 1_000_000},
        ]
        assert total_reads_by_flowcell(rows) == [
            {"flowcell": "FC_A", "read_count": 1_000_000},
            {"flowcell": "FC_B", "read_count": 10_000_000},
        ]

    def test_built_from_the_per_group_table(self):
        other_dir = "testproj/ORD01/G2/raw"
        keys = (
            flowcell_files()
            + fastqs()
            + flowcell_files(raw_dir=other_dir)
            + fastqs(raw_dir=other_dir, group="G2")
        )
        meta = {
            name: {
                "read_count": n,
                "errors": [],
                "__source_key": f"{d}/{name}-metadata.json",
            }
            for d, name, n in (
                (RAW, f"{FC}_G1_GEX_S1_L001_R1_001.fastq.gz", 3_000_000),
                (other_dir, f"{FC}_G2_GEX_S1_L001_R1_001.fastq.gz", 5_000_000),
            )
        }
        rows = summarize_illumina_flowcells(keys, {}, meta)
        assert total_reads_by_flowcell(rows) == [
            {"flowcell": FC, "read_count": 8_000_000}
        ]

    def test_empty_table(self):
        assert total_reads_by_flowcell([]) == []
