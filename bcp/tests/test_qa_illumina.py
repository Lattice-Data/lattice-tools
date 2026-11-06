"""
Tests for Illumina delivery support (``sequencing_type="Illumina"``).

The QA pipeline grew up around Ultima (UG100): wafers, the Ultima trimmer, and
fastqs the trimmer emits under the single literal name ``_S1_L001_R1_001``.
An Illumina delivery is a different shape, and this module pins it:

    {proj}/{Order#}/{GroupID}/raw/        -- flat, one directory
        {FC}_CopyComplete.txt             -- eight per-flowcell run files
        ...
        {FC}_RunParameters.xml            -- says which reads to expect
        {FC}-{GroupID}_{Assay}_S{n}_L00#_{READ}_001.fastq.gz
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
    check_extra_raw_files,
    summarize_illumina_flowcells,
    validate_fastq_counts,
    validate_read_metadata,
)
from qa_constants import ALLOWED_RAW_ASSAYS, ILLUMINA_FLOWCELL_FILE_SUFFIXES
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
            f = f"{raw_dir}/{flowcell}-{group}_{assay}_S{sample}_L{lane}_{read}_001.fastq.gz"
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
    NAME = f"{RAW}/{FC}-G1_GEX_S3_L002_R1_001.fastq.gz"

    def test_parses_into_the_five_tuple(self):
        assert parse_raw_filename(self.NAME, "10x", sequencing_type="Illumina") == (
            FC,
            "G1",
            "GEX",
            "S3",
            "L002",
        )

    def test_group_id_containing_a_hyphen(self):
        name = f"{RAW}/{FC}-LAT-123_GEX_S1_L001_R2_001.fastq.gz"
        run, group, assay, _, _ = parse_raw_filename(
            name, "10x", sequencing_type="Illumina"
        )
        # The group must survive intact: _process_raw_file compares it against
        # the folder name and raises WRONG GROUP on a mismatch.
        assert (run, group, assay) == (FC, "LAT-123", "GEX")

    def test_combined_assay_token_is_not_split(self):
        name = f"{RAW}/{FC}-G1_GEX_hash_oligo_S1_L001_R1_001.fastq.gz"
        assert parse_raw_filename(name, "10x", sequencing_type="Illumina")[2] == (
            "GEX_hash_oligo"
        )

    def test_flowcell_run_file_is_not_a_library(self):
        assert (
            parse_raw_filename(
                f"{RAW}/{FC}_RunInfo.xml", "10x", sequencing_type="Illumina"
            )
            is None
        )

    def test_no_ultima_fallback_under_illumina(self):
        # The Ultima fallback splitter would read a hyphenated name as a
        # library and invent an assay token, firing WRONG ASSAY on every
        # flowcell file in the delivery.
        assert (
            parse_raw_filename(
                f"{RAW}/{FC}-Logs.gzip", "10x", sequencing_type="Illumina"
            )
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
        # The one hyphen-separated run file, so it gets its own case.
        keys = flowcell_files(omit=("-Logs.gzip",)) + fastqs()
        _, raw_lost, _ = self._check(keys, planned={FC: ["R1", "R2"]})
        assert any(r.get("-Logs.gzip") == f"{RAW}/{FC}-Logs.gzip" for r in raw_lost)

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
        keys.remove(f"{RAW}/{FC}-G1_GEX_S1_L001_R2_001.fastq.gz-metadata.json")
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
        ultima_stray = f"{RAW}/{FC}-G1_GEX_S1_L001_trimmer-stats.csv"
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
        assert rows[0]["lane_count"] == 4
        assert rows[0]["lanes"] == "L001, L002, L003, L004"
        assert rows[0]["library_count"] == 2

    def test_unknown_planned_reads_are_named_as_such(self):
        rows = summarize_illumina_flowcells(flowcell_files() + fastqs(), {})
        assert rows[0]["planned_reads"].startswith("UNKNOWN")


class TestIlluminaFastqCounts:
    """Lane counts differ legitimately between modalities; that is not a defect."""

    def _log(self, gex_lanes, cri_lanes):
        return {
            "G1": {
                "GEX": [
                    f"{FC}-G1_GEX_S1_L{lane}_{r}_001.fastq.gz"
                    for lane in gex_lanes
                    for r in ("R1", "I1", "I2", "R2")
                ],
                "CRI": [
                    f"{FC}-G1_CRI_S2_L{lane}_{r}_001.fastq.gz"
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
            f"{FC}-G1_GEX_S9_L001_{r}_001.fastq.gz" for r in ("R1", "R2")
        ]
        errors = validate_fastq_counts(log, "10x", sequencing_type="Illumina")
        assert len(errors) == 1
        assert "MISMATCH FQ COUNTS" in errors[0]


class TestIndexReadsAreNotLibraryReads:
    """An index read is a sequencing artifact, not a library read."""

    def _metadata(self, reads):
        return {
            f"{RAW}/{FC}-G1_GEX_S1_L001_{r}_001.fastq.gz": {"read_count": 5_000_000}
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
        meta[f"{RAW}/{FC}-G1_GEX_S1_L001_R2_001.fastq.gz"]["read_count"] = 4_999_999
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

    def test_q30_is_skipped_but_not_silently(self):
        data = self._gather(*self._delivery())
        assert data.pct_q30_values == {}
        # Cell 2 writes warnings into {order}_errors.txt; a skipped check that
        # leaves a clean log is exactly what that rule exists to prevent.
        assert any("Q30 SKIPPED" in w for w in data.gathering_warnings)
