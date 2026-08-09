"""
Tests for SeaHub lab raw QA support (seahub_sci assay).
"""

from __future__ import annotations

import os
import threading

import pytest

from qa_checks import (
    check_expected_raw_files,
    check_extra_raw_files,
    validate_fastq_counts,
    validate_read_metadata,
)
from qa_gather import gather_qa_data
from qa_mods import (
    QARunContext,
    apply_seahub_trim_fail_blocks,
    normalize_raw_assay,
    parse_raw_filename,
    parse_seahub_raw_path,
    parse_seahub_trim_fail_csv,
    resolve_qa_run_context,
    seahub_file_stem,
    seahub_trimmer_failure_storage_key,
    seahub_trimmer_group_storage_key,
)

from tests.test_qa_gather import MockS3Client, _make_ctx

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
QA_FIXTURES_DIR = os.path.join(FIXTURES_DIR, "qa")

SEAHUB_BASE = (
    "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
    "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
)
SEAHUB_KEY_CRAM = f"{SEAHUB_BASE}.trim.cram"
SEAHUB_KEY_FAIL = f"{SEAHUB_BASE}.trim_fail.csv"
SEAHUB_KEY_TRIM = f"{SEAHUB_BASE}.trim.csv"
SEAHUB_TRIM_FAIL = os.path.join(QA_FIXTURES_DIR, "seahub_trim_fail_sample.csv")
SEAHUB_TRIM = os.path.join(QA_FIXTURES_DIR, "seahub_trim_sample.csv")


class TestSeahubConstantsAndContext:
    def test_normalize_seahub_sci(self):
        assert normalize_raw_assay("seahub_sci") == "seahub_sci"

    def test_resolve_context_labalpha(self):
        ctx = resolve_qa_run_context(
            data_source="s3",
            raw_assay="seahub_sci",
            s3_path="s3://czi-labalpha/labalpha-seahub-bcp/REF3",
        )
        assert ctx.bucket == "czi-labalpha"
        assert ctx.provider == "labalpha"
        assert ctx.proj == "labalpha-seahub-bcp"
        assert ctx.order == "REF3"
        assert ctx.listing_prefix == "labalpha-seahub-bcp/REF3/"


class TestSeahubParsing:
    def test_parse_seahub_raw_path(self):
        info = parse_seahub_raw_path(SEAHUB_KEY_CRAM)
        assert info == {
            "experiment_id": "REF3",
            "sublibrary": "P05_1",
            "wafer": "430479",
        }

    def test_seahub_file_stem(self):
        assert (
            seahub_file_stem(
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim.cram"
            )
            == "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        )
        assert (
            seahub_file_stem(
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv"
            )
            == "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        )

    def test_parse_raw_filename_seahub_cram(self):
        # Group is the folder sublibrary, not the ExperimentID: the folder value
        # is SOP-anchored and stable, whereas returning the ExperimentID lumped
        # every sublibrary of an experiment into one row.
        assert parse_raw_filename(SEAHUB_KEY_CRAM, "seahub_sci") == (
            "430479",
            "P05_1",
            "GEX_hash_oligo",
            "Z0097",
            "CAGTCAGTTGCAGAT",
        )

    def test_parse_raw_filename_seahub_trim_fail(self):
        assert parse_raw_filename(SEAHUB_KEY_FAIL, "seahub_sci") == (
            "430479",
            "P05_1",
            "GEX_hash_oligo",
            "Z0097",
            "CAGTCAGTTGCAGAT",
        )

    def test_storage_keys(self):
        assert seahub_trimmer_group_storage_key(SEAHUB_KEY_FAIL) == "REF3/P05_1"
        assert seahub_trimmer_failure_storage_key(SEAHUB_KEY_FAIL) == (
            "430479",
            "430479",
        )

    def test_parse_seahub_raw_path_rejects_wrong_segment_count(self):
        # Missing the wafer level -> only 5 segments.
        assert parse_seahub_raw_path("proj/REF3/raw/P05_1/file.trim.cram") is None
        # "raw" not at the documented position.
        assert (
            parse_seahub_raw_path("proj/REF3/notraw/P05_1/430479/file.trim.cram")
            is None
        )

    def test_parse_seahub_raw_path_anchors_by_position_not_name(self):
        # A segment other than index 2 is literally "raw" -> still parses
        # correctly because the lookup is positional, not a name search.
        info = parse_seahub_raw_path("proj/raw/raw/P05_1/430479/file.trim.cram")
        assert info == {
            "experiment_id": "raw",
            "sublibrary": "P05_1",
            "wafer": "430479",
        }


class TestSeahubTrimAdapters:
    def test_grab_seahub_trim_fail_csv_is_per_format(self):
        """One trimmer-fail value per format block, each < 100% (regression:
        the shared single-denominator parser summed both modalities' failures
        over one modality's total and produced >100%)."""
        stats: dict = {}
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL), stats, "430479"
        )
        assert "430479" in stats
        # SeaHub inputs start from RSQ-passing reads -> no RSQ rows.
        assert stats["430479"]["rsq"] == []
        # Two format blocks (JumboSciHash, JumboSciGEX) -> two values.
        fails = stats["430479"]["trimmer_fail"]
        assert len(fails) == 2
        assert all(0 <= pct < 100 for pct in fails)

    def test_grab_seahub_trim_fail_csv_denominator_per_group(self):
        """Each modality's failures divide by its own total read count."""
        stats: dict = {}
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL), stats, "430479"
        )
        hash_pct, gex_pct = stats["430479"]["trimmer_fail"]
        assert hash_pct == 100 * 157722381 / 260527531
        assert gex_pct == 100 * 100735571 / 158233602

    def test_grab_seahub_trim_fail_csv_appends(self):
        """Repeated wells accumulate into the same distribution."""
        stats: dict = {}
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL), stats, "430479"
        )
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL), stats, "430479"
        )
        assert len(stats["430479"]["trimmer_fail"]) == 4

    def _write_inconsistent_totals_csv(self, tmp_path) -> str:
        path = os.path.join(tmp_path, "inconsistent_totals.csv")
        with open(path, "w") as f:
            f.write("format,reason,failed_read_count,total_read_count\n")
            f.write("JumboSciHash,no match,100,1000\n")
            f.write("JumboSciHash,no match,50,2000\n")
        return path

    def test_grab_seahub_trim_fail_csv_warns_on_inconsistent_totals(self, tmp_path):
        """A format block with a non-uniform total_read_count still falls
        back to the first value, but surfaces the anomaly via `warnings`
        (the list the QA notebook prints to the user) instead of silently
        dropping the divergent totals."""
        path = self._write_inconsistent_totals_csv(tmp_path)
        stats: dict = {}
        warnings: list[str] = []
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(path, warnings=warnings), stats, "exp1"
        )
        assert len(warnings) == 1
        assert "JumboSciHash" in warnings[0]
        assert "1000" in warnings[0] and "2000" in warnings[0]
        # Fallback behavior: first total_read_count still used.
        assert stats["exp1"]["trimmer_fail"] == [100 * 150 / 1000]

    def test_grab_seahub_trim_fail_csv_no_warnings_list_is_fine(self, tmp_path):
        """Omitting `warnings` (the default for every pre-existing call) must
        not raise, even when totals are inconsistent."""
        path = self._write_inconsistent_totals_csv(tmp_path)
        stats: dict = {}
        apply_seahub_trim_fail_blocks(parse_seahub_trim_fail_csv(path), stats, "exp1")
        assert stats["exp1"]["trimmer_fail"] == [100 * 150 / 1000]


class TestSeahubChecks:
    def _full_stem_files(self, stem: str, raw_dir: str) -> list[str]:
        suffixes = [
            ".trim.cram",
            ".trim.csv",
            ".trim.stderr",
            ".trim.stdout",
            ".trim_fail.csv",
        ]
        return [f"{raw_dir}/{stem}{suf}" for suf in suffixes]

    def test_expected_raw_files_complete(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir)
        good, lost, found = check_expected_raw_files(all_raw, "seahub_sci")
        assert good == 1
        assert lost == []
        assert len(found) == 5

    def test_expected_raw_files_missing_suffix(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir)[:-1]
        _good, lost, _found = check_expected_raw_files(all_raw, "seahub_sci")
        assert len(lost) == 1
        assert ".trim_fail.csv" in lost[0]

    def test_extra_raw_allows_metadata(self):
        raw_dir = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        stem = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
        all_raw = self._full_stem_files(stem, raw_dir) + [
            f"{raw_dir}/{stem}.trim.cram-metadata.json"
        ]
        _good, _lost, found = check_expected_raw_files(all_raw, "seahub_sci")
        extra = check_extra_raw_files(all_raw, found, "seahub_sci")
        assert extra == []

    def test_fastq_counts_gex_hash_only(self):
        fastq_log = {
            "REF3": {
                "GEX_hash_oligo": ["a.trim.cram", "b.trim.cram"],
            }
        }
        assert validate_fastq_counts(fastq_log, "seahub_sci") == []

    def test_read_metadata_skips_r1_r2_pairing(self):
        meta = {
            SEAHUB_KEY_CRAM: {
                "read_count": 2_000_000,
                "__reported_filename": SEAHUB_KEY_CRAM.split("/")[-1],
            }
        }
        counts, errors, pairing = validate_read_metadata(meta, "seahub_sci")
        assert pairing["r1_without_r2_metadata"] == []
        assert pairing["r2_without_r1_metadata"] == []


class TestSeahubGather:
    def _seahub_keys(self) -> list[str]:
        """One well at SOP depth, as a key list.

        A key list rather than a hand-written page map: a page map states the
        listing's *answers*, so it can only describe objects the walk already
        looked for. Every SeaHub S3 test used one, which is why no test could
        express an object above wafer depth and why the gatherer silently
        dropped them.
        """
        return [
            SEAHUB_KEY_CRAM,
            SEAHUB_KEY_FAIL,
            SEAHUB_KEY_TRIM,
            f"{SEAHUB_BASE}.trim.stderr",
            f"{SEAHUB_BASE}.trim.stdout",
        ]

    def test_gather_from_s3_seahub_layout(self):
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            SEAHUB_KEY_TRIM: open(SEAHUB_TRIM).read(),
        }
        s3 = MockS3Client(keys=self._seahub_keys(), file_contents=file_contents)
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        data = gather_qa_data(ctx, s3)
        assert data.has_raw is True
        assert len(data.all_raw_files) == 5
        assert "P05_1" in data.fastq_log
        assert "GEX_hash_oligo" in data.fastq_log["P05_1"]
        assert "430479" in data.trimmer_failure_stats
        assert "REF3/P05_1" in data.group_failure_stats
        assert data.discovered_wafers == {"430479"}

    def test_gather_manifest_enriches_seahub(self):
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        ctx = self._manifest_ctx(manifest)
        trim_fail_key = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv"
        )
        trim_key = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim.csv"
        )
        trim_fail_a2 = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim_fail.csv"
        )
        trim_a2 = (
            "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
            "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim.csv"
        )
        file_contents = {
            trim_fail_key: open(SEAHUB_TRIM_FAIL).read(),
            trim_key: open(SEAHUB_TRIM).read(),
            trim_fail_a2: open(SEAHUB_TRIM_FAIL).read(),
            trim_a2: open(SEAHUB_TRIM).read(),
        }
        s3 = MockS3Client(file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        assert len(data.all_raw_files) == 10
        assert "430479" in data.trimmer_failure_stats

    def _manifest_ctx(self, manifest: str, **kwargs) -> QARunContext:
        """Resolve a manifest context the way the notebook does.

        Deliberately not hand-built: a hand-built context can assert an ``order``
        the resolver is incapable of producing, which is exactly how manifest mode
        shipped with an empty ExperimentID under a green suite.
        """
        return resolve_qa_run_context(
            data_source="manifest",
            raw_assay="seahub_sci",
            manifest_path=manifest,
            manifest_delimiter="\t",
            manifest_s3_column=6,
            manifest_has_header=True,
            run_label="REF3",
            **kwargs,
        )

    def _missing_metadata_warnings(self, ctx) -> list[str]:
        file_contents = {
            k: open(SEAHUB_TRIM_FAIL).read()
            for k in (
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv",
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim_fail.csv",
            )
        }
        data = gather_qa_data(ctx, MockS3Client(file_contents=file_contents))
        return [w for w in data.gathering_warnings if w.startswith("METADATA MISSING")]

    def test_missing_metadata_json_warns_once_with_a_count(self):
        """One warning for the whole upload, not one per CRAM.

        CZI generates these sidecars for the upload as a whole, so absence is an
        upload-wide fact: the per-well form printed 336 identical lines on REF3
        and 864 on GENE7, burying every other warning.
        """
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        warnings = self._missing_metadata_warnings(self._manifest_ctx(manifest))

        assert len(warnings) == 1
        assert "2 CRAM(s) have no matching metadata sidecar" in warnings[0]

    def test_missing_metadata_warning_names_the_delivered_spelling(self):
        """A bare-family well must not be told about .trim.* files it never had."""
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        warnings = self._missing_metadata_warnings(self._manifest_ctx(manifest))

        assert "2x .trim.cram-metadata.json" in warnings[0]
        assert ".trim.cram" in warnings[0]

    def test_missing_metadata_warning_carries_examples(self):
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        warnings = self._missing_metadata_warnings(self._manifest_ctx(manifest))

        assert (
            "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT" in (warnings[0])
        )
        # Two examples is the cap; with only two missing there is no "more".
        assert "more" not in warnings[0]

    def test_many_missing_sidecars_are_still_one_warning(self, tmp_path):
        """The case that motivated this: a whole plate with no sidecars."""
        base = "labalpha-seahub-bcp/REF3/raw/P05_1/430479"
        keys = [
            f"{base}/430479-REF3_P05_1_A{i}_GEX_hash_oligo-Z{i:04d}-CAGTCAGTTGCAGAT{s}"
            for i in range(1, 21)
            for s in (".trim.cram", ".trim.csv", ".trim.stderr", ".trim.stdout")
        ]
        manifest = tmp_path / "listing.tsv"
        manifest.write_text(
            "S3_Full_Path\n" + "".join(f"s3://czi-labalpha/{k}\n" for k in keys)
        )
        ctx = resolve_qa_run_context(
            data_source="manifest",
            raw_assay="seahub_sci",
            manifest_path=str(manifest),
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=True,
            run_label="REF3",
        )
        data = gather_qa_data(ctx, MockS3Client())
        warnings = [
            w for w in data.gathering_warnings if w.startswith("METADATA MISSING")
        ]

        assert len(warnings) == 1
        assert "20 CRAM(s)" in warnings[0]
        assert "and 18 more" in warnings[0]

    def test_present_metadata_json_no_warn(self):
        """A .trim.cram with its .trim.cram-metadata.json sidecar does not warn."""
        meta_key = f"{SEAHUB_BASE}.trim.cram-metadata.json"
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            meta_key: '{"read_count": 1000, "filename": "%s"}'
            % SEAHUB_KEY_CRAM.split("/")[-1],
        }
        s3 = MockS3Client(
            keys=self._seahub_keys() + [meta_key], file_contents=file_contents
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
        missing = [
            w for w in data.gathering_warnings if w.startswith("METADATA MISSING")
        ]
        assert missing == []


class TestSeahubGatherShallowObjects:
    """Objects above wafer depth must reach the SOP table.

    The S3 walk read only ``CommonPrefixes`` from the ``raw/`` and sublibrary
    listings and threw their ``Contents`` away, so an object sitting directly in
    ``{exp}/raw/`` or in ``{exp}/raw/{sublibrary}/`` never entered
    ``all_raw_files``. ``bad_path_depth`` could therefore only fire for keys that
    were too *deep* -- never too shallow, which is the commoner human error --
    and S3 mode and manifest mode disagreed about what the bucket held. Measured
    on a real GENE7 listing: a 5.3 GB extensionless object at
    ``GENE7/raw/P10/436516`` is flagged from a manifest and invisible from S3.
    """

    PROJ = "labalpha-seahub-bcp"
    PREFIX = "labalpha-seahub-bcp/REF3/"
    DEEP = [SEAHUB_KEY_CRAM]

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj=self.PROJ,
            order="REF3",
            listing_prefix=self.PREFIX,
        )

    def _gather(self, keys: list[str]):
        return gather_qa_data(self._ctx(), MockS3Client(keys=keys))

    def test_object_at_sublibrary_level_is_gathered(self):
        stray = f"{self.PROJ}/REF3/raw/P05_1/436516"
        data = self._gather(self.DEEP + [stray])
        assert stray in data.all_raw_files
        assert stray in data.raw_file_sizes
        assert "bad_path_depth" in {v["type"] for v in data.sop_violations}

    def test_object_directly_in_raw_is_gathered(self):
        stray = f"{self.PROJ}/REF3/raw/loose_notes.txt"
        data = self._gather(self.DEEP + [stray])
        assert stray in data.all_raw_files
        assert "bad_path_depth" in {v["type"] for v in data.sop_violations}

    def test_raw_holding_only_loose_objects_is_not_called_missing(self):
        """Objects present but no wafer folder is a depth problem, not an absence."""
        stray = f"{self.PROJ}/REF3/raw/loose_notes.txt"
        data = self._gather([stray])
        assert data.has_raw is True
        assert data.all_raw_files == [stray]
        assert [w for w in data.gathering_warnings if "MISSING or empty" in w] == []

    def test_genuinely_empty_raw_still_warns(self):
        data = self._gather([f"{self.PROJ}/REF3/processed/x.bam"])
        assert data.has_raw is False
        assert [w for w in data.gathering_warnings if "MISSING or empty" in w]

    def test_wafers_are_still_discovered_from_the_folder_walk(self):
        data = self._gather(self.DEEP + [f"{self.PROJ}/REF3/raw/P05_1/436516"])
        assert data.discovered_wafers == {"430479"}

    def test_s3_and_manifest_agree_on_one_key_set(self, tmp_path):
        """The two modes are the same question asked twice; they must not differ."""
        keys = self.DEEP + [
            f"{self.PROJ}/REF3/raw/P05_1/436516",
            f"{self.PROJ}/REF3/raw/loose_notes.txt",
        ]
        manifest = tmp_path / "listing.tsv"
        manifest.write_text(
            "S3_Full_Path\n" + "".join(f"s3://czi-labalpha/{k}\n" for k in keys)
        )
        from_s3 = self._gather(keys)
        from_manifest = gather_qa_data(
            resolve_qa_run_context(
                data_source="manifest",
                raw_assay="seahub_sci",
                manifest_path=str(manifest),
                manifest_delimiter="\t",
                manifest_s3_column=0,
                manifest_has_header=True,
                run_label="REF3",
            ),
            MockS3Client(),
        )
        assert sorted(from_s3.all_raw_files) == sorted(from_manifest.all_raw_files)
        assert sorted(v["type"] for v in from_s3.sop_violations) == sorted(
            v["type"] for v in from_manifest.sop_violations
        )


class TestSeahubTrimFailIsParsedOnce:
    """Each per-well fail CSV is read from disk once, not once per distribution.

    The counts feed two distributions -- wafer-level and sublibrary-level -- and
    the gatherer used to call the combined parse-and-apply helper for each,
    reading the same downloaded file twice. On an upload with one of these per
    well that is one wasted parse per well.
    """

    def test_one_read_csv_per_downloaded_file(self, monkeypatch):
        import qa_mods

        reads: list[str] = []
        real_read_csv = qa_mods.pd.read_csv

        def counting_read_csv(path, *args, **kwargs):
            reads.append(str(path))
            return real_read_csv(path, *args, **kwargs)

        monkeypatch.setattr(qa_mods.pd, "read_csv", counting_read_csv)

        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )
        s3 = MockS3Client(
            keys=[SEAHUB_KEY_CRAM, SEAHUB_KEY_FAIL],
            file_contents={SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read()},
        )
        data = gather_qa_data(ctx, s3)

        assert len(reads) == 1
        # and both distributions still got the numbers
        assert data.trimmer_failure_stats["430479"]["trimmer_fail"]
        assert data.group_failure_stats["REF3/P05_1"]["trimmer_fail"]

    def test_fail_counts_survive_the_split(self):
        from qa_mods import apply_seahub_trim_fail_blocks, parse_seahub_trim_fail_csv

        counts: dict = {}
        apply_seahub_trim_fail_blocks(
            parse_seahub_trim_fail_csv(SEAHUB_TRIM_FAIL),
            {},
            "w1",
            fail_counts=counts,
            stem_key="stem1",
        )

        assert set(counts["stem1"]) == {"JumboSciHash", "JumboSciGEX"}
        assert counts["stem1"]["JumboSciGEX"]["total"] == 158233602


class TestSeahubIngestGates:
    """What reaches the trimmer-fail histograms, and what fastq_log is keyed by."""

    PROJ = "labalpha-seahub-bcp"
    DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479"
    STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj=self.PROJ,
            order="REF3",
            listing_prefix=f"{self.PROJ}/REF3/",
        )

    def _csv(self, total: int) -> str:
        return (
            "format,reason,failed read count,total read count\n"
            f"JumboSciGEX,barcode,{total // 10},{total}\n"
        )

    def _run(self, extra: dict[str, str]):
        keys = [
            f"{self.DIR}/{self.STEM}.trim.cram",
            f"{self.DIR}/{self.STEM}.trim_fail.csv",
        ]
        contents = {f"{self.DIR}/{self.STEM}.trim_fail.csv": self._csv(1000)}
        keys += list(extra)
        contents.update(extra)
        return gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

    def test_the_experiment_id_is_not_a_fastq_log_key(self):
        """It was seeded and never filled: fastq_log is keyed by sublibrary."""
        data = self._run({})

        assert "REF3" not in data.fastq_log
        assert "REF3_P05_1" in data.fastq_log

    def test_fastq_log_keys_match_between_the_two_modes(self, tmp_path):
        """The seed was the only thing that differed; manifest mode never made one."""
        keys = [f"{self.DIR}/{self.STEM}{s}" for s in (".trim.cram", ".trim_fail.csv")]
        contents = {f"{self.DIR}/{self.STEM}.trim_fail.csv": self._csv(1000)}
        manifest = tmp_path / "listing.tsv"
        manifest.write_text(
            "S3_Full_Path\n" + "".join(f"s3://czi-labalpha/{k}\n" for k in keys)
        )
        from_manifest = gather_qa_data(
            resolve_qa_run_context(
                data_source="manifest",
                raw_assay="seahub_sci",
                manifest_path=str(manifest),
                manifest_delimiter="\t",
                manifest_s3_column=0,
                manifest_has_header=True,
                run_label="REF3",
            ),
            MockS3Client(file_contents=contents),
        )
        from_s3 = gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

        assert sorted(from_s3.fastq_log) == sorted(from_manifest.fastq_log)

    def test_a_bare_fail_csv_with_an_unparseable_stem_is_not_ingested(self):
        """`_fail.csv` is generic, so the bare family needs its stem to parse.

        Without the gate a file the SOP table reports as `unexpected_suffix`
        still contributed its numbers to the wafer and sublibrary histograms.
        """
        junk = f"{self.DIR}/README_fail.csv"

        data = self._run({junk: self._csv(999999)})

        assert data.trimmer_failure_stats["430479"]["trimmer_fail"] == [100 / 10]

    def test_a_trim_fail_csv_with_an_odd_stem_is_still_ingested(self):
        """The asymmetry is deliberate: `.trim.*` is distinctive enough that a
        malformed stem is a stem defect, not a different kind of file."""
        odd = f"{self.DIR}/README.trim_fail.csv"

        data = self._run({odd: self._csv(2000)})

        assert len(data.trimmer_failure_stats["430479"]["trimmer_fail"]) == 2


class TestFailCountsAreKeyedByFolderAndStem:
    """A bare stem is not unique across sublibrary folders.

    The reconciliation reads these to decide whether the trimmer consumed the
    whole delivered file, and it looks them up per well. Keyed on the stem
    alone, two wells in different folders that share one merged their
    per-format counts, so one well answered with the other's totals.
    """

    PROJ = "labalpha-seahub-bcp"
    STEM = "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT"
    FOLDERS = ("REF3_P05_1", "REF3_P05_2")

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj=self.PROJ,
            order="REF3",
            listing_prefix=f"{self.PROJ}/REF3/",
        )

    def _csv(self, total: int) -> str:
        return (
            "format,reason,failed read count,total read count\n"
            f"JumboSciGEX,barcode,{total // 10},{total}\n"
        )

    def _dirs(self) -> list[str]:
        return [f"{self.PROJ}/REF3/raw/{f}/430479" for f in self.FOLDERS]

    def test_one_entry_per_folder_for_a_shared_stem(self):
        """The same wafer and stem under two sublibrary folders."""
        dirs = self._dirs()
        keys = [
            f"{d}/{self.STEM}{s}"
            for d in dirs
            for s in (".trim.cram", ".trim_fail.csv")
        ]
        contents = {
            f"{d}/{self.STEM}.trim_fail.csv": self._csv(1000 * (n + 1))
            for n, d in enumerate(dirs)
        }

        data = gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

        assert sorted(data.seahub_fail_counts) == [(d, self.STEM) for d in sorted(dirs)]

    def test_the_counts_are_not_merged_across_folders(self):
        dirs = self._dirs()
        keys = [
            f"{d}/{self.STEM}{s}"
            for d in dirs
            for s in (".trim.cram", ".trim_fail.csv")
        ]
        contents = {
            f"{d}/{self.STEM}.trim_fail.csv": self._csv(1000 * (n + 1))
            for n, d in enumerate(dirs)
        }

        data = gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

        totals = {
            key: entry["JumboSciGEX"]["total"]
            for key, entry in data.seahub_fail_counts.items()
        }
        assert sorted(totals.values()) == [1000, 2000]

    def test_the_reconciliation_finds_them_under_the_same_key(self):
        """Writer and reader have to agree, or every well reads as unavailable."""
        from qa_seahub_recon import reconcile_trimming
        from qa_seahub_source import SourceEntry, index_trimmed_upload

        d = self._dirs()[0]
        keys = [f"{d}/{self.STEM}{s}" for s in (".trim.cram", ".trim_fail.csv")]
        contents = {f"{d}/{self.STEM}.trim_fail.csv": self._csv(1000)}
        data = gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )
        source = {
            ("430479", "Z0097"): SourceEntry(
                wafer="430479",
                ug="Z0097",
                barcode="CAGTCAGTTGCAGAT",
                group="REF3_P05_1_A1",
                assay="GEX_hash_oligo",
                cram_key="v/REF3/raw/430479/x.cram",
                read_count=1000,
            )
        }

        report = reconcile_trimming(
            source, index_trimmed_upload(keys), data.seahub_fail_counts
        )

        assert [r["category"] for r in report.rows] == []


class TestOneBadMetadataSidecarDoesNotKillTheRun:
    """864 sidecars go through the pool for a seahub_sci run.

    ``_should_download_metadata_json`` gained ``seahub_sci``, so this path
    carries one object per well now. One truncated or non-JSON sidecar used to
    raise out of ``future.result()`` and abort the whole gather. A sidecar that
    cannot be read leaves the well with no read count -- the same state a
    missing one leaves it in -- so it degrades instead, with the failures
    collapsed into one warning.
    """

    DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P05_1/430479"
    STEMS = (
        "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT",
        "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT",
    )

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
        )

    def _run(self, bad: str | None):
        keys = [
            f"{self.DIR}/{stem}{s}"
            for stem in self.STEMS
            for s in (".trim.cram", ".trim.cram-metadata.json")
        ]
        good = '{"read_count": 1000, "filename": "x.trim.cram"}'
        contents = {
            f"{self.DIR}/{self.STEMS[0]}.trim.cram-metadata.json": good,
        }
        if bad is not None:
            contents[f"{self.DIR}/{self.STEMS[1]}.trim.cram-metadata.json"] = bad
        return gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

    @pytest.mark.parametrize(
        "payload,label",
        [
            (None, "object absent"),
            ("", "zero bytes"),
            ("{not json at all", "truncated"),
            ("[1, 2, 3]", "JSON, but not an object"),
        ],
    )
    def test_the_good_sidecar_still_lands(self, payload, label):
        data = self._run(payload)

        assert len(data.read_metadata) == 1, label

    @pytest.mark.parametrize("payload", [None, "", "{not json at all", "[1, 2, 3]"])
    def test_the_failure_is_reported_once(self, payload):
        data = self._run(payload)

        warnings = [
            w for w in data.gathering_warnings if w.startswith("METADATA UNREADABLE")
        ]
        assert len(warnings) == 1
        assert "1 object(s) could not be read" in warnings[0]


class TestSeahubTrimFailIsFetchedInParallel:
    """The per-well trim failure CSVs share the metadata path's thread pool.

    There is one per well, so fetching them inside the walk was 336 sequential
    round-trips on REF3 and 864 on GENE7 -- while the ``.cram-metadata.json``
    sidecars sitting beside them already went through a 16-worker pool.

    Downloading and parsing run in the pool; applying runs single-threaded in
    listing order. That is what keeps the concurrency lock-free, and what keeps
    the output identical rather than merely equivalent: the structures these
    feed are appended to, so a parallel apply would reorder them run to run.
    """

    PROJ = "labalpha-seahub-bcp"
    WAFERS = ("430479", "430480", "430481", "430482")

    def _keys(self) -> list[str]:
        """One well per wafer, so every well has its own fail CSV."""
        return [
            f"{self.PROJ}/REF3/raw/REF3_P05_1/{wafer}/"
            f"{wafer}-REF3_P05_1_A{i}_GEX_hash_oligo-Z{i:04d}-CAGTCAGTTGCAGAT{suffix}"
            for i, wafer in enumerate(self.WAFERS, start=1)
            for suffix in (".trim.cram", ".trim_fail.csv")
        ]

    def _fail_keys(self, keys: list[str]) -> list[str]:
        return [k for k in keys if k.endswith(".trim_fail.csv")]

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj=self.PROJ,
            order="REF3",
            listing_prefix=f"{self.PROJ}/REF3/",
        )

    def _fail_csv(self, total: int) -> str:
        return (
            "format,reason,failed read count,total read count\n"
            f"JumboSciGEX,barcode,{total // 7},{total}\n"
        )

    def test_the_downloads_overlap(self):
        """A barrier, not a stopwatch: serial code cannot reach it and hangs.

        Two fail-CSV downloads must be in flight at once. If they are fetched
        one at a time the second never arrives, the barrier times out, and this
        fails deterministically rather than flakily.
        """
        # The barrier releases in pairs, so an odd number of files would leave
        # the last thread waiting alone until the timeout.
        assert len(self.WAFERS) % 2 == 0
        keys = self._keys()
        barrier = threading.Barrier(2, timeout=10)
        reached = []

        class BarrierS3(MockS3Client):
            def download_file(self, bucket, key, local):
                if key.endswith(".trim_fail.csv"):
                    barrier.wait()
                    reached.append(key)
                return super().download_file(bucket, key, local)

        s3 = BarrierS3(
            keys=keys,
            file_contents={k: self._fail_csv(1000) for k in self._fail_keys(keys)},
        )
        gather_qa_data(self._ctx(), s3)

        assert len(reached) == len(self.WAFERS)

    def test_every_fail_csv_is_fetched_exactly_once(self):
        keys = self._keys()
        fetched: list[str] = []
        lock = threading.Lock()

        class CountingS3(MockS3Client):
            def download_file(self, bucket, key, local):
                with lock:
                    fetched.append(key)
                return super().download_file(bucket, key, local)

        s3 = CountingS3(
            keys=keys,
            file_contents={k: self._fail_csv(1000) for k in self._fail_keys(keys)},
        )
        gather_qa_data(self._ctx(), s3)

        assert sorted(fetched) == sorted(self._fail_keys(keys))

    def test_the_distributions_follow_listing_order(self):
        """Applied serially in listing order, so the appended values do not
        reorder between runs even though the fetches finish out of order."""
        keys = self._keys()
        totals = {k: 1000 * (n + 1) for n, k in enumerate(self._fail_keys(keys))}
        s3 = MockS3Client(
            keys=keys,
            file_contents={k: self._fail_csv(t) for k, t in totals.items()},
        )

        data = gather_qa_data(self._ctx(), s3)

        # One sublibrary, so every well lands in one distribution in order.
        by_group = data.group_failure_stats["REF3/REF3_P05_1"]["trimmer_fail"]
        expected = [100 * (t // 7) / t for t in totals.values()]
        assert by_group == expected

    def test_warnings_are_merged_in_listing_order(self):
        """Each worker collects its own, so they cannot interleave by timing."""
        keys = self._keys()
        inconsistent = (
            "format,reason,failed read count,total read count\n"
            "JumboSciGEX,barcode,10,{a}\n"
            "JumboSciGEX,barcode,20,{b}\n"
        )
        contents = {
            k: inconsistent.format(a=1000 * (n + 1), b=2000 * (n + 1))
            for n, k in enumerate(self._fail_keys(keys))
        }
        s3 = MockS3Client(keys=keys, file_contents=contents)

        data = gather_qa_data(self._ctx(), s3)

        totals = [
            w.split("values [")[1].split(",")[0]
            for w in data.gathering_warnings
            if w.startswith("INCONSISTENT TOTAL")
        ]
        assert totals == ["1000", "2000", "3000", "4000"]

    def test_a_single_fail_csv_needs_no_pool(self):
        keys = self._keys()[:2]
        s3 = MockS3Client(
            keys=keys,
            file_contents={k: self._fail_csv(1000) for k in self._fail_keys(keys)},
        )

        data = gather_qa_data(self._ctx(), s3)

        assert data.trimmer_failure_stats["430479"]["trimmer_fail"]

    @pytest.mark.parametrize(
        "payload,label",
        [
            (None, "object absent"),
            # The two shapes an interrupted upload actually leaves behind.
            ("", "zero bytes -> EmptyDataError"),
            ('format,reason\n"half-writ', "truncated mid-write -> ParserError"),
        ],
    )
    def test_one_unreadable_fail_csv_does_not_abort_the_run(self, payload, label):
        """There is one of these per well, so aborting loses the other 863.

        This reverses the earlier call to let it raise. The reasoning then was
        that no graceful category existed -- but a well with no trimmer counts
        already reconciles as ``metadata_unavailable``, and an interrupted
        upload is exactly what leaves a zero-byte or half-written CSV behind.
        """
        keys = self._keys()
        fail_keys = self._fail_keys(keys)
        contents = {k: self._fail_csv(1000) for k in fail_keys}
        if payload is None:
            del contents[fail_keys[2]]
        else:
            contents[fail_keys[2]] = payload

        data = gather_qa_data(
            self._ctx(), MockS3Client(keys=keys, file_contents=contents)
        )

        # the three readable wells still landed
        assert len(data.trimmer_failure_stats) == len(self.WAFERS) - 1, label
        warnings = [
            w for w in data.gathering_warnings if w.startswith("TRIM FAIL UNREADABLE")
        ]
        assert len(warnings) == 1
        assert "1 object(s) could not be read" in warnings[0]

    def test_many_unreadable_fail_csvs_are_one_warning(self):
        keys = self._keys()
        s3 = MockS3Client(keys=keys, file_contents={})

        data = gather_qa_data(self._ctx(), s3)

        warnings = [
            w for w in data.gathering_warnings if w.startswith("TRIM FAIL UNREADABLE")
        ]
        assert len(warnings) == 1
        assert f"{len(self.WAFERS)} object(s) could not be read" in warnings[0]
        assert "and 2 more" in warnings[0]


class TestSeahubListingTruncation:
    """A listing cut short at 1000 entries must say so.

    ``list_objects`` reports the fact only in ``IsTruncated``, which nothing
    read: past the cap the remaining folders did not exist as far as QA was
    concerned and the run still looked clean. Unreachable with anything observed
    so far -- REF3 has 7 sublibraries, GENE7 has 9 -- but the failure mode is
    silence, which is the wrong one for a limit nobody thinks to check.
    """

    PROJ = "labalpha-seahub-bcp"

    def _ctx(self):
        return _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj=self.PROJ,
            order="REF3",
            listing_prefix=f"{self.PROJ}/REF3/",
        )

    def _gather(self, keys: list[str], limit: int | None):
        return gather_qa_data(self._ctx(), MockS3Client(keys=keys, list_limit=limit))

    def _wells(self, sublibraries: int, wafers: int) -> list[str]:
        return [
            f"{self.PROJ}/REF3/raw/P{s:02d}/{430000 + w}/"
            f"{430000 + w}-REF3_P{s:02d}_A1_GEX_hash_oligo-Z{w:04d}-CAGTCAGTTGCAGAT"
            ".trim.cram"
            for s in range(1, sublibraries + 1)
            for w in range(1, wafers + 1)
        ]

    def _truncation_warnings(self, data) -> list[str]:
        return [w for w in data.gathering_warnings if w.startswith("LISTING TRUNCATED")]

    def test_truncated_sublibrary_listing_warns(self):
        data = self._gather(self._wells(sublibraries=3, wafers=1), limit=2)

        assert self._truncation_warnings(data)
        assert "sublibrary folders" in self._truncation_warnings(data)[0]

    def test_truncated_wafer_listing_warns(self):
        data = self._gather(self._wells(sublibraries=1, wafers=3), limit=2)

        assert any("wafer folders" in w for w in self._truncation_warnings(data))

    def test_truncated_top_listing_warns(self):
        """The site the review did not name: this one hides `processed/`."""
        keys = self._wells(sublibraries=1, wafers=1) + [
            f"{self.PROJ}/REF3/processed/x.bam",
            f"{self.PROJ}/REF3/other/y.txt",
        ]
        data = self._gather(keys, limit=1)

        assert any(
            "the experiment folder" in w for w in self._truncation_warnings(data)
        )

    def test_an_untruncated_listing_is_silent(self):
        data = self._gather(self._wells(sublibraries=3, wafers=3), limit=None)

        assert self._truncation_warnings(data) == []
        assert len(data.discovered_wafers) == 3

    def test_truncation_no_longer_costs_raw_objects(self):
        """Objects are paginated since the recursive walk, so only folders go."""
        keys = self._wells(sublibraries=3, wafers=1)
        data = self._gather(keys, limit=1)

        assert sorted(data.all_raw_files) == sorted(keys)
        assert self._truncation_warnings(data)


class TestSeahubManifestExperimentId:
    """Manifest mode must name the ExperimentID it is checking.

    ``resolve_qa_run_context`` used to hard-code ``order=""`` for every manifest
    run and drop the ``order=`` argument on the floor. The empty value then
    reached the cross-experiment check as the *expected* ExperimentID, making it
    unequal to every object's, and disabled ``index_untrimmed_sources``' own
    ``if experiment_id:`` filter. Both failures were invisible because the tests
    hand-built a context the resolver could not produce.
    """

    LISTING = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")

    def _resolve(self, manifest: str, **kwargs) -> QARunContext:
        return resolve_qa_run_context(
            data_source="manifest",
            raw_assay="seahub_sci",
            manifest_path=manifest,
            manifest_delimiter="\t",
            manifest_s3_column=6,
            manifest_has_header=True,
            run_label="run1",
            **kwargs,
        )

    def _manifest(self, tmp_path, keys: list[str]) -> str:
        path = tmp_path / "listing.tsv"
        path.write_text(
            "S3_Full_Path\n" + "".join(f"s3://czi-labalpha/{k}\n" for k in keys)
        )
        return str(path)

    def _resolve_simple(self, manifest: str, **kwargs) -> QARunContext:
        return resolve_qa_run_context(
            data_source="manifest",
            raw_assay="seahub_sci",
            manifest_path=manifest,
            manifest_delimiter="\t",
            manifest_s3_column=0,
            manifest_has_header=True,
            run_label="run1",
            **kwargs,
        )

    def test_experiment_id_derived_from_manifest_keys(self):
        assert self._resolve(self.LISTING).order == "REF3"

    def test_explicit_order_is_not_discarded(self):
        assert self._resolve(self.LISTING, order="GENE7").order == "GENE7"

    def test_run_label_still_names_the_outputs(self):
        ctx = self._resolve(self.LISTING)
        assert ctx.output_label == "run1"
        assert ctx.order == "REF3"

    def test_mixed_experiments_raise_once_at_resolve_time(self, tmp_path):
        manifest = self._manifest(
            tmp_path,
            [
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/a.trim.cram",
                "labalpha-seahub-bcp/GENE7/raw/P02/437685/b.trim.cram",
            ],
        )
        with pytest.raises(ValueError, match="mixes SeaHub ExperimentIDs"):
            self._resolve_simple(manifest)

    def test_non_seahub_manifest_keeps_empty_order(self):
        ctx = resolve_qa_run_context(
            data_source="manifest",
            raw_assay="sci_plex",
            manifest_path=self.LISTING,
            manifest_delimiter="\t",
            manifest_s3_column=6,
            manifest_has_header=True,
            run_label="run1",
        )
        assert ctx.order == ""

    def test_gather_reports_no_false_wrong_experiment(self):
        """The regression: one error per object, every object, on a clean upload."""
        ctx = self._resolve(self.LISTING)
        fail_csv = open(SEAHUB_TRIM_FAIL).read()
        s3 = MockS3Client(
            file_contents={
                f"labalpha-seahub-bcp/REF3/raw/P05_1/430479/430479-REF3_P05_1_"
                f"{well}_GEX_hash_oligo-{ug}-{bc}.trim_fail.csv": fail_csv
                for well, ug, bc in (
                    ("A1", "Z0097", "CAGTCAGTTGCAGAT"),
                    ("A2", "Z0105", "CATGGCGCAGTGCTGAT"),
                )
            }
        )
        data = gather_qa_data(ctx, s3)
        assert data.all_raw_files
        assert [e for e in data.gathering_errors if "WRONG EXPERIMENT" in e] == []

    def test_wrong_experiment_is_one_error_not_one_per_object(self, tmp_path):
        """A genuine mixup is still reported -- once, with a count and examples."""
        keys = [
            f"labalpha-seahub-bcp/REF3/raw/P05_1/430479/w{i}.trim.cram"
            for i in range(6)
        ]
        ctx = self._resolve_simple(self._manifest(tmp_path, keys), order="GENE7")
        data = gather_qa_data(ctx, MockS3Client())
        errors = [e for e in data.gathering_errors if "WRONG EXPERIMENT" in e]
        assert len(errors) == 1
        assert "6 object(s) belong to REF3, not GENE7" in errors[0]
        assert "and 4 more" in errors[0]

    def test_unresolvable_experiment_warns_once_and_skips_the_check(self, tmp_path):
        """No ExperimentID means the check cannot run -- say so once, not per object."""
        keys = [f"labalpha-seahub-bcp/REF3/raw/loose{i}.trim.cram" for i in range(4)]
        ctx = self._resolve_simple(self._manifest(tmp_path, keys))
        assert ctx.order == ""
        data = gather_qa_data(ctx, MockS3Client())
        assert [e for e in data.gathering_errors if "WRONG EXPERIMENT" in e] == []
        unknown = [
            w for w in data.gathering_warnings if w.startswith("EXPERIMENT UNKNOWN")
        ]
        assert len(unknown) == 1


class TestSeahubRawFileSizes:
    """Sizes feed the trimmed-vs-vendor size check, so the gatherer records them."""

    STEM = "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT"
    DIR = "labalpha-seahub-bcp/REF3/raw/REF3_P04_1/437120"

    def _keys(self) -> list[str]:
        return [
            f"{self.DIR}/{self.STEM}{suffix}"
            for suffix in (".trim.cram", ".trim.csv", ".trim.stderr", ".trim.stdout")
        ]

    def test_sizes_are_recorded_for_a_seahub_s3_run(self):
        keys = self._keys()
        cram = f"{self.DIR}/{self.STEM}.trim.cram"
        s3 = MockS3Client(keys=keys, sizes={cram: 14_400_000_000})
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
            output_label="REF3",
        )

        data = gather_qa_data(ctx, s3)

        assert data.raw_file_sizes[cram] == 14_400_000_000
        assert set(data.raw_file_sizes) == set(keys)

    def test_sizes_stay_empty_for_a_10x_run(self):
        """Only the SeaHub listing path collects them; 10x must be untouched."""
        s3 = MockS3Client(keys=["testproj/ORD01/G1/raw/x.fastq.gz"])

        data = gather_qa_data(_make_ctx(raw_assay="10x"), s3)

        assert data.raw_file_sizes == {}

    def test_a_missing_size_reads_as_unknown(self):
        s3 = MockS3Client(keys=self._keys())
        ctx = _make_ctx(
            raw_assay="seahub_sci",
            bucket="czi-labalpha",
            provider="labalpha",
            proj="labalpha-seahub-bcp",
            order="REF3",
            listing_prefix="labalpha-seahub-bcp/REF3/",
            output_label="REF3",
        )

        data = gather_qa_data(ctx, s3)

        assert set(data.raw_file_sizes.values()) == {0}
