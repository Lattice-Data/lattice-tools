"""
Tests for SeaHub lab raw QA support (seahub_sci assay).
"""

from __future__ import annotations

import os

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
    grab_seahub_trim_fail_csv,
    normalize_raw_assay,
    parse_raw_filename,
    parse_seahub_raw_path,
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
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
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
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        hash_pct, gex_pct = stats["430479"]["trimmer_fail"]
        assert hash_pct == 100 * 157722381 / 260527531
        assert gex_pct == 100 * 100735571 / 158233602

    def test_grab_seahub_trim_fail_csv_appends(self):
        """Repeated wells accumulate into the same distribution."""
        stats: dict = {}
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
        grab_seahub_trim_fail_csv(stats, "430479", SEAHUB_TRIM_FAIL)
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
        grab_seahub_trim_fail_csv(stats, "exp1", path, warnings=warnings)
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
        grab_seahub_trim_fail_csv(stats, "exp1", path)
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
    def _seahub_paginated_pages(self) -> dict:
        o = "labalpha-seahub-bcp/REF3/"
        raw = f"{o}raw/"
        sublib = f"{raw}P05_1/"
        wafer = f"{sublib}430479/"
        return {
            (o, "/"): [
                {"CommonPrefixes": [{"Prefix": raw}]},
            ],
            (raw, "/"): [
                {"CommonPrefixes": [{"Prefix": sublib}]},
            ],
            (sublib, "/"): [
                {"CommonPrefixes": [{"Prefix": wafer}]},
            ],
            (wafer, ""): [
                {
                    "Contents": [
                        {"Key": SEAHUB_KEY_CRAM},
                        {"Key": SEAHUB_KEY_FAIL},
                        {"Key": SEAHUB_KEY_TRIM},
                        {"Key": f"{SEAHUB_BASE}.trim.stderr"},
                        {"Key": f"{SEAHUB_BASE}.trim.stdout"},
                    ]
                }
            ],
        }

    def test_gather_from_s3_seahub_layout(self):
        pages = self._seahub_paginated_pages()
        paginated = dict(pages)
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            SEAHUB_KEY_TRIM: open(SEAHUB_TRIM).read(),
        }
        s3 = MockS3Client(paginated_pages=paginated, file_contents=file_contents)
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

    def test_missing_metadata_json_warns(self):
        """Each .trim.cram lacking a .trim.cram-metadata.json sidecar warns."""
        manifest = os.path.join(QA_FIXTURES_DIR, "seahub_s3_listing.tsv")
        ctx = self._manifest_ctx(manifest)
        file_contents = {
            k: open(SEAHUB_TRIM_FAIL).read()
            for k in (
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT.trim_fail.csv",
                "labalpha-seahub-bcp/REF3/raw/P05_1/430479/"
                "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT.trim_fail.csv",
            )
        }
        s3 = MockS3Client(file_contents=file_contents)
        data = gather_qa_data(ctx, s3)
        missing = [
            w for w in data.gathering_warnings if w.startswith("METADATA MISSING")
        ]
        assert len(missing) == 2
        assert all(".trim.cram-metadata.json" in w for w in missing)

    def test_present_metadata_json_no_warn(self):
        """A .trim.cram with its .trim.cram-metadata.json sidecar does not warn."""
        pages = dict(self._seahub_paginated_pages())
        wafer_key = ("labalpha-seahub-bcp/REF3/raw/P05_1/430479/", "")
        meta_key = f"{SEAHUB_BASE}.trim.cram-metadata.json"
        pages[wafer_key] = [
            {"Contents": pages[wafer_key][0]["Contents"] + [{"Key": meta_key}]}
        ]
        file_contents = {
            SEAHUB_KEY_FAIL: open(SEAHUB_TRIM_FAIL).read(),
            meta_key: '{"read_count": 1000, "filename": "%s"}'
            % SEAHUB_KEY_CRAM.split("/")[-1],
        }
        s3 = MockS3Client(paginated_pages=pages, file_contents=file_contents)
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
