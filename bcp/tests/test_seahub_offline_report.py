"""Tests for the offline SeaHub report tool (tests/seahub_offline_report.py).

The tool exists so a claim about a real upload can be re-derived by someone who
has the listing. That is only worth anything if the tool still runs, and it is
not exercised by any other test, so these keep it honest against the production
code it drives.
"""

from __future__ import annotations

import json

import pytest

from tests.seahub_offline_report import (
    PARITY_FIELDS,
    load_listing,
    main,
    parity_diff,
    report,
)

BUCKET = "czi-labalpha"
PROJ = "labalpha-seahub-bcp"
RAW = f"{PROJ}/REF3/raw"
# The SOP folder carries the ExperimentID prefix; "P05_1" alone is the truncated
# form, which test_a_defect_shows_up uses deliberately.
WELLS = (
    (
        "REF3_P05_1",
        "430479",
        "430479-REF3_P05_1_A1_GEX_hash_oligo-Z0097-CAGTCAGTTGCAGAT",
    ),
    (
        "REF3_P05_1",
        "430479",
        "430479-REF3_P05_1_A2_GEX_hash_oligo-Z0105-CATGGCGCAGTGCTGAT",
    ),
)
SUFFIXES = (".trim.cram", ".trim.csv", ".trim.stderr", ".trim.stdout", ".trim_fail.csv")


def _keys() -> list[str]:
    return [
        f"{RAW}/{folder}/{wafer}/{stem}{suffix}"
        for folder, wafer, stem in WELLS
        for suffix in SUFFIXES
    ]


def _listing(tmp_path, keys=None, header=True, delimiter="\t"):
    """A listing in the shape the console export produces: URI plus a size."""
    keys = _keys() if keys is None else keys
    suffix = "tsv" if delimiter == "\t" else "csv"
    path = tmp_path / f"REF3.{suffix}"
    rows = []
    if header:
        rows.append(delimiter.join(["File_Name", "Size_Bytes", "S3_Full_Path"]))
    for key in keys:
        rows.append(
            delimiter.join([key.split("/")[-1], "12345", f"s3://{BUCKET}/{key}"])
        )
    path.write_text("\n".join(rows) + "\n")
    return path


class TestLoadListing:
    def test_finds_the_uri_column_by_content(self, tmp_path):
        """Listings come from several tools and agree only on holding URIs."""
        bucket, keys, sizes = load_listing(_listing(tmp_path))

        assert bucket == BUCKET
        assert sorted(keys) == sorted(_keys())
        assert set(sizes.values()) == {12345}

    def test_a_headerless_csv_works_too(self, tmp_path):
        bucket, keys, _ = load_listing(_listing(tmp_path, header=False, delimiter=","))

        assert bucket == BUCKET
        assert len(keys) == len(_keys())

    def test_a_file_with_no_uris_is_a_clear_error(self, tmp_path):
        path = tmp_path / "empty.tsv"
        path.write_text("File_Name\tSize_Bytes\nfoo.cram\t1\n")

        with pytest.raises(SystemExit, match="no s3:// URIs"):
            load_listing(path)


class TestReport:
    def test_it_recovers_the_experiment_id_from_the_keys(self, tmp_path):
        """Nothing names REF3 on the command line; it is read off the listing."""
        assert report(_listing(tmp_path), scratch=tmp_path)["experiment_id"] == "REF3"

    def test_it_reports_the_wells_the_listing_describes(self, tmp_path):
        result = report(_listing(tmp_path), scratch=tmp_path)

        assert result["listing_objects"] == len(_keys())
        assert result["all_raw_files"] == len(_keys())
        assert result["wells"] == len(WELLS)

    def test_a_clean_upload_reports_no_sop_rows(self, tmp_path):
        result = report(_listing(tmp_path), scratch=tmp_path)

        assert result["sop_rows"] == 0
        assert result["well_verdicts"]["COMPLIANT"] == len(WELLS)

    def test_a_defect_shows_up(self, tmp_path):
        """A truncated sublibrary folder, the commonest real defect."""
        keys = [k.replace("/raw/REF3_P05_1/", "/raw/P05_1/") for k in _keys()]
        result = report(_listing(tmp_path, keys=keys), scratch=tmp_path)

        assert "sublibrary_folder_truncated" in result["sop_rules"]

    def test_the_two_modes_agree(self, tmp_path):
        listing = _listing(tmp_path)

        from_s3 = report(listing, mode="s3", scratch=tmp_path)
        from_manifest = report(listing, mode="manifest", scratch=tmp_path)

        assert parity_diff(from_s3, from_manifest) == {}

    def test_parity_excludes_the_field_that_cannot_match(self, tmp_path):
        """discovered_wafers is folder-walk only, so s3 mode alone populates it."""
        listing = _listing(tmp_path)

        from_s3 = report(listing, mode="s3", scratch=tmp_path)
        from_manifest = report(listing, mode="manifest", scratch=tmp_path)

        assert "discovered_wafers" not in PARITY_FIELDS
        assert from_s3["discovered_wafers"] != from_manifest["discovered_wafers"]

    def test_parity_diff_names_a_real_disagreement(self, tmp_path):
        assert parity_diff({**_zeroed(), "wells": 1}, _zeroed()) == {"wells": (1, 0)}


def _zeroed() -> dict:
    return dict.fromkeys(PARITY_FIELDS, 0)


class TestCommandLine:
    def test_json_output_is_parseable_and_free_of_progress_noise(
        self, tmp_path, capsys
    ):
        """--json is meant to be diffed, so stdout must carry nothing else."""
        assert (
            main([str(_listing(tmp_path)), "--json", "--scratch", str(tmp_path)]) == 0
        )

        payload = json.loads(capsys.readouterr().out)
        assert payload["REF3/s3"]["experiment_id"] == "REF3"

    def test_both_modes_prints_a_parity_verdict(self, tmp_path, capsys):
        main([str(_listing(tmp_path)), "--mode", "both", "--scratch", str(tmp_path)])

        assert "s3 and manifest agree" in capsys.readouterr().out
