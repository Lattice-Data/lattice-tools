"""Tests for file_extract.geo_md5 and the streaming md5 in s3_utils."""

from __future__ import annotations

import csv
import hashlib
import io
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pytest
from botocore.exceptions import IncompleteReadError

from file_extract.geo_md5 import (
    GeoFile,
    GeoMd5Error,
    character_warning,
    default_geo_md5_output_name,
    extension_warning,
    extract_geo_md5,
    geo_file_name,
    name_warnings,
    plan_geo_files,
    read_table,
)
from file_extract.s3_utils import parse_s3_object_uri, stream_md5
from tests.file_extract_helpers import MockS3Client

BUCKET = "example-bucket"
CR_BASENAME = "sample_filtered_feature_bc_matrix.h5"
H5_A = f"proj/libA/processed/outs/per_sample_outs/s1/{CR_BASENAME}"
H5_B = f"proj/libB/processed/outs/per_sample_outs/s2/{CR_BASENAME}"
H5AD_A = "proj/order/processed/run_date/samples/sampleA.QSR-1_anndata.h5ad"
H5AD_B = "proj/order/processed/run_date/samples/sampleB.QSR-2_anndata.h5ad"
BODY_A = b"a" * 1000 + b"tail"
BODY_B = bytes(range(256)) * 3


def _uri(key: str) -> str:
    return f"s3://{BUCKET}/{key}"


def _md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


def _file(row: int, key: str, name: str | None = None) -> GeoFile:
    return GeoFile(row=row, bucket=BUCKET, key=key, name=name or key.rsplit("/", 1)[-1])


def _write_table(
    path: Path,
    header: list[str],
    rows: list[list[str]],
    *,
    delimiter: str = "\t",
    encoding: str = "utf-8",
) -> Path:
    with path.open("w", newline="", encoding=encoding) as fh:
        writer = csv.writer(fh, delimiter=delimiter)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _read_csv(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.reader(fh))


# --- parse_s3_object_uri -----------------------------------------------------


def test_parse_s3_object_uri() -> None:
    assert parse_s3_object_uri(_uri(H5_A)) == (BUCKET, H5_A)
    assert parse_s3_object_uri(f"  {_uri(H5_A)}\n") == (BUCKET, H5_A)


def test_parse_s3_object_uri_keeps_query_and_fragment_characters() -> None:
    key = "proj/odd?name#1.h5"
    assert parse_s3_object_uri(_uri(key)) == (BUCKET, key)


@pytest.mark.parametrize(
    "uri, message",
    [
        ("https://example.com/a.h5", "Invalid S3 URI"),
        ("s3://", "Invalid S3 URI"),
        (f"s3://{BUCKET}", "does not name an object"),
        (f"s3://{BUCKET}/", "does not name an object"),
        (f"s3://{BUCKET}/proj/", "names a prefix"),
        (f"s3://{BUCKET}//proj/a.h5", "starts with '/'"),
        (f"s3://{BUCKET}/proj/a\nb.h5", "line break"),
    ],
)
def test_parse_s3_object_uri_rejects(uri: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_s3_object_uri(uri)


# --- geo_file_name -----------------------------------------------------------


@pytest.mark.parametrize(
    "cell, expected",
    [
        ("libA-s1.h5", "libA-s1.h5"),
        ("  libA-s1.h5 ", "libA-s1.h5"),
        ("example-lab:libA-s1.h5", "libA-s1.h5"),
        ('["example-lab:libA-s1.h5"]', "libA-s1.h5"),
        ('[ "example-lab:libA-s1.h5" ]', "libA-s1.h5"),
        # Everything after the first colon, as db2-flattener strips it.
        ("example-lab:odd:name.h5", "odd:name.h5"),
    ],
)
def test_geo_file_name(cell: str, expected: str) -> None:
    assert geo_file_name(cell) == expected


@pytest.mark.parametrize(
    "cell, message",
    [
        ('["example-lab:a.h5", "example-lab:b.h5"]', "lists 2 aliases"),
        ("[]", "lists 0 aliases"),
        ("[not json", "not a JSON alias array"),
        ("[1]", "not a list of alias strings"),
        ("", "empty"),
        ("   ", "empty"),
        ("example-lab:", "empty"),
        (".", "not a file name"),
        ('["example-lab:.."]', "not a file name"),
    ],
)
def test_geo_file_name_rejects(cell: str, message: str) -> None:
    with pytest.raises(GeoMd5Error, match=message):
        geo_file_name(cell)


# --- plan_geo_files ----------------------------------------------------------


def test_plan_defaults_to_basename_in_row_order() -> None:
    rows = [{"s3_uri": _uri(H5AD_B)}, {"s3_uri": _uri(H5AD_A)}]
    files = plan_geo_files(rows)
    assert [f.name for f in files] == [
        "sampleB.QSR-2_anndata.h5ad",
        "sampleA.QSR-1_anndata.h5ad",
    ]
    assert [f.row for f in files] == [1, 2]
    assert files[0].s3_uri == _uri(H5AD_B)
    assert files[0].basename == "sampleB.QSR-2_anndata.h5ad"


def test_plan_name_column_strips_namespace() -> None:
    rows = [
        {"aliases": '["example-lab:libA-s1.h5"]', "s3_uri": _uri(H5_A)},
        {"aliases": "example-lab:libB-s2.h5", "s3_uri": _uri(H5_B)},
    ]
    files = plan_geo_files(rows, name_column="aliases")
    assert [f.name for f in files] == ["libA-s1.h5", "libB-s2.h5"]


def test_plan_honours_uri_column() -> None:
    rows = [{"uri": _uri(H5AD_A), "s3_uri": "junk"}]
    files = plan_geo_files(rows, uri_column="uri")
    assert files[0].key == H5AD_A


def test_plan_missing_columns_name_the_columns_present() -> None:
    rows = [{"#response": "", "uri": _uri(H5AD_A), "": "extra"}]
    with pytest.raises(GeoMd5Error, match=r"'s3_uri' not in table.*'#response', 'uri'"):
        plan_geo_files(rows)
    with pytest.raises(GeoMd5Error, match=r"'aliases' not in table"):
        plan_geo_files(rows, uri_column="uri", name_column="aliases")


def test_plan_no_rows() -> None:
    with pytest.raises(GeoMd5Error, match="no data rows"):
        plan_geo_files([])


def test_plan_empty_uri_is_an_error_not_a_skip() -> None:
    rows = [{"s3_uri": _uri(H5AD_A)}, {"s3_uri": "  "}]
    with pytest.raises(GeoMd5Error, match="row 2: empty 's3_uri'"):
        plan_geo_files(rows)


def test_plan_invalid_uri_names_the_row() -> None:
    rows = [{"s3_uri": _uri(H5AD_A)}, {"s3_uri": f"s3://{BUCKET}/proj/"}]
    with pytest.raises(GeoMd5Error, match="row 2: S3 URI names a prefix"):
        plan_geo_files(rows)


def test_plan_bad_name_cell_names_the_row() -> None:
    rows = [{"aliases": '["a", "b"]', "s3_uri": _uri(H5_A)}]
    with pytest.raises(GeoMd5Error, match="row 1: name cell .* lists 2 aliases"):
        plan_geo_files(rows, name_column="aliases")


def test_plan_duplicate_basenames_suggest_name_column() -> None:
    rows = [{"s3_uri": _uri(H5_A)}, {"s3_uri": _uri(H5_B)}]
    with pytest.raises(GeoMd5Error) as exc_info:
        plan_geo_files(rows)
    message = str(exc_info.value)
    assert f"'{CR_BASENAME}'" in message
    assert "GEO rejects duplicate file names" in message
    assert "--name-column" in message


def test_plan_duplicate_names_from_name_column_do_not_suggest_it() -> None:
    rows = [
        {"aliases": "example-lab:same.h5", "s3_uri": _uri(H5_A)},
        {"aliases": "example-lab:same.h5", "s3_uri": _uri(H5_B)},
    ]
    with pytest.raises(GeoMd5Error) as exc_info:
        plan_geo_files(rows, name_column="aliases")
    message = str(exc_info.value)
    assert "'same.h5'" in message
    assert "--name-column" not in message


def test_plan_same_object_twice() -> None:
    rows = [
        {"aliases": "example-lab:a.h5", "s3_uri": _uri(H5_A)},
        {"aliases": "example-lab:b.h5", "s3_uri": _uri(H5_A)},
    ]
    with pytest.raises(GeoMd5Error, match="listed more than once"):
        plan_geo_files(rows, name_column="aliases")


@pytest.mark.parametrize("name", ["dir/file.h5", "dir\\file.h5"])
def test_plan_rejects_path_separators(name: str) -> None:
    rows = [{"aliases": f"example-lab:{name}", "s3_uri": _uri(H5_A)}]
    with pytest.raises(GeoMd5Error, match="path separator"):
        plan_geo_files(rows, name_column="aliases")


# --- extension_warning -------------------------------------------------------


def test_extension_warning_none_when_extensions_match() -> None:
    assert extension_warning(_file(1, H5_A, "libA-s1.h5")) is None
    assert extension_warning(_file(1, H5AD_A)) is None
    mtx = "proj/scaleplex/s.QSR-1-SCALEPLEX.filtered.matrix/matrix.mtx.gz"
    assert extension_warning(_file(1, mtx, "s.QSR-1.matrix.mtx.gz")) is None
    assert extension_warning(_file(1, "proj/noext", "renamed")) is None


def test_extension_warning_alias_without_extension() -> None:
    warning = extension_warning(_file(3, H5_A, "libA-s1"))
    assert warning == (
        "row 3: GEO file name 'libA-s1' ends in no extension but the S3 object "
        "ends in .h5"
    )


def test_extension_warning_different_extension() -> None:
    warning = extension_warning(_file(2, H5AD_A, "sampleA.h5"))
    assert warning is not None
    assert "ends in .h5 but the S3 object ends in .h5ad" in warning
    # A dotted alias has an "extension" too, which is the point of the warning.
    assert extension_warning(_file(1, H5_A, "lib-003.sampleA")) is not None


def test_extension_warning_sees_past_a_compression_suffix() -> None:
    # .tsv.gz and .mtx.gz both end in .gz; the comparison must pair the
    # compression suffix with the one before it or it cannot tell them apart.
    mtx = "proj/scaleplex/s.QSR-1.filtered.matrix/matrix.mtx.gz"
    warning = extension_warning(_file(1, mtx, "s.QSR-1.features.tsv.gz"))
    assert warning is not None
    assert "ends in .tsv.gz but the S3 object ends in .mtx.gz" in warning
    warning = extension_warning(_file(2, mtx, "s.QSR-1.matrix.mtx"))
    assert warning is not None
    assert "ends in .mtx but the S3 object ends in .mtx.gz" in warning


def test_character_warning() -> None:
    assert character_warning(_file(1, H5_A, "lib-003_s1.v2.h5")) is None
    warning = character_warning(_file(4, H5_A, "lib 003 s1.h5"))
    assert warning is not None
    assert warning.startswith("row 4: GEO file name 'lib 003 s1.h5' has characters")
    assert character_warning(_file(1, H5_A, "a,b.h5")) is not None
    assert character_warning(_file(1, H5_A, "odd:name.h5")) is not None


def test_name_warnings_in_table_order() -> None:
    files = [
        _file(1, H5_A, "fine.h5"),
        _file(2, H5_A, "no ext"),
        _file(3, H5_A, "x.h5"),
    ]
    warnings = name_warnings(files)
    assert [w.split(":")[0] for w in warnings] == ["row 2", "row 2"]
    assert "no extension" in warnings[0]
    assert "characters" in warnings[1]


# --- read_table / output name -----------------------------------------------


def test_read_table_tsv_and_csv(tmp_path: Path) -> None:
    tsv = _write_table(tmp_path / "t.tsv", ["s3_uri", "x"], [[_uri(H5AD_A), "a,b"]])
    assert read_table(tsv) == [{"s3_uri": _uri(H5AD_A), "x": "a,b"}]
    csv_path = _write_table(
        tmp_path / "t.csv", ["s3_uri", "x"], [[_uri(H5AD_A), "a\tb"]], delimiter=","
    )
    assert read_table(csv_path) == [{"s3_uri": _uri(H5AD_A), "x": "a\tb"}]


def test_read_table_strips_excel_bom(tmp_path: Path) -> None:
    table = _write_table(
        tmp_path / "t.tsv", ["s3_uri"], [[_uri(H5AD_A)]], encoding="utf-8-sig"
    )
    assert list(read_table(table)[0]) == ["s3_uri"]


def test_read_table_rejects_non_utf8(tmp_path: Path) -> None:
    table = tmp_path / "utf16.tsv"
    table.write_bytes("s3_uri\n".encode("utf-16"))
    with pytest.raises(GeoMd5Error, match="not UTF-8 text"):
        read_table(table)


def test_read_table_empty_file(tmp_path: Path) -> None:
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    with pytest.raises(GeoMd5Error, match="expected a header row"):
        read_table(empty)


def test_default_geo_md5_output_name() -> None:
    assert default_geo_md5_output_name("x/order_h5_info.tsv") == (
        "order_h5_info_geo_md5.csv"
    )
    assert default_geo_md5_output_name("RawMatrixFile.csv") == (
        "RawMatrixFile_geo_md5.csv"
    )


# --- stream_md5 --------------------------------------------------------------


def test_stream_md5_matches_hashlib_and_counts_bytes() -> None:
    client = MockS3Client(object_bodies={H5_A: BODY_A})
    assert stream_md5(client, BUCKET, H5_A) == (_md5(BODY_A), len(BODY_A))
    assert client.get_object_calls == [
        {"Bucket": BUCKET, "Key": H5_A, "ChecksumMode": "ENABLED"}
    ]


def test_stream_md5_chunked_read_is_the_same_digest() -> None:
    client = MockS3Client(object_bodies={H5_B: BODY_B})
    assert stream_md5(client, BUCKET, H5_B, chunk_bytes=7) == (
        _md5(BODY_B),
        len(BODY_B),
    )


def test_stream_md5_empty_object() -> None:
    client = MockS3Client(object_bodies={H5_A: b""})
    assert stream_md5(client, BUCKET, H5_A) == (_md5(b""), 0)


def test_stream_md5_short_read_is_an_error() -> None:
    class TruncatingClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {"Body": io.BytesIO(BODY_A[:-3]), "ContentLength": len(BODY_A)}

    client = TruncatingClient(object_bodies={H5_A: BODY_A})
    # botocore's own type, so retry_with_backoff treats it as transient too.
    with pytest.raises(
        IncompleteReadError, match=r"1001 read, but total bytes expected is 1004"
    ):
        stream_md5(client, BUCKET, H5_A)


def test_stream_md5_over_long_read_is_an_error() -> None:
    class PaddingClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {"Body": io.BytesIO(BODY_A + b"xyz"), "ContentLength": len(BODY_A)}

    with pytest.raises(
        IncompleteReadError, match=r"1007 read, but total bytes expected is 1004"
    ):
        stream_md5(PaddingClient(), BUCKET, H5_A)


def test_stream_md5_without_content_length_trusts_the_stream() -> None:
    class NoLengthClient(MockS3Client):
        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            return {"Body": io.BytesIO(BODY_A)}

    assert stream_md5(NoLengthClient(), BUCKET, H5_A) == (_md5(BODY_A), len(BODY_A))


def test_stream_md5_with_a_real_boto3_client(monkeypatch: pytest.MonkeyPatch) -> None:
    """The mock above cannot vouch for ChecksumMode or StreamingBody; moto can."""
    pytest.importorskip("moto")
    import boto3
    from moto import mock_aws

    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-east-1")
    with mock_aws():
        client = boto3.client("s3", region_name="us-east-1")
        client.create_bucket(Bucket=BUCKET)
        # Stored with a checksum botocore computes natively, so the response
        # carries x-amz-checksum-crc32 and the validating body is exercised.
        client.put_object(
            Bucket=BUCKET, Key=H5_B, Body=BODY_B, ChecksumAlgorithm="CRC32"
        )
        assert stream_md5(client, BUCKET, H5_B, chunk_bytes=100) == (
            _md5(BODY_B),
            len(BODY_B),
        )


# --- extract_geo_md5 ---------------------------------------------------------


def test_extract_writes_geo_header_and_rows_in_table_order(tmp_path: Path) -> None:
    client = MockS3Client(object_bodies={H5_A: BODY_A, H5_B: BODY_B})
    files = [_file(1, H5_B, "libB-s2.h5"), _file(2, H5_A, "libA-s1.h5")]
    out = tmp_path / "geo_md5.csv"

    summary = extract_geo_md5(client, files, out, show_progress=False)

    assert _read_csv(out) == [
        ["file name", "file checksum"],
        ["libB-s2.h5", _md5(BODY_B)],
        ["libA-s1.h5", _md5(BODY_A)],
    ]
    assert summary.total == 2
    assert summary.md5_ok == 2
    assert summary.bytes_total == len(BODY_A) + len(BODY_B)
    assert summary.failures == []
    assert not summary.has_failures


def test_extract_failed_object_leaves_cell_empty_and_is_listed(tmp_path: Path) -> None:
    client = MockS3Client(object_bodies={H5AD_A: BODY_A})
    files = [_file(1, H5AD_A), _file(2, H5AD_B)]
    out = tmp_path / "geo_md5.csv"

    summary = extract_geo_md5(client, files, out, show_progress=False)

    assert _read_csv(out) == [
        ["file name", "file checksum"],
        ["sampleA.QSR-1_anndata.h5ad", _md5(BODY_A)],
        ["sampleB.QSR-2_anndata.h5ad", ""],
    ]
    assert summary.md5_ok == 1
    assert summary.bytes_total == len(BODY_A)
    assert len(summary.failures) == 1
    uri, err, detail = summary.failures[0]
    assert uri == _uri(H5AD_B)
    assert "NoSuchKey" in err
    assert detail == ""
    assert summary.has_failures


def test_extract_quotes_a_name_holding_a_comma(tmp_path: Path) -> None:
    client = MockS3Client(object_bodies={H5_A: BODY_A})
    out = tmp_path / "geo_md5.csv"
    extract_geo_md5(client, [_file(1, H5_A, "a,b.h5")], out, show_progress=False)
    assert _read_csv(out)[1] == ["a,b.h5", _md5(BODY_A)]
    assert out.read_text(encoding="utf-8").splitlines()[1].startswith('"a,b.h5",')


def test_extract_pool_size_is_min_of_workers_and_file_count(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import file_extract.geo_md5 as geo_md5_module

    seen: list[int | None] = []

    class RecordingExecutor(ThreadPoolExecutor):
        def __init__(self, max_workers: int | None = None, **kwargs: Any) -> None:
            seen.append(max_workers)
            super().__init__(max_workers=max_workers, **kwargs)

    monkeypatch.setattr(geo_md5_module, "ThreadPoolExecutor", RecordingExecutor)
    client = MockS3Client(object_bodies={H5_A: BODY_A, H5_B: BODY_B, H5AD_A: b"x"})
    one = [_file(1, H5_A)]
    three = [_file(1, H5_A), _file(2, H5_B), _file(3, H5AD_A)]
    out = tmp_path / "geo_md5.csv"

    assert (
        extract_geo_md5(client, one, out, workers=64, show_progress=False).md5_ok == 1
    )
    assert (
        extract_geo_md5(client, three, out, workers=2, show_progress=False).md5_ok == 3
    )
    assert extract_geo_md5(client, three, out, show_progress=False).md5_ok == 3
    assert seen == [1, 2, 3]


def test_extract_retries_a_transient_get(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from botocore.exceptions import ClientError

    class FlakyClient(MockS3Client):
        """Fails the first GET of every key with a throttle, then serves it."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.throttled: set[str] = set()

        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            if Key not in self.throttled:
                self.throttled.add(Key)
                raise ClientError(
                    {"Error": {"Code": "SlowDown", "Message": "slow"}}, "GetObject"
                )
            return super().get_object(Bucket=Bucket, Key=Key, **kwargs)

    monkeypatch.setattr("file_extract.retry.time.sleep", lambda _seconds: None)
    out = tmp_path / "geo_md5.csv"

    summary = extract_geo_md5(
        FlakyClient(object_bodies={H5_A: BODY_A}),
        [_file(1, H5_A)],
        out,
        retries=2,
        show_progress=False,
    )
    assert summary.md5_ok == 1 and not summary.failures
    assert _read_csv(out)[1] == [CR_BASENAME, _md5(BODY_A)]

    summary = extract_geo_md5(
        FlakyClient(object_bodies={H5_A: BODY_A}),
        [_file(1, H5_A)],
        out,
        retries=1,
        show_progress=False,
    )
    assert summary.md5_ok == 0
    assert len(summary.failures) == 1 and "SlowDown" in summary.failures[0][1]
    assert _read_csv(out)[1] == [CR_BASENAME, ""]


# --- review follow-ups --------------------------------------------------------


def test_case_warnings_flag_names_that_differ_only_in_case() -> None:
    rows = [
        {"s3_uri": _uri(H5_A), "aliases": "example-lab:Sample1.h5"},
        {"s3_uri": _uri(H5AD_A), "aliases": "example-lab:other.h5ad"},
        {"s3_uri": _uri(H5_B), "aliases": "example-lab:sample1.h5"},
    ]
    # Not a duplicate for GEO, so the plan goes through; the warning is advisory.
    files = plan_geo_files(rows, name_column="aliases")
    assert len(files) == 3
    assert name_warnings(files) == [
        "rows 1, 3: GEO file names 'Sample1.h5', 'sample1.h5' differ only in "
        "case; a case-insensitive filesystem merges them during the renaming step"
    ]


def test_read_table_skips_rows_of_only_delimiters(tmp_path: Path) -> None:
    table = tmp_path / "t.tsv"
    table.write_text(
        f"s3_uri\tnote\n{_uri(H5AD_A)}\tx\n\t\n{_uri(H5AD_B)}\t\n\t\t\n",
        encoding="utf-8",
    )
    assert [row["s3_uri"] for row in read_table(table)] == [_uri(H5AD_A), _uri(H5AD_B)]

    # A row with any content is still a row, and a missing URI still an error.
    table.write_text(f"s3_uri\tnote\n{_uri(H5AD_A)}\tx\n\tkept\n", encoding="utf-8")
    with pytest.raises(GeoMd5Error, match="row 2: empty 's3_uri'"):
        plan_geo_files(read_table(table))


def test_read_table_rejects_repeated_column_names(tmp_path: Path) -> None:
    # Trailing empty header cells are tolerated; a repeated real name is not,
    # because csv.DictReader would silently keep only the last one.
    table = _write_table(
        tmp_path / "t.tsv",
        ["s3_uri", "x", "s3_uri", "", ""],
        [[_uri(H5AD_A), "a", _uri(H5AD_B), "", ""]],
    )
    with pytest.raises(GeoMd5Error, match="repeats column names: 's3_uri'"):
        read_table(table)


def test_extract_fails_before_any_download_when_output_cannot_be_written(
    tmp_path: Path,
) -> None:
    client = MockS3Client(object_bodies={H5_A: BODY_A})
    with pytest.raises(OSError):
        extract_geo_md5(
            client,
            [_file(1, H5_A)],
            tmp_path / "missing" / "out.csv",
            show_progress=False,
        )
    assert client.get_object_calls == []


def test_extract_interrupt_cancels_queued_downloads(tmp_path: Path) -> None:
    class InterruptingClient(MockS3Client):
        """The first GET is interrupted; later ones are slow enough to observe."""

        def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
            first = not self.get_object_calls
            result = super().get_object(Bucket=Bucket, Key=Key, **kwargs)
            if first:
                raise KeyboardInterrupt
            time.sleep(0.2)
            return result

    keys = [f"proj/lib{i}/matrix.h5" for i in range(6)]
    client = InterruptingClient(object_bodies={key: BODY_A for key in keys})
    files = [_file(i + 1, key, f"s{i}.h5") for i, key in enumerate(keys)]
    out = tmp_path / "geo_md5.csv"
    with pytest.raises(KeyboardInterrupt):
        extract_geo_md5(client, files, out, workers=1, show_progress=False)
    # The interrupted stream, plus at most the one the worker had already
    # picked up; the rest were cancelled before they started. Without the
    # cancel, all six run to completion first.
    assert len(client.get_object_calls) <= 2
    assert _read_csv(out) == [["file name", "file checksum"]]
