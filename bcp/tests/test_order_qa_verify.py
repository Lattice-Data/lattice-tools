"""Tests for order_qa.verify and order_qa.integrity.

No network: a small stub client stands in for S3, the same way test_qa_gather
does, so the listing, quiescence, capability and version paths are all exercised
including their denial branches.
"""

from datetime import datetime, timedelta, timezone

from botocore.exceptions import ClientError

from order_qa.integrity import (
    collect_version_state,
    collect_versions_via_head,
    probe_capabilities,
)
from order_qa.verify import (
    Listing,
    S3Object,
    assess_quiescence,
    compare_to_manifest,
    list_order_objects,
    read_manifest_keys,
)

NOW = datetime(2026, 8, 18, 12, 0, tzinfo=timezone.utc)


def obj(key, size=1000, minutes_old=60, etag='"abc123"'):
    return S3Object(
        key=key,
        size_bytes=size,
        last_modified=NOW - timedelta(minutes=minutes_old),
        etag=etag,
    )


def denial(operation):
    return ClientError(
        {"Error": {"Code": "AccessDenied", "Message": f"{operation} denied"}},
        operation,
    )


class StubS3:
    """Enough of the S3 API for these paths, with configurable denials."""

    def __init__(self, pages=None, deny=(), parts_count=12):
        self._pages = pages or []
        self._deny = set(deny)
        self._parts_count = parts_count
        self.head_calls = []

    def get_paginator(self, operation):
        if operation in self._deny:
            raise denial(operation)
        stub = self

        class _Paginator:
            def paginate(self, **kwargs):
                if operation in stub._deny:
                    raise denial(operation)
                return iter(stub._pages)

        return _Paginator()

    def head_object(self, Bucket, Key, PartNumber=None):
        if "head_object" in self._deny:
            raise denial("HeadObject")
        self.head_calls.append(Key)
        if PartNumber is not None:
            if "head_object_part_number" in self._deny:
                raise denial("HeadObject")
            return {"PartsCount": self._parts_count, "ContentLength": 8 << 20}
        return {
            "ContentLength": 1000,
            "VersionId": f"ver-{Key}",
            "ServerSideEncryption": "aws:kms",
        }

    def get_object_attributes(self, **kwargs):
        if "get_object_attributes" in self._deny:
            raise denial("GetObjectAttributes")
        return {"ObjectSize": 1000, "ObjectParts": {"TotalPartsCount": 12}}

    def list_object_versions(self, **kwargs):
        if "list_object_versions" in self._deny:
            raise denial("ListObjectVersions")
        return {"Versions": [], "DeleteMarkers": []}


class TestListing:
    def test_folder_markers_are_not_deliverables(self):
        """A zero-byte 'dir/' key is not a file the vendor uploaded."""
        listing = Listing("b", "p/", [obj("p/raw/", size=0), obj("p/raw/a.fastq.gz")])
        assert len(listing.objects) == 2
        assert len(listing.files) == 1
        assert listing.total_bytes == 1000

    def test_paginated_listing_collects_every_page(self):
        pages = [
            {
                "Contents": [
                    {"Key": "p/a", "Size": 1, "LastModified": NOW, "ETag": '"x"'}
                ]
            },
            {
                "Contents": [
                    {"Key": "p/b", "Size": 2, "LastModified": NOW, "ETag": '"y"'}
                ]
            },
        ]
        listing = list_order_objects(StubS3(pages=pages), "b", "p/")
        assert listing.complete
        assert listing.keys == {"p/a", "p/b"}

    def test_failed_listing_is_reported_not_raised(self):
        """A partial listing looks exactly like a partial upload.

        So it must never come back as a clean, complete listing.
        """
        listing = list_order_objects(StubS3(deny={"list_objects_v2"}), "b", "p/")
        assert listing.complete is False
        assert "AccessDenied" in listing.error

    def test_naive_timestamps_do_not_crash_quiescence(self):
        """boto3 returns aware datetimes; fixtures may not."""
        pages = [
            {
                "Contents": [
                    {
                        "Key": "p/a",
                        "Size": 1,
                        "LastModified": datetime(2026, 8, 18, 10, 0),
                        "ETag": '"x"',
                    }
                ]
            }
        ]
        listing = list_order_objects(StubS3(pages=pages), "b", "p/")
        assert assess_quiescence(listing, threshold_minutes=15, now=NOW).quiet


class TestMultipartEtag:
    def test_multipart_etag_is_not_md5(self):
        """FASTQs arrive multipart, where the ETag is a digest of part digests."""
        multipart = obj("p/a.fastq.gz", etag='"d41d8cd98f00b204e9800998ecf8427e-12"')
        assert multipart.etag_is_md5 is False
        assert multipart.multipart_parts == 12

    def test_single_part_etag_is_md5(self):
        single = obj("p/a.csv", etag='"d41d8cd98f00b204e9800998ecf8427e"')
        assert single.etag_is_md5 is True
        assert single.multipart_parts is None


class TestQuiescence:
    def test_quiet_when_newest_is_older_than_threshold(self):
        q = assess_quiescence(
            Listing("b", "p/", [obj("p/a", minutes_old=60)]),
            threshold_minutes=15,
            now=NOW,
        )
        assert q.quiet and q.checked

    def test_not_quiet_while_objects_are_still_arriving(self):
        q = assess_quiescence(
            Listing("b", "p/", [obj("p/a", minutes_old=2)]),
            threshold_minutes=15,
            now=NOW,
        )
        assert q.quiet is False
        assert "NOT quiet" in q.summary

    def test_newest_object_decides_not_the_oldest(self):
        listing = Listing(
            "b", "p/", [obj("p/a", minutes_old=600), obj("p/b", minutes_old=1)]
        )
        assert assess_quiescence(listing, threshold_minutes=15, now=NOW).quiet is False

    def test_force_records_that_the_gate_was_skipped(self):
        q = assess_quiescence(
            Listing("b", "p/", [obj("p/a", minutes_old=1)]),
            threshold_minutes=15,
            force=True,
            now=NOW,
        )
        assert q.quiet is True
        assert q.forced is True
        assert "forced" in q.summary

    def test_empty_prefix_is_unchecked_and_not_quiet(self):
        """Nothing there is a verification failure, not a settled upload."""
        q = assess_quiescence(Listing("b", "p/", []), now=NOW)
        assert q.checked is False
        assert q.quiet is False


class TestManifest:
    def test_no_manifest_is_reported_unchecked_never_passed(self):
        comparison = compare_to_manifest(Listing("b", "p/", [obj("p/a")]), None)
        assert comparison.checked is False
        assert comparison.ok is False
        assert "naming convention" in comparison.reason

    def test_missing_and_extra_keys(self, tmp_path):
        manifest = tmp_path / "m.tsv"
        manifest.write_text("s3://b/p/a\ts3://b/p/gone\n")
        listing = Listing("b", "p/", [obj("p/a"), obj("p/surprise")])
        comparison = compare_to_manifest(listing, manifest)
        assert comparison.checked
        assert comparison.missing == ["p/gone"]
        assert comparison.extra == ["p/surprise"]
        assert comparison.ok is False

    def test_size_mismatch_is_a_failure(self, tmp_path):
        manifest = tmp_path / "m.csv"
        manifest.write_text("s3://b/p/a,2048\n")
        comparison = compare_to_manifest(
            Listing("b", "p/", [obj("p/a", size=1000)]), manifest
        )
        assert comparison.size_mismatches == [("p/a", 2048, 1000)]
        assert comparison.ok is False

    def test_uris_found_regardless_of_column(self, tmp_path):
        """Manifest column order varies; a key in the wrong column is still a key."""
        manifest = tmp_path / "m.tsv"
        manifest.write_text("sample1\tGEX\ts3://b/p/a\t1000\n")
        keys = read_manifest_keys(manifest)
        assert keys == {"p/a": 1000}

    def test_manifest_with_no_uris_is_unchecked(self, tmp_path):
        manifest = tmp_path / "m.csv"
        manifest.write_text("just,some,headers\n")
        comparison = compare_to_manifest(Listing("b", "p/", [obj("p/a")]), manifest)
        assert comparison.checked is False
        assert "no s3:// URIs" in comparison.reason


class TestCapabilityProbe:
    def test_denials_are_labelled_denied_not_broken(self):
        report = probe_capabilities(
            StubS3(deny={"get_object_attributes", "list_object_versions"}),
            "czi-psomagen",
            probe_key="p/a.fastq.gz",
            prefix="p/",
        )
        assert report.get("head_object").status == "available"
        assert report.get("get_object_attributes").status == "denied"
        assert report.has("get_object_attributes") is False
        assert report.has("head_object") is True

    def test_no_probe_key_gives_a_reason_not_a_pass(self):
        report = probe_capabilities(StubS3(), "b", probe_key="", prefix="p/")
        assert report.capabilities == []
        assert "no object available" in report.reason


class TestVersionState:
    def test_falls_back_to_head_when_listing_versions_is_denied(self):
        stub = StubS3(deny={"list_object_versions"})
        capabilities = probe_capabilities(stub, "b", probe_key="p/a", prefix="p/")
        state = collect_version_state(
            stub, "b", "p/", ["p/a", "p/b"], capabilities=capabilities, head_limit=10
        )
        assert state.method == "head_object"
        assert state.version_by_key == {"p/a": "ver-p/a", "p/b": "ver-p/b"}

    def test_bounded_sweep_reports_itself_as_partial(self):
        """A cap reported as a complete check is how a 4000-file order gets
        released on the strength of its first 200 files."""
        state = collect_versions_via_head(
            StubS3(), "b", [f"p/{i}" for i in range(10)], limit=3
        )
        assert state.covered == 3
        assert state.total == 10
        assert state.partial is True

    def test_no_capability_at_all_means_no_claim(self):
        stub = StubS3(deny={"head_object", "list_object_versions"})
        capabilities = probe_capabilities(stub, "b", probe_key="p/a", prefix="p/")
        state = collect_version_state(
            stub, "b", "p/", ["p/a"], capabilities=capabilities, head_limit=10
        )
        assert state.checked is False
        assert "cannot be detected" in state.reason
