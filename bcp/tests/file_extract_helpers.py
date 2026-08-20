"""Shared mocks and helpers for file_extract tests."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

FIXTURES = Path(__file__).parent / "fixtures" / "file_extract"


class MockPaginator:
    """Minimal paginator wrapping MockS3Client.list_objects_v2 pages."""

    def __init__(self, s3_client: "MockS3Client") -> None:
        self._s3 = s3_client

    def paginate(self, **kwargs: Any):
        bucket = kwargs.get("Bucket", "")
        prefix = kwargs.get("Prefix", "")
        self._s3.paginate_calls += 1
        pages = self._s3._paginated_pages.get((bucket, prefix))
        if pages is not None:
            yield from pages
            return
        yield self._s3.list_objects_v2(Bucket=bucket, Prefix=prefix)


class MockS3Client:
    """S3 mock for file_extract tests."""

    def __init__(
        self,
        keys: list[str] | None = None,
        sizes: dict[str, int] | None = None,
        object_bodies: dict[str, bytes | str] | None = None,
        crc_by_key: dict[str, str] | None = None,
        paginated_pages: dict[tuple[str, str], list[dict]] | None = None,
    ) -> None:
        self._keys = keys or []
        self._sizes = sizes or {}
        self._object_bodies = object_bodies or {}
        self._crc_by_key = crc_by_key or {}
        self._paginated_pages = paginated_pages or {}
        # Lets tests assert how many times a prefix was walked.
        self.paginate_calls = 0

    def get_paginator(self, operation: str) -> MockPaginator:
        assert operation == "list_objects_v2"
        return MockPaginator(self)

    def list_objects_v2(self, Bucket: str = "", Prefix: str = "") -> dict:
        pages = self._paginated_pages.get((Bucket, Prefix))
        if pages is not None:
            return pages[0] if pages else {}

        contents = []
        for key in self._keys:
            if key.startswith(Prefix):
                contents.append({"Key": key, "Size": self._sizes.get(key, 0)})
        return {"Contents": contents} if contents else {}

    def get_object_attributes(
        self,
        Bucket: str = "",
        Key: str = "",
        ObjectAttributes: list[str] | None = None,
    ) -> dict:
        if Key not in self._crc_by_key:
            raise RuntimeError("No ChecksumCRC64NVME in object attributes")
        return {"Checksum": {"ChecksumCRC64NVME": self._crc_by_key[Key]}}

    def get_object(self, Bucket: str = "", Key: str = "", **kwargs: Any) -> dict:
        if Key not in self._object_bodies:
            from botocore.exceptions import ClientError

            raise ClientError(
                {"Error": {"Code": "NoSuchKey", "Message": "Not found"}},
                "GetObject",
            )
        body = self._object_bodies[Key]
        if isinstance(body, str):
            body = body.encode("utf-8")
        range_header = kwargs.get("Range")
        if isinstance(range_header, str) and range_header.startswith("bytes="):
            spec = range_header[len("bytes=") :]
            start_s, _, end_s = spec.partition("-")
            start = int(start_s or 0)
            end = int(end_s) + 1 if end_s else len(body)
            body = body[start:end]
        return {"Body": io.BytesIO(body)}
