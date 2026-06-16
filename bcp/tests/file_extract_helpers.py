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

    def get_object(self, Bucket: str = "", Key: str = "") -> dict:
        if Key not in self._object_bodies:
            from botocore.exceptions import ClientError

            raise ClientError(
                {"Error": {"Code": "NoSuchKey", "Message": "Not found"}},
                "GetObject",
            )
        body = self._object_bodies[Key]
        if isinstance(body, str):
            body = body.encode("utf-8")
        return {"Body": io.BytesIO(body)}
