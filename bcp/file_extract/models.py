from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class S3Location:
    """Parsed S3 URI."""

    bucket: str
    prefix: str


@dataclass(frozen=True)
class ListedObject:
    """S3 object from a listing."""

    key: str
    size_bytes: int


@dataclass
class RunSummary:
    """Aggregate counts from an extraction run."""

    total: int = 0
    crc_ok: int = 0
    enrichment_ok: int = 0
    failures: list[tuple[str, str, str]] = field(default_factory=list)
    read_tally: dict[str, int] = field(default_factory=dict)
    set_count: int = 0
    warnings: list[str] = field(default_factory=list)
    empty_raw_subdirs: list[str] = field(default_factory=list)
    md5_ok: int = 0

    @property
    def has_failures(self) -> bool:
        return bool(self.failures)


# An exception whose str() is empty still failed; the row must say so.
NO_ERROR_MESSAGE = "failed without an error message"


def failure_detail(other_err: str, md5_err: str) -> str:
    """A failure's detail: the other step's error, then the md5 one, labelled."""
    parts = [other_err] if other_err else []
    if md5_err:
        parts.append(f"md5: {md5_err}")
    return "; ".join(parts)
