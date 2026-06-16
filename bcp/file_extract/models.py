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

    @property
    def has_failures(self) -> bool:
        return bool(self.failures)
