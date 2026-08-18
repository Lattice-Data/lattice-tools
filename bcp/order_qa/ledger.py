"""Prior-run state, so a re-run can say what changed.

Two facts make this necessary rather than a nicety:

* Resequencing of low-Q30 samples is routine, so an order that has already been
  reported on will legitimately be QA'd again. "Did anything change since last
  time" is the question a re-run is actually asking.
* The buckets are versioned. Re-uploading the same key creates a new version
  instead of overwriting, so a silent re-upload leaves the key count and the byte
  total identical -- it is invisible to any check that counts and sums. What does
  change is the key's current ``VersionId``, and that is only detectable against a
  record of what it was before.

The ledger is one small JSON file per order, replaced atomically, sitting beside
the timestamped run directories it summarizes. Run directories are never touched
after they are written.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

LEDGER_NAME = "ledger.json"
LEDGER_VERSION = 1


@dataclass
class Ledger:
    """What the last run of this order saw."""

    order_key: str = ""
    run_at: str = ""
    verdict: str = ""
    file_count: int = 0
    total_bytes: int = 0
    newest_object: str = ""
    version_by_key: dict[str, str] = field(default_factory=dict)
    version_method: str = ""
    version_coverage: str = ""
    report_dir: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "ledger_version": LEDGER_VERSION,
            "order_key": self.order_key,
            "run_at": self.run_at,
            "verdict": self.verdict,
            "file_count": self.file_count,
            "total_bytes": self.total_bytes,
            "newest_object": self.newest_object,
            "version_method": self.version_method,
            "version_coverage": self.version_coverage,
            "report_dir": self.report_dir,
            "version_by_key": self.version_by_key,
        }


def ledger_path(order_root: Path) -> Path:
    return order_root / LEDGER_NAME


def load_ledger(order_root: Path) -> Ledger | None:
    """Read the previous ledger, or None if this is the first run.

    A corrupt or unreadable ledger returns None with a warning rather than
    raising: losing the comparison against last time is a degraded re-run, but
    refusing to QA the order at all because a bookkeeping file went bad is worse.
    """
    path = ledger_path(order_root)
    if not path.exists():
        return None
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        log.warning(
            "Previous ledger at %s is unreadable (%s); "
            "this run cannot be compared against the last one.",
            path,
            exc,
        )
        return None
    if not isinstance(raw, dict):
        log.warning("Previous ledger at %s is not an object; ignoring.", path)
        return None
    return Ledger(
        order_key=str(raw.get("order_key", "")),
        run_at=str(raw.get("run_at", "")),
        verdict=str(raw.get("verdict", "")),
        file_count=int(raw.get("file_count", 0) or 0),
        total_bytes=int(raw.get("total_bytes", 0) or 0),
        newest_object=str(raw.get("newest_object", "")),
        version_by_key=dict(raw.get("version_by_key") or {}),
        version_method=str(raw.get("version_method", "")),
        version_coverage=str(raw.get("version_coverage", "")),
        report_dir=str(raw.get("report_dir", "")),
    )


def save_ledger(order_root: Path, ledger: Ledger) -> Path:
    """Write the ledger atomically.

    Written to a temp file in the same directory and then ``os.replace``d, so a
    run killed mid-write leaves the previous ledger intact rather than a truncated
    file that the next run would have to discard.
    """
    order_root.mkdir(parents=True, exist_ok=True)
    path = ledger_path(order_root)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(order_root), prefix=".ledger-", suffix=".tmp"
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(ledger.as_dict(), handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, path)
    except BaseException:
        # Includes KeyboardInterrupt/SystemExit: a signal during the write should
        # not leave a .ledger-*.tmp behind on a shared EFS home.
        with contextlib.suppress(OSError):
            os.unlink(tmp_name)
        raise
    return path


@dataclass
class LedgerDiff:
    """What changed between the previous run and this one."""

    compared: bool = False
    reason: str = ""
    previous_run_at: str = ""
    previous_verdict: str = ""
    new_keys: list[str] = field(default_factory=list)
    removed_keys: list[str] = field(default_factory=list)
    reuploaded_keys: list[str] = field(default_factory=list)
    unchanged: int = 0
    uncomparable_keys: int = 0

    @property
    def any_change(self) -> bool:
        return bool(self.new_keys or self.removed_keys or self.reuploaded_keys)

    @property
    def summary(self) -> str:
        if not self.compared:
            return (
                f"first run for this order ({self.reason})"
                if self.reason
                else "first run for this order"
            )
        if not self.any_change:
            note = (
                f", {self.uncomparable_keys} key(s) not version-checked"
                if self.uncomparable_keys
                else ""
            )
            return (
                f"unchanged since {self.previous_run_at} "
                f"({self.unchanged} key(s) identical{note})"
            )
        parts = []
        if self.new_keys:
            parts.append(f"{len(self.new_keys)} new")
        if self.removed_keys:
            parts.append(f"{len(self.removed_keys)} gone")
        if self.reuploaded_keys:
            parts.append(f"{len(self.reuploaded_keys)} re-uploaded")
        return f"changed since {self.previous_run_at}: {', '.join(parts)}"


def diff_against_ledger(
    previous: Ledger | None,
    *,
    current_keys: set[str],
    current_versions: dict[str, str],
) -> LedgerDiff:
    """Compare this run's keys and version identities against the last run's.

    ``reuploaded_keys`` is the payload: a key present both times whose VersionId
    changed was re-uploaded between runs. That is the case a count-and-size
    comparison cannot see, and on these buckets it is also the case that most
    often means "resequenced", so it has to be surfaced rather than inferred.

    Keys that either run could not version-check are counted, not guessed at. A
    bounded HeadObject sweep leaves most of a large order unchecked, and reporting
    those as unchanged would turn a sampling limit into a false assurance.
    """
    if previous is None:
        return LedgerDiff(compared=False)
    if not previous.version_by_key:
        return LedgerDiff(
            compared=False,
            reason=(
                "the previous run recorded no version identities "
                f"({previous.version_method or 'method unknown'}), so re-uploads "
                "cannot be detected"
            ),
            previous_run_at=previous.run_at,
            previous_verdict=previous.verdict,
        )

    previous_keys = set(previous.version_by_key)
    comparable = previous_keys & current_keys
    reuploaded = [
        key
        for key in sorted(comparable)
        if key in current_versions
        and previous.version_by_key[key]
        and current_versions[key] != previous.version_by_key[key]
    ]
    unchanged = sum(
        1
        for key in comparable
        if key in current_versions
        and current_versions[key] == previous.version_by_key[key]
    )
    return LedgerDiff(
        compared=True,
        previous_run_at=previous.run_at,
        previous_verdict=previous.verdict,
        new_keys=sorted(current_keys - previous_keys),
        removed_keys=sorted(previous_keys - current_keys),
        reuploaded_keys=reuploaded,
        unchanged=unchanged,
        uncomparable_keys=len(comparable) - unchanged - len(reuploaded),
    )
