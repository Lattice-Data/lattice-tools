"""The run manifest: what was judged, against what, with which decisions.

Without this, a validation result from three months ago cannot be explained, and
ChEBI asking "why did you say this in September" has no answer. With it a run is
replayable: the manifest names every input and reference table by SHA256, the
decisions data by hash and row count, which checks ran, which could not, and how
many findings came from evidence of an independent origin.

The run identifier is derived from the inputs, not from the clock. That is what
makes the SDF outputs byte-identical across runs -- handoff 4.4 -- while still
letting a held-back record say which configuration produced it. In the prototype
the identifier was a timestamp written into every annotated record, so no two runs
could agree byte for byte and the invariant could not even be tested. The
wall-clock time is still recorded, once, here.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

TOOL = "chebi_gate"
TOOL_VERSION = "1.0"

MANIFEST_FILENAME = "run_manifest.json"

# Keys whose value names a filesystem location or a moment, not a thing. They are
# recorded in the manifest and excluded from the run identifier: "path" for every
# input, table and cache, "release_dir" and "generated" from the ChEBI index
# manifest, which is embedded whole.
RUN_ID_IGNORED_KEYS = frozenset({"path", "release_dir", "generated"})


def run_id_identity(value):
    """``value`` with every :data:`RUN_ID_IGNORED_KEYS` key removed, recursively."""
    if isinstance(value, dict):
        return {
            key: run_id_identity(item)
            for key, item in value.items()
            if key not in RUN_ID_IGNORED_KEYS
        }
    if isinstance(value, list):
        return [run_id_identity(item) for item in value]
    return value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def directory_digest(directory: str | Path, pattern: str = "*") -> dict:
    """Pin a cache directory by what is in it, not by how many files it has.

    The PubChem and CAS Common Chemistry caches were recorded as a path and a
    ``*.json`` count, so re-fetching a cache in place -- the ordinary thing to do
    when a record's entry is wrong -- left the manifest and the run identifier
    byte for byte unchanged. Two runs over completely different evidence claimed
    to be the same run.

    The digest is over sorted ``name\0sha256\0`` pairs, so it does not depend on
    directory order and does notice a file being renamed rather than edited.
    """
    directory = Path(directory)
    files = sorted(p for p in directory.glob(pattern) if p.is_file())
    rolling = hashlib.sha256()
    for path in files:
        rolling.update(path.name.encode("utf-8"))
        rolling.update(b"\0")
        rolling.update(sha256_file(path).encode("ascii"))
        rolling.update(b"\0")
    return {
        "path": str(directory),
        "records": len(files),
        "sha256": rolling.hexdigest(),
    }


def _git(repo: str | Path | None, *args: str) -> str | None:
    """One git command's stdout, or None when git could not answer.

    Never raises: a run from an exported tarball is legitimate and must still
    produce a manifest, just one that says the commit is unknown.
    """
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(repo) if repo else None,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def git_commit(repo: str | Path | None = None) -> str:
    """The commit the gate ran from, or "unknown" outside a work tree."""
    out = _git(repo, "rev-parse", "HEAD")
    return (out or "").strip() or "unknown"


def git_worktree_state(repo: str | Path | None = None) -> dict:
    """Whether the tracked files differ from the commit, and how.

    The commit alone does not say what ran. An uncommitted edit to a check --
    which is the normal state of a working session -- left the manifest naming a
    commit whose code produced different findings from the ones recorded beside
    it, with nothing to indicate the discrepancy. The manifest exists to explain a
    verdict months later; a verdict attributed to the wrong code is worse than one
    attributed to none.

    Untracked files are excluded deliberately: a run directory or an editor
    scratch file sitting in the tree does not change what the gate does, and
    counting them would make every manifest dirty and the flag worth nothing.

    ``dirty`` is None when git could not answer, which is a different thing from
    a clean tree and has to read differently.
    """
    status = _git(repo, "status", "--porcelain", "--untracked-files=no")
    if status is None:
        return {"git_dirty": None, "git_diff_sha256": "unknown"}
    if not status.strip():
        return {"git_dirty": False, "git_diff_sha256": ""}
    diff = _git(repo, "diff", "HEAD") or ""
    return {
        "git_dirty": True,
        # So two runs from the same uncommitted edit are recognisably the same
        # code, and two different edits are not.
        "git_diff_sha256": hashlib.sha256(diff.encode("utf-8")).hexdigest(),
    }


@dataclass
class Manifest:
    """Everything needed to explain or replay one run."""

    inputs: list[dict] = field(default_factory=list)
    reference: dict = field(default_factory=dict)
    decisions: dict = field(default_factory=dict)
    checks: dict = field(default_factory=dict)
    results: dict = field(default_factory=dict)
    environment: dict = field(default_factory=dict)
    timestamp: str = ""
    run_id: str = ""

    def add_input(self, path: str | Path, records: int) -> None:
        path = Path(path)
        self.inputs.append(
            {
                "path": str(path),
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "records": records,
            }
        )

    def compute_run_id(self) -> str:
        """A stable identifier for this configuration of inputs and evidence.

        Derived from the inputs, the reference hashes, the decisions hashes and the
        tool version -- deliberately not from the time, so two runs over the same
        material produce the same identifier and therefore the same output bytes.

        "The same material" has to mean the same *content*. The payload used to be
        the reference and decisions blocks verbatim, and those carry the absolute
        path of every table, cache and decisions file, so the same evidence checked
        out at a different place gave a different identifier -- the manifest could
        not be replayed on another machine, which is most of what it is for. Worse,
        the ChEBI index manifest is embedded whole and carries ``generated``, the
        wall-clock string from when the index was distilled. A clock value reached
        the run identifier despite the docstring above saying it could not.

        :data:`RUN_ID_IGNORED_KEYS` is the list of keys that say where something is
        or when it was built rather than what it is. They stay in the manifest,
        which is the right place for them.
        """
        payload = json.dumps(
            {
                "tool": TOOL,
                "version": TOOL_VERSION,
                "inputs": [
                    {"name": i["name"], "sha256": i["sha256"]} for i in self.inputs
                ],
                "reference": run_id_identity(self.reference),
                "decisions": run_id_identity(self.decisions),
                "checks": run_id_identity(self.checks),
            },
            sort_keys=True,
        )
        self.run_id = hashlib.sha256(payload.encode()).hexdigest()[:12]
        return self.run_id

    def record_environment(self, repo: str | Path | None = None) -> None:
        from rdkit import rdBase

        self.environment = {
            "tool": TOOL,
            "tool_version": TOOL_VERSION,
            "git_commit": git_commit(repo),
            **git_worktree_state(repo),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "rdkit": rdBase.rdkitVersion,
        }

    def to_dict(self) -> dict:
        return {
            "run_id": self.run_id,
            "timestamp": self.timestamp,
            "environment": self.environment,
            "inputs": self.inputs,
            "reference": self.reference,
            "decisions": self.decisions,
            "checks": self.checks,
            "results": self.results,
        }

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
        return path
