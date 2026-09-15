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


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(repo: str | Path | None = None) -> str:
    """The commit the gate ran from, or "unknown" outside a work tree.

    Never raises: a run from an exported tarball is legitimate and must still
    produce a manifest, just one that says the commit is unknown.
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo) if repo else None,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return result.stdout.strip() or "unknown"


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
        """
        payload = json.dumps(
            {
                "tool": TOOL,
                "version": TOOL_VERSION,
                "inputs": [
                    {"name": i["name"], "sha256": i["sha256"]} for i in self.inputs
                ],
                "reference": self.reference,
                "decisions": self.decisions,
                "checks": self.checks,
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
        path.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=False) + "\n")
        return path
