"""No collaborator identifier may reach the tree, which is a public repo.

Scope is every *tracked* text file, not just ``bcp/tests``. The guard started
life covering test assets only, and the rest of the tree quietly accumulated the
same identifiers it exists to keep out -- lab names in module docstrings and the
notebook's config comments, real vendor order numbers in ``--cro-order`` help
text, and a real-uploads table in ``SEAHUB_QA.md`` that tied named experiments to
their defect counts. Read together on a public repo, that reads as naming a
collaborator for a naming mistake.

Tracked files only, deliberately: what gets pushed is what is exposed, and a
working tree routinely holds real S3 listings, QA outputs and scratch CSVs that
must not fail somebody's commit. That does mean this guard says nothing about
untracked files -- ``git add`` is where it takes effect.

Sanitized stand-ins, so a replacement stays recognisable as an example:
``labalpha`` / ``labbeta`` for labs, ``NVUS0000000000-NN`` for vendor orders,
``AN0000000N`` for accessions, and ``{lab}`` / ``{ExperimentID}`` where the text
is a path template rather than an example.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]

SCANNED_SUFFIXES = frozenset(
    {".py", ".ipynb", ".md", ".csv", ".tsv", ".json", ".html", ".txt", ".yml", ".yaml"}
)

BLOCKED_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:weissman|trapnell|hamazaki|marson|lange|ucsf|califano)\b", re.IGNORECASE
    ),
    re.compile(r"\bNVUS(?!0000000000-)\d{10}-\d+\b"),
    re.compile(r"\bAN(?!0000000[0-9])\d{8}\b"),
    re.compile(r"/ORPROJ1/"),
    re.compile(r"\b(?:liguo|pennyyang)\b", re.IGNORECASE),
)

SELF = Path(__file__).name


def _tracked_text_files() -> list[Path]:
    """Every tracked file this guard can read, minus itself.

    ``git ls-files`` rather than ``rglob``: an untracked real listing sitting in
    the working tree is not published and must not fail a commit that never
    touches it.
    """
    listed = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return sorted(
        REPO_ROOT / name
        for name in listed.split("\0")
        if name
        and Path(name).suffix.lower() in SCANNED_SUFFIXES
        and Path(name).name != SELF
    )


def test_tracked_files_are_sanitized() -> None:
    paths = _tracked_text_files()
    assert paths, "no tracked text files found -- is this a git checkout?"

    violations: list[str] = []
    for path in paths:
        if not path.is_file():
            continue
        content = path.read_text(encoding="utf-8", errors="replace")
        for pattern in BLOCKED_PATTERNS:
            for match in pattern.finditer(content):
                violations.append(
                    f"{path.relative_to(REPO_ROOT)} :: {pattern.pattern} "
                    f":: {match.group(0)}"
                )

    assert not violations, "Unsanitized identifiers found:\n" + "\n".join(violations)
