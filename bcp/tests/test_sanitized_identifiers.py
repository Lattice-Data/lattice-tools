"""Keep collaborator identifiers out of the tree, which is a public repo.

Two layers, because neither is sufficient alone.

**Known names, anywhere** (:data:`BLOCKED_PATTERNS`). A denylist is the only
thing that can catch a name in *prose* -- "the trapnell example" in a docstring
sits in no structured position -- and it is the only layer that knows a
particular string is real. It cannot catch a name nobody has added to it.

**Structured slots, allowlisted** (:data:`ALLOWED_SLOT_TOKENS`). That ceiling is
not theoretical: the commit that first widened this guard beyond ``bcp/tests``
still shipped ``wang-tetrapod-atlas``, on the same notebook line whose order
number it rewrote, because "wang" was not in the denylist. Nor would a
``czi-<lab>`` or ``*-seahub-bcp`` rule have found it -- it sat in the project
segment after a *vendor* bucket, a third slot. So the two slots where a
collaborator name actually lands are allowlisted instead, and an unrecognised
token fails whether or not anyone has seen it before.

Scope is every *tracked* text file: what gets pushed is what is exposed. Tracked
only, deliberately -- a working tree routinely holds real S3 listings, QA outputs
and scratch CSVs, and those must not fail somebody's unrelated commit. ``git
add`` is where this takes effect.

Sanitized stand-ins, so a replacement stays recognisable as an example:
``labalpha`` / ``labbeta`` for labs, ``NVUS0000000000-NN`` for vendor orders,
``AN0000000N`` for accessions, ``LOCALPROJ`` for the local project, and ``{lab}``
/ ``{ExperimentID}`` where the text is a path template rather than an example.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SELF = Path(__file__).resolve().relative_to(REPO_ROOT).as_posix()

SCANNED_SUFFIXES = frozenset(
    {".py", ".ipynb", ".md", ".csv", ".tsv", ".json", ".html", ".txt", ".yml", ".yaml"}
)

# Layer 1. The digit runs are ``(?!0{6,})\d{6,}`` rather than a fixed width: the
# stand-in convention is a long zero pad, so exempting that is what the lookahead
# is for, while pinning an exact length let a real identifier through on nothing
# but a typo. ``NVUS202101701-66`` has nine digits where every other observed
# order has ten, and the previous ``\d{10}`` could not see it -- a string this
# guard's own commit had to delete by hand.
BLOCKED_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:weissman|trapnell|hamazaki|marson|lange|ucsf|califano)\b", re.IGNORECASE
    ),
    re.compile(r"\bNVUS(?!0{6,})\d{6,}-\d+\b"),
    re.compile(r"\bAN(?!0{6,})\d{6,}\b"),
    re.compile(r"/ORPROJ1/"),
    re.compile(r"\b(?:liguo|pennyyang)\b", re.IGNORECASE),
)

# Layer 2. The lab slot of a bucket, and the project slot beneath it -- which by
# the documented convention is ``{lastname}-{projectname}``, so it is precisely
# where a surname lands. Both require a letter first, which exempts the ``{lab}``
# and ``<lab>`` placeholders without listing them.
LAB_SLOT = re.compile(r"czi-([A-Za-z][A-Za-z0-9-]*)")
PROJECT_SLOT = re.compile(r"czi-[A-Za-z0-9{}<>-]+/([A-Za-z][A-Za-z0-9_.-]*)")

# Every lab, project and ExperimentID token the tree is allowed to name. Adding
# one is meant to be a deliberate act: if a real name needs an example, invent a
# synthetic stand-in instead and add that.
ALLOWED_SLOT_TOKENS = frozenset(
    {
        # Sequencing vendors. Not collaborators, and functional -- `provider`
        # values drive vendor-specific parsing, so these cannot be placeholders.
        "novogene",
        "psomagen",
        # Synthetic labs.
        "lab",
        "labalpha",
        "labbeta",
        "synthetic",
        "other",
        "wrong",
        "x",
        "never-heard-of-it",
        # A deliberately hyphenated lab, for the exact-prefix project rule.
        "van-der-berg",
        # Synthetic projects.
        "Lab-Seahub-Beta",
        "Project_Embryo_Alpha",
        "Project_Scaling_Alpha",
        "lab-alpha",
        "lab-seahub-alpha",
        "lab-seahub-beta",
        "labalpha-seahub-bcp",
        "labalpha-tetrapod-atlas",
        "myproj",
        "p",
        "proj",
        "project",
        "project-alpha",
        "project-devdelta",
        "project-embryo-alpha",
        "project-immune-gamma",
        "project-killifish",
        "project-persona",
        "project-perturb-alpha",
        "project-scaling-alpha",
        # Starts with "raw" on purpose: proves the `raw` folder check is
        # segment-wise, not substring.
        "rawdata-seahub-bcp",
        "synthetic-project-alpha",
        "test-project",
        "unrelated",
        # ExperimentIDs, which share the project slot in bucket+order paths.
        "REF3",
        "REF3_P05_1",
    }
)


def _tracked_text_files() -> list[Path]:
    """Every tracked file this guard can read, minus itself.

    ``git ls-files`` rather than ``rglob``: an untracked real listing in the
    working tree is not published and must not fail a commit that never touches
    it. A git failure fails the test rather than skipping it -- a guard that
    passes without having checked anything is worse than one that errors.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:  # pragma: no cover - git absent
        pytest.fail("no git binary, so tracked files cannot be listed")
    if result.returncode != 0:
        pytest.fail(
            f"git ls-files exited {result.returncode} in {REPO_ROOT}, so this "
            f"guard cannot see what is tracked: {result.stderr.strip()}"
        )

    paths = sorted(
        REPO_ROOT / name
        for name in result.stdout.split("\0")
        if name and Path(name).suffix.lower() in SCANNED_SUFFIXES and name != SELF
    )
    assert paths, f"git listed no tracked text files under {REPO_ROOT}"
    return paths


def _readable_files() -> list[tuple[Path, str]]:
    return [
        (p, p.read_text(encoding="utf-8", errors="replace"))
        for p in _tracked_text_files()
        if p.is_file()
    ]


def known_identifiers(text: str) -> list[str]:
    """Layer 1 hits in ``text``."""
    return [m.group(0) for pattern in BLOCKED_PATTERNS for m in pattern.finditer(text)]


def unrecognised_slot_tokens(text: str) -> list[tuple[str, str]]:
    """Layer 2 hits in ``text``, as ``(slot, token)``."""
    return [
        (slot, m.group(1))
        for slot, pattern in (("lab", LAB_SLOT), ("project", PROJECT_SLOT))
        for m in pattern.finditer(text)
        if m.group(1) not in ALLOWED_SLOT_TOKENS
    ]


def test_no_known_identifier_appears_anywhere() -> None:
    violations = [
        f"{path.relative_to(REPO_ROOT)} :: {hit}"
        for path, content in _readable_files()
        for hit in known_identifiers(content)
    ]

    assert not violations, "Unsanitized identifiers found:\n" + "\n".join(violations)


def test_path_slots_use_synthetic_tokens() -> None:
    """The layer that can catch a name nobody thought to add to the denylist."""
    violations = [
        f"{path.relative_to(REPO_ROOT)} :: {slot} slot :: {token!r}"
        for path, content in _readable_files()
        for slot, token in unrecognised_slot_tokens(content)
    ]

    assert not violations, (
        "Unrecognised token in a lab/project path slot:\n"
        + "\n".join(sorted(set(violations)))
        + "\n\nThat slot is where a collaborator name lands. Replace it with a "
        "synthetic stand-in, or add it to ALLOWED_SLOT_TOKENS if it already is one."
    )


@pytest.mark.parametrize(
    "text,flagged",
    [
        # Layer 1, both holes this file's own history produced. A nine-digit
        # order where every observed one has ten, and a short accession: the
        # previous fixed-width patterns could see neither.
        ("NVUS202101701-66", True),
        ("AN1234567", True),
        ("the trapnell example", True),
        # ... while the stand-ins stay exempt, however long the zero pad.
        ("NVUS0000000000-11", False),
        ("AN0000000", False),
        ("AN00000001", False),
        # Layer 2, the case no denylist can reach: a lab nobody has named here.
        ("s3://czi-novogene/wang-tetrapod-atlas/x", True),
        ("s3://czi-schmidt/schmidt-seahub-bcp/EXP1", True),
        # ... and the shapes the tree legitimately uses.
        ("s3://czi-novogene/labalpha-tetrapod-atlas/x", False),
        ("s3://czi-{lab}/{lab}-seahub-bcp/REF3", False),
        ("bucket is not of the form czi-<lab>", False),
    ],
)
def test_the_patterns_actually_fire(text: str, flagged: bool) -> None:
    """A guard that only ever passes is indistinguishable from no guard.

    The tree is clean, so neither test above can tell a working pattern from a
    broken one. These probes can.
    """
    hits = known_identifiers(text) + [t for _slot, t in unrecognised_slot_tokens(text)]

    assert bool(hits) is flagged, f"{text!r} -> {hits}"
