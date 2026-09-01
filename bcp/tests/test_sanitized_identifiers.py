"""Keep collaborator identifiers out of the tree, which is a public repo.

Two layers, because neither is sufficient alone.

**Known names, anywhere** (:data:`BLOCKED_PATTERNS`). A denylist is the only
thing that can catch a name in *prose* -- a docstring reading "the <lab> example"
sits in no structured position -- and it is the only layer that knows a
particular string is real. It cannot catch a name nobody has added to it.

This file names no real identifier outside that list, deliberately: it excludes
itself from the scan (see :data:`SELF`), so anything quoted here as an example
could never be flagged. Every narrative example below is invented, and so is
every probe bar one: a positive probe for the name alternation has nowhere to
draw from but the alternation itself, so that single line quotes a real surname
the list already carries.

**Structured slots, allowlisted** (:data:`ALLOWED_SLOT_TOKENS`). That ceiling is
not theoretical, twice over. The commit that first widened this guard beyond
``bcp/tests`` still shipped a surname in a notebook example -- on the same line
whose order number it rewrote -- because that surname was not in the denylist.
And the first version of *this* layer covered only the project segment written
directly after a ``czi-…`` bucket, which is not the form the tree mostly uses --
so ``smith-seahub-bcp/REF3/raw/…`` passed both layers, and a real person-shaped
lab name sat in ``--lab`` help text reading "Lab name or path, e.g. <name>" until
review found it. :data:`SLOT_PATTERNS` now covers all five ways one gets written.

An unrecognised token in any of those slots fails whether or not anyone has seen
it before, which is the whole point: adding to the denylist requires knowing a
name is real, and adding to the allowlist requires only knowing yours is not.

Scope is every *tracked* text file: what gets pushed is what is exposed. Tracked
only, deliberately -- a working tree routinely holds real S3 listings, QA outputs
and scratch CSVs, and those must not fail somebody's unrelated commit. ``git
add`` is where this takes effect.

Sanitized stand-ins, so a replacement stays recognisable as an example:
``labalpha`` / ``labbeta`` for labs, ``example-lab`` for a lab namespace,
``NVUS0000000000-NN`` for vendor orders, ``AN0000000N`` for accessions,
``LOCALPROJ`` for the local project, and ``{lab}`` / ``{ExperimentID}`` where the
text is a path template rather than an example. The zero pads are what the
lookaheads exempt, so any length of pad works.
"""

from __future__ import annotations

from collections.abc import Iterator
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
# but a typo. One order this guard's own commit had to delete by hand has nine
# digits where every other observed one has ten -- shaped like the invented
# ``NVUS123456789-66``, which the previous ``\d{10}`` could not see either. The
# real string is not quoted here: this file is excluded from its own scan, so a
# leak parked in it is the one leak nothing can catch.
#
# Names and identifiers are each one alternation rather than one pattern apiece.
# The scan covers 77 MB across ~600 tracked files and is dominated by regex, not
# I/O -- reading every file costs 0.05s, while a single case-insensitive name
# alternation costs 1.11s. Five patterns took 2.90s per pass and three take
# 1.96s, which is why they are merged; caching the reads would have bought 0.05s.
BLOCKED_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:weissman|trapnell|hamazaki|marson|lange|ucsf|califano"
        r"|liguo|pennyyang|wang|marlow)\b",
        re.IGNORECASE,
    ),
    re.compile(r"\b(?:NVUS(?!0{6,})\d{6,}-\d+|AN(?!0{6,})\d{6,})\b"),
    re.compile(r"/ORPROJ1/"),
)

# Layer 2. Five slots over six patterns, because that is how many ways a lab or
# project name gets written and only the first was covered at first:
# ``smith-seahub-bcp/REF3/raw/...`` passed both layers while
# ``labalpha-mapping-grns-perturb-seq`` was not even in the allowlist,
# which is how a ``--lab`` help string kept a real person-shaped name. Bare keys
# are the bigger surface -- of ~146 ``*-seahub-bcp`` occurrences, almost all are
# written without the bucket.
#
# Every pattern requires a letter first, which exempts the ``{lab}`` and
# ``<lab>`` placeholders without listing them. A general "hyphenated path
# segment" rule was tried and abandoned: it surfaces 189 tokens, mostly minified
# JavaScript from the Cell Ranger HTML fixtures.
#
# Each slot carries the literal it cannot match without, and the scan skips the
# regex when that literal is absent. Not premature: the big Cell Ranger HTML and
# barcode-whitelist fixtures are most of the 77 MB and contain none of these
# markers, and ``\b([A-Za-z][A-Za-z0-9-]*)-seahub-bcp\b`` has to be retried at
# every word boundary without it. Measured over the tracked tree, same match
# count either way: 1.16s to 0.21s.
SLOT_PATTERNS: tuple[tuple[str, tuple[str, ...], re.Pattern[str]], ...] = (
    ("lab", ("czi-",), re.compile(r"czi-([A-Za-z][A-Za-z0-9-]*)")),
    (
        "project",
        ("czi-",),
        re.compile(r"czi-[A-Za-z0-9{}<>-]+/([A-Za-z][A-Za-z0-9_.-]*)"),
    ),
    # The ``{lastname}`` of a SeaHub project, wherever it is written.
    (
        "project lastname",
        ("-seahub-bcp",),
        re.compile(r"\b([A-Za-z][A-Za-z0-9-]*)-seahub-bcp\b"),
    ),
    # CLI example values, which are prose to every path pattern above. Split by
    # flag, because the two slots have different shapes and therefore different
    # safe patterns.
    #
    # A project is ``{lastname}-{projectname}``, always hyphenated, so requiring a
    # hyphen or dot costs no coverage and drops the false positive the loose form
    # had: a sentence like "--project is required" reported the token ``is``.
    # Not hypothetical -- a review reply discussing this very pattern tripped it.
    (
        "cli flag",
        ("--project",),
        re.compile(r"--project[= ]+([A-Za-z][A-Za-z0-9_.]*[-.][A-Za-z0-9_.-]*)"),
    ),
    # A lab namespace is ``{lastname}`` and may be a single word, so this one
    # cannot require a hyphen -- ``--lab smith`` is the case that matters most,
    # and the hyphen rule above would wave it straight through. The cost is that
    # it also captures the next word of any sentence about the flag, so prose
    # reading "--lab is optional" reports ``is``. That is why the prose stopwords
    # are in :data:`ALLOWED_SLOT_TOKENS`; without them the sentence you are
    # reading would fail the build, and only self-exclusion hid the last one.
    (
        "cli flag",
        ("--lab",),
        re.compile(r"--lab[= ]+([A-Za-z][A-Za-z0-9_.-]*)"),
    ),
    # The lab namespace as a path, which is where a real person-shaped name was
    # found: ``--lab`` help text reading "Lab name or path, e.g. <name>".
    ("labs path", ("/labs/",), re.compile(r"/labs/([A-Za-z][A-Za-z0-9_.-]*)")),
)

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
        # Synthetic lab namespaces, as written in CLI examples and /labs/ paths.
        "example-lab",
        "other-lab",
        "labgamma",
        "notlabalpha",
        "otherlab",
        # A lab whose name is a prefix of another lab's, for the exact-prefix rule.
        "van",
        # Prose stopwords. In effect these serve the ``--lab`` slot alone: it is
        # the one pattern that cannot require a hyphen, so it also captures the
        # next word of any sentence explaining the flag. Twice now such a
        # sentence has tripped this guard -- once in a review reply, and once on
        # the very line of :data:`SLOT_PATTERNS` that asserted no such line
        # existed. No lab namespace is an English function word, so these cost
        # nothing against the real threat.
        "and",
        "are",
        "defaults",
        "is",
        "may",
        "must",
        "name",
        "not",
        "or",
        "should",
        "value",
        "was",
        # Starts with "raw" on purpose: proves the `raw` folder check is
        # segment-wise, not substring.
        "rawdata",
        # Synthetic projects.
        "labalpha-mapping-grns-perturb-seq",
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


def _readable_files() -> Iterator[tuple[Path, str]]:
    """Stream ``(path, text)``, so the whole tree is never resident at once.

    A list held ~77 MB of file text per test. Reading is not the cost here --
    every tracked file reads in 0.05s from cache -- so there is nothing to gain
    by keeping it.
    """
    for path in _tracked_text_files():
        if path.is_file():
            yield path, path.read_text(encoding="utf-8", errors="replace")


def known_identifiers(text: str) -> list[str]:
    """Layer 1 hits in ``text``."""
    return [m.group(0) for pattern in BLOCKED_PATTERNS for m in pattern.finditer(text)]


def unrecognised_slot_tokens(text: str) -> list[tuple[str, str]]:
    """Layer 2 hits in ``text``, as ``(slot, token)``."""
    return [
        (slot, m.group(1))
        for slot, literals, pattern in SLOT_PATTERNS
        if any(literal in text for literal in literals)
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
        ("NVUS123456789-66", True),
        ("AN1234567", True),
        ("the trapnell example", True),
        # ... while the stand-ins stay exempt, however long the zero pad.
        ("NVUS0000000000-11", False),
        ("AN0000000", False),
        ("AN00000001", False),
        # Layer 2, the case no denylist can reach: a lab nobody has named here.
        ("s3://czi-novogene/smith-tetrapod-atlas/x", True),
        ("s3://czi-schmidt/schmidt-seahub-bcp/EXP1", True),
        # The three slots the bucket-prefixed pattern could not see. A bare key
        # is the tree's commonest form, and a CLI help string is prose to every
        # path pattern -- which is how a person-shaped lab name survived in
        # `--lab` help until review.
        ("smith-seahub-bcp/REF3/raw/P03/432640/x.trim.cram", True),
        ("--project smith-mapping-grns-perturb-seq", True),
        ("--lab smith-jones", True),
        ("Lab name or path, e.g. smith-jones or /labs/smith-jones/", True),
        # ... and the shapes the tree legitimately uses.
        ("s3://czi-novogene/labalpha-tetrapod-atlas/x", False),
        ("labalpha-seahub-bcp/REF3/raw/P03/432640/x.trim.cram", False),
        ("--lab example-lab or /labs/example-lab/", False),
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
