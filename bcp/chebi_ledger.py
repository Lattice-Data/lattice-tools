"""
The cumulative record of every compound this project has already identified.

Each batch of experimental perturbations is mostly the previous batch again. In the
second real batch, **987 of 1106 distinct compounds — 2340 of 2612 source rows —
needed no network call at all**, purely because the first batch's verdicts were still
on disk to compare against. That is the largest single saving in the method and the
only one that gets bigger every time it is used.

It is also the one that cannot be regenerated. A network cache is slow to rebuild but
rebuildable; these verdicts embody a round of questions to the experimental team, a
set of hand-adjudicated promotions, and a stereochemistry correction that withdrew
eight identifiers. Losing the file means asking the lab the same 113 questions again.

So the ledger is committed, and it is deliberately **not** the per-run `master.tsv`:
that file carries 50-plus columns of intermediate state for one spreadsheet and is
derived data. This is the answer and the evidence for it, one row per compound, which
is cumulative curation knowledge and belongs in history.

**A later batch may change an answer, but never silently.** `merge()` returns every
key whose verdict or identifier moved, so a run reports the change rather than
absorbing it — the eight withdrawn identifiers above are exactly the case this exists
for.
"""

from __future__ import annotations

import csv
import re
import unicodedata
from pathlib import Path
from typing import Iterable

HERE = Path(__file__).resolve().parent
LEDGER_PATH = HERE / "chebi_answer_ledger.tsv"

FIELDS = [
    "experimental_perturbation",
    "cas",
    "verdict",
    "term_id",
    "chebi_name",
    "reason_code",
    "inchikey",
    "identity_evidence",
    "first_seen_batch",
    "last_updated_batch",
]

VERIFIED = "verified"
UNVERIFIED = "unverified"

_DASHES = re.compile(r"[‐-―−－]")


def answer_key(name: str, cas: str) -> tuple[str, str]:
    """
    The identity of a row for ledger purposes: `(normalised name, CAS)`.

    The name is case-, space-, underscore-, hyphen- and comma-insensitive, and every
    dash variant folds to `-`. That collapses `PHA 767491`, `PHA-767491` and
    `PHA_767491` onto one key, which is the whole point.

    What it deliberately does **not** fold is parentheses, stereo markers and the `?`
    a failed encoding leaves behind. Collapsing `(±)-Blebbistatin` onto
    `Blebbistatin` would hand a racemate the answer found for an unspecified-stereo
    record, which is the exact substitution the identification policy forbids; and
    letting a mojibake-damaged name inherit its repaired twin's ChEBI ID would
    launder a guess into a verdict. Both stay distinct and pay for their own lookup.

    The CAS is included in the key because the same name legitimately appears with
    two different registry numbers — a free base and its salt — and they are not the
    same answer.
    """
    folded = unicodedata.normalize("NFKC", name or "")
    folded = re.sub(r"[\s_\-,]+", "", _DASHES.sub("-", folded)).casefold()
    return folded, (cas or "").strip()


def load(path: Path | str = LEDGER_PATH) -> dict[tuple[str, str], dict[str, str]]:
    """Read the ledger, keyed by `answer_key`. Missing file reads as empty."""
    path = Path(path)
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {answer_key(r["experimental_perturbation"], r["cas"]): r for r in rows}


def lookup(
    ledger: dict[tuple[str, str], dict[str, str]], name: str, cas: str
) -> tuple[dict[str, str] | None, str]:
    """
    Find a previous answer for a row. Returns `(entry, how)`.

    `how` is `"name+cas"` for a full match, `"name_only"` or `"cas_only"` for a half
    match, and `""` when nothing matched.

    **A half match is not an answer.** Same name with a different CAS, or the same
    CAS under a different name, is the signal the method exists to interrogate: in
    the second batch nine such rows produced two real findings, one of them a
    doxorubicin row that had swapped to the free base's registry number. The caller
    must send a half match through the lookup phases with this entry attached for
    comparison, never adopt it. This function reports the distinction and refuses to
    make the judgement.
    """
    key = answer_key(name, cas)
    if key in ledger:
        return ledger[key], "name+cas"
    name_key = key[0]
    for (other_name, other_cas), entry in ledger.items():
        if other_name == name_key:
            return entry, "name_only"
    if key[1]:
        for (_, other_cas), entry in ledger.items():
            if other_cas == key[1]:
                return entry, "cas_only"
    return None, ""


def merge(
    ledger: dict[tuple[str, str], dict[str, str]],
    rows: Iterable[dict[str, str]],
    batch: str,
) -> list[tuple[tuple[str, str], dict[str, str], dict[str, str]]]:
    """
    Fold a finished batch into the ledger in place. Returns the changes it made.

    Each change is `(key, before, after)`, and a caller is expected to report them.
    The later batch wins, because a re-run under corrected logic is exactly why an
    answer would move — but an answer moving is a finding about the previous batch's
    published identifiers, not bookkeeping, so it is returned rather than swallowed.
    """
    changed = []
    for row in rows:
        key = answer_key(row["experimental_perturbation"], row["cas"])
        entry = dict(row)
        prior = ledger.get(key)
        entry["first_seen_batch"] = prior["first_seen_batch"] if prior else batch
        entry["last_updated_batch"] = batch
        if prior and (
            prior.get("verdict") != entry.get("verdict")
            or prior.get("term_id") != entry.get("term_id")
        ):
            changed.append((key, prior, entry))
        ledger[key] = entry
    return changed


def save(
    ledger: dict[tuple[str, str], dict[str, str]], path: Path | str = LEDGER_PATH
) -> int:
    """Write the ledger sorted by compound name. Returns the row count."""
    path = Path(path)
    rows = sorted(
        ledger.values(),
        key=lambda r: (r["experimental_perturbation"].casefold(), r["cas"]),
    )
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)
