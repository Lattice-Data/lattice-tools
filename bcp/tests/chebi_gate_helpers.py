"""Fixtures and record builders for the ChEBI submission gate tests.

The molfiles live as data files under ``fixtures/chebi_gate/mol/`` rather than as
string literals here. They were generated once with RDKit and are never rebuilt at
test time: a test that builds its own fixture with the library under test cannot
fail when that library changes, and the gate's contract is about *bytes*, so its
fixtures have to be bytes that do not move.

The record framing copies the real submission files exactly, which the round-trip
tests depend on: ``M  END`` is followed directly by ``> <NAME>`` with no blank line
between them, each data field is ``> <TAG>`` then its value then a blank line, so a
record ends with two newlines before its ``$$$$``.

One fixture deserves a note, because it encodes a limitation rather than a
convenience. ``decylamine_maleate`` uses a ten-carbon amine, not the ethylamine
that would read more naturally, because the gate calls the fragment with the most
heavy atoms the parent and everything else a counterion. Maleic acid has eight
heavy atoms, so with a small base the acid becomes the "parent", the base becomes
the "counterion", and the maleate/fumarate geometry check finds nothing to look
at. Real submission records are drugs, always larger than their counterions, so
the assumption holds there -- but a fixture has to respect it to exercise the
check at all.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path

MOL_DIR = Path(__file__).parent / "fixtures" / "chebi_gate" / "mol"

# The five field tags the ChEBI submission template uses, in the order the real
# files write them.
TEMPLATE_TAGS = ("NAME", "SYNONYM", "IUPAC_NAME", "CAS_NO", "RELATIONSHIP")

# ChEBI relationship codes used by the fixtures, from the class table the gate
# validates against.
ISA_HYDROCHLORIDE = "ISA36807"
ISA_MALEATE = "ISA50221"
ISA_FUMARATE = "ISA50921"
ISA_HYDROBROMIDE = "ISA48367"


@lru_cache(maxsize=None)
def molblock(name: str) -> str:
    """Load a molfile fixture by stem, e.g. ``molblock("amine_hcl")``.

    Cached because several tests build many records from the same block, and read
    as text rather than bytes so a test can retitle it with string slicing.
    """
    path = MOL_DIR / f"{name}.mol"
    if not path.exists():
        available = ", ".join(sorted(p.stem for p in MOL_DIR.glob("*.mol")))
        raise FileNotFoundError(f"no molfile fixture {name!r}; have: {available}")
    return path.read_text()


def retitle(block: str, title: str) -> str:
    """Replace a molfile's title line, which is the first line of the block.

    A rename has to update this line as well as the NAME field, or the gate's
    title-equals-NAME check fires on the record that was just fixed.
    """
    return title + block[block.index("\n") :]


def record_text(block: str, fields: Mapping[str, str]) -> str:
    """One SDF record as text, framed the way the real submission files are."""
    body = block if block.endswith("\n") else block + "\n"
    for tag, value in fields.items():
        body += f"> <{tag}>\n{value}\n\n"
    return body


def sdf_bytes(*records: str, terminate_last: bool = True) -> bytes:
    """An SDF file's bytes from record texts.

    ``terminate_last=False`` omits the final ``$$$$``, which is how a truncated
    download or a hand-edited file arrives.
    """
    out = "".join(f"{rec}$$$$\n" for rec in records)
    if not terminate_last and out.endswith("$$$$\n"):
        out = out[: -len("$$$$\n")]
    return out.encode("ascii")


def salt_record(
    name: str = "ethylamine hydrochloride",
    *,
    mol: str = "amine_hcl",
    cas: str = "557-66-4",
    synonym: str | None = "ethanamine hydrochloride",
    iupac: str | None = "ethanamine;hydrochloride",
    relationship: str | None = ISA_HYDROCHLORIDE,
    title: str | None = None,
) -> str:
    """A well-formed salt record, to be perturbed one field at a time by a test.

    Passing ``None`` for an optional field omits the tag entirely, which is a
    different condition from an empty value and the checks treat it as such.
    ``title`` defaults to ``name`` so the record is self-consistent.
    """
    block = retitle(molblock(mol), name if title is None else title)
    fields: dict[str, str] = {"NAME": name}
    if synonym is not None:
        fields["SYNONYM"] = synonym
    if iupac is not None:
        fields["IUPAC_NAME"] = iupac
    fields["CAS_NO"] = cas
    if relationship is not None:
        fields["RELATIONSHIP"] = relationship
    return record_text(block, fields)


def neutral_record(
    name: str = "ethanol",
    *,
    mol: str = "ethanol",
    cas: str = "64-17-5",
    synonym: str | None = "ethyl alcohol",
    iupac: str | None = "ethanol",
    relationship: str | None = None,
    title: str | None = None,
) -> str:
    """A well-formed non-salt record."""
    return salt_record(
        name,
        mol=mol,
        cas=cas,
        synonym=synonym,
        iupac=iupac,
        relationship=relationship,
        title=title,
    )
