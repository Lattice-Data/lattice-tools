"""The CAS registry export as a versioned input artifact.

There is no CAS batch download. The registry is a licensed product, scripted
retrieval and bulk download are forbidden by the academic terms, and API access is
not in the site licence. What exists is a human procedure: paste the CAS numbers
into SciFinder's Advanced Search "Substance RN" field in blocks of 60, export each
result set as PDF, parse the PDFs into one table.

That is not a defect to engineer away. It makes the CAS data a human-supplied,
versioned artifact, which is *better* for reproducibility than a live API call: the
table and the PDFs behind it can be hashed, kept, and pointed at months later when
ChEBI asks why we said something in September. So this module only ever reads a
parsed table, and it never fetches anything.

The consequences are designed for rather than papered over:

- A CAS number with no row cannot be checked on this axis, and the gate says so
  instead of treating silence as agreement.
- Stoichiometry comes from ``molecular_formula`` and never from ``inchi`` or
  ``canonical_smiles``. The registry gives zolantidine dimaleate a formula of
  ``C22H27N3OS.2C4H4O4`` and an InChI carrying a single maleate. The formula is
  authoritative about composition; the structure string is not, and reading the
  ratio off it silently halves the salt.
- The table is not committed. It is hashed into the run manifest, because it is
  ~160 KB derived from PDFs that are themselves the evidence.
"""

from __future__ import annotations

import csv
import io
import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

from cas_registry import CAS_VALID, classify_cas

log = logging.getLogger(__name__)

# The columns parse_scifinder_pdf.py writes. Asserted in full, because a table
# missing a column the checks read would otherwise look like a table of records
# that merely have nothing to say.
REGISTRY_COLUMNS = (
    "cas",
    "molecular_formula",
    "registry_name",
    "index_name",
    "name_ratio",
    "n_components",
    "component_cas",
    "component_formulas",
    "stereo_note",
    "molecular_weight",
    "canonical_smiles",
    "isomeric_smiles",
    "inchi",
    "inchikey",
    "n_references",
    "n_suppliers",
    "other_names",
    "source_pdf",
)

# Columns that describe composition, and are therefore allowed to drive a
# stoichiometry comparison.
STOICHIOMETRY_COLUMNS = ("molecular_formula", "name_ratio", "n_components")

# Columns that describe a structure. Useful as evidence for a human, never as a
# source of component ratios.
STRUCTURE_COLUMNS = ("canonical_smiles", "isomeric_smiles", "inchi", "inchikey")


class RegistryError(Exception):
    """Raised when the CAS registry table is missing or malformed."""


@dataclass(frozen=True)
class RegistryRow:
    """One substance as the CAS registry describes it."""

    cas: str
    molecular_formula: str
    registry_name: str
    name_ratio: str
    stereo_note: str
    molecular_weight: str
    inchikey: str
    source_pdf: str
    raw: dict[str, str]

    @property
    def structure_strings(self) -> dict[str, str]:
        """The structure columns, for a human reading the evidence.

        Deliberately grouped behind a name that says what they are, so that using
        one of them for stoichiometry is a visible decision rather than an
        accident of reaching into ``raw``.
        """
        return {name: self.raw.get(name, "") for name in STRUCTURE_COLUMNS}


def normalise_cas(raw: str) -> str:
    """The join key for a CAS number: normalised when it is recoverable, else raw.

    The table was indexed on whatever the CSV spelled while every lookup passes
    ``RecordContext.cas_key``, which is normalised -- so a row written
    ``0557-66-4`` could never be found, and EXT-01 and EXT-04 reported "no
    registry row" for a record whose row was sitting in the table. The registry is
    assembled by hand from PDF exports, so a leading zero or a stray space is the
    expected kind of defect, not an exotic one.

    Normalising on both sides rather than rejecting the row: decisions.py refuses
    a repairable CAS because a decision that cannot match anything is worse than
    no decision, but a registry row that cannot match is evidence thrown away, and
    the evidence is what the whole module is short of. Same reasoning as
    :attr:`RecordContext.cas_key`, and deliberately the same rule.
    """
    value = (raw or "").strip()
    if not value:
        return ""
    normalised, verdict, _ = classify_cas(value)
    return normalised if verdict == CAS_VALID else value


@dataclass(frozen=True)
class Registry:
    """A loaded CAS registry export, indexed by CAS number."""

    rows: dict[str, RegistryRow]
    path: Path | None = None
    sha256: str = ""
    source_pdfs: tuple[str, ...] = ()

    def get(self, cas: str) -> RegistryRow | None:
        return self.rows.get(normalise_cas(cas)) if cas else None

    def coverage(self, cas_values: set[str]) -> dict[str, int]:
        """How much of a batch this table can speak to.

        Reported in the run summary rather than inferred: clearance on a record
        with no registry row means "internally consistent", not "confirmed".
        """
        present = {c for c in cas_values if normalise_cas(c) in self.rows}
        return {
            "records": len(cas_values),
            "with_registry_row": len(present),
            "without_registry_row": len(cas_values - present),
            "registry_rows": len(self.rows),
        }

    def __len__(self) -> int:
        return len(self.rows)


EMPTY = Registry(rows={})


def load(path: str | Path) -> Registry:
    """Load a parsed SciFinder export. Raises rather than returning a partial table."""
    path = Path(path)
    if not path.exists():
        raise RegistryError(
            f"CAS registry table not found: {path}. It is produced by parsing the "
            "SciFinder PDF exports; the gate never fetches CAS data."
        )
    # Read once. The hash and the rows came from two separate reads of the same
    # path, so a table rewritten between them would be hashed as one file and
    # parsed as another -- and the manifest pins the run by that hash.
    data = path.read_bytes()
    with io.StringIO(data.decode("utf-8")) as handle:
        reader = csv.DictReader(handle)
        header = tuple(reader.fieldnames or ())
        missing = [c for c in REGISTRY_COLUMNS if c not in header]
        if missing:
            raise RegistryError(f"{path.name} is missing columns {missing}")

        rows: dict[str, RegistryRow] = {}
        pdfs: set[str] = set()
        for line_no, raw in enumerate(reader, start=2):
            cas = (raw.get("cas") or "").strip()
            if not cas:
                raise RegistryError(f"{path.name} row {line_no}: empty cas")
            # Indexed on the join key, not on the spelling. `cas` below keeps the
            # spelling, because that is what the table actually says and a report
            # quoting it should quote it verbatim.
            key = normalise_cas(cas)
            if key in rows:
                raise RegistryError(
                    f"{path.name} row {line_no}: duplicate cas {cas}; a registry "
                    "table must have one row per substance"
                )
            source_pdf = (raw.get("source_pdf") or "").strip()
            if source_pdf:
                pdfs.add(source_pdf)
            rows[key] = RegistryRow(
                cas=cas,
                molecular_formula=(raw.get("molecular_formula") or "").strip(),
                registry_name=(raw.get("registry_name") or "").strip(),
                name_ratio=(raw.get("name_ratio") or "").strip(),
                stereo_note=(raw.get("stereo_note") or "").strip(),
                molecular_weight=(raw.get("molecular_weight") or "").strip(),
                inchikey=(raw.get("inchikey") or "").strip(),
                source_pdf=source_pdf,
                raw=dict(raw),
            )

    log.info("loaded %d CAS registry rows from %s", len(rows), path)
    return Registry(
        rows=rows,
        path=path,
        sha256=hashlib.sha256(data).hexdigest(),
        source_pdfs=tuple(sorted(pdfs)),
    )
