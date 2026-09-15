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
import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

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


@dataclass(frozen=True)
class Registry:
    """A loaded CAS registry export, indexed by CAS number."""

    rows: dict[str, RegistryRow]
    path: Path | None = None
    sha256: str = ""
    source_pdfs: tuple[str, ...] = ()

    def get(self, cas: str) -> RegistryRow | None:
        return self.rows.get(cas.strip()) if cas else None

    def coverage(self, cas_values: set[str]) -> dict[str, int]:
        """How much of a batch this table can speak to.

        Reported in the run summary rather than inferred: clearance on a record
        with no registry row means "internally consistent", not "confirmed".
        """
        present = {c for c in cas_values if c in self.rows}
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
    data = path.read_bytes()
    with path.open(newline="", encoding="utf-8") as handle:
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
            if cas in rows:
                raise RegistryError(
                    f"{path.name} row {line_no}: duplicate cas {cas}; a registry "
                    "table must have one row per substance"
                )
            source_pdf = (raw.get("source_pdf") or "").strip()
            if source_pdf:
                pdfs.add(source_pdf)
            rows[cas] = RegistryRow(
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
