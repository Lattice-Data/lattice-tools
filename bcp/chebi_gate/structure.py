"""RDKit-derived properties of one SDF record: the gate's single source of truth.

Every check that compares a *name* or an *asserted class* against what was actually
drawn reads its facts from here, and so do the external checks that compare the
drawn structure against the CAS registry or against ChEBI. Sharing one derivation
is what stops the validator and the cross-checks from disagreeing: in the
prototype this module was ported from, the gate read its structure facts from
``validate_chebi_sdf.analyse_structure`` while the cross-check scripts read theirs
from ``scripts/common.structure_keys``, and the two disagreed about whether a
record even had a molecular formula.

That disagreement was not theoretical. ``analyse_structure`` never set ``formula``,
the gate's EXT-04 passed ``s.get("formula")`` -- ``None`` -- to the registry-formula
comparator, and the comparator's ``if not a or not b: return ""`` turned every
comparison into silence. Measured on the 195-record salt file with a full evidence
cache, the gate printed ``checks run: local, EXT-01, EXT-04, EXT-02`` and produced
zero EXT-04 findings. The check announced itself and could not fire. Hence
:attr:`Structure.formula`, and hence one module.

RDKit is a hard dependency of this package and is pinned in bcp/requirements.txt.
There is no fallback: an InChIKey cannot be derived from a molfile without it, and
a gate that silently skipped structure checks when RDKit was absent would clear
records on the strength of having checked nothing.
"""

from __future__ import annotations

import logging
import re
from collections import Counter
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import inchi as rdkit_inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

from .sdf import SdfRecord

log = logging.getLogger(__name__)

# RDKit logs parse diagnostics to stderr at import-time severity. The gate reports
# a parse failure as an INT-07 finding with the record's identity attached, which
# is more useful than an unattributed warning, so the raw stream is silenced.
RDLogger.DisableLog("rdApp.*")

# Counterion formulae treated as the salt component when counting stoichiometry.
# RDKit writes these as CalcMolFormula strings, so they carry charges the way RDKit
# spells them ("Cl-", not "[Cl-]").
HALIDE_COUNTERIONS = {
    "hydrochloride": ("Cl", "HCl", "Cl-"),
    "hydrobromide": ("Br", "HBr", "Br-"),
    "iodide": ("I", "HI", "I-"),
}

WATER_FORMULA = "H2O"
ETHANOL_FORMULA = "C2H6O"


@dataclass(frozen=True)
class Structure:
    """What the drawn structure of one record actually is.

    ``parse`` is False when RDKit could not read the molfile, in which case every
    other field is at its default and the caller must not draw conclusions from
    them -- a zero ``stereo_total`` on an unparseable record means "unknown", not
    "no stereocentres".
    """

    parse: bool = False
    formula: str | None = None
    parent_formula: str | None = None
    inchikey: str | None = None
    base_inchikey: str | None = None
    n_frag: int = 0
    frags: dict[str, int] = field(default_factory=dict)
    counter: list[str] = field(default_factory=list)
    water: int = 0
    ethanol: int = 0
    net_charge: int | None = None
    zero_coords: bool = False
    cation_n_noh: int = 0
    geoms: list[str] = field(default_factory=list)
    stereo_total: int = 0
    stereo_unspec: int = 0
    tetra_specified: int = 0
    chiral_flag: int = 0

    @property
    def skeleton(self) -> str | None:
        """The connectivity block of the InChIKey: stereo and protonation stripped.

        The first 14 characters. Two records sharing a skeleton are the same
        constitution drawn with different stereochemistry, charge or salt form.
        """
        return self.inchikey[:14] if self.inchikey else None

    @property
    def parent_skeleton(self) -> str | None:
        """The connectivity block of the neutralised largest fragment's InChIKey."""
        return self.base_inchikey[:14] if self.base_inchikey else None


def analyse(record: SdfRecord) -> Structure:
    """Derive every structural fact the checks need from one record's molfile."""
    mol = Chem.MolFromMolBlock(record.mol_block, removeHs=False)
    if mol is None:
        log.debug("record %d (%s) does not parse", record.index, record.name)
        return Structure(parse=False, chiral_flag=_chiral_flag(record.counts_line))

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    frags = Chem.GetMolFrags(mol, asMols=True)
    # Ties go to the first fragment, matching the prototype. Which fragment counts
    # as the parent decides the base InChIKey and therefore every parent-skeleton
    # relationship, so the tie-break has to be stable rather than merely defined.
    biggest = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    # Identity comparison, not formula comparison: a record drawn as two copies of
    # the same fragment must keep the second copy in `counter`, or its
    # stoichiometry reads as 1:0.
    others = [f for f in frags if f is not biggest]

    return Structure(
        parse=True,
        formula=rdMolDescriptors.CalcMolFormula(mol),
        parent_formula=rdMolDescriptors.CalcMolFormula(biggest),
        inchikey=_inchikey(mol),
        base_inchikey=_base_inchikey(biggest),
        n_frag=len(frags),
        frags=dict(Counter(rdMolDescriptors.CalcMolFormula(f) for f in frags)),
        counter=[rdMolDescriptors.CalcMolFormula(f) for f in others],
        water=sum(
            1 for f in others if rdMolDescriptors.CalcMolFormula(f) == WATER_FORMULA
        ),
        ethanol=sum(
            1 for f in others if rdMolDescriptors.CalcMolFormula(f) == ETHANOL_FORMULA
        ),
        net_charge=Chem.GetFormalCharge(mol),
        zero_coords=_zero_coords(mol),
        cation_n_noh=_cation_n_noh(biggest),
        geoms=_double_bond_geometries(others),
        **_stereo_counts(mol),
        chiral_flag=_chiral_flag(record.counts_line),
    )


def _inchikey(mol: Chem.Mol) -> str | None:
    try:
        return rdkit_inchi.MolToInchiKey(mol) or None
    except Exception as exc:  # pragma: no cover - RDKit raises bare Exception
        log.debug("InChIKey generation failed: %s", exc)
        return None


def _base_inchikey(biggest: Chem.Mol) -> str | None:
    """InChIKey of the neutralised largest fragment.

    Neutralising is what makes a hydrochloride and its free base share a key, which
    is how the gate finds the parent of a salt and how it spots two records that
    are the same compound in different salt forms.
    """
    try:
        base = rdMolStandardize.Uncharger().uncharge(Chem.Mol(biggest))
        return rdkit_inchi.MolToInchiKey(base) or None
    except Exception as exc:  # pragma: no cover - RDKit raises bare Exception
        log.debug("base InChIKey generation failed: %s", exc)
        return None


def _zero_coords(mol: Chem.Mol) -> bool:
    """Whether every atom sits at the origin, i.e. the record carries no drawing."""
    if not mol.GetNumConformers():
        return False
    conf = mol.GetConformer()
    return all(
        abs(conf.GetAtomPosition(i).x) < 1e-6 and abs(conf.GetAtomPosition(i).y) < 1e-6
        for i in range(mol.GetNumAtoms())
    )


def _cation_n_noh(biggest: Chem.Mol) -> int:
    """Count quaternary-style cationic nitrogens: charged, no H, not an N-oxide.

    A hydrohalide donates its proton to a basic nitrogen, so the cation keeps an H.
    A quaternary ammonium or a pyridinium has no H to have received, which means a
    record asserting the hydrochloride class while drawing one of these is
    misclassified however well its formula adds up. Excluding N-oxides keeps that
    from firing on a neutral functional group that merely looks charged.

    ``includeNeighbors=True`` is not optional. Bare ``GetTotalNumHs()`` counts only
    implicit and property hydrogens, so a nitrogen whose hydrogens are drawn as
    explicit atoms in the molfile reports zero -- and a correctly drawn ionic
    ammonium chloride then looks like a quaternary salt, failing CON-01 for the
    hydrochloride class it legitimately holds. The prototype had the bare call.
    It never showed, because all 290 records in the batch it was written against
    draw no explicit hydrogens at all; a future batch generated by a different
    path would trip it silently.
    """
    return sum(
        1
        for atom in biggest.GetAtoms()
        if atom.GetSymbol() == "N"
        and atom.GetFormalCharge() == 1
        and atom.GetTotalNumHs(includeNeighbors=True) == 0
        and not any(
            n.GetSymbol() == "O" and n.GetFormalCharge() == -1
            for n in atom.GetNeighbors()
        )
    )


def _double_bond_geometries(others: list[Chem.Mol]) -> list[str]:
    """E/Z of the C=C in each butenedioate counterion, or "undefined".

    Maleate and fumarate are the same formula and differ only here, so this is the
    only evidence that separates the two classes. "undefined" is reported rather
    than guessed: a maleate drawn with an unspecified double bond is a defect that
    a chemist has to resolve, not a maleate.
    """
    geoms: list[str] = []
    for frag in others:
        if not re.fullmatch(
            r"C4H[234]O4[+-]?\d?", rdMolDescriptors.CalcMolFormula(frag)
        ):
            continue
        for bond in frag.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetBeginAtom().GetSymbol() == "C"
                and bond.GetEndAtom().GetSymbol() == "C"
            ):
                geoms.append(
                    {
                        Chem.BondStereo.STEREOZ: "Z",
                        Chem.BondStereo.STEREOE: "E",
                    }.get(bond.GetStereo(), "undefined")
                )
    return geoms


def _stereo_counts(mol: Chem.Mol) -> dict[str, int]:
    """Total, unspecified and specified-tetrahedral stereo element counts.

    Double bonds are counted only when both ends are carbon. RDKit reports
    potential stereo on C=N and N=N too, and those are routinely left unspecified
    in a correctly drawn record, so including them would make almost every record
    look stereochemically incomplete.
    """
    elements = []
    for element in Chem.FindPotentialStereo(mol):
        if str(element.type).endswith("Bond_Double"):
            bond = mol.GetBondWithIdx(element.centeredOn)
            if not (
                bond.GetBeginAtom().GetSymbol() == "C"
                and bond.GetEndAtom().GetSymbol() == "C"
            ):
                continue
        elements.append(element)

    return {
        "stereo_total": len(elements),
        "stereo_unspec": sum(
            1
            for e in elements
            if str(e.specified).endswith("Unspecified")
            or str(e.specified).endswith("Unknown")
        ),
        "tetra_specified": sum(
            1
            for e in elements
            if str(e.type).endswith("Atom_Tetrahedral")
            and str(e.specified).endswith("Specified")
        ),
    }


def _chiral_flag(counts_line: str) -> int:
    """The molfile chiral flag, columns 13-15 of the counts line.

    0 means the drawing asserts relative stereochemistry only, 1 means absolute.
    Parsed defensively: the prototype's ``int(counts[12:15] or 0)`` raises
    ValueError on a counts line whose chiral-flag field is blank-padded rather
    than absent, because ``" "`` is truthy and ``int(" ")`` is an error.
    """
    raw = counts_line[12:15].strip()
    if not raw:
        return 0
    try:
        return int(raw)
    except ValueError:
        log.warning("unparseable chiral flag %r in counts line", counts_line[12:15])
        return 0
