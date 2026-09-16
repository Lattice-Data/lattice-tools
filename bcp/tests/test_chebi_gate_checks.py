"""One test per check: it fires when it should, and stays quiet when it should not.

A check that never fires and a check that fires on everything both produce a report
that looks plausible, so each check is exercised in both directions. Several of
these encode a defect that actually happened; those say so.
"""

from __future__ import annotations

import pytest

from chebi_gate import checks, classes, sdf, structure
from chebi_gate.checks import HIGH, INFO, LOW, MEDIUM
from tests.chebi_gate_helpers import (
    ISA_FUMARATE,
    ISA_HYDROCHLORIDE,
    ISA_MALEATE,
    ISA_SULFONATE,
    molblock,
    neutral_record,
    record_text,
    retitle,
    salt_record,
    sdf_bytes,
)


def contexts(*records: str, role: str = checks.ROLE_AUTO) -> list[checks.RecordContext]:
    parsed = sdf.parse_bytes(sdf_bytes(*records))
    return [
        checks.RecordContext(record=r, structure=structure.analyse(r), role=role)
        for r in parsed.records
    ]


def run(*records: str, role: str = checks.ROLE_AUTO) -> list[checks.Finding]:
    """Every finding, record-scoped and file-scoped, for a whole file."""
    parsed = sdf.parse_bytes(sdf_bytes(*records))
    ctxs = [
        checks.RecordContext(record=r, structure=structure.analyse(r), role=role)
        for r in parsed.records
    ]
    out: list[checks.Finding] = []
    for ctx in ctxs:
        out.extend(checks.run_record_checks(ctx))
    out.extend(
        checks.run_file_checks(checks.FileContext(contexts=ctxs, raw=parsed.raw))
    )
    return out


def ids(findings: list[checks.Finding]) -> set[tuple[str, str]]:
    return {(f.check, f.severity) for f in findings}


def detail(findings: list[checks.Finding], check: str) -> str:
    matches = [f.detail for f in findings if f.check == check]
    assert matches, f"no {check} finding in {ids(findings)}"
    return matches[0]


# --------------------------------------------------------------- the registry


def test_every_registered_check_declares_a_title_and_severities():
    assert checks.CHECKS
    for check in checks.CHECKS.values():
        assert check.title and check.title[0].isupper()
        assert check.severities
        assert set(check.severities) <= set(checks.SEVERITY_RANK)


def test_a_finding_reports_whether_its_comparison_was_independent():
    """The summary must never add circular and independent findings together."""
    finding = checks.Finding(check="CON-02", severity=HIGH, detail="x")
    assert finding.is_independent is True

    circular = checks.Finding(
        check="CON-02", severity=HIGH, detail="x", independent=False
    )
    assert circular.is_independent is False


def test_a_check_emitting_an_undeclared_severity_is_a_hard_error():
    """Otherwise a record is held or cleared for a reason nobody can look up."""
    bad = checks.Finding(check="CON-09", severity=HIGH, detail="x")
    with pytest.raises(ValueError, match="declared"):
        checks.validate_finding(bad, "CON-09")


def test_a_finding_for_an_unregistered_check_is_a_hard_error():
    bad = checks.Finding(check="NOPE-99", severity=HIGH, detail="x")
    with pytest.raises(ValueError, match="unregistered"):
        checks.validate_finding(bad, "CON-09")


def test_worst_severity_orders_high_above_info():
    findings = [
        checks.Finding(check="REL-01", severity=INFO, detail="x"),
        checks.Finding(check="CON-09", severity=LOW, detail="x"),
        checks.Finding(check="CON-02", severity=HIGH, detail="x"),
    ]
    assert checks.worst_severity(findings) == HIGH
    assert checks.worst_severity([]) is None


def test_a_waiver_downgrades_a_finding_without_hiding_it():
    finding = checks.Finding(
        check="CON-04", severity=MEDIUM, detail="stereo unspecified"
    )
    waived = checks.waive(finding, "racemate, reviewed", "casreg.csv")

    assert waived.severity == INFO
    assert "stereo unspecified" in waived.detail
    assert "WAIVED: racemate, reviewed" in waived.detail
    assert waived.waiver == "racemate, reviewed"
    assert waived.evidence == "casreg.csv"


def test_findings_carry_the_cas_number_as_their_join_key():
    """Handoff 4.1: renaming a record must not break any lookup keyed on it."""
    findings = run(salt_record(name="foo dihydrate", cas="557-66-4", iupac="foo"))
    assert all(f.cas == "557-66-4" for f in findings)
    assert all(f.name for f in findings)


# ------------------------------------------------------------ record integrity


def test_int02_reports_a_title_that_does_not_match_the_name():
    """Handoff case 15: a rename has to update the molfile title line as well."""
    renamed = salt_record(name="new name", title="old name")
    assert ("INT-02", MEDIUM) in ids(run(renamed))
    assert "mol title differs from NAME" in detail(run(renamed), "INT-02")


def test_int02_accepts_a_rename_that_also_updated_the_title():
    block = retitle(molblock("amine_hcl"), "new name")
    fixed = record_text(
        block,
        {
            "NAME": "new name",
            "SYNONYM": "x",
            "IUPAC_NAME": "x",
            "CAS_NO": "557-66-4",
            "RELATIONSHIP": ISA_HYDROCHLORIDE,
        },
    )
    assert ("INT-02", MEDIUM) not in ids(run(fixed))


def test_int02_reports_a_tag_outside_the_template():
    record = record_text(
        retitle(molblock("ethanol"), "ethanol"),
        {"NAME": "ethanol", "CAS_NO": "64-17-5", "COMMENT": "note to self"},
    )
    assert "unexpected tags ['COMMENT']" in detail(run(record), "INT-02")


def test_int03_reports_missing_required_fields():
    record = record_text(retitle(molblock("ethanol"), "ethanol"), {"NAME": "ethanol"})
    found = run(record)
    assert ("INT-03", HIGH) in ids(found)
    assert ("INT-03", MEDIUM) in ids(found)
    details = [f.detail for f in found if f.check == "INT-03"]
    assert "CAS_NO missing or empty" in details
    assert "SYNONYM absent" in details
    assert "IUPAC_NAME absent" in details


def test_int03_reports_a_salt_with_no_asserted_class():
    assert ("INT-03", HIGH) in ids(run(salt_record(relationship=None)))


def test_int03_accepts_a_complete_record():
    assert ("INT-03", HIGH) not in ids(run(salt_record()))


@pytest.mark.parametrize(
    "cas,fires",
    [
        ("557-66-4", False),
        ("64-17-5", False),
        ("64-17-4", True),
        ("5-17-64", True),
        ("not-a-cas", True),
    ],
)
def test_int04_checks_the_cas_check_digit(cas, fires):
    """Handoff case 8: valid and invalid check digits, both directions."""
    found = ids(run(salt_record(cas=cas)))
    assert (("INT-04", HIGH) in found) is fires


@pytest.mark.parametrize("cas", ["42971-09-05", "0557-66-4", "557-66-4 00:00:00"])
def test_int04_refuses_a_cas_that_only_verifies_after_repair(cas):
    """The string in the file is what ChEBI stores, so a repairable value is wrong.

    bcp/cas_registry.classify_cas calls these ``valid`` with a repair code, because
    the intended number is recoverable. That is right for ingest and wrong for a
    submission gate, so the gate reports the repair instead of applying it.
    """
    found = run(salt_record(cas=cas))
    assert ("INT-04", HIGH) in ids(found)
    assert "malformed" in detail(found, "INT-04")


def test_int04_stays_quiet_when_the_cas_is_merely_absent():
    """INT-03 already reports that; saying it twice doubles the record's count."""
    record = record_text(
        retitle(molblock("ethanol"), "ethanol"),
        {"NAME": "ethanol", "SYNONYM": "x", "IUPAC_NAME": "x", "CAS_NO": ""},
    )
    found = ids(run(record))
    assert ("INT-03", HIGH) in found
    assert ("INT-04", HIGH) not in found


def test_int07_reports_a_molfile_that_does_not_parse():
    record = record_text(
        molblock("unparseable"),
        {
            "NAME": "broken record",
            "SYNONYM": "x",
            "IUPAC_NAME": "x",
            "CAS_NO": "64-17-5",
        },
    )
    assert ("INT-07", HIGH) in ids(run(record))


def test_int07_reports_a_counts_line_that_is_not_v2000():
    block = molblock("ethanol").replace("V2000", "V3000")
    record = record_text(
        retitle(block, "ethanol"),
        {"NAME": "ethanol", "SYNONYM": "x", "IUPAC_NAME": "x", "CAS_NO": "64-17-5"},
    )
    assert ("INT-07", HIGH) in ids(run(record))


def test_int08_reports_a_structure_with_no_drawing():
    record = record_text(
        retitle(molblock("zero_coords"), "no drawing"),
        {"NAME": "no drawing", "SYNONYM": "x", "IUPAC_NAME": "x", "CAS_NO": "64-17-5"},
    )
    assert ("INT-08", MEDIUM) in ids(run(record))
    assert ("INT-08", MEDIUM) not in ids(run(neutral_record(iupac="x")))


def test_int09_reports_a_net_charge():
    """A record drawn as an ion is either missing a counterion or misnamed.

    This is the defect that held D-AP5 back: drawn as the monoanion, so it had to
    be neutralised and then re-checked against ChEBI, because neutralising it
    might make it match an entry that already exists.
    """
    record = neutral_record(name="ethanolate anion", mol="ethanolate_anion", iupac="x")
    found = run(record)
    assert ("INT-09", HIGH) in ids(found)
    assert "net charge -1" in detail(found, "INT-09")


def test_int09_accepts_an_ionically_drawn_but_neutral_salt():
    record = salt_record(mol="ionic_amine_chloride", iupac="x")
    assert ("INT-09", HIGH) not in ids(run(record))


# ---------------------------------------------------------- class and its proof


def test_int03_does_not_demand_a_class_from_a_record_that_is_only_solvated():
    """A parent plus water is not a salt, and there is no class it could assert.

    `is_salt` read "more than one fragment", so a pure hydrate was required to name
    a RELATIONSHIP class -- at `high`, which always holds the record -- while the
    class table contains no hydrate or solvate code. There was no value that would
    have released it.
    """
    record = salt_record(
        name="decylamine hydrate",
        mol="amine_hydrate",
        iupac="decan-1-amine;hydrate",
        relationship=None,
    )
    assert ("INT-03", HIGH) not in ids(run(record))
    assert not [name for name in classes.ALL_CLASSES.values() if "hydrate" in name]


def test_int03_still_demands_a_class_from_a_hydrated_salt():
    """The solvate exclusion must not let a real counterion through unclassed."""
    found = run(
        salt_record(
            name="decane-1,10-diamine dihydrochloride dihydrate",
            mol="diamine_2hcl_2h2o",
            iupac="x",
            relationship=None,
        )
    )
    assert ("INT-03", HIGH) in ids(found)
    assert "no RELATIONSHIP class asserted" in detail(found, "INT-03")


def test_int03_still_demands_a_class_from_an_unsolvated_salt():
    record = salt_record(mol="amine_hcl", iupac="x", relationship=None)
    assert ("INT-03", HIGH) in ids(run(record))


def test_con01_accepts_a_hydrochloride_drawn_as_one():
    assert ("CON-01", HIGH) not in ids(run(salt_record(iupac="x")))


def test_con01_reports_a_class_the_structure_contradicts():
    found = run(salt_record(name="foo maleate", relationship=ISA_MALEATE, iupac="x"))
    assert ("CON-01", HIGH) in ids(found)
    assert "class maleate not supported" in detail(found, "CON-01")


def test_con01_reports_a_relationship_code_the_gate_does_not_know():
    """Handoff case 17: ISA48369 broke two records until the gate learned it.

    An unknown code is a finding, not a pass: the right failure is "confirm this
    code" rather than silent acceptance of a class nobody has checked.
    """
    found = run(salt_record(relationship="ISA99999", iupac="x"))
    assert ("CON-01", HIGH) in ids(found)
    assert "unknown or unverified relationship ISA99999" in detail(found, "CON-01")


def test_con01_separates_maleate_from_fumarate_by_geometry():
    maleate = salt_record(
        name="decylamine maleate",
        mol="decylamine_maleate",
        relationship=ISA_MALEATE,
        iupac="x",
    )
    fumarate = salt_record(
        name="decylamine fumarate",
        mol="decylamine_fumarate",
        relationship=ISA_FUMARATE,
        iupac="x",
    )
    assert ("CON-01", HIGH) not in ids(run(maleate))
    assert ("CON-01", HIGH) not in ids(run(fumarate))

    swapped = salt_record(
        name="decylamine maleate",
        mol="decylamine_fumarate",
        relationship=ISA_MALEATE,
        iupac="x",
    )
    assert ("CON-01", HIGH) in ids(run(swapped))


def test_con01_reads_butenedioate_geometry_when_the_counterion_is_the_parent():
    """Maleic acid is 8 heavy atoms, so a smaller base makes it the parent.

    The geometry was read from every fragment *except* the largest, so GABA
    maleate -- 7 heavy atoms against 8 -- had no geometry to look at, `geoms` came
    back empty and CON-01 reported "class maleate not supported" at high on a
    correctly drawn record. Same shape as the sulfonate rule, which scans every
    fragment for the same reason.
    """
    record = salt_record(
        name="GABA maleate", mol="gaba_maleate", relationship=ISA_MALEATE, iupac="x"
    )
    assert ("CON-01", HIGH) not in ids(run(record))


def test_con01_still_separates_the_two_geometries_when_the_base_is_smaller():
    """Scanning every fragment must not make the check stop discriminating."""
    record = salt_record(
        name="GABA fumarate", mol="gaba_maleate", relationship=ISA_FUMARATE, iupac="x"
    )
    assert ("CON-01", HIGH) in ids(run(record))


def test_con01_does_not_let_a_lone_butenedioate_satisfy_its_own_class():
    """Reading geometry from every fragment includes the record's own, if alone."""
    record = neutral_record(
        name="maleic acid", mol="maleic_acid", relationship=ISA_MALEATE, iupac="x"
    )
    assert ("CON-01", HIGH) in ids(run(record))


def test_con01_rejects_the_sulfonate_class_when_only_the_parent_carries_s_and_o3():
    """The counterion here is HCl. The parent's own sulfonamide is not a sulfonate.

    The rule tested the substrings "S" and "O3" against every fragment, parent
    included, so any drug with a sulfur and three oxygens satisfied ISA64382 with
    no sulfonate drawn at all -- six records of the reference batch have that
    shape, and ISA64382 is one keystroke from the ISA36807 they actually carry.
    """
    found = run(
        salt_record(
            name="N-(2-methoxyphenyl)benzenesulfonamide hydrochloride",
            mol="sulfonamide_hcl",
            relationship=ISA_SULFONATE,
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) in ids(found)
    assert "class sulfonate not supported" in detail(found, "CON-01")


def test_con01_accepts_a_tosylate_as_a_sulfonate():
    record = salt_record(
        name="decylamine tosylate",
        mol="amine_tosylate",
        relationship=ISA_SULFONATE,
        iupac="x",
    )
    assert ("CON-01", HIGH) not in ids(run(record))


def test_con01_accepts_a_sulfonate_that_is_larger_than_its_own_base():
    """Why the fix cannot simply be "look at the counterions".

    3-Methyl-GABA napadisylate draws two small bases against one
    naphthalenedisulfonate, so the counterion is the fragment with the most heavy
    atoms and `Structure.counter` holds two copies of the base instead.
    """
    record = salt_record(
        name="3-methyl-GABA napadisylate",
        mol="gaba_napadisylate",
        relationship=ISA_SULFONATE,
        iupac="x",
    )
    assert ("CON-01", HIGH) not in ids(run(record))


def test_con01_reports_the_sulfonate_class_on_a_record_with_no_counterion():
    """One fragment is no counterion, whatever that fragment is made of."""
    found = run(
        neutral_record(
            name="4-methylbenzenesulfonic acid",
            mol="tosylic_acid",
            relationship=ISA_SULFONATE,
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) in ids(found)


# One row per entry in ALL_CLASSES: a fixture whose drawn counterion supports the
# class, and one that does not. Ten of the fourteen rules had no test of any kind,
# and each of them emits a HIGH that holds a record.
CLASS_FIXTURES = {
    "hydrochloride": ("amine_hcl", "amine_hi"),
    "hydrobromide": ("amine_hbr", "amine_hcl"),
    "iodide": ("amine_hi", "amine_hcl"),
    "maleate": ("decylamine_maleate", "decylamine_fumarate"),
    "fumarate": ("decylamine_fumarate", "decylamine_maleate"),
    "oxalate": ("amine_oxalate", "amine_hcl"),
    "tartrate": ("amine_tartrate", "amine_oxalate"),
    "sulfate": ("amine_sulfate", "amine_mesylate"),
    "methanesulfonate": ("amine_mesylate", "amine_sulfate"),
    "sulfonate": ("amine_tosylate", "sulfonamide_hcl"),
    "sodium": ("carboxylate_na", "carboxylate_k"),
    "potassium": ("carboxylate_k", "carboxylate_na"),
    "organic-bromide": ("quat_ammonium_br", "quat_ammonium_br_cl"),
    "quaternary-ammonium": ("quat_ammonium_br", "amine_hcl"),
}


def _structure(mol: str):
    parsed = sdf.parse_bytes(sdf_bytes(salt_record("x", mol=mol, iupac="y")))
    return structure.analyse(parsed.records[0])


def test_every_relationship_code_the_gate_knows_has_a_fixture():
    """A class added without a test is a HIGH nobody has exercised."""
    assert set(CLASS_FIXTURES) == set(classes.ALL_CLASSES.values())


@pytest.mark.parametrize("name", sorted(CLASS_FIXTURES))
def test_each_class_rule_accepts_its_own_counterion_and_rejects_another(name):
    supporting, contradicting = CLASS_FIXTURES[name]
    assert classes.class_supported(name, _structure(supporting)), name
    assert not classes.class_supported(name, _structure(contradicting)), name


# Every class whose rule requires *all* its counterions to be a given ion, with a
# fixture of that salt carrying a water of crystallisation.
SOLVATED_FIXTURES = {
    "hydrochloride": "diamine_2hcl_2h2o",
    "hydrobromide": "amine_hbr_hydrate",
    "iodide": "amine_hi_hydrate",
    "organic-bromide": "quat_ammonium_br_hydrate",
}


@pytest.mark.parametrize("name", sorted(SOLVATED_FIXTURES))
def test_a_solvate_does_not_stop_a_class_being_supported(name):
    """The fix that reached the hydrohalides and not the other two rules.

    `all(x in allowed for x in counter)` is False the moment a water of
    crystallisation is drawn, so a correctly drawn hydrate was held at HIGH -- and
    ISA48369 exists for quaternary bromides, which are commonly hydrates. The
    message did not even name a reason: `unsupported_reason` explains only the
    hydrohalide branch, so it read "counterions ['Br-', 'H2O'], geometry []".
    """
    drawn = _structure(SOLVATED_FIXTURES[name])
    assert "H2O" in drawn.counter, name
    assert classes.class_supported(name, drawn), name


def test_con01_accepts_a_quaternary_ammonium_under_its_own_class():
    """classes.py devotes a paragraph to ISA48369 not getting the hydrohalide test.

    Nothing pinned it: the stricter rule "would look like the safer choice and
    would reject every correctly reclassified record".
    """
    record = salt_record(
        name="tetramethylammonium bromide",
        mol="quat_ammonium_br",
        iupac="x",
        relationship="ISA35273",
    )
    assert ("CON-01", HIGH) not in ids(run(record))


def test_a_fragment_that_merely_begins_like_a_tartrate_is_not_one():
    """`startswith` is not formula equality: C4H6O6 is a prefix of C4H6O6S."""
    from chebi_gate.structure import Structure

    assert not classes.class_supported(
        "tartrate", Structure(parse=True, frags={"C4H6O6S": 1}, counter=["C4H6O6S"])
    )
    assert classes.class_supported(
        "tartrate", Structure(parse=True, frags={"C4H5O6-": 1}, counter=["C4H5O6-"])
    )


@pytest.mark.parametrize(
    ("name", "formula"),
    [
        ("oxalate", "C2H2O4"),
        ("tartrate", "C4H6O6"),
        ("sulfate", "H2O4S"),
        ("methanesulfonate", "CH4O3S"),
        ("maleate", "C4H4O4"),
        ("sulfonate", "C7H8O3S"),
    ],
)
def test_a_lone_counterion_does_not_satisfy_its_own_class(name, formula):
    """Four of these scanned every fragment with no "is there a counterion" test.

    So a single-fragment oxalic acid record asserting ISA64148 was "supported",
    and nothing else caught it: `is_salt` is True because a RELATIONSHIP is set,
    `has_counterion` is False so INT-03 asks for nothing, and CON-07 only fires on
    a non-salt that carries a class. The record cleared.
    """
    from chebi_gate.structure import Structure

    lone = Structure(parse=True, frags={formula: 1}, counter=[], parent_formula=formula)
    assert not classes.class_supported(name, lone), name


def test_int02_reports_a_repeated_data_tag():
    """The parser keeps the first and drops the rest, so nothing could see it.

    A duplicate never enters `fields`, so it is absent from `tag_order` and the
    unexpected-tag test cannot find it either. The record cleared -- and a cleared
    record is byte-identical, so both NAME fields went to ChEBI. A warning on
    stderr is not a verdict.
    """
    record = salt_record(iupac="x").replace(
        "> <SYNONYM>", "> <NAME>\nA SECOND NAME\n\n> <SYNONYM>", 1
    )
    found = run(record)
    assert ("INT-02", MEDIUM) in ids(found)
    assert "repeated tags ['NAME']" in detail(found, "INT-02")


def test_int02_stays_quiet_when_every_tag_appears_once():
    assert ("INT-02", MEDIUM) not in ids(run(salt_record(iupac="x")))


def test_con01_sees_a_quaternary_cation_smaller_than_its_counterion():
    """`cation_n_noh` was counted over the largest fragment, often the counterion.

    Tetramethylammonium is 5 heavy atoms against tosylate's 11, so the scan ran
    over the anion and returned 0: a record correctly asserting ISA35273 was held
    at high with no edit that would release it.
    """
    record = salt_record(
        name="tetramethylammonium tosylate",
        mol="quat_ammonium_tosylate",
        relationship="ISA35273",
        iupac="x",
    )
    assert ("CON-01", HIGH) not in ids(run(record))


def test_the_hydrohalide_guard_still_sees_a_cation_smaller_than_its_counterion():
    """The mirror case: the guard stopped catching a quaternary cation."""
    found = run(
        salt_record(
            name="tetramethylammonium tosylate",
            mol="quat_ammonium_tosylate",
            relationship=ISA_HYDROCHLORIDE,
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) in ids(found)


def test_int02_does_not_say_the_title_differs_from_a_name_that_is_absent():
    """INT-03 already reports the missing NAME; saying it twice adds nothing."""
    record = salt_record(iupac="x").replace(
        "> <NAME>\nethylamine hydrochloride\n\n", "", 1
    )
    found = run(record)
    assert ("INT-03", HIGH) in ids(found)
    assert not [f for f in found if f.check == "INT-02" and "title" in f.detail]


def test_con01_reports_a_quaternary_cation_asserted_as_a_hydrohalide():
    found = run(
        salt_record(
            name="tetramethylammonium bromide",
            mol="quat_ammonium_br",
            relationship="ISA48367",
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) in ids(found)
    assert "not a hydrohalide" in detail(found, "CON-01")


def test_con01_accepts_a_hydrated_quaternary_bromide_as_an_organic_bromide():
    """The solvate exemption belongs to every rule that demands *all* counterions
    match, not only the hydrohalides.

    It first existed inside the hydrohalide branch alone, so ``iodide`` and
    ``organic-bromide`` still tested the raw counterion list: a quaternary bromide
    monohydrate gave ``['Br-', 'H2O']``, failed ``all()``, and was held at HIGH.
    Quaternary bromides are commonly hydrates, and ISA48369 exists precisely
    because records were reclassified into it -- so this held a sound record, and
    ``unsupported_reason`` could not even name why: ``cation_n_noh`` is consulted
    only for the hydrohalide classes.
    """
    found = run(
        salt_record(
            name="tetramethylammonium bromide hydrate",
            mol="quat_ammonium_br_hydrate",
            relationship="ISA48369",
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) not in ids(found)


def test_con01_still_reports_an_organic_bromide_carrying_a_foreign_counterion():
    """The exemption must skip solvates only, never another salt-forming ion.

    Without this, widening the halide rules to ignore solvent would be
    indistinguishable from widening them to ignore everything.
    """
    found = run(
        salt_record(
            name="tetramethylammonium bromide chloride",
            mol="quat_ammonium_br_cl",
            relationship="ISA48369",
            iupac="x",
        )
    )
    assert ("CON-01", HIGH) in ids(found)
    assert "class organic-bromide not supported" in detail(found, "CON-01")


def test_con05_is_used_when_the_geometry_is_merely_undefined():
    """Missing evidence and contradicted evidence are different repairs.

    The real case was (+/-)-1-(1,2-diphenylethyl)piperidine maleate, whose C=C
    carried MDL stereo flag 3, "either". The fix is to set the bond to Z, not to
    change the asserted class, so it must not be reported as a wrong class.
    """
    record = salt_record(
        name="decylamine maleate",
        mol="maleate_bond_either",
        relationship=ISA_MALEATE,
        iupac="x",
    )
    found = ids(run(record))
    assert ("CON-05", HIGH) in found
    assert ("CON-01", HIGH) not in found


def test_int05_sees_two_spellings_of_one_cas_as_a_duplicate():
    """Every other join uses the normalised key; this one used the raw string."""
    found = run(
        salt_record(name="one", cas="557-66-4", iupac="x"),
        salt_record(name="two", cas="0557-66-4", iupac="x"),
    )
    duplicates = [f for f in found if f.check == "INT-05" and "CAS" in f.detail]
    assert duplicates, ids(found)


def test_int05_does_not_call_two_different_cas_numbers_duplicates():
    found = run(
        salt_record(name="one", cas="557-66-4", iupac="x"),
        salt_record(name="two", cas="64-17-5", mol="ethanol", iupac="x"),
    )
    assert not [f for f in found if f.check == "INT-05" and "CAS" in f.detail]


# -------------------------------------------------------------- stoichiometry


def test_con02_accepts_a_monohydrochloride_with_one_counterion():
    assert ("CON-02", HIGH) not in ids(run(salt_record(iupac="x")))


def test_con02_reports_a_di_prefix_with_one_counterion():
    found = run(
        salt_record(name="ethylamine dihydrochloride", mol="amine_hcl", iupac="x")
    )
    assert ("CON-02", HIGH) in ids(found)
    assert 'name prefix "di" implies 2 counterion(s)' in detail(found, "CON-02")


def test_con02_accepts_a_di_prefix_with_two_counterions():
    record = salt_record(name="ethylamine dihydrochloride", mol="amine_2hcl", iupac="x")
    assert ("CON-02", HIGH) not in ids(run(record))


def test_con02_counts_an_ionically_drawn_counterion_as_the_salt():
    """A counterion drawn as its anion is still a counterion, not more base.

    CON-02 kept its own counterion table, which listed the neutral acid of every
    organic counterion and almost none of the anions. So a tosylate drawn
    ionically -- cation plus ``C7H7O3S-``, net charge 0, which INT-09 is happy
    with -- had its counterion counted as *base*: ``salt_n`` 0 against ``base_n``
    2, ratio 0.0 against a wanted 1, and a high finding reading "implies 1
    counterion(s) per base, structure has 0" on a record whose stoichiometry is
    exactly right. Hydrogen maleate, oxalate, tartrate and mesylate all had it.

    Latent on the reference batch -- every ionic drawing there is a halide, an
    alkali metal or a charged parent -- so no baseline number moves, and the only
    thing standing between this and a wrongly held deposit was the drawing
    convention the depositor happened to choose.
    """
    found = run(
        salt_record(
            name="decylamine mono tosylate",
            mol="amine_tosylate_ionic",
            relationship=ISA_SULFONATE,
            iupac="x",
        )
    )
    assert ("CON-02", HIGH) not in ids(found)


def test_con02_still_reports_a_wrong_ratio_on_an_ionic_drawing():
    """Recognising the anion must not stop the check counting.

    The same fixture with a "di" prefix has to fail: one counterion against a
    prefix asserting two. Without this, counting anionic counterions would be
    indistinguishable from not counting at all.
    """
    found = run(
        salt_record(
            name="decylamine ditosylate",
            mol="amine_tosylate_ionic",
            relationship=ISA_SULFONATE,
            iupac="x",
        )
    )
    assert ("CON-02", HIGH) in ids(found)
    assert 'name prefix "di" implies 2 counterion(s)' in detail(found, "CON-02")


def test_con02_counts_an_ionic_halide_whose_charge_is_spelled_out():
    """The charge suffix is stripped before the lookup, and that is load-bearing.

    COUNTERIONS holds charge-stripped spellings, so ``"Cl-"`` is found as ``"Cl"``.
    Drop the stripping and an ionically drawn chloride counts as base -- ratio 0
    against a wanted 1 -- and unlike a tosylate there is no composition test to
    rescue it, because a halide is not a sulfonate. This is the spelling most
    likely to actually arrive: the reference batch draws Na+, I-, Br- and K+ this
    way, and those escaped the old table only because it happened to list "Cl-",
    "Br-", "I-" and "Na+" verbatim -- the very coincidence this replaces.
    """
    found = run(
        salt_record(
            name="tetramethylammonium monochloride",
            mol="ionic_amine_chloride",
            iupac="x",
        )
    )
    assert ("CON-02", HIGH) not in ids(found)


DIHYDRATE_NAME = "decane-1,10-diamine dihydrochloride dihydrate"


def _dihydrate_record(mol: str) -> str:
    return salt_record(
        name=DIHYDRATE_NAME,
        mol=mol,
        cas="6055-52-3",
        synonym="decamethylenediamine dihydrochloride dihydrate",
        iupac="x",
        relationship=ISA_HYDROCHLORIDE,
    )


def test_con02_reads_the_same_ratio_whichever_order_the_fragments_are_drawn_in():
    """The two fixtures are the same substance; only the molfile atom order differs.

    `base_counts[0]` took the first entry of a dict keyed in RDKit fragment order,
    so the denominator was a fact about the drawing. Drawn amine-first this record
    passed and drawn water-first it reported 'di implies 2, structure has 1' -- a
    high finding, which holds the record back, on a correct dihydrochloride.
    """
    parent_first = ids(run(_dihydrate_record("diamine_2hcl_2h2o")))
    water_first = ids(run(_dihydrate_record("diamine_2hcl_2h2o_water_first")))

    assert parent_first == water_first
    assert ("CON-02", HIGH) not in water_first


def test_con02_leaves_solvate_out_of_the_counterion_ratio():
    """Whether a hydrate is present is CON-03's question, not a change of ratio."""
    anhydrous = ids(
        run(salt_record(name="ethylamine dihydrochloride", mol="amine_2hcl", iupac="x"))
    )
    hydrated = ids(run(_dihydrate_record("diamine_2hcl_2h2o")))

    assert ("CON-02", HIGH) not in anhydrous
    assert ("CON-02", HIGH) not in hydrated


def test_con02_still_reports_a_wrong_ratio_on_a_solvated_record():
    """The solvate exclusion must not turn the check off for hydrates."""
    found = run(
        salt_record(
            name="decane-1,10-diamine trihydrochloride dihydrate",
            mol="diamine_2hcl_2h2o",
            cas="6055-52-3",
            iupac="x",
            relationship=ISA_HYDROCHLORIDE,
        )
    )
    assert ("CON-02", HIGH) in ids(found)
    assert 'name prefix "tri" implies 3 counterion(s)' in detail(found, "CON-02")


def test_syn04_does_not_call_an_unreadable_record_anhydrous():
    """Structure's defaults mean "unknown", not "absent".

    A record RDKit cannot read has `water == 0` because nothing was measured, so
    a `monohydrate` synonym drew a medium saying the entry is anhydrous. It is
    held by INT-07 regardless, but its report carried a claim the gate never
    established -- and the report is what a chemist works from.
    """
    found = run(
        salt_record(
            name="bad record",
            mol="unparseable",
            synonym="something monohydrate",
            iupac="x",
        )
    )
    assert ("INT-07", HIGH) in ids(found)
    assert not [f for f in found if f.check == "SYN-04" and "hydrate" in f.detail]


def test_syn04_still_reports_a_hydrate_synonym_on_a_readable_anhydrous_record():
    found = run(
        salt_record(
            name="ethylamine hydrochloride",
            mol="amine_hcl",
            synonym="ethanamine hydrochloride monohydrate",
            iupac="x",
        )
    )
    assert "hydrate synonyms on an anhydrous entry" in detail(found, "SYN-04")


def test_syn04_still_compares_two_names_when_the_structure_is_unreadable():
    """The racemic clause needs no structure, so `parse` must not gate it too."""
    found = run(
        salt_record(
            name="(R)-something hydrochloride",
            mol="unparseable",
            synonym="(+/-)-something hydrochloride",
            iupac="x",
        )
    )
    assert "racemic/rel synonyms on an enantiopure entry" in detail(found, "SYN-04")


def _con02(frags: dict[str, int], name: str) -> list[checks.Finding]:
    """CON-02 against a stated fragment composition, with no molfile in the way."""
    from chebi_gate.structure import Structure

    parsed = sdf.parse_bytes(sdf_bytes(salt_record(name, iupac="x")))
    ctx = checks.RecordContext(
        record=parsed.records[0],
        structure=Structure(parse=True, frags=frags, parent_formula="C20H24N2"),
        role=checks.ROLE_SALT,
    )
    return list(checks.con02_stoichiometry_word(ctx))


def test_con02_does_not_count_a_fragment_that_merely_begins_like_a_counterion():
    """A near-miss formula is base, on both sides of the ratio.

    CON-02 used to carry its own counterion pattern, matched as a *prefix*, so
    C4H2O4S read as a fumarate. That table is gone -- `classes.is_counterion` is
    the single one now -- but the property it got wrong is the same property the
    shared table has to keep, so this asserts through the check rather than
    against whatever the table is spelled as today.
    """
    spurious = _con02({"C20H24N2": 1, "C4H2O4S": 2}, "foo dihydrochloride")
    assert spurious, "two C4H2O4S fragments are not two counterions"
    # "structure has 0 (" and not "structure has 0": counting the look-alike in
    # the numerator alone gives 2/3, and "structure has 0.666667" contains the
    # shorter string.
    assert "structure has 0 (" in spurious[0].detail, spurious[0].detail

    # And in the denominator: a look-alike that is the base must be counted as
    # base. Excluded from both sides, base_n is 0 and the check silently declines
    # to judge a record whose ratio is plainly wrong.
    as_base = _con02({"C4H2O4S": 2, "HCl": 2}, "foo dihydrochloride")
    assert as_base, "C4H2O4S is the base here, so the ratio is 2 HCl to 2 base"
    assert "structure has 1 (" in as_base[0].detail, as_base[0].detail

    assert not _con02({"C20H24N2": 1, "C4H4O4": 2}, "foo dihydrochloride")
    assert not _con02({"C20H24N2": 1, "C4H2O4-2": 2}, "foo dihydrochloride"), (
        "the trailing charge is still a counterion"
    )


def test_con02_finds_a_base_when_every_fragment_reads_as_a_counterion():
    """`is_counterion` recognises a sulfonate by composition, so a parent can match.

    C, H, O and S with three oxygens per sulfur is a shape a drug parent can have.
    With no fragment left as base, base_n was 0, `ratio` was None and the check
    returned in silence on a record it was looking straight at -- the worst
    available failure. None of the reference batch's 290 records has that shape;
    this is the latent case.
    """
    assert not _con02_parent({"C8H8O3S": 1, "HCl": 1}, "x hydrochloride")

    wrong = _con02_parent({"C8H8O3S": 1, "HCl": 2}, "x hydrochloride")
    assert wrong, "two HCl against one base is not a monohydrochloride"
    assert "structure has 2" in wrong[0].detail


def _con02_parent(frags: dict[str, int], name: str) -> list[checks.Finding]:
    from chebi_gate.structure import Structure

    parsed = sdf.parse_bytes(sdf_bytes(salt_record(name, iupac="x")))
    ctx = checks.RecordContext(
        record=parsed.records[0],
        structure=Structure(parse=True, frags=frags, parent_formula=next(iter(frags))),
        role=checks.ROLE_SALT,
    )
    return list(checks.con02_stoichiometry_word(ctx))


def test_con08_compares_the_systematic_name_against_the_structure():
    found = run(salt_record(mol="amine_hcl", iupac="ethanamine;dihydrochloride"))
    assert ("CON-08", HIGH) in ids(found)
    assert "vs 1 HCl in structure" in detail(found, "CON-08")


def test_con08_accepts_agreeing_stoichiometry():
    record = salt_record(mol="amine_2hcl", iupac="ethanamine;dihydrochloride")
    assert ("CON-08", HIGH) not in ids(run(record))


# --------------------------------------------------------------------- solvates


def test_con03_reports_a_hydrate_with_no_water():
    found = run(salt_record(name="ethylamine hydrochloride hydrate", iupac="x"))
    assert ("CON-03", HIGH) in ids(found)
    assert "no water in structure" in detail(found, "CON-03")


def test_con03_accepts_a_hydrate_that_has_water():
    record = salt_record(
        name="decylamine hydrate", mol="amine_hydrate", iupac="x", relationship=None
    )
    assert ("CON-03", HIGH) not in ids(run(record))


def test_con03_reports_an_ethanolate_with_no_ethanol():
    record = salt_record(name="ethylamine hydrochloride monoethanolate", iupac="x")
    assert ("CON-03", HIGH) in ids(run(record))


@pytest.mark.parametrize(
    "name,iupac",
    [
        ("ethanol", "ethanol"),
        ("2-phenylethanol", "2-phenylethanol"),
        ("2-aminoethanol", "2-aminoethanol"),
        ("sodium salt dehydrate process", "x"),
    ],
)
def test_con03_does_not_fire_on_a_compound_merely_named_like_a_solvate(name, iupac):
    """Regression: the prototype used ``"ethanol" in iupac``, a substring test.

    Every alcohol whose parent chain is ethane contains the substring, so each one
    would be reported as an ethanolate with no ethanol drawn. The batch this was
    ported from contained no such compound, so it never misfired -- luck, not
    correctness.
    """
    record = neutral_record(name=name, iupac=iupac, mol="ethanol")
    assert ("CON-03", HIGH) not in ids(run(record))


def test_con03_reads_a_solvate_written_as_a_machine_name_component():
    record = salt_record(name="foo", iupac="foo;ethanol")
    assert ("CON-03", HIGH) in ids(run(record))


# ------------------------------------------------------------- stereochemistry


def test_con04_reports_a_racemate_drawn_fully_specified():
    record = salt_record(
        name="(+/-)-decylamine maleate",
        mol="decylamine_maleate",
        relationship=ISA_MALEATE,
        iupac="x",
    )
    found = run(record)
    assert ("CON-04", HIGH) in ids(found)
    assert "structure fully specified" in detail(found, "CON-04")


def test_con04_reports_a_stereospecific_name_with_an_unspecified_centre():
    record = salt_record(
        name="(Z)-decylamine maleate",
        mol="maleate_bond_either",
        relationship=ISA_MALEATE,
        iupac="x",
    )
    found = run(record)
    assert ("CON-04", MEDIUM) in ids(found)
    assert "stereo elements unspecified" in detail(found, "CON-04")


def test_con04_ignores_a_vendor_code_that_looks_like_a_stereo_descriptor():
    """ "RS 56812" must not read as an (RS) racemate, or every such record fires."""
    ctx = contexts(salt_record(name="RS 56812 hydrochloride", iupac="x"))[0]
    assert ctx.racemic is False
    assert ctx.specific is False


def test_con09_is_a_low_finding_that_never_holds_a_record_back():
    """MDL convention: chiral flag 0 on an enantiopure drawing means relative stereo.

    Worth telling a curator about and not worth refusing a submission over, which
    is exactly what the low severity encodes. 24 records in the batch carried it.
    """
    assert checks.CHECKS["CON-09"].severities == (LOW,)


# ------------------------------------------------------------- name formatting


@pytest.mark.parametrize(
    "name,severity,fragment",
    [
        ("Foo (kinase inhibitor)", MEDIUM, "annotation in parentheses"),
        ("Foo, sodium salt", LOW, "comma-style salt suffix"),
        ("PMPA", HIGH, "PMPA is ambiguous"),
    ],
)
def test_con10_reports_names_that_are_not_chemical_names(name, severity, fragment):
    found = run(neutral_record(name=name, iupac="x"))
    assert ("CON-10", severity) in ids(found)
    assert fragment in detail(found, "CON-10")


def test_con10_reports_a_shortened_greek_letter():
    found = run(neutral_record(name="(S)-(+)-a-Methylhistamine", iupac="x"))
    assert ("CON-10", LOW) in ids(found)
    assert '"a-" used where "alpha-" intended' in detail(found, "CON-10")


def test_con10_accepts_an_ordinary_chemical_name():
    assert ("CON-10", HIGH) not in ids(run(neutral_record(iupac="x")))


def test_con07_reports_a_non_salt_carrying_a_class():
    record = neutral_record(iupac="x", relationship=ISA_HYDROCHLORIDE)
    found = run(record, role=checks.ROLE_NEUTRAL)
    assert ("CON-07", MEDIUM) in ids(found)


def test_con07_reports_a_multi_fragment_record_in_a_neutral_file():
    found = run(salt_record(iupac="x"), role=checks.ROLE_NEUTRAL)
    assert ("CON-07", HIGH) in ids(found)


# ------------------------------------------------------------ synonym hygiene


def test_syn01_reports_catalogue_identifiers():
    found = run(
        salt_record(
            synonym="ethanamine hydrochloride;CHEMBL12345;SCHEMBL999", iupac="x"
        )
    )
    assert ("SYN-01", MEDIUM) in ids(found)
    assert "2 catalogue/database identifiers" in detail(found, "SYN-01")


def test_syn02_reports_a_bare_component_token():
    """The delimiter is ";" with no space, so a bare token cannot be told from a part."""
    found = run(salt_record(synonym="ethylamine hcl;hydrochloride", iupac="x"))
    assert ("SYN-02", MEDIUM) in ids(found)
    assert "hydrochloride" in detail(found, "SYN-02")


def test_syn03_reports_case_and_spacing_duplicates():
    found = run(salt_record(synonym="Ethylamine HCl;ethylamine hcl", iupac="x"))
    assert ("SYN-03", LOW) in ids(found)


def test_syn04_reports_a_racemic_synonym_on_an_enantiopure_entry():
    record = salt_record(
        name="(R)-ethylamine hydrochloride",
        synonym="(+/-)-ethylamine hydrochloride",
        iupac="x",
    )
    assert ("SYN-04", MEDIUM) in ids(run(record))


def test_syn04_reports_a_hydrate_synonym_on_an_anhydrous_entry():
    found = run(salt_record(synonym="ethylamine hydrochloride hydrate", iupac="x"))
    assert ("SYN-04", MEDIUM) in ids(found)
    assert "hydrate synonyms on an anhydrous entry" in detail(found, "SYN-04")


def test_syn05_reports_the_delimiter_written_with_a_space():
    """194 of 195 originals use ';' with no space, and '; ' trips this check."""
    found = run(salt_record(synonym="ethylamine hcl; ethanamine hcl", iupac="x"))
    assert ("SYN-05", LOW) in ids(found)


def test_syn05_accepts_the_bare_delimiter():
    record = salt_record(synonym="ethylamine hcl;ethanamine hcl", iupac="x")
    assert ("SYN-05", LOW) not in ids(run(record))


def test_syn06_reports_an_inchikey_written_as_a_synonym():
    """Handoff case 14: an InChIKey has no digits, so a word-shaped test accepts it."""
    found = run(
        salt_record(synonym="ethylamine hcl;XWBDWHCCBGMXKG-UHFFFAOYSA-N", iupac="x")
    )
    assert ("SYN-06", MEDIUM) in ids(found)
    assert "XWBDWHCCBGMXKG-UHFFFAOYSA-N" in detail(found, "SYN-06")


def test_syn06_reports_an_inchi_string_written_as_a_synonym():
    found = run(salt_record(synonym="ethylamine hcl;InChI=1S/C2H7N/c1-2-3", iupac="x"))
    assert ("SYN-06", MEDIUM) in ids(found)


def test_synonym_checks_stay_quiet_when_the_field_is_absent():
    """Absent is INT-03's business; the hygiene checks have nothing to say."""
    found = ids(run(salt_record(synonym=None, iupac="x")))
    assert ("INT-03", MEDIUM) in found
    for check in ("SYN-01", "SYN-02", "SYN-03", "SYN-04", "SYN-05", "SYN-06"):
        assert not any(c == check for c, _ in found)


def test_iup01_reports_a_machine_joined_systematic_name():
    assert ("IUP-01", LOW) in ids(run(salt_record()))
    assert ("IUP-01", LOW) not in ids(run(salt_record(iupac="ethanamine")))


# ------------------------------------------------------------- whole-file checks


def test_int05_reports_a_duplicate_cas_against_every_record_involved():
    found = [
        f
        for f in run(salt_record(), salt_record(name="other name"))
        if f.check == "INT-05"
    ]
    labels = {(f.record_index, f.severity) for f in found}
    assert labels == {(1, HIGH), (2, HIGH)}
    assert any("duplicate CAS" in f.detail for f in found)


def test_int05_identifies_the_other_record_by_cas_not_only_by_name():
    """A name-keyed message goes stale the moment one of the pair is renamed."""
    found = run(salt_record(), salt_record(name="other name"))
    message = detail([f for f in found if f.check == "INT-05"], "INT-05")
    assert "557-66-4" in message


def test_int05_reports_a_duplicate_structure_under_different_names():
    found = run(
        salt_record(name="name one", cas="557-66-4"),
        salt_record(name="name two", cas="64-17-5"),
    )
    assert any("duplicate InChIKey" in f.detail for f in found if f.check == "INT-05")


def test_int05_is_quiet_on_distinct_records():
    found = run(salt_record(), neutral_record(iupac="x"))
    assert ("INT-05", HIGH) not in ids(found)


def test_int06_reports_a_non_ascii_file_once_not_once_per_record():
    """The prototype attributed this to every record, holding back a whole file."""
    raw = sdf_bytes(salt_record(), neutral_record(synonym="cafe_ alcohol", iupac="x"))
    raw = raw.replace(b"cafe_", b"caf\xe9")
    parsed = sdf.parse_bytes(raw)
    ctxs = [
        checks.RecordContext(record=r, structure=structure.analyse(r))
        for r in parsed.records
    ]
    found = checks.run_file_checks(checks.FileContext(contexts=ctxs, raw=parsed.raw))
    encoding = [f for f in found if f.check == "INT-06"]
    assert len(encoding) == 1
    assert encoding[0].severity == HIGH
    assert encoding[0].record_index == 0


def test_int06_reports_crlf_at_low_severity():
    raw = sdf_bytes(neutral_record(iupac="x")).replace(b"\n", b"\r\n")
    parsed = sdf.parse_bytes(raw)
    ctxs = [
        checks.RecordContext(record=r, structure=structure.analyse(r))
        for r in parsed.records
    ]
    found = checks.run_file_checks(checks.FileContext(contexts=ctxs, raw=parsed.raw))
    assert ("INT-06", LOW) in {(f.check, f.severity) for f in found}


def test_rel01_reports_shared_parent_skeletons_per_record():
    """Attributed per record, so holding one record carries the information."""
    found = [
        f
        for f in run(
            salt_record(name="mono", mol="amine_hcl", cas="557-66-4"),
            salt_record(name="di", mol="amine_2hcl", cas="64-17-5"),
        )
        if f.check == "REL-01"
    ]
    assert len(found) == 2
    assert all(f.severity == INFO for f in found)
    assert {f.record_index for f in found} == {1, 2}


def test_rel01_never_holds_a_record_back():
    assert checks.CHECKS["REL-01"].severities == (INFO,)


@pytest.mark.parametrize(
    "synonym,fires",
    [
        ("ethylamine hydrochloride hydrate", True),
        ("compound monohydrate", True),
        ("carbohydrate derivative", False),
        ("sodium dehydrate", False),
    ],
)
def test_syn04_matches_the_solvate_word_not_a_substring(synonym, fires):
    """SYN-04 had the same bare-substring bug CON-03 did, and kept it longer.

    It fired on a synonym reading "carbohydrate derivative" or "sodium dehydrate",
    reporting a hydrate synonym on an anhydrous entry where there was none.
    """
    found = run(salt_record(synonym=f"ethylamine hcl;{synonym}", iupac="x"))
    assert (("SYN-04", MEDIUM) in ids(found)) is fires
