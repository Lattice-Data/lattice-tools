"""Discovering, distilling and querying a ChEBI flat-file release.

"Is this already in ChEBI?" needs a reverse lookup from a structure, which the
REST API does not offer, so the gate reads a bulk release. The release is ~870 MB
and is not in git, so these tests build miniature tables that reproduce the shapes
that actually broke: a molfile with embedded newlines inside a quoted column, a
definition field with embedded newlines, markup in `name` but not in `ascii_name`,
an accession table where most rows are not CAS, and a secondary identifier whose
primary is missing.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from chebi_gate import chebi_release

LISTING = """
<html><body>
<a href="LICENSE">LICENSE</a>
<a href="compounds.tsv.gz">compounds.tsv.gz</a>
<a href="database_accession.tsv.gz">database_accession.tsv.gz</a>
<a href="names.tsv.gz">names.tsv.gz</a>
<a href="reference.tsv.gz">reference.tsv.gz</a>
<a href="secondary_ids.tsv.gz">secondary_ids.tsv.gz</a>
<a href="status.tsv.gz">status.tsv.gz</a>
<a href="structure_registry.tsv.gz">structure_registry.tsv.gz</a>
<a href="structures.tsv.gz">structures.tsv.gz</a>
</body></html>
"""

# A molfile with real newlines, quoted, as structures.tsv carries it.
MOLFILE = "ethanol\n  RDKit\n\n  1  0  0  0\nM  END\n"

STRUCTURES = (
    "id\tcompound_id\tstatus_id\tmolfile\tsmiles\tstandard_inchi\t"
    "standard_inchi_key\tdimension\tdefault_structure\n"
    f'1\t16236\t1\t"{MOLFILE}"\tCCO\tInChI=1S/C2H6O\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N\t2D\tY\n'
    f'2\t99999\t1\t"{MOLFILE}"\tCCO\tInChI=1S/C2H6O\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N\t2D\tY\n'
    f'3\t30089\t1\t"{MOLFILE}"\tCC(O)=O\tInChI=1S/C2H4O2\tQTBSBXVTEAMEQO-UHFFFAOYSA-N\t2D\tY\n'
    f'4\t77777\t1\t"{MOLFILE}"\tCCO\tInChI=1S/C2H6O\tLFQSCWFLJHTTHZ-ZZZZZZZZZZ-N\t2D\tY\n'
    f'5\t55555\t1\t"{MOLFILE}"\t\t\t\t2D\tY\n'
)

# `definition` carries embedded newlines, which is what makes a line-split parse
# invent rows that do not exist.
COMPOUNDS = (
    "id\tname\tstatus_id\tsource\tparent_id\tmerge_type\tchebi_accession\t"
    "definition\tascii_name\tstars\tmodified_on\trelease_date\n"
    "16236\t<i>ethanol</i>\t3\tChEBI\t\t\tCHEBI:16236\t"
    '"A primary alcohol.\nIt is also a solvent."\tethanol\t3\t2026-01-01\t255\n'
    "99999\tsecond ethanol entry\t3\tChEBI\t\t\tCHEBI:99999\t"
    "\tsecond ethanol entry\t2\t2026-01-01\t255\n"
    "30089\tacetic acid\t3\tChEBI\t\t\tCHEBI:30089\t\tacetic acid\t3\t2026-01-01\t255\n"
    "77777\tethanol relative\t3\tChEBI\t\t\tCHEBI:77777\t\tethanol relative\t1\t"
    "2026-01-01\t255\n"
    "55555\tno structure\t3\tChEBI\t\t\tCHEBI:55555\t\tno structure\t1\t2026-01-01\t255\n"
    "90\t(-)-epicatechin\t3\tChEBI\t\t\tCHEBI:90\t\t(-)-epicatechin\t3\t2026-01-01\t255\n"
)

ACCESSIONS = (
    "id\tcompound_id\taccession_number\ttype\tstatus_id\tsource_id\n"
    "1\t16236\t64-17-5\tCAS\t1\t1\n"
    "2\t30089\t64-19-7\tCAS\t1\t1\n"
    "3\t16236\tsome-ref\tMANUAL_X_REF\t1\t1\n"
    "4\t16236\t12345\tCITATION\t1\t1\n"
    "5\t99999\t64-17-5\tCAS\t1\t1\n"
    "6\t16236\t200-578-6\tREGISTRY_NUMBER\t1\t1\n"
)

SECONDARY = (
    "compound_id\tsecondary_id\n"
    "90\t18484\n"
    "16236\t44594\n"
    "31415\t99998\n"  # a primary that is absent from compounds.tsv
)

STATUS = "id\tname\n1\tCHECKED\n3\tOK\n9\tSUBMITTED\n"


@pytest.fixture
def release(tmp_path: Path) -> Path:
    directory = tmp_path / "release"
    directory.mkdir()
    (directory / "structures.tsv").write_text(STRUCTURES)
    (directory / "compounds.tsv").write_text(COMPOUNDS)
    (directory / "database_accession.tsv").write_text(ACCESSIONS)
    (directory / "secondary_ids.tsv").write_text(SECONDARY)
    (directory / "status.tsv").write_text(STATUS)
    return directory


@pytest.fixture
def index(release: Path, tmp_path: Path) -> chebi_release.Index:
    return chebi_release.distil(
        release, tmp_path / "index", generated="2026-09-15T00:00:00Z"
    )


# ------------------------------------------------------------------ discovery


def test_discovery_matches_exact_filenames():
    """Handoff case 13: all four hard-coded paths in the original plan 404 today."""
    found = chebi_release.discover(LISTING)
    assert found["structures"] == "structures.tsv.gz"
    assert found["compounds"] == "compounds.tsv.gz"
    assert found["secondary_ids"] == "secondary_ids.tsv.gz"


def test_discovery_does_not_confuse_structure_registry_for_structures():
    """The prototype pattern-matched and broke ties on filename length.

    Its patterns matched both structure_registry.tsv.gz and structures.tsv.gz, and
    it picked the right one only because 17 characters is fewer than 25.
    """
    found = chebi_release.discover(LISTING)
    assert found["structures"] == "structures.tsv.gz"
    assert "structure_registry" not in found["structures"]


def test_discovery_accepts_an_uncompressed_release():
    listing = LISTING.replace(".tsv.gz", ".tsv")
    found = chebi_release.discover(listing)
    assert found["structures"] == "structures.tsv"


def test_a_missing_required_table_fails_loudly_with_the_listing():
    """The prototype printed "NOT FOUND automatically" and exited 0.

    A run that retrieved 0 of 6 tables therefore looked like a success to CI.
    """
    listing = LISTING.replace('<a href="structures.tsv.gz">structures.tsv.gz</a>', "")
    with pytest.raises(chebi_release.ReleaseError) as exc:
        chebi_release.discover(listing)
    assert "structures.tsv" in str(exc.value)
    assert "compounds.tsv.gz" in str(exc.value), "the listing must be in the message"


def test_a_missing_optional_table_only_warns(caplog):
    listing = LISTING.replace(
        '<a href="secondary_ids.tsv.gz">secondary_ids.tsv.gz</a>', ""
    )
    found = chebi_release.discover(listing)
    assert "secondary_ids" not in found
    assert "optional table" in caplog.text


# --------------------------------------------------------------- distillation


def test_distillation_indexes_structures_by_inchikey(index):
    assert index.exact("QTBSBXVTEAMEQO-UHFFFAOYSA-N") == ("CHEBI:30089",)
    assert index.exact("NOSUCHKEYXXXXX-UHFFFAOYSA-N") == ()
    assert index.exact(None) == ()


def test_a_key_owned_by_several_entries_keeps_all_of_them(index):
    """1,453 of 188,378 distinct keys in release 255 are held by more than one entry.

    Returning one silently hides the rest, and "already in ChEBI as X" naming only
    one of three is a report a curator cannot act on.
    """
    assert index.exact("LFQSCWFLJHTTHZ-UHFFFAOYSA-N") == ("CHEBI:16236", "CHEBI:99999")


def test_skeleton_lookup_ignores_stereo_and_charge(index):
    """A relative, not a duplicate: it stays in the submission but must be named."""
    relatives = index.skeleton("LFQSCWFLJHTTHZ-QQQQQQQQQQ-N")
    assert set(relatives) == {"CHEBI:16236", "CHEBI:99999", "CHEBI:77777"}


def test_a_structure_row_with_no_inchikey_is_skipped(index):
    assert "CHEBI:55555" not in [i for ids in index.by_inchikey.values() for i in ids]
    assert index.manifest["sources"]["structures"]["without_key"] == 1


def test_a_quoted_multi_line_molfile_does_not_corrupt_the_parse(index):
    """structures.tsv embeds molfiles in a quoted column and is 755 MB.

    Handoff case 11: it must be parsed as real delimited text with the field-size
    limit raised, never by splitting lines on tabs.
    """
    assert index.manifest["sources"]["structures"]["rows"] == 5
    assert index.manifest["measured"]["structures_with_inchikey"] == 4


def test_an_embedded_newline_in_a_definition_does_not_invent_rows(index):
    """Measured on the real release: a line-split pass reports 3 phantom rows."""
    assert index.manifest["sources"]["compounds"]["rows"] == 6


def test_labels_prefer_ascii_name_over_the_marked_up_name(index):
    """Handoff case 12, verified: `name` carries markup in 28,243 of 218,533 rows.

    `ascii_name` carries it in 1,349 and is never blank.
    """
    assert index.label("CHEBI:16236") == "ethanol"
    assert "<i>" not in index.label("CHEBI:16236")
    assert index.label("16236") == "ethanol"


def test_stars_are_carried_through(index):
    assert index.stars["16236"] == 3
    assert index.stars["77777"] == 1


# ----------------------------------------------------------------- CAS index


def test_only_cas_typed_accessions_are_indexed(index):
    """The type column holds codes, and most rows are not CAS.

    In release 255, 40,203 of 422,942 accession rows are CAS.
    """
    assert index.cas("64-17-5") == "CHEBI:16236"
    assert index.cas("64-19-7") == "CHEBI:30089"
    assert index.cas("some-ref") is None
    assert index.cas("200-578-6") is None
    assert index.cas(None) is None


def test_the_first_writer_wins_for_a_cas_claimed_twice(index):
    """438 CAS values map to more than one compound, so CAS is evidence, not identity."""
    assert index.cas("64-17-5") == "CHEBI:16236"


# ------------------------------------------------------- secondary identifiers


def test_a_secondary_identifier_resolves_to_its_primary(index):
    """The flat files carry only primary ids, so an older identifier resolves to
    nothing unless it is mapped. Handoff case 5 is this failure: one ChEBI
    identifier quoted in the original planning document does not exist.
    """
    assert index.resolve("CHEBI:18484") == "CHEBI:90"
    assert index.resolve("18484") == "CHEBI:90"
    assert index.label(index.resolve("CHEBI:18484")) == "(-)-epicatechin"


def test_a_primary_identifier_resolves_to_itself(index):
    assert index.resolve("CHEBI:90") == "CHEBI:90"


def test_an_unknown_identifier_passes_through_qualified(index):
    assert index.resolve("CHEBI:1234567") == "CHEBI:1234567"


def test_a_secondary_whose_primary_is_absent_is_tolerated(index):
    """288 of the referenced primaries are absent from compounds.tsv on release 255."""
    assert index.resolve("CHEBI:99998") == "CHEBI:31415"
    assert index.label("CHEBI:31415") == ""


def test_no_structure_bearing_id_is_ever_a_secondary(index):
    """Measured on release 255: 0 of 189,896 structure ids are secondary ids.

    This is why there is no liveness test. A structure match always lands on a
    primary entry, so the prototype's "matched entry is merged, not a duplicate"
    branch was unreachable by construction. The handoff's stated mechanism --
    liveness from parent_id -- does not work either: parent_id is non-empty in 0
    of 218,533 csv-parsed rows.
    """
    structure_ids = {
        i.split(":")[-1] for ids in index.by_inchikey.values() for i in ids
    }
    assert structure_ids.isdisjoint(index.secondary_to_primary)


# -------------------------------------------------------- pinning and reuse


def test_the_manifest_pins_every_source_table(index):
    sources = index.manifest["sources"]
    for kind in ("structures", "compounds", "database_accession", "secondary_ids"):
        assert len(sources[kind]["sha256"]) == 64
        assert sources[kind]["bytes"] > 0
        assert sources[kind]["rows"] > 0


def test_the_manifest_records_what_was_measured(index):
    assert index.manifest["measured"] == {
        "structures_with_inchikey": 4,
        "distinct_cas": 2,
        "secondary_ids": 3,
        "compounds": 6,
    }
    assert index.manifest["generated"] == "2026-09-15T00:00:00Z"


def test_distilling_twice_gives_byte_identical_index_files(release, tmp_path):
    """Handoff 4.4: the same inputs must produce the same outputs, every time."""
    first = chebi_release.distil(release, tmp_path / "a", generated="2026-09-15")
    second = chebi_release.distil(release, tmp_path / "b", generated="2026-09-15")
    assert len(first) == len(second)
    for filename in chebi_release.INDEX_FILES.values():
        assert (tmp_path / "a" / filename).read_bytes() == (
            tmp_path / "b" / filename
        ).read_bytes()


def test_an_index_reloads_to_the_same_lookups(release, tmp_path):
    built = chebi_release.distil(release, tmp_path / "index", generated="x")
    loaded = chebi_release.load_index(tmp_path / "index")
    assert loaded.by_inchikey == built.by_inchikey
    assert loaded.by_cas == built.by_cas
    assert loaded.secondary_to_primary == built.secondary_to_primary
    assert loaded.labels == built.labels


def test_a_gzipped_release_reads_the_same_as_a_plain_one(release, tmp_path):
    gzipped = tmp_path / "gz"
    gzipped.mkdir()
    for path in release.iterdir():
        with gzip.open(gzipped / f"{path.name}.gz", "wt", encoding="utf-8") as out:
            out.write(path.read_text())
    plain = chebi_release.distil(release, tmp_path / "i1", generated="x")
    compressed = chebi_release.distil(gzipped, tmp_path / "i2", generated="x")
    assert plain.by_inchikey == compressed.by_inchikey
    assert plain.labels == compressed.labels


def test_loading_a_missing_index_says_how_to_build_one(tmp_path):
    with pytest.raises(chebi_release.ReleaseError, match="flat_files"):
        chebi_release.load_index(tmp_path / "nowhere")


def _restate(index_dir: Path) -> None:
    """Rewrite the recorded stats to match the files, as a forger would have to."""
    manifest_path = index_dir / chebi_release.INDEX_MANIFEST
    manifest = json.loads(manifest_path.read_text())
    manifest["index"] = {
        name: chebi_release._file_stats(index_dir / filename)
        for name, filename in chebi_release.INDEX_FILES.items()
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def test_a_tampered_index_header_fails_the_load(release, tmp_path):
    """The column check, reached by restating the hash so it is not caught first."""
    chebi_release.distil(release, tmp_path / "index", generated="x")
    path = tmp_path / "index" / chebi_release.INDEX_FILES["cas"]
    path.write_text("cas\twrong_column\n64-17-5\t16236\n")
    _restate(tmp_path / "index")
    with pytest.raises(chebi_release.ReleaseError, match="columns are"):
        chebi_release.load_index(tmp_path / "index")


def test_a_repeated_structure_row_is_indexed_once(release, tmp_path):
    """Nothing guarantees one structure row per compound across releases.

    A repeated (compound_id, inchikey) pair made Index.exact return the same id
    twice, so EXT-02 named it twice in a single finding.
    """
    # A second structure row for compound 16236, same key. The molfile column of
    # the existing rows carries embedded newlines, so the duplicate is written as a
    # whole row rather than by copying a line.
    structures = release / "structures.tsv"
    duplicate = (
        '6\t16236\t1\t""\tCCO\tInChI=1S/C2H6O\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N\t2D\tY\n'
    )
    structures.write_text(structures.read_text() + duplicate)

    index = chebi_release.distil(release, tmp_path / "index", generated="x")
    for ids in index.by_inchikey.values():
        assert len(ids) == len(set(ids)), ids


def test_an_index_edited_after_distillation_is_refused(release, tmp_path):
    """The hash in the run manifest was the one recorded at distil time.

    So an index edited afterwards reported its original hash, and the run manifest
    -- whose whole purpose is pinning what was judged against what -- pinned a file
    that no longer existed. Every EXT-02 and EXT-03 verdict drawn from it was
    attributed to evidence that was not consulted.
    """
    index_dir = tmp_path / "index"
    chebi_release.distil(release, index_dir, generated="x")
    path = index_dir / chebi_release.INDEX_FILES["cas"]
    path.write_text(path.read_text() + "999999\t1-1-1\n")

    with pytest.raises(chebi_release.ReleaseError, match="does not match"):
        chebi_release.load_index(index_dir)


def test_an_untouched_index_loads_and_keeps_its_recorded_hashes(release, tmp_path):
    index_dir = tmp_path / "index"
    chebi_release.distil(release, index_dir, generated="x")
    loaded = chebi_release.load_index(index_dir)
    for name, filename in chebi_release.INDEX_FILES.items():
        stats = chebi_release._file_stats(index_dir / filename)
        assert loaded.manifest["index"][name]["sha256"] == stats["sha256"]


def test_an_index_with_no_recorded_stats_is_measured_rather_than_refused(
    release, tmp_path
):
    """Nothing was claimed about it, so there is nothing to contradict."""
    index_dir = tmp_path / "index"
    chebi_release.distil(release, index_dir, generated="x")
    manifest_path = index_dir / chebi_release.INDEX_MANIFEST
    manifest = json.loads(manifest_path.read_text())
    del manifest["index"]
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    loaded = chebi_release.load_index(index_dir)
    assert set(loaded.manifest["index"]) == set(chebi_release.INDEX_FILES)
    assert len(loaded.manifest["index"]["cas"]["sha256"]) == 64


def test_a_release_missing_a_table_on_disk_fails_the_distillation(release, tmp_path):
    (release / "structures.tsv").unlink()
    with pytest.raises(chebi_release.ReleaseError, match="structures.tsv not found"):
        chebi_release.distil(release, tmp_path / "index")


def test_a_release_without_secondary_ids_still_distils(release, tmp_path, caplog):
    (release / "secondary_ids.tsv").unlink()
    index = chebi_release.distil(release, tmp_path / "index", generated="x")
    assert index.secondary_to_primary == {}
    assert index.resolve("CHEBI:18484") == "CHEBI:18484"
    assert "secondary_ids table absent" in caplog.text


def test_the_empty_index_is_usable_as_a_no_op():
    """So a run with no ChEBI evidence takes the same path and reports it unchecked."""
    assert len(chebi_release.EMPTY) == 0
    assert chebi_release.EMPTY.exact("LFQSCWFLJHTTHZ-UHFFFAOYSA-N") == ()
    assert chebi_release.EMPTY.cas("64-17-5") is None
    assert chebi_release.EMPTY.resolve("CHEBI:18484") == "CHEBI:18484"


def test_the_written_manifest_is_json_on_disk(release, tmp_path):
    chebi_release.distil(release, tmp_path / "index", generated="x")
    path = tmp_path / "index" / chebi_release.INDEX_MANIFEST
    assert json.loads(path.read_text())["release_url"] == chebi_release.FLAT_FILES_URL
