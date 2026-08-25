"""
Corrected SOP names for a non-compliant SeaHub upload, plus a per-well verdict.

The SOP rules in :mod:`qa_seahub_sop` each report one defect and the repair for
that defect alone.  This module composes them into a single corrected key per
object, and rolls the result up to one verdict per well.

Where the corrected sublibrary comes from
-----------------------------------------
Not from the vendor path.  The vendor layout is
``{project}/{order}/{ExperimentID}/raw/{wafer}/`` -- measured on a real delivery,
the segment before ``raw`` is the ExperimentID alone (``REF3``) and carries no
sublibrary at all.  It only looks like a sublibrary in older fixtures that use
``{ExperimentID}_{sublibrary}``.

Instead it comes from :func:`qa_seahub_sop.seahub_group_parts` applied to the
*current* folder and the *current* filename, which reconciles the two without
consulting the vendor at all.  The vendor index is needed for exactly one thing
the trimmed upload cannot supply: a missing sublibrary *type* token -- plus, as a
last resort, a filename the folder cannot explain either way.

Two names, not one
------------------
The sublibrary name that goes in the *filename* and the folder segment that goes
in the *path* are tracked separately, because they are allowed to differ: the
folder may elide the redundant ``{ExperimentID}_`` prefix the filename carries
(``P03/…-CHEM16_P03_A1_…``), and both spellings are SOP-clean.  So the proposal
keeps the folder exactly as delivered and rebuilds the filename around the name
the filename itself gives.  Collapsing the two into one variable, as this did, is
what made an accepted upload propose a move of every object into
``{ExperimentID}_{folder}/`` -- or, had the folder been preferred instead,
rewrite every filename.  They diverge only when the vendor has to arbitrate a
``sublibrary_mismatch``, where its name becomes both.

Safety properties
-----------------
* **Self-validating.** Every proposed key is fed back through
  :func:`qa_seahub_sop.validate_seahub_key`; if the proposal is not itself
  SOP-clean it is withheld and reported as unresolved. A rename plan that
  produces another non-compliant name is worse than no plan.
* **Idempotent.** A compliant object proposes itself and contributes no row.
* **Collision-safe.** Two objects can never be handed the same destination, and a
  destination that already exists is blocked rather than overwritten.
* **Advisory only.** Nothing here moves, deletes or rewrites an object. The
  ``.cram-metadata.json`` sidecar carries the old name in its own ``filename``
  field; a rename leaves that stale, and QA does not edit object contents.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Mapping

from qa_constants import (
    SEAHUB_BARE_TO_TRIM_SUFFIX,
    SEAHUB_DOUBLED_WAFER_RE,
    SEAHUB_RENAME_FIELD_SEP,
    SEAHUB_WELL_RE,
    SEAHUB_WELL_VERDICTS,
    raw_expected,
)
from qa_mods import parse_seahub_raw_path, seahub_stem_and_family
from qa_seahub_sop import (
    is_non_sequencing_artifact,
    normalize_doubled_wafer,
    seahub_group_parts,
    unrecognized_suffix,
    validate_seahub_key,
)
from qa_seahub_source import parse_seahub_stem_fields

__all__ = [
    "RenameMapping",
    "RenameProposal",
    "WellRollup",
    "build_rename_mapping",
    "expected_trimmed_key",
    "roll_up_wells",
    "rollup_summary",
]

RENAME_COLUMNS = (
    "current_s3_uri",
    "expected_s3_uri",
    "defects",
    "name_source",
    "status",
    "unresolved",
    "wafer",
    "ug",
    "sublibrary",
    "well",
    "artifact",
)

WELL_COLUMNS = (
    "verdict",
    "wafer",
    "ug",
    "sublibrary",
    "well",
    "objects",
    "missing_artifacts",
    "defects",
    "name_source",
    "detail",
)


@dataclass(frozen=True)
class RenameProposal:
    """The corrected location for one S3 object, or why there isn't one."""

    current_s3_uri: str
    expected_s3_uri: str = ""
    defects: tuple[str, ...] = ()
    unresolved: tuple[str, ...] = ()
    name_source: str = "inferred"
    wafer: str = ""
    ug: str = ""
    sublibrary: str = ""
    well: str = ""
    artifact: str = ""

    @property
    def compliant(self) -> bool:
        return not self.defects and not self.unresolved

    @property
    def renameable(self) -> bool:
        return bool(self.expected_s3_uri) and not self.unresolved and bool(self.defects)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _vendor_group_sans_wafer(group: str, wafer: str) -> str:
    """Drop a leading wafer token a vendor group still carries.

    ``index_untrimmed_sources`` does not normalize a doubled vendor wafer, so a
    stem like ``438514-438514-REF3_P07_1_A3-...`` leaves the second wafer inside
    the group. Both readers of a vendor group have to strip it, and for a while
    only one did: the gap row reached the CSV as sublibrary
    ``438514-REF3_P07_1``, and the ``sublibrary_mismatch`` repair proposed a key
    that re-validated as ``duplicated_wafer_token`` and was therefore withheld --
    turning a well that could have been renamed into ``UNKNOWN``. The self-check
    meant no bad key was ever proposed, so this cost a repair rather than
    producing a wrong one.
    """
    return group[len(wafer) + 1 :] if wafer and group.startswith(f"{wafer}-") else group


def _vendor_sublibrary(group: str) -> str:
    """Strip a trailing well token off a vendor filename group.

    ``REF3_P07_1_A3`` gives ``REF3_P07_1``. The vendor stem is the only place the
    sublibrary appears on that side, since the path does not carry it.
    """
    head, _sep, tail = group.rpartition("_")
    return head if head and SEAHUB_WELL_RE.match(tail) else group


def _vendor_sublibrary_and_well(group: str) -> tuple[str, str]:
    """Split a vendor filename group into ``(sublibrary, well)``.

    ``REF3_P07_1_A3`` gives ``("REF3_P07_1", "A3")``; a group with no trailing
    well token gives ``(group, "")``. The vendor stem is the only place either
    appears on that side, since the vendor path carries neither.

    No second well-token check on the way out: :func:`_vendor_sublibrary` only
    returns a shorter string when the tail already matched, so the trailing
    segment is a well exactly when there is one.
    """
    sublibrary = _vendor_sublibrary(group)
    return sublibrary, group[len(sublibrary) + 1 :] if group != sublibrary else ""


def expected_trimmed_key(
    bucket: str, s3_key: str, source: Any | None = None
) -> RenameProposal:
    """Compose the SOP location for one trimmed object.

    ``source`` is the matching vendor :class:`~qa_seahub_source.SourceEntry`, or
    None. It is consulted for a missing sublibrary type token, and to arbitrate a
    filename neither folder spelling explains.
    """
    current = f"s3://{bucket}/{s3_key}"
    basename = s3_key.split("/")[-1]

    path_info = parse_seahub_raw_path(s3_key)
    if path_info is None:
        return RenameProposal(current, unresolved=("bad_path_depth",))

    parsed = seahub_stem_and_family(basename)
    if parsed is None:
        if is_non_sequencing_artifact(basename):
            return RenameProposal(
                current, defects=("non_sequencing_artifact",), artifact=basename
            )
        return RenameProposal(current, unresolved=("unexpected_suffix",))

    raw_stem, suffix, family = parsed
    trim_suffix = (
        suffix if family == "trim" else SEAHUB_BARE_TO_TRIM_SUFFIX.get(suffix, "")
    )
    if not trim_suffix:
        return RenameProposal(current, unresolved=("unknown_artifact_suffix",))

    defects: list[str] = []
    if family == "bare":
        defects.append("missing_trim_infix")

    doubled = SEAHUB_DOUBLED_WAFER_RE.match(raw_stem)
    if doubled is not None and doubled.group("first") != doubled.group("second"):
        # Two different wafers is not a repair QA may guess at.
        return RenameProposal(
            current, unresolved=("conflicting_wafer_tokens",), artifact=trim_suffix
        )
    stem, was_doubled = normalize_doubled_wafer(raw_stem)
    if was_doubled:
        defects.append("duplicated_wafer_token")

    fields = parse_seahub_stem_fields(stem)
    if fields is None:
        return RenameProposal(
            current, unresolved=("unparseable_stem",), artifact=trim_suffix
        )

    wafer = str(fields["wafer"])
    ug = str(fields["ug"])
    barcode = str(fields["barcode"])
    group = str(fields["group"])
    if wafer != path_info["wafer"]:
        return RenameProposal(
            current,
            unresolved=("wafer_mismatch",),
            wafer=wafer,
            ug=ug,
            artifact=trim_suffix,
        )

    sublibrary, well, matched = seahub_group_parts(
        group, path_info["sublibrary"], path_info["experiment_id"]
    )
    # Kept as delivered: both spellings of the folder are clean, so a proposal
    # that moved the object between them would be a move for no reason.
    folder = path_info["sublibrary"]
    if not matched:
        # The trimmed side is self-contradictory; the vendor stem is the only
        # remaining authority for the sublibrary name.
        vendor_group = getattr(source, "group", "") if source is not None else ""
        if not vendor_group:
            return RenameProposal(
                current,
                unresolved=("sublibrary_mismatch",),
                wafer=wafer,
                ug=ug,
                artifact=trim_suffix,
            )
        sublibrary, well = _vendor_sublibrary_and_well(
            _vendor_group_sans_wafer(vendor_group, getattr(source, "wafer", ""))
        )
        folder = sublibrary
        defects.append("sublibrary_mismatch")

    if well and not SEAHUB_WELL_RE.match(well):
        return RenameProposal(
            current,
            unresolved=("bad_well",),
            wafer=wafer,
            ug=ug,
            sublibrary=sublibrary,
            artifact=trim_suffix,
        )

    assay = fields["assay"]
    vendor_assay = getattr(source, "assay", None) if source is not None else None
    name_source = "inferred"
    if assay is None:
        if not vendor_assay:
            return RenameProposal(
                current,
                unresolved=("invalid_sublibrary_type",),
                wafer=wafer,
                ug=ug,
                sublibrary=sublibrary,
                well=well,
                artifact=trim_suffix,
            )
        assay = vendor_assay
        name_source = "vendor"
        defects.append("invalid_sublibrary_type")
    elif vendor_assay and vendor_assay != assay:
        return RenameProposal(
            current,
            unresolved=("conflicting_sublibrary_type",),
            wafer=wafer,
            ug=ug,
            sublibrary=sublibrary,
            well=well,
            artifact=trim_suffix,
        )

    well_part = f"_{well}" if well else ""
    # `sublibrary` in the name, `folder` in the path: see "Two names, not one".
    expected_stem = f"{wafer}-{sublibrary}{well_part}_{assay}-{ug}-{barcode}"
    project = s3_key.split("/")[0]
    expected_key = (
        f"{project}/{path_info['experiment_id']}/raw/"
        f"{folder}/{wafer}/{expected_stem}{trim_suffix}"
    )

    # A proposal that is not itself SOP-clean is worse than no proposal -- but
    # only for facts the move could change. `bad_bucket` and
    # `lab_project_mismatch` are upload-scope: the proposal keeps the same bucket
    # and the same project, so they hold before and after, and withholding on
    # them rejected proposals that strictly reduce the violation set. One wrong
    # bucket then blacked out every rename in the upload -- measured on the
    # shared fixture, 10 moveable objects to 0 and the wells from
    # {COMPLIANT:1, RENAMEABLE:2, DATA_GAP:1, UNKNOWN:2} to {DATA_GAP:1, UNKNOWN:5}.
    # Filtered here rather than in validate_seahub_key, whose own behaviour is
    # correct and is what the SOP table reports.
    residual = [
        v for v in validate_seahub_key(bucket, expected_key) if v.scope != "upload"
    ]
    if residual:
        return RenameProposal(
            current,
            unresolved=tuple(f"proposal_not_sop:{v.type}" for v in residual),
            wafer=wafer,
            ug=ug,
            sublibrary=sublibrary,
            well=well,
            artifact=trim_suffix,
        )

    return RenameProposal(
        current_s3_uri=current,
        expected_s3_uri=f"s3://{bucket}/{expected_key}",
        defects=tuple(sorted(defects)),
        name_source=name_source,
        wafer=wafer,
        ug=ug,
        sublibrary=sublibrary,
        well=well,
        artifact=trim_suffix,
    )


@dataclass
class RenameMapping:
    """One row per object that needs moving, plus what could not be resolved."""

    rows: list[dict[str, Any]] = field(default_factory=list)
    proposals: list[RenameProposal] = field(default_factory=list)
    compliant_objects: int = 0
    collisions: list[dict[str, Any]] = field(default_factory=list)
    total_objects: int = 0

    @property
    def counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for row in self.rows:
            counts[row["status"]] = counts.get(row["status"], 0) + 1
        return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))

    def moveable(self) -> list[dict[str, Any]]:
        """Rows that are safe to act on."""
        return [r for r in self.rows if r["status"] == "rename"]


def _status_for(proposal: RenameProposal) -> str:
    if "non_sequencing_artifact" in proposal.defects:
        return "not_data"
    if proposal.unresolved:
        return "unresolved"
    return "rename"


def build_rename_mapping(
    bucket: str,
    all_raw_files: list[str],
    untrimmed_index: Mapping[tuple[str, str], Any] | None = None,
) -> RenameMapping:
    """Propose a corrected location for every object that needs one.

    One row per *object*: an S3 move is an object operation, the five artifacts of
    a well map to five different suffixes, and a well may be missing one.
    Compliant objects are counted but contribute no row, so a compliant upload
    yields an empty mapping.
    """
    index = untrimmed_index or {}
    proposals: list[RenameProposal] = []
    for key in sorted(all_raw_files):
        parsed = seahub_stem_and_family(key.split("/")[-1])
        source = None
        if parsed is not None:
            stem, _suffix, _family = parsed
            normalized, _doubled = normalize_doubled_wafer(stem)
            fields = parse_seahub_stem_fields(normalized)
            if fields is not None:
                source = index.get((str(fields["wafer"]), str(fields["ug"])))
        proposals.append(expected_trimmed_key(bucket, key, source))

    mapping = RenameMapping(total_objects=len(all_raw_files), proposals=proposals)

    rows: list[dict[str, Any]] = []
    for proposal in proposals:
        if proposal.compliant:
            mapping.compliant_objects += 1
            continue
        row = {
            "current_s3_uri": proposal.current_s3_uri,
            "expected_s3_uri": proposal.expected_s3_uri,
            "defects": SEAHUB_RENAME_FIELD_SEP.join(proposal.defects),
            "name_source": proposal.name_source,
            "status": _status_for(proposal),
            "unresolved": SEAHUB_RENAME_FIELD_SEP.join(proposal.unresolved),
            "wafer": proposal.wafer,
            "ug": proposal.ug,
            "sublibrary": proposal.sublibrary,
            "well": proposal.well,
            "artifact": proposal.artifact,
        }
        rows.append(row)

    # Collision safety. Blocked rows keep expected_s3_uri for diagnosis but must
    # never be moved.
    existing = {f"s3://{bucket}/{key}" for key in all_raw_files}
    claimants: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if row["status"] == "rename" and row["expected_s3_uri"]:
            claimants.setdefault(row["expected_s3_uri"], []).append(row)

    for destination, contenders in sorted(claimants.items()):
        if len(contenders) > 1:
            for row in contenders:
                row["status"] = "blocked"
            mapping.collisions.append(
                {
                    "kind": "many_to_one",
                    "expected_s3_uri": destination,
                    "detail": (
                        f"{len(contenders)} objects would move to the same key: "
                        + ", ".join(r["current_s3_uri"] for r in contenders)
                    ),
                }
            )
            continue
        row = contenders[0]
        if destination in existing and destination != row["current_s3_uri"]:
            row["status"] = "blocked"
            mapping.collisions.append(
                {
                    "kind": "destination_exists",
                    "expected_s3_uri": destination,
                    "detail": (
                        f"{row['current_s3_uri']} would overwrite an object that "
                        "already exists"
                    ),
                }
            )

    # The rows + compliant == total invariant is asserted by
    # test_every_object_is_accounted_for, not here: an assert in a
    # notebook library can only turn a reported finding into a dead kernel.
    mapping.rows = sorted(rows, key=lambda r: r["current_s3_uri"])
    return mapping


@dataclass
class WellRollup:
    """One verdict per well, which is the headline the notebook prints."""

    rows: list[dict[str, Any]]
    unaccounted: int = 0


def _delivered_artifacts(suffixes: set[str]) -> set[str]:
    """The artifact *kinds* a well delivered, whatever they were named.

    Each suffix maps to its SOP counterpart, so ``.cram`` and ``.trim.cram``
    count as the same kind. Completeness then asks only whether every kind
    arrived; the spelling is a separate axis, reported by the SOP validator as
    ``missing_trim_infix``. Judging by family instead let one optional sidecar
    carrying the other family's name decide the whole requirement set, which is
    how this and :func:`qa_checks.check_expected_raw_files` came to disagree.
    """
    return {SEAHUB_BARE_TO_TRIM_SUFFIX.get(s, s) for s in suffixes}


# The data-bearing artifact, under either spelling. An empty .stdout or
# .stderr is ordinary -- a step that printed nothing -- so emptiness is only a
# finding for this one.
SEAHUB_CRAM_SUFFIXES = (".cram", ".trim.cram")


def roll_up_wells(
    bucket: str,
    all_raw_files: list[str],
    source_index: Mapping[tuple[str, str], Any] | None = None,
    sizes: Mapping[str, int] | None = None,
) -> WellRollup:
    """Reduce an upload to one verdict per well.

    Precedence, first match wins:

    1. ``DATA_GAP`` -- the vendor delivered this well and nothing was uploaded
       for it. The largest gap there is, so it is named as one: as ``UNKNOWN`` it
       was the only kind of gap the notebook did not write to ``errors.txt``,
       which writes ``DATA_GAP`` rows alone. Nothing was *unidentifiable* about
       it either -- the vendor names it exactly.
    2. ``UNKNOWN`` -- the well cannot be *identified* (unparseable stem, or one
       identity claimed by two different names).
    3. ``DATA_GAP`` -- identified, but a required artifact of the delivered family
       is absent. This outranks the un-nameable form of UNKNOWN so that a
       genuinely missing CRAM stays visible even when a vendor order is missing
       from the source list.
    4. ``UNKNOWN`` -- identified and complete, but no corrected key is derivable.
    5. ``RENAMEABLE`` -- every object is either already clean or carries a
       corrected key.
    6. ``UNKNOWN`` -- defects a rename cannot repair. Unreachable as the code
       stands, since every proposal that reaches here has a corrected key, and
       kept only so that a future defect carrying no proposal cannot fall
       through and read as ``COMPLIANT``.
    7. ``COMPLIANT``.

    Branch 5 asks :attr:`RenameProposal.renameable` rather than testing the
    defect names against :data:`SEAHUB_RENAMEABLE_SOP_TYPES`. The two headline
    CSVs then decide moveability the same way and cannot contradict each other:
    with the membership test, a ``sublibrary_mismatch`` the vendor had already
    repaired produced an object the rename CSV said to move and a well the status
    CSV called ``UNKNOWN``.
    """
    index = source_index or {}
    wells: dict[tuple[str, str], dict[str, Any]] = {}
    unaccounted = 0

    for key in sorted(all_raw_files):
        basename = key.split("/")[-1]
        parsed = seahub_stem_and_family(basename)
        path_info = parse_seahub_raw_path(key)
        if path_info is None:
            unaccounted += 1
            continue
        if parsed is None:
            # A well path with a misspelled artifact suffix: the stem still
            # parses, but seahub_stem_and_family rejects the name. Without this
            # branch the well vanished from the headline CSV -- no row, not
            # counted in unaccounted -- while the SOP table reported the suffix.
            # Browser junk under a well folder (login.html) also has an
            # unrecognized suffix but no parseable stem; skip those here.
            suffix = unrecognized_suffix(basename)
            if not suffix:
                continue
            raw_stem = basename[: -len(suffix)]
            family = ""
        else:
            raw_stem, suffix, family = parsed
        normalized, _doubled = normalize_doubled_wafer(raw_stem)
        # Fall back to the raw stem when normalizing makes it unparseable, so
        # this keys wells the same way index_trimmed_upload does. Without it the
        # two disagreed on exactly one shape -- a doubled wafer with no group,
        # 123456-123456-Z0001-ACGT, which parses only *before* normalization --
        # and a well that was uploaded got a phantom "nothing was uploaded" row
        # beside its real one.
        fields = parse_seahub_stem_fields(normalized) or parse_seahub_stem_fields(
            raw_stem
        )
        if parsed is None and fields is None:
            # Browser junk under a well folder (login.html) has an unrecognized
            # suffix and no parseable stem; only the misspelled-artifact case
            # above belongs in the roll-up.
            continue
        identity = (
            (str(fields["wafer"]), str(fields["ug"]))
            if fields is not None
            else ("", f"{path_info['sublibrary']}/{raw_stem}")
        )
        well = wells.setdefault(
            identity,
            {
                "identity": identity,
                "names": set(),
                "suffixes": set(),
                "sublibrary": path_info["sublibrary"],
                "keys": [],
                "parseable": fields is not None,
                "empty_cram": False,
            },
        )
        well["names"].add(("/".join(key.split("/")[:-1]), raw_stem))
        well["suffixes"].add(suffix)
        well["keys"].append(key)
        # A CRAM that is present and empty carries no data, so it is not a
        # delivered artifact. Only a *known* zero counts: a key absent from the
        # mapping means the size was never collected (manifest mode collects
        # none), and treating unknown as empty would fail every well.
        if suffix in SEAHUB_CRAM_SUFFIXES and (sizes or {}).get(key) == 0:
            well["empty_cram"] = True

    for identity, entry in index.items():
        wells.setdefault(
            identity,
            {
                "identity": identity,
                "names": set(),
                "suffixes": set(),
                "sublibrary": "",
                "keys": [],
                "parseable": True,
                "vendor_only": True,
            },
        )

    rows: list[dict[str, Any]] = []
    for identity in sorted(wells, key=lambda k: (k[0], k[1])):
        well = wells[identity]
        wafer, ug = identity
        source = index.get(identity)
        proposals = [expected_trimmed_key(bucket, k, source) for k in well["keys"]]
        defects = sorted({d for p in proposals for d in p.defects})
        unresolved = sorted({u for p in proposals for u in p.unresolved})
        name_source = (
            "vendor"
            if any(p.name_source == "vendor" for p in proposals)
            else "inferred"
        )

        # Unconditional: a vendor well with nothing uploaded is missing all five,
        # and zeroing the requirement made missing_artifacts print blank on the
        # one row where the whole set is absent.
        required = tuple(raw_expected["seahub_sci"])
        delivered = _delivered_artifacts(well["suffixes"])
        missing = [s for s in required if s not in delivered]
        # An uploaded-but-empty CRAM is the shape an interrupted copy leaves
        # behind, and it used to read as COMPLIANT: the completeness check asks
        # only whether the artifact arrived, and the trimmed-vs-vendor size
        # check skips a zero because it cannot tell "empty" from "size not
        # collected". Reported as a gap in the data, which is what it is, and
        # without needing a vendor delivery to compare against.
        empty_cram = bool(well.get("empty_cram"))

        sublibrary = well["sublibrary"]
        well_id = next((p.well for p in proposals if p.well), "")

        if not well["keys"]:
            verdict = "DATA_GAP"
            detail = "delivered by the vendor but nothing was uploaded for this well"
            if source is not None:
                detail += f"; the vendor delivered it as {source.cram_key}"
                # Both come off the uploaded objects everywhere else, and this is
                # the row that has none. Leaving them blank made the highest
                # severity verdict the least locatable one: the notebook printed
                # "wafer 438514 Z0305 ():" and the CSV carried empty cells for
                # the well an operator most needs to find. The vendor stem names
                # both, which is the same authority the mismatch repair uses.
                #
                # This is the vendor's SOP name for the sublibrary, which is not
                # always the folder the upload used -- every real upload elides
                # the ExperimentID prefix, which is an accepted spelling, so an
                # uploaded row says "P07_1" where this one says "REF3_P07_1" for
                # the same sublibrary. That difference is now permanent rather
                # than something a rename resolves. Documented in SEAHUB_QA.md and
                # pinned by test_a_gap_row_names_the_sublibrary_the_vendor_used.
                group = source.group
                sublibrary, well_id = _vendor_sublibrary_and_well(
                    _vendor_group_sans_wafer(group, source.wafer)
                )
        elif not well["parseable"] or len(well["names"]) > 1:
            verdict = "UNKNOWN"
            detail = (
                f"one well appears under {len(well['names'])} different names"
                if len(well["names"]) > 1
                else "stem does not decompose into the SOP fields"
            )
        elif missing or empty_cram:
            verdict = "DATA_GAP"
            if empty_cram:
                detail = "the uploaded CRAM is 0 bytes, so this well carries no data"
                if missing:
                    detail += f"; also absent: {', '.join(missing)}"
            else:
                detail = f"absent from the delivered set: {', '.join(missing)}"
            if source is not None:
                detail += f"; the vendor delivered it as {source.cram_key}"
        elif unresolved:
            verdict = "UNKNOWN"
            detail = f"no corrected name derivable: {', '.join(unresolved)}"
        elif defects and all(p.compliant or p.renameable for p in proposals):
            verdict = "RENAMEABLE"
            detail = "complete; every defect is repairable by renaming"
        elif defects:
            verdict = "UNKNOWN"
            detail = f"defects a rename cannot repair: {', '.join(defects)}"
        else:
            verdict = "COMPLIANT"
            detail = ""

        rows.append(
            {
                "verdict": verdict,
                "wafer": wafer,
                "ug": ug,
                "sublibrary": sublibrary,
                "well": well_id,
                "objects": len(well["keys"]),
                "missing_artifacts": SEAHUB_RENAME_FIELD_SEP.join(missing),
                "defects": SEAHUB_RENAME_FIELD_SEP.join(defects),
                "name_source": name_source,
                "detail": detail,
            }
        )
    return WellRollup(rows=rows, unaccounted=unaccounted)


def rollup_summary(rows: list[dict[str, Any]]) -> dict[str, int]:
    """Count verdicts in the documented display order."""
    counts = {verdict: 0 for verdict in SEAHUB_WELL_VERDICTS}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    return counts
