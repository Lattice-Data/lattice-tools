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

Instead the folder decision comes from :func:`qa_seahub_sop.seahub_group_parts`
applied to the *current* folder and the *current* filename: when the filename
says ``{ExperimentID}_{folder}``, the folder is the truncated side and the
filename already carries the authoritative name.  That covers every REF3 defect
without consulting the vendor at all.  The vendor index is needed for exactly one
thing the trimmed upload cannot supply: a missing sublibrary *type* token.

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
    "source_sublibrary_segment",
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


def source_sublibrary_segment(cram_key: str) -> str:
    """The vendor path segment before ``raw``.

    Returns the ExperimentID for the real vendor layout, and the full sublibrary
    only for the older ``{ExperimentID}_{sublibrary}`` shape. Used purely as a
    last-resort corroborator when the trimmed folder and filename cannot be
    reconciled with each other -- never as the primary source.
    """
    parts = cram_key.split("/")
    for i in range(len(parts) - 1, 1, -1):
        if parts[i] == "raw" and len(parts) - i == 3:
            return parts[i - 1]
    return ""


def _vendor_sublibrary(group: str) -> str:
    """Strip a trailing well token off a vendor filename group.

    ``REF3_P07_1_A3`` gives ``REF3_P07_1``. The vendor stem is the only place the
    sublibrary appears on that side, since the path does not carry it.
    """
    head, _sep, tail = group.rpartition("_")
    return head if head and SEAHUB_WELL_RE.match(tail) else group


def expected_trimmed_key(
    bucket: str, s3_key: str, source: Any | None = None
) -> RenameProposal:
    """Compose the SOP location for one trimmed object.

    ``source`` is the matching vendor :class:`~qa_seahub_source.SourceEntry`, or
    None. It is consulted only for a missing sublibrary type token.
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

    folder, well, state = seahub_group_parts(
        group, path_info["sublibrary"], path_info["experiment_id"]
    )
    if state == "truncated":
        defects.append("sublibrary_folder_truncated")
    elif state == "mismatch":
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
        folder = _vendor_sublibrary(vendor_group)
        well = ""
        vendor_well = vendor_group[len(folder) + 1 :] if vendor_group != folder else ""
        if SEAHUB_WELL_RE.match(vendor_well):
            well = vendor_well
        defects.append("sublibrary_mismatch")

    if well and not SEAHUB_WELL_RE.match(well):
        return RenameProposal(
            current,
            unresolved=("bad_well",),
            wafer=wafer,
            ug=ug,
            sublibrary=folder,
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
                sublibrary=folder,
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
            sublibrary=folder,
            well=well,
            artifact=trim_suffix,
        )

    well_part = f"_{well}" if well else ""
    expected_stem = f"{wafer}-{folder}{well_part}_{assay}-{ug}-{barcode}"
    project = s3_key.split("/")[0]
    expected_key = (
        f"{project}/{path_info['experiment_id']}/raw/"
        f"{folder}/{wafer}/{expected_stem}{trim_suffix}"
    )

    # A proposal that is not itself SOP-clean is worse than no proposal.
    residual = validate_seahub_key(bucket, expected_key)
    if residual:
        return RenameProposal(
            current,
            unresolved=tuple(f"proposal_not_sop:{v.type}" for v in residual),
            wafer=wafer,
            ug=ug,
            sublibrary=folder,
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
        sublibrary=folder,
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

    mapping.rows = sorted(rows, key=lambda r: r["current_s3_uri"])
    assert len(mapping.rows) + mapping.compliant_objects == mapping.total_objects
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


def roll_up_wells(
    bucket: str,
    all_raw_files: list[str],
    source_index: Mapping[tuple[str, str], Any] | None = None,
) -> WellRollup:
    """Reduce an upload to one verdict per well.

    Precedence, first match wins:

    1. ``UNKNOWN`` -- the well cannot be *identified* (unparseable stem, one
       identity claimed by two different names, or a vendor well with nothing
       delivered at all).
    2. ``DATA_GAP`` -- identified, but a required artifact of the delivered family
       is absent. This outranks the un-nameable form of UNKNOWN so that a
       genuinely missing CRAM stays visible even when a vendor order is missing
       from the source list.
    3. ``UNKNOWN`` -- identified and complete, but no corrected key is derivable.
    4. ``RENAMEABLE`` -- every object is either already clean or carries a
       corrected key.
    5. ``UNKNOWN`` -- defects a rename cannot repair. Unreachable as the code
       stands, since every proposal that reaches here has a corrected key, and
       kept only so that a future defect carrying no proposal cannot fall
       through and read as ``COMPLIANT``.
    6. ``COMPLIANT``.

    Branch 4 asks :attr:`RenameProposal.renameable` rather than testing the
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
            # Junk has no well; it is reported by the SOP table and the mapping.
            continue
        raw_stem, suffix, family = parsed
        normalized, _doubled = normalize_doubled_wafer(raw_stem)
        fields = parse_seahub_stem_fields(normalized)
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
                "families": set(),
                "sublibrary": path_info["sublibrary"],
                "keys": [],
                "parseable": fields is not None,
            },
        )
        well["names"].add(("/".join(key.split("/")[:-1]), raw_stem))
        well["suffixes"].add(suffix)
        well["families"].add(family)
        well["keys"].append(key)

    for identity, entry in index.items():
        wells.setdefault(
            identity,
            {
                "identity": identity,
                "names": set(),
                "suffixes": set(),
                "families": set(),
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

        required = tuple(raw_expected["seahub_sci"]) if well["keys"] else ()
        delivered = _delivered_artifacts(well["suffixes"])
        missing = [s for s in required if s not in delivered]

        if not well["keys"]:
            verdict = "UNKNOWN"
            detail = "delivered by the vendor but nothing was uploaded for this well"
        elif not well["parseable"] or len(well["names"]) > 1:
            verdict = "UNKNOWN"
            detail = (
                f"one well appears under {len(well['names'])} different names"
                if len(well["names"]) > 1
                else "stem does not decompose into the SOP fields"
            )
        elif missing:
            verdict = "DATA_GAP"
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
                "sublibrary": well["sublibrary"],
                "well": next((p.well for p in proposals if p.well), ""),
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
