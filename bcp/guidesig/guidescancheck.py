"""Report how many targeting guides in a sheet have an exact genome match.

Runs ``bcp/guidescan_pipeline.py``'s ``GuidescanPipeline`` with its own fixed
settings against the species' guidescan index, downloaded from S3 into
``bcp/.guidescan_cache/<species>/`` on first use. The exact-match rate is
reported, with a warning below ``WARN_BELOW_PERCENT``; it never fails the run.

Before running, the sheet's Ensembl gene IDs are cross-checked against
``--species``, since a run against the wrong genome would report a misleadingly
low rate. That cross-check does fail the run.

This module is imported only for ``--guidescan`` runs.
"""

from __future__ import annotations

import csv
import logging
import re
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .core import (
    CONTROL_ROLE,
    DEFAULT_COLUMN,
    FILE_FORMATS,
    ID_COLUMN,
    ROLE_COLUMN,
    TARGET_COLUMN,
    GuideSigError,
)

PAM_COLUMN = "guide_PAM"
SPCAS9_PAM = "NGG"
WARN_BELOW_PERCENT = 90.0
EXAMPLE_LIMIT = 10

REFERENCE_BUCKET = "submissions-lattice"
REFERENCE_PREFIX = "BCP/references/guidescan_references/"
# FASTA basename per species; the index files are named after it.
REFERENCE_FASTA = {
    "mouse": "Mus_musculus.GRCm39.dna.primary_assembly.fa",
    "human": "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
}
# Files downloaded per species, as suffixes of the FASTA basename. The FASTA
# itself is not downloaded: guidescan_pipeline.py only requires the index.
REFERENCE_SUFFIXES = (
    ".index.forward",
    ".index.reverse",
    ".index.gs",
    ".forward.dna",
    ".reverse.dna",
)
CACHE_DIR = Path(__file__).resolve().parent.parent / ".guidescan_cache"

# Ensembl gene ID prefixes, optionally versioned (ENSG00000141510.18). Anything
# else, including placeholders such as ENSGFAKE..., says nothing about species.
_SPECIES_ID_PATTERNS = {
    "human": re.compile(r"ENSG\d+(\.\d+)?"),
    "mouse": re.compile(r"ENSMUSG\d+(\.\d+)?"),
}


class GuidescanCheckError(GuideSigError):
    """Raised when the guidescan check cannot run or the species disagrees."""


@dataclass(frozen=True)
class TargetingGuide:
    guide_id: str
    protospacer: str
    pam: str


@dataclass
class GuidescanReport:
    total: int = 0
    matched: int = 0
    unmatched: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def percent(self) -> float:
        return 100.0 * self.matched / self.total if self.total else 0.0

    @property
    def below_threshold(self) -> bool:
        return self.total > 0 and self.percent < WARN_BELOW_PERCENT


# --- reading the sheet --------------------------------------------------------


def _data_rows(path: Path, file_format: str):
    """Yield (physical_row, cells-by-header) for the sheet's data rows.

    Only called after the core checks have passed, so the file is known to
    parse; blank and comment rows are skipped by the same rules as core.
    """
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=FILE_FORMATS[file_format])
        header = next(reader)
        record_line = reader.line_num + 1
        for row in reader:
            physical_row = record_line
            record_line = reader.line_num + 1
            if not row or all(cell == "" for cell in row):
                continue
            if row[0].lstrip().startswith("#"):
                continue
            yield physical_row, dict(zip(header, row))


def parse_pam(raw: str) -> str:
    """Turn a sheet PAM such as ``3' NGG`` into the bare ``NGG`` guidescan takes.

    Any length is accepted, as guidescan_pipeline.py reads the length from the
    PAM itself. A 5' PAM is rejected: the pipeline never passes guidescan's
    ``--start``, so it only searches for a PAM after the protospacer.
    """
    text = raw.strip().replace("’", "'")
    match = re.fullmatch(r"(?:(3|5)'\s*)?([ACGTNRYKMSWBDHV]+)", text, re.IGNORECASE)
    if not match:
        raise ValueError(f"unrecognized PAM {raw!r}")
    side, pam = match.group(1), match.group(2).upper()
    if side == "5":
        raise ValueError(f"5' PAM {raw!r} is not supported by guidescan_pipeline.py")
    return pam


def pipeline_pam(pam: str) -> str:
    """Map a parsed sheet PAM to the pattern guidescan is run with.

    Sheets record either the pattern (``NGG``) or the actual bases at each
    guide's site (``AGG``, ``TGG``, ...). Every SpCas9 site is an instance of
    NGG, so any PAM fitting NGG becomes NGG; the protospacer must still match
    exactly, and only the sheet's specific first PAM base goes unchecked. Any
    other PAM is passed through unchanged; ``read_targeting_guides`` then
    requires one pattern across the sheet, as guidescan_pipeline.py does.
    """
    if re.fullmatch(r"[ACGTN]GG", pam):
        return SPCAS9_PAM
    return pam


def read_targeting_guides(
    path: str | Path, *, file_format: str, column: str = DEFAULT_COLUMN
) -> tuple[list[TargetingGuide], set[str]]:
    """Return the sheet's targeting guides and all its guide_target_gene_id values.

    A targeting guide is any row whose role is not ``control`` and that has a
    protospacer. Its PAM comes from ``guide_PAM``, mapped by ``pipeline_pam``.

    Raises:
        GuidescanCheckError: If the guides do not share one PAM pattern, which
            guidescan_pipeline.py requires, among the other row failures.
    """
    file_path = Path(path)
    guides: dict[str, TargetingGuide] = {}
    target_ids: set[str] = set()
    for physical_row, cells in _data_rows(file_path, file_format):
        target = cells.get(TARGET_COLUMN, "").strip()
        if target:
            target_ids.add(target)
        if cells.get(ROLE_COLUMN, "").strip().lower() == CONTROL_ROLE:
            continue
        protospacer = cells.get(column, "").strip().upper()
        if not protospacer:
            continue
        where = f"{file_path}: row {physical_row}"
        if PAM_COLUMN not in cells:
            raise GuidescanCheckError(
                f"{file_path}: required column {PAM_COLUMN!r} absent from header"
            )
        try:
            pam = pipeline_pam(parse_pam(cells[PAM_COLUMN]))
        except ValueError as exc:
            raise GuidescanCheckError(f"{where}: {exc}") from exc
        guide_id = cells.get(ID_COLUMN, "").strip()
        if not guide_id:
            raise GuidescanCheckError(f"{where}: targeting guide has no {ID_COLUMN}")
        guide = TargetingGuide(guide_id, protospacer, pam)
        previous = guides.setdefault(guide_id, guide)
        if previous != guide:
            raise GuidescanCheckError(
                f"{where}: {ID_COLUMN} {guide_id!r} appears with two different "
                "protospacers or PAMs"
            )
    patterns = sorted({guide.pam for guide in guides.values()})
    if len(patterns) > 1:
        raise GuidescanCheckError(
            f"{file_path}: targeting guides use more than one PAM pattern "
            f"({', '.join(patterns)}); guidescan_pipeline.py needs one PAM for "
            "all guides"
        )
    return list(guides.values()), target_ids


def check_species(target_ids: set[str], species: str) -> str | None:
    """Fail if the sheet's Ensembl IDs point to another species.

    Returns a note when no ID says anything about species, else None.
    """
    found = {
        name
        for name, pattern in _SPECIES_ID_PATTERNS.items()
        for gid in target_ids
        if pattern.fullmatch(gid)
    }
    if not found:
        return (
            "species cross-check skipped: no Ensembl human (ENSG) or mouse "
            "(ENSMUSG) gene IDs in guide_target_gene_id"
        )
    if len(found) > 1:
        raise GuidescanCheckError(
            "guide_target_gene_id mixes human (ENSG) and mouse (ENSMUSG) gene IDs"
        )
    (detected,) = found
    if detected != species:
        raise GuidescanCheckError(
            f"--species {species} but guide_target_gene_id holds {detected} "
            "Ensembl gene IDs"
        )
    return None


# --- references ---------------------------------------------------------------


def _make_s3_client() -> Any:
    try:
        import boto3
    except ImportError as exc:
        raise GuidescanCheckError(
            "--guidescan needs boto3 to download references; install it from "
            "requirements.txt"
        ) from exc
    return boto3.client("s3")


def ensure_references(
    species: str,
    *,
    cache_dir: Path = CACHE_DIR,
    s3_client: Any = None,
    log: Callable[[str], None] = lambda message: None,
) -> Path:
    """Download any missing or size-mismatched index file; return the index prefix.

    Each file is downloaded to a temporary name and renamed when complete, so
    an interrupted download is never mistaken for a cached file.
    """
    fasta = REFERENCE_FASTA[species]
    target_dir = cache_dir / species
    target_dir.mkdir(parents=True, exist_ok=True)
    client = s3_client if s3_client is not None else _make_s3_client()

    for suffix in REFERENCE_SUFFIXES:
        name = fasta + suffix
        key = REFERENCE_PREFIX + name
        local = target_dir / name
        try:
            size = client.head_object(Bucket=REFERENCE_BUCKET, Key=key)["ContentLength"]
        except Exception as exc:
            raise GuidescanCheckError(
                f"could not find reference s3://{REFERENCE_BUCKET}/{key}: {exc}"
            ) from exc
        if local.exists() and local.stat().st_size == size:
            continue
        log(f"downloading {name} ({size / 1e9:.2f} GB) to {target_dir}")
        partial = local.with_name(local.name + ".partial")
        try:
            client.download_file(REFERENCE_BUCKET, key, str(partial))
            partial.replace(local)
        except Exception as exc:
            partial.unlink(missing_ok=True)
            raise GuidescanCheckError(
                f"could not download s3://{REFERENCE_BUCKET}/{key}: {exc}"
            ) from exc

    return target_dir / (fasta + ".index")


# --- running guidescan --------------------------------------------------------


def _run_pipeline(index_prefix: Path, guides_file: Path, output_dir: Path) -> Any:
    """Run GuidescanPipeline with its fixed settings; return its result frame."""
    try:
        from guidescan_pipeline import GuidescanPipeline
    except ImportError as exc:
        raise GuidescanCheckError(
            f"--guidescan could not load guidescan_pipeline (missing: {exc.name}); "
            "run from bcp/ with requirements.txt installed"
        ) from exc

    # CRITICAL, not WARNING: the pipeline logs its own failures at ERROR, and
    # guidesig reports every failure itself as one `error:` line, so letting
    # ERROR through would print each failure twice. guidescan_pipeline still
    # calls logging.basicConfig at import; that side effect is left as is.
    pipeline_logger = logging.getLogger("guidescan_pipeline")
    previous_level = pipeline_logger.level
    pipeline_logger.setLevel(logging.CRITICAL)
    try:
        pipeline = GuidescanPipeline(
            genome_index=str(index_prefix),
            guides_file=str(guides_file),
            output_dir=str(output_dir),
        )
        return pipeline.run()
    except Exception as exc:
        raise GuidescanCheckError(f"guidescan pipeline failed: {exc}") from exc
    finally:
        pipeline_logger.setLevel(previous_level)


def run_guidescan_check(
    guides: list[TargetingGuide],
    index_prefix: Path,
    *,
    runner: Callable[[Path, Path, Path], Any] | None = None,
) -> GuidescanReport:
    """Run guidescan on ``guides`` and count those with an exact match.

    Guides are passed to guidescan under synthetic IDs (g000001, ...) and
    mapped back, so sheet IDs that pandas would reinterpret (e.g. ``0033``
    read back as 33) cannot be miscounted.
    """
    if runner is None:
        runner = _run_pipeline
    report = GuidescanReport(total=len(guides))
    if not guides:
        report.notes.append("no targeting guides; guidescan not run")
        return report

    by_token = {f"g{i:06d}": guide for i, guide in enumerate(guides, start=1)}
    with tempfile.TemporaryDirectory(prefix="guidesig_guidescan_") as tmp:
        tmp_path = Path(tmp)
        guides_file = tmp_path / "guides.csv"
        with guides_file.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["id", "sequence", "pam", "chromosome", "position", "sense"]
            )
            for token, guide in by_token.items():
                writer.writerow([token, guide.protospacer, guide.pam, "", "", ""])
        result = runner(index_prefix, guides_file, tmp_path / "out")

    matched_tokens = {
        str(token)
        for token, chromosome in zip(result["id"], result["chromosome"])
        if str(chromosome) != "NA"
    }
    report.matched = sum(token in matched_tokens for token in by_token)
    report.unmatched = sorted(
        guide.guide_id
        for token, guide in by_token.items()
        if token not in matched_tokens
    )
    return report


def stderr_log(message: str) -> None:
    print(f"guidescan: {message}", file=sys.stderr)


def format_report(report: GuidescanReport) -> list[str]:
    """Lines for stderr: the rate, unmatched examples, and any warning."""
    lines = [f"note: {note}" for note in report.notes]
    if report.total == 0:
        return lines
    lines.append(
        f"guidescan: {report.matched} of {report.total} targeting guides "
        f"({report.percent:.1f}%) have an exact match"
    )
    if report.unmatched:
        examples = ", ".join(report.unmatched[:EXAMPLE_LIMIT])
        lines.append(f"guidescan: no exact match: {examples}")
    if report.below_threshold:
        lines.append(
            f"warning: guidescan exact-match rate {report.percent:.1f}% is below "
            f"{WARN_BELOW_PERCENT:.0f}%"
        )
    return lines
