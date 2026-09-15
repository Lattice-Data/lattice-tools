"""Command line for the ChEBI submission gate."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import casreg, chebi_release, checks, client, decisions, external
from . import io as gate_io

log = logging.getLogger(__name__)


EXIT_OK = 0
EXIT_HELD = 1
EXIT_USAGE = 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chebi_gate",
        description=(
            "Split an SDF into records cleared for ChEBI submission and records "
            "held back with written reasons."
        ),
        epilog=(
            "Exit status: 0 if every record cleared, 1 if any was held, 2 on a "
            "usage or input error. Suitable for a build pipeline."
        ),
    )
    parser.add_argument(
        "sdf", nargs="?", help="the SDF you intend to submit (omit with --distil)"
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="where to write the outputs (default: alongside the input)",
    )
    parser.add_argument(
        "--cas-registry",
        default=None,
        metavar="CSV",
        help=(
            "parsed SciFinder export. The only independent source with full "
            "coverage; enables EXT-04 and improves EXT-01."
        ),
    )
    parser.add_argument(
        "--chebi-index",
        default=None,
        metavar="DIR",
        help="distilled ChEBI release index, built by --distil; enables EXT-02/03",
    )
    parser.add_argument(
        "--pubchem-cache",
        default=None,
        metavar="DIR",
        help=(
            "PubChem JSON cache. Circular evidence: the SDFs were generated from "
            "it, so its agreement is reported separately and never counted as "
            "independent."
        ),
    )
    parser.add_argument(
        "--cas-common-chemistry",
        default=None,
        metavar="DIR",
        help="CAS Common Chemistry JSON cache (independent, partial coverage)",
    )
    parser.add_argument(
        "--decisions",
        default=None,
        metavar="DIR",
        help=(
            "directory holding waivers.csv and quarantine.csv "
            f"(default: {decisions.DECISIONS_DIR})"
        ),
    )
    parser.add_argument(
        "--no-decisions",
        action="store_true",
        help="ignore the decisions data entirely (for auditing what it suppresses)",
    )
    parser.add_argument(
        "--role",
        choices=list(checks.ROLES),
        default=checks.ROLE_AUTO,
        help=(
            "treat every record as a salt, as a neutral compound, or decide per "
            "record from its fragments and asserted class (default: auto)"
        ),
    )
    parser.add_argument(
        "--allow-medium",
        action="store_true",
        help="clear records whose worst finding is medium",
    )
    parser.add_argument(
        "--distil",
        nargs=2,
        metavar=("RELEASE_DIR", "INDEX_DIR"),
        default=None,
        help=(
            "instead of running the gate, distil a downloaded ChEBI flat-file "
            f"release into a queryable index, without gating anything. Releases "
            f"live at {chebi_release.FLAT_FILES_URL}"
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    return parser


def _cache_dir(value: str | None, option: str) -> Path | None:
    """A cache directory the caller named, or a usage error saying it is not there.

    One mechanism, not two: this and a separate loop in ``_evidence`` were written
    against the same defect from different directions, and two guards that agree
    cannot be told apart by a test. This one is kept because it says *which* way
    the path is wrong, which is the difference between a typo and a file passed
    where a directory belongs.

    Both JSON caches used to be wrapped in ``Path()`` and handed straight to
    ``Evidence``, which tests ``is_dir()`` and quietly treats a false as "no such
    source configured". A mistyped path therefore ran the whole gate with the
    source missing and reported the reduced coverage as the truth: records cleared,
    exit 0, and a summary reading "with no independent source, clearance means
    internally consistent, not confirmed" -- on a run where the caller had asked
    for one. --cas-registry and --chebi-index already fail loudly on a bad path,
    through casreg.load and load_index. These two were the exception.
    """
    if not value:
        return None
    path = Path(value)
    if path.is_dir():
        return path
    problem = "is not a directory" if path.exists() else "does not exist"
    raise client.GateError(f"{option} {path} {problem}")


def _evidence(args: argparse.Namespace) -> external.Evidence:
    registry = casreg.EMPTY
    if args.cas_registry:
        registry = casreg.load(args.cas_registry)
    index = chebi_release.EMPTY
    if args.chebi_index:
        index = chebi_release.load_index(args.chebi_index)
    return external.Evidence(
        registry=registry,
        chebi=index,
        pubchem_dir=_cache_dir(args.pubchem_cache, "--pubchem-cache"),
        common_chemistry_dir=_cache_dir(
            args.cas_common_chemistry, "--cas-common-chemistry"
        ),
    )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )

    if args.distil:
        release_dir, index_dir = args.distil
        index = chebi_release.distil(release_dir, index_dir)
        print(f"indexed {len(index)} structures into {index_dir}")
        return EXIT_OK

    if not args.sdf:
        print("no SDF given; pass one, or use --distil", file=sys.stderr)
        return EXIT_USAGE
    sdf_path = Path(args.sdf)
    if not sdf_path.exists():
        print(f"no such file: {sdf_path}", file=sys.stderr)
        return EXIT_USAGE

    try:
        evidence = _evidence(args)
        loaded = (
            None
            if args.no_decisions
            else decisions.load(args.decisions or decisions.DECISIONS_DIR)
        )
        gate_run = client.run(
            sdf_path,
            evidence=evidence,
            decisions=loaded,
            role=args.role,
            allow_medium=args.allow_medium,
            repo=Path(__file__).resolve().parent.parent,
        )
    except (
        client.GateError,
        decisions.DecisionsError,
        casreg.RegistryError,
        chebi_release.ReleaseError,
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE

    out_dir = Path(args.out_dir) if args.out_dir else sdf_path.parent
    outputs = gate_io.write(gate_run, out_dir, stem=sdf_path.stem)
    print(gate_io.summary(gate_run, outputs))
    return EXIT_HELD if gate_run.held else EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
