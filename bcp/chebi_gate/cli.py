"""Command line for the ChEBI submission gate."""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime, timezone
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


# Every failure the gate reports as a usage or input error rather than a crash.
# OSError is in the list because the documented contract is about outcomes, not
# about which layer raised: an unwritable --out-dir and a mistyped --cas-registry
# are the same thing to whoever typed the command.
INPUT_ERRORS = (
    client.GateError,
    decisions.DecisionsError,
    casreg.RegistryError,
    chebi_release.ReleaseError,
    OSError,
    # JSONDecodeError and UnicodeDecodeError both subclass it, which is how a
    # corrupt index_manifest.json and a registry CSV saved as cp1252 -- Excel's
    # default on Windows, and that table is assembled by hand -- escaped as a
    # traceback with exit 1. Exit 1 is the status that means "records were held",
    # which is the one confusion a pipeline must not have.
    ValueError,
)


def _distil(args: argparse.Namespace) -> int:
    release_dir, index_dir = args.distil
    index = chebi_release.distil(
        release_dir,
        index_dir,
        # The wall clock, once, so index_manifest.json records when it was built.
        # Nothing passed this but the tests, so every real manifest recorded
        # `"generated": ""` while the docs described that field as the reason the
        # manifest is not byte-identical between builds. It cannot affect a run
        # identifier: `generated` is in manifest.RUN_ID_IGNORED_KEYS.
        generated=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )
    print(f"indexed {len(index)} structures into {index_dir}")
    return EXIT_OK


def _gate(args: argparse.Namespace, sdf_path: Path) -> int:
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
    out_dir = Path(args.out_dir) if args.out_dir else sdf_path.parent
    outputs = gate_io.write(gate_run, out_dir, stem=sdf_path.stem)
    print(gate_io.summary(gate_run, outputs))
    return EXIT_HELD if gate_run.held else EXIT_OK


def main(argv: list[str] | None = None) -> int:
    """Run the gate, or distil a release. Never lets an input error reach a traceback.

    Both branches sit inside the handler. They did not: --distil ran before it, so
    a mistyped release directory came out as a ReleaseError traceback, and
    gate_io.write ran after it, so an unwritable --out-dir did the same -- both
    against a documented contract of "2 on a usage or input error" that the
    parser's own epilog repeats.
    """
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )

    try:
        if args.distil:
            if args.sdf:
                # Silently gating nothing is the one outcome a pipeline cannot
                # notice: --distil returns 0, which is also "every record cleared".
                print(
                    f"--distil does not gate; drop {args.sdf} or drop --distil",
                    file=sys.stderr,
                )
                return EXIT_USAGE
            return _distil(args)
        if not args.sdf:
            print("no SDF given; pass one, or use --distil", file=sys.stderr)
            return EXIT_USAGE
        sdf_path = Path(args.sdf)
        if not sdf_path.exists():
            print(f"no such file: {sdf_path}", file=sys.stderr)
            return EXIT_USAGE
        return _gate(args, sdf_path)
    except INPUT_ERRORS as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE


if __name__ == "__main__":
    sys.exit(main())
