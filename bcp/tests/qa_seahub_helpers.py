"""
Shared SeaHub test fixtures: a listing reproducing every real REF3 defect.

Derived from two real S3 listings reviewed alongside these changes -- a trimmed
upload of 2035 objects and a 2597-object vendor delivery -- with lab and order
identifiers sanitized. The vendor side uses the measured layout, where the
segment before ``raw`` is the ExperimentID alone and carries no sublibrary.

Wells, one per defect combination:

===== ============================================================= ========
well  defect                                                        verdict
===== ============================================================= ========
A     truncated folder + doubled wafer + no type tag + no ``.trim``  RENAMEABLE
B     correct filename, truncated folder only                        RENAMEABLE
C     four sidecars and no CRAM anywhere                             DATA_GAP
D     fully SOP-compliant                                            COMPLIANT
G1    no vendor match, but the type tag is present                   RENAMEABLE
G2    no vendor match and no type tag, so not inferable              UNKNOWN
===== ============================================================= ========
"""

from __future__ import annotations

BUCKET = "czi-labalpha"
PROJECT = "labalpha-seahub-bcp"
EXPERIMENT = "REF3"
RAW = f"{PROJECT}/{EXPERIMENT}/raw"

VENDOR_BUCKET = "czi-novogene"
VENDOR_ORDER = "NVUS0000000000-11"
VENDOR_ORDER_2 = "NVUS0000000000-12"

BARE_SUFFIXES = (".cram", ".csv", ".stderr", ".stdout", "_fail.csv")
TRIM_SUFFIXES = (
    ".trim.cram",
    ".trim.csv",
    ".trim.stderr",
    ".trim.stdout",
    ".trim_fail.csv",
)
SIDECAR_BARE = ".cram-metadata.json"
SIDECAR_TRIM = ".trim.cram-metadata.json"

JUNK_NAMES = (
    "login.html",
    "login.1.html",
    "ug-icon.png",
    "urls.txt",
    "objects_list_CRO-929.txt",
)

# (folder, wafer, stem, family, omit_cram, sidecar)
WELLS = (
    (
        "P04_1",
        "437120",
        "437120-437120-REF3_P04_1_A1-Z0001-CAGCTCGAATGCGAT",
        "bare",
        False,
        True,
    ),
    (
        "P05_1",
        "436830",
        "436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT",
        "trim",
        False,
        False,
    ),
    (
        "P07_1",
        "438514",
        "438514-438514-REF3_P07_1_A3-Z0305-CACACACAACATGAT",
        "bare",
        True,
        False,
    ),
    (
        "REF3_P05_2",
        "436831",
        "436831-REF3_P05_2_A11_GEX_hash_oligo-Z0170-CTCGCAATAGATGAC",
        "trim",
        False,
        True,
    ),
    (
        "P06_1",
        "439000",
        "439000-439000-REF3_P06_1_B2_GEX_hash_oligo-Z0400-ACGTACGTACGTACG",
        "bare",
        False,
        False,
    ),
    (
        "P06_1",
        "439000",
        "439000-439000-REF3_P06_1_B3-Z0401-ACGTACGTACGTACA",
        "bare",
        False,
        False,
    ),
)

# Wells the vendor delivery accounts for. G1/G2 (wafer 439000) are deliberately
# absent, mirroring REF3_P05_1 missing from the real order.
VENDOR_WELLS = (
    (
        VENDOR_ORDER,
        "437120",
        "437120-REF3_P04_1_A1_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT",
    ),
    (
        VENDOR_ORDER,
        "436830",
        "436830-REF3_P05_1_A10_GEX_hash_oligo-Z0169-CTCGCAATAGATGAT",
    ),
    (
        VENDOR_ORDER_2,
        "438514",
        "438514-REF3_P07_1_A3_GEX_hash_oligo-Z0305-CACACACAACATGAT",
    ),
    (
        VENDOR_ORDER_2,
        "436831",
        "436831-REF3_P05_2_A11_GEX_hash_oligo-Z0170-CTCGCAATAGATGAC",
    ),
)


# A vendor wafer the trimmed upload never mentions at all -- the inverse of
# wafer 439000 above, which is trimmed but undelivered. Opt-in via
# ``ref3_vendor_keys(extra_wafer=True)``, because every caller of the default
# fixture pins a well count or a verdict summary.
#
# Two wells rather than one, so a test can tell "this whole wafer is absent"
# from "one of its wells is absent": the defect this exists to pin is that
# locating vendor prefixes from the trimmed side makes an untouched wafer
# unreachable, and a single-well version of that is indistinguishable from an
# ordinary missing well.
UNTRIMMED_WAFER = "440000"
UNTRIMMED_WELLS = (
    (
        VENDOR_ORDER,
        UNTRIMMED_WAFER,
        f"{UNTRIMMED_WAFER}-REF3_P08_1_C1_GEX_hash_oligo-Z0500-ACGTACGTACGTACC",
    ),
    (
        VENDOR_ORDER,
        UNTRIMMED_WAFER,
        f"{UNTRIMMED_WAFER}-REF3_P08_1_C2_GEX_hash_oligo-Z0501-ACGTACGTACGTACT",
    ),
)


def ref3_trimmed_keys(include_junk: bool = True) -> list[str]:
    """Every object of the modelled trimmed upload."""
    keys: list[str] = []
    for folder, wafer, stem, family, omit_cram, sidecar in WELLS:
        suffixes = list(TRIM_SUFFIXES if family == "trim" else BARE_SUFFIXES)
        if omit_cram:
            suffixes = [s for s in suffixes if not s.endswith(".cram")]
        if sidecar:
            suffixes.append(SIDECAR_TRIM if family == "trim" else SIDECAR_BARE)
        keys.extend(f"{RAW}/{folder}/{wafer}/{stem}{s}" for s in suffixes)
    if include_junk:
        keys.extend(f"{RAW}/P04_1/437120/{name}" for name in JUNK_NAMES)
    return sorted(keys)


def ref3_vendor_keys(
    orders: tuple[str, ...] = (VENDOR_ORDER, VENDOR_ORDER_2),
    extra_wafer: bool = False,
) -> list[str]:
    """Vendor CRAMs and sidecars, in the measured vendor layout.

    ``extra_wafer`` adds :data:`UNTRIMMED_WELLS` -- a delivered wafer with no
    trimmed counterpart of any kind.
    """
    keys: list[str] = []
    wells = VENDOR_WELLS + (UNTRIMMED_WELLS if extra_wafer else ())
    for order, wafer, stem in wells:
        if order not in orders:
            continue
        base = f"{PROJECT}/{order}/{EXPERIMENT}/raw/{wafer}/{stem}"
        keys.extend([f"{base}.cram", f"{base}{SIDECAR_BARE}"])
        keys.append(
            f"{PROJECT}/{order}/{EXPERIMENT}/raw/{wafer}/{wafer}_LibraryInfo.xml"
        )
    return sorted(set(keys))


def ref3_sizes() -> dict[str, int]:
    """A plausible size per object: trimmed CRAMs well below their source."""
    sizes: dict[str, int] = {}
    for key in ref3_trimmed_keys():
        if key.endswith(".cram"):
            sizes[key] = 14_400_000_000
        elif key.endswith(".stdout"):
            sizes[key] = 0
        else:
            sizes[key] = 100
    return sizes


def vendor_uri(order: str) -> str:
    return f"s3://{VENDOR_BUCKET}/{PROJECT}/{order}"
