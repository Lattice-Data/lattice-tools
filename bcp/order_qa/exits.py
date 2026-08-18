"""Exit codes.

A caller has to be able to tell "the vendor is still uploading" (wait and retry)
from "the upload is wrong" (tell the vendor) from "the data is wrong" (tell the
lab) without parsing prose, so each gets its own code. ``DEGRADED`` exists for
the same reason ``structure_check`` exits 1 on an upstream outage: a run that
could not ask its questions is not a run that got clean answers, and a wrapper
doing ``order-qa ... || alert`` must not stay silent through one.

Findings are the product of a QA run, not a failure of it -- but unlike
``structure_check``, where a finding is a note on one row, a finding here means
the order is not releasable. So ``QA_FINDINGS`` is non-zero.
"""

from __future__ import annotations

from typing import Final

OK: Final = 0
"""Verified complete, quiescent, and every applicable check passed."""

INTERNAL: Final = 1
"""The run could not be completed or recorded -- e.g. the report could not be
written. Named rather than left to an unhandled traceback's implicit 1, so a
caller can tell "this tool broke" from any verdict about the order."""

USAGE: Final = 2
"""Bad arguments. argparse's own code, kept so the two cannot disagree."""

STILL_UPLOADING: Final = 3
"""Objects were still arriving inside the quiescence window. Retry later."""

VERIFICATION_FAILED: Final = 4
"""The upload itself is wrong: missing files, or unexpected extras."""

QA_FINDINGS: Final = 5
"""Upload verified, but QA found defects in the data."""

DEGRADED: Final = 6
"""At least one applicable check could not run. The order has no verdict."""

RESOLUTION_FAILED: Final = 7
"""The spec named nothing QA-able: unknown prefix, or a non-order path."""

_NAMES: Final = {
    OK: "ok",
    INTERNAL: "internal_error",
    USAGE: "usage",
    STILL_UPLOADING: "still_uploading",
    VERIFICATION_FAILED: "verification_failed",
    QA_FINDINGS: "qa_findings",
    DEGRADED: "degraded",
    RESOLUTION_FAILED: "resolution_failed",
}


def name_for(code: int) -> str:
    """Human name for an exit code, for the status line and the report."""
    return _NAMES.get(code, "unknown")
