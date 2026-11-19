"""Planning and enrichment shared by every output set of a fastq or cram run.

A run writes one set of outputs (an info TSV and the two submission sheets) per
order that owns files. A single-order run has one set; a run with
``--resequencing-orders`` can have several. The sets are planned, checked and
enriched together, and only written separately:

- alias collisions are checked across *all* sets at once, because a filename
  repeated in two sets is still one Lattice alias claimed twice by one lab;
- every file is fetched in one pass, so the worker pool and progress bar cover
  the whole run rather than restarting per set.

Row shaping stays in fastq.py and cram.py, which differ there.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Callable, Sequence

from .enrich import fetch_results
from .models import ListedObject
from .s3_utils import s3_uri_for
from .sheets import SequenceFileRecord, SheetOptions, validate_plan

RecordBuilder = Callable[..., SequenceFileRecord]


@dataclass(frozen=True)
class OutputTarget:
    """One set of outputs and the files written into it.

    Each file carries the order it was sequenced in, which lands in the set's
    ``CRO_order`` column and the info TSV's ``source_order``. A single-order
    run passes ``""``, meaning the run's ``--cro-order``.
    """

    output_path: str
    sheets: SheetOptions | None
    files: Sequence[tuple[ListedObject, str]]


def plan_and_fetch(
    s3_client: Any,
    bucket: str,
    targets: Sequence[OutputTarget],
    build_record: RecordBuilder,
    *,
    retries: int = 5,
    workers: int | None = None,
    show_progress: bool = True,
    inline: bool = False,
) -> tuple[list[list[SequenceFileRecord]], dict[str, dict[str, object]]]:
    """Plan each target in S3-key order, validate all of them, then fetch once.

    Returns one sorted record list per target, and the fetch results keyed by
    S3 key. Validation runs before any per-file request, so an unsubmittable
    run fails in the first second rather than after a request per file.
    """
    plans = [
        sorted(
            (
                replace(
                    build_record(
                        key=obj.key,
                        s3_uri=s3_uri_for(bucket, obj.key),
                        file_size=obj.size_bytes,
                    ),
                    cro_order=cro_order,
                )
                for obj, cro_order in target.files
            ),
            key=lambda record: record.sort_key,
        )
        for target in targets
    ]
    namespaces = {t.sheets.lab.namespace for t in targets if t.sheets is not None}
    if len(namespaces) > 1:
        raise ValueError("every output set of one run must share a lab namespace")
    if namespaces:
        validate_plan(
            [record for plan in plans for record in plan],
            namespace=namespaces.pop(),
        )

    keys = [record.key for plan in plans for record in plan]
    if not keys:
        # The pool is sized to the key count and refuses zero workers.
        return plans, {}
    results = fetch_results(
        s3_client,
        bucket,
        keys,
        retries=retries,
        workers=workers,
        show_progress=show_progress,
        inline=inline,
    )
    return plans, results
