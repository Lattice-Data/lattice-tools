"""Non-interactive sequencing-order QA.

Wraps the checks that ``qa.ipynb`` runs cell by cell so an order can be QA'd in
one command, with every check accounted for in the output. The notebook now
covers seven assays across raw and processed data, and its gating conditions
live in the cells themselves -- so a cell that prints "skipped" scrolls past as
easily as one that passes. Here the set of checks is declared data (see
``registry``), every declared check lands in the report with a status, and
"skipped" is a reported outcome rather than an absence.

Reuses ``resolve_qa_run_context``, ``gather_qa_data`` and the ``qa_checks``
validators unchanged; nothing about what a check means is re-implemented here.
"""

from __future__ import annotations
