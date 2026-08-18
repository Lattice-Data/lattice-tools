# Handoff: finish the adversarial review of the ChEBI logic promotion

Paste this whole file as the opening message of a new Claude Code session, with the
working directory at the repository root
(`/Users/gabdank/github/Lattice-Data/lattice-tools`).

---

## Situation

Reusable chemical-informatics logic has been promoted out of two disposable per-batch run
directories (`bcp/mashroom/`, `bcp/mashroom2/` — both gitignored, both slated for
deletion) into the permanent `bcp/` packages, so the reasoning survives the data. The
implementation is **done and green**. What is unfinished is the adversarial review of it.

Everything is committed on the branch **`chebi-promote-reasoning`**. Nothing is on `main`.

```bash
cd /Users/gabdank/github/Lattice-Data/lattice-tools
git checkout chebi-promote-reasoning
cd bcp && pytest tests/ -q     # expect 1878 passed, 13 skipped, 21 deselected
```

**Read `bcp/REVIEW_PROGRESS.md` first.** It is the authoritative state: what was built,
every one of the 44 review findings with a `DONE`/`OPEN` status, and the nine open
important items spelled out. `bcp/.chebi_review_findings.json` is the same list, machine
readable.

## Your job, in order

1. **Get verdicts for the 17 findings that never got one.** A review workflow
   (`wf_b7f1006e-d80`) spawned 49 agents; 17 verify agents died on a session limit, so
   those findings were never refuted or confirmed. Run `bcp/resume_review.sh`, which
   re-invokes that workflow with `resumeFromRunId` so completed agents replay from cache
   and only the dead ones re-run. **If you hit a session limit again, the script is
   idempotent — just run it again.** Keep going until every finding has a verdict.

2. **Act on whatever survives refutation.** The nine open `important` items are listed in
   `REVIEW_PROGRESS.md` with enough detail to act without this conversation. The single
   most substantive one: `bcp/CHEBI_IDENTIFICATION.md` credits `structure_check` with a
   name-resolution ladder whose rules mostly still live in
   `bcp/mashroom2/phase3_name_structure.py` and were **not** promoted. Either promote them
   or stop claiming they survive.

3. **Do not trust a refuter over a reproduction.** Several refuters in the first pass were
   wrong: they refuted findings that direct execution confirmed, including a blocking
   `cas.py` bug. Run the code yourself before dismissing a claim, and fix what reproduces.

4. **Keep the gate green** after every change: `ruff check bcp/`,
   `ruff format --check bcp/`, `pytest tests/` from `bcp`.

5. **Delete the process files before any PR** — `REVIEW_PROGRESS.md`, `REVIEW_HANDOFF.md`,
   `resume_review.sh`, `.chebi_review_findings.json`. They are committed only so a
   restarted session can resume from the branch.

## Standing constraints

- **Commit on `chebi-promote-reasoning`, never on `main`.**
- Two new modules are pure and offline by design: no network calls in
  `structure_check/inchi.py` or `structure_check/cas.py`, and no module-global caches
  (`test_structure_check.py` has an autouse fixture defusing `_parent_cache` precisely
  because cross-test bleed made results order-dependent).
- House documentation style in `bcp/` is decision-record prose: `##` H2 top level, never
  `#` or `####`, a `---` before every H2 but the first, rationale under its own headings,
  rejected alternatives recorded with the measured damage they did.
- Tests are plain functions annotated `-> None`, offline by default, with live API calls
  behind the `pubchem` and `chebi` markers. Fixtures are **real compounds**, never invented
  strings.
- The load-bearing principle throughout: **an answer from upstream and a silence from
  upstream must never render as the same output.** Most defects found here were violations
  of it in one direction or another.

## One thing to watch

Mid-task, two files changed underneath the working tree in ways that were not the
assistant's edits: `chebi_lookup/client.py` had `props.get("ConnectivitySMILES", "")`
replaced with `props.get("GoneAgain", "")`, and the `properties.json` unit fixture reverted
to PubChem's retired key names. A test caught it instantly. If a previously green test
starts failing, diff against the branch before assuming you broke it.
