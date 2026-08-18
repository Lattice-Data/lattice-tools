# ChEBI promotion — review progress and resume state

Working notes for an in-flight task. **Delete this file (and `REVIEW_HANDOFF.md`,
`resume_review.sh`, `.chebi_review_findings.json`) before opening a PR** — they are
process state, not repo artifacts. They are committed only so a restarted session can
pick the work up from a branch.

Branch: `chebi-promote-reasoning`. Run from the `bcp` directory throughout.

## What the task is

Promote the reusable chemical-informatics logic out of the disposable per-batch run
directories (`bcp/mashroom/`, `bcp/mashroom2/`, both gitignored and slated for deletion)
into the permanent packages, so the reasoning survives the data. Agreed scope:

1. `bcp/CHEBI_IDENTIFICATION.md` — the method record. **done**
2. `bcp/structure_check/inchi.py` — InChI layer comparison, no network. **done**
3. `bcp/structure_check/cas.py` — CAS validation and repair, no network. **done**
4. Offline regression tests carrying the real compounds that failed. **done**
5. Four confirmed bugs in the permanent packages, each with a test. **done**
6. `bcp/.gitignore` for the per-run directories. **done**

## Gate state

| check | state |
|---|---|
| `pytest tests/` | **1883 passed**, 13 skipped, 21 deselected |
| `ruff check bcp/` | clean |
| `ruff format --check bcp/` | clean |
| new offline tests | 95 (38 inchi + 41 cas + 16 chebi_lookup) |

## Adversarial review: complete

The review workflow (`wf_b7f1006e-d80`) raised **44 findings** across five dimensions.
The first pass lost **17 of 49 verify agents to a session limit**; the run was resumed
with `resumeFromRunId`, which replayed the cached agents and re-ran only the dead ones.
It now reports **49 of 49 agents done, 0 errors** — every finding has a verdict.

- **57 refuted**, **14 survived** refutation.
- Every survivor has been acted on: fixed, or (for the one irreducible ambiguity and the
  one scope boundary) documented explicitly.

**Nothing is open.** The gate is green at 1883 tests.

Where a reviewer's claim was reproduced by direct execution, it was fixed regardless of
what a refuter said — several refuters were wrong. Where it could not be reproduced, it
was left alone.

### The four blocking findings, all fixed

| file | defect | fix |
|---|---|---|
| `cas.py` | `classify_cas("3-4-07689")` → `('07689-03-4','valid')`; rotation output never re-normalised, and leading zeros are invisible to the checksum | rotated result re-runs the leading-zero repair |
| `inchi.py` | `/c` and `/h` never compared, so ethanol vs dimethyl ether → `layers_identical` | tiered comparison incl. `CONSTITUTION_LAYERS` |
| `inchi.py` | layers keyed by assignment, so `/t /m /s` re-emitted after `/i` overwrote the real ones | `setdefault` — first occurrence wins |
| tests | (same as the `/c`/`/h` finding, reported from the test dimension) | 3 tests added |

### What the second pass found that the first missed

- **`comparison_verdict_from_inchi()` did not exist**, though `inchi.py`'s docstring told
  callers to use it. Written, and it asserts at import time that it maps every verdict.
- **Three doc numbers contradicted the run data**: the batch-1 alternate-registry count
  (six → **17**), the batch-2 xref coverage (94% → **91%** of resolved CAS numbers), and
  the regression count (five of six → **12 of 14**).
- **The name ladder was credited to `structure_check`** when five of its six rules live
  only in the per-run script being deleted. The section now says so in its first sentence
  and reads as a specification to re-implement from. Promoting the ladder is the clear
  next piece of work and is *not* done.
- **Four tests passed for the wrong reason**, each demonstrated by a mutation that stayed
  green: `"@" not in` satisfied by an empty string; a numeric synonym silently dropped;
  a short-CSV-row test that patched nothing and would have reached live PubChem; and the
  multiplier stripping in `classify_pair` with no coverage at all — the ABT-702 fixture
  carries its multiplier on the *counterion*, so it could never pin it. A real
  apomorphine hydrochloride / hemihydrate pair does, and the mutation now fails.
- **The `properties.json` unit fixture was genuinely dead.** Live routing keys off the
  CID and every mocked lookup yields CID 702, so the recorded fixture always wins. It is
  now re-keyed to PubChem's current shape *and* pinned against the recorded one, so
  neither can drift.

### Deliberately not fixed, with the reason

- **A genuine date can still classify as a valid CAS.** `2020-01-01 00:00:00` → `2020-01-1`
  passes the checksum. No string rule separates it from `2113-05-05`, whose real number is
  year-shaped. Documented, and pinned by a test so it is known rather than surprising.
- **`FORM_DIFFERS` returns before the stereo comparison**, so a salt pair that also
  mismatches on stereo reports only the form difference. It still refuses the row, which is
  the safe direction; a combined verdict would be more informative but changes the
  vocabulary.
- **Two committed CSVs** (`cas_batch_input*_mapped.csv`) carry the empty SMILES columns the
  fixed bug produced. Regenerating them needs live PubChem calls — a separate job.

## An injected fault, and why it matters

Mid-task, `chebi_lookup/client.py` line 265 changed to `props.get("GoneAgain", "")` and
`tests/fixtures/chebi_lookup/properties.json` reverted to the retired key names — neither
was my edit. `test_lookup_cas_populates_both_smiles_columns` failed immediately with
`assert '' == 'CCO'`. Both restored. Keep that test: it is the only thing that objected.

## How to resume

```bash
cd /Users/gabdank/github/Lattice-Data/lattice-tools/bcp
git checkout chebi-promote-reasoning
./resume_review.sh          # re-runs only the verify agents that died, from cache
pytest tests/ -q           # must stay green
```

`resume_review.sh` re-invokes the review workflow with `resumeFromRunId`, so every agent
that already completed replays from cache and only the dead ones re-run. It loops with
backoff and is safe to run repeatedly.
