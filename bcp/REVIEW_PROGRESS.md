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
| `pytest tests/` | **1878 passed**, 13 skipped, 21 deselected |
| `ruff check bcp/` | clean |
| `ruff format --check bcp/` | clean |
| new offline tests | 90 (35 inchi + 40 cas + 15 chebi_lookup) |

## Adversarial review: 44 findings, 27 verified, 17 never verified

The review workflow (`wf_b7f1006e-d80`) raised 44 findings across five dimensions.
**17 of 49 verify agents died on a session limit**, so those findings have no
refutation verdict. `.chebi_review_findings.json` holds all 44 with a `status` field.

- **20 addressed** in the working tree (4 blocking, 8 important, 8 nit).
- **24 open** (0 blocking, 9 important, 15 nit).

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

### Open, important — these need a decision or a fix

1. **`inchi.py`: `FORM_DIFFERS` returns before the stereo comparison.** A salt pair that
   *also* has a stereo mismatch reports only `form_differs`, so the reason code is
   incomplete. Safe direction (it still refuses), but possibly worth a combined verdict.
2. **`CHEBI_IDENTIFICATION.md`: "six apparent disagreements in the first batch"** is a
   pre-correction snapshot; re-running the audit gave 17. Check and correct.
3. **`CHEBI_IDENTIFICATION.md`: "Five of six regression assertions fail against the
   implementation they replaced"** — was measured against the *stereo-only* predecessor,
   before the `/c`/`/h` fix. Re-measure or reword.
4. **`CHEBI_IDENTIFICATION.md`: the name ladder is attributed to `structure_check`**, but
   five of its six rules live only in `mashroom2/phase3_name_structure.py` and were NOT
   promoted. Either promote them or stop crediting the package for them. **This is the
   most substantive open item: the doc currently overstates what survives deletion.**
5. **`tests`: the "never identical" parametrize only covers junk with no `/`** — add a
   truncated-InChI case (a real one now exists elsewhere in the file).
6. **`tests`: `test_map_cas_file_tolerates_a_short_row` may reach live PubChem** — it
   patches nothing. Verify it is genuinely offline; if not, patch `requests.get`.
7. **`tests`: `"@" not in canonical_smiles` is satisfied by an empty string** — assert the
   exact expected value instead.
8. **`tests`: `test_numeric_synonyms...` passes if numeric synonyms are silently dropped**
   — assert the numeric one survives, not just that "ethanol" is present.
9. **`properties.json`: the fixture edit may be dead** — `route_pubchem_get` prefers the
   recorded `pubchem_live/64-17-5/` fixture, which already had the new keys. Confirm which
   fixture the SMILES tests actually exercise; if the unit fixture is unused by them, the
   test does not pin what it claims.

### Open, nit — lower value, listed for completeness

`inchi.py` docstring points at `structure_check.client.comparison_verdict_from_inchi()`,
**which does not exist** (either write it or drop the reference); no whitespace
normalisation on InChI input; `defined_stereo` counts rather than checks containment;
batch-2 xref "94%" not reproducible from run data (89% or 95% depending on denominator);
`test_every_verdict_returned_is_declared` is close to a tautology; `defined_stereo`'s last
assertion is self-referential; two committed CSVs (`cas_batch_input*_mapped.csv`) carry the
empty SMILES columns the fixed bug produced; `mashroom/HANDOFF_NEXT_BATCH.md` still lists
the SMILES rename as open (that file is gitignored and being deleted — ignore).

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
