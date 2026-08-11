## Structure check

Cross-check a curation sheet **by structure instead of by name**. Part of BCP tooling under `bcp/structure_check`.

Comparing chemical names as strings does not work. ChEBI calls `CHEBI:92401` *"6-(1,3-dioxo-2-benzo[de]isoquinolinyl)-N-hydroxyhexanamide"*; a sheet calls it *"Scriptaid"*. Both are correct, and no string rule separates that harmless case from a genuinely wrong row.

InChIKeys do work, because they are derived from the structure itself. So this tool resolves each identifier in a row to a structure **independently** — the CAS number and the compound name through PubChem, the ChEBI ID through ChEBI — and compares those.

---

## The two questions it answers separately

| Column | Question |
|---|---|
| `id_cas_verdict` | Does the ChEBI ID hold the same molecule as the CAS number? — *"has ChEBI got something else under this ID?"* |
| `name_cas_verdict` | Is the CAS number right for the compound the row claims to be? |

Keeping them apart is the point. A row can have a perfectly correct ChEBI ID **for a CAS number that was itself wrong** — the ID check passes, the name check fails, and only then is the actual defect visible.

---

## CLI entry point

Run from the **`bcp`** directory:

```bash
cd bcp
python -m structure_check --help
```

### Batch mode (CSV)

`--cas-column` is the pivot for both checks. Pass `--chebi-column`, `--name-column`, or both.

| Flag | Description |
|------|-------------|
| `--input`, `-i` | Input CSV |
| `--cas-column` | CAS number column (default: `CAS`) |
| `--chebi-column` | Optional: ChEBI ID column to check against the CAS |
| `--name-column` | Optional: compound name column to check against the CAS |
| `--output`, `-o` | Output CSV (default: `<input_stem>_structure_checked.csv`) |
| `-v`, `--verbose` | Debug logging |

```bash
python -m structure_check --input curated.csv \
  --cas-column "CAS number" --chebi-column "ChEBI ID" --name-column "Name" \
  --output reviewed.csv
```

### Single row

```bash
python -m structure_check --cas 22573-88-2 \
  --name "Alexidine_dihydrochloride" --chebi CHEBI:27391
```

---

## The rule: sort by `review`

One column tells a curator what to do, ranked by what it would cost to be wrong:

| Value | Meaning |
|-------|---------|
| `investigate` | A **different molecular skeleton**. The row points at another compound. |
| `check` | Same skeleton with **different stereochemistry**, or the **same compound in a different salt form**. A form question, not a wrong compound. |
| `ok` | At least one comparison succeeded and every comparison made agreed. |
| `unverified` | No comparison could be made at all. Not a finding — see the verdict columns for why. |

Sort by the **`review_rank`** column (1–4, worst first) — everything needing a human lands at the top. `review` carries the same answer in words, but alphabetically `check` sorts before `investigate`, so the rank column is there to make "sort by it" literally true.

### Why the skeleton/stereo split matters

An InChIKey is `SKELETON(14)-STEREO+FLAGS(10)-PROTONATION(1)`. The first block hashes molecular connectivity, so a difference there means a genuinely different molecule. A difference only *after* it means the same skeleton with different stereochemistry.

That split is what makes this a real signal rather than a fuzzy one — and it earns its keep: `UCN-01` and the CAS its row records are **epimers of 7-hydroxystaurosporine**, which surfaces as `check`, not `investigate`. A name comparison could never have made that distinction.

### Salt forms are not wrong molecules

A counterion is part of an InChIKey's connectivity block, so "Doxorubicin" and doxorubicin **hydrochloride** hash to different skeletons — and a sheet naming one while citing the other's CAS is ordinary looseness, not an error. Left unhandled this was the tool's largest false-positive class: 10 of 16 initial `investigate` rows on a real 117-row sheet.

So when a skeleton difference is found, both sides are re-resolved to their **desalted parent** via PubChem and compared again. Agreeing parents demote the row to `salt_differs` / `check`.

This runs in one direction only: it can *downgrade* a difference already found, never claim a match. PubChem's parent assignment is unreliable for multi-component salts — for pyrvinium pamoate it returns the **pamoic acid counterion** — so a wrong answer there leaves a row flagged rather than waving one through. Pyrvinium is exactly that case and stays `investigate`.

**One caveat, stated precisely.** Because the demotion runs on both sides, two genuinely unrelated compounds cited as salts of the *same* counterion could each resolve to that counterion, agree, and have a real skeleton difference demoted from `investigate` to `check`. The row is still flagged and still needs a human, so nothing is waved through — but "false positives only" would be too strong a claim. In practice the demotion needs both `inchikey → CID → parent CID → InChIKey` chains to succeed, and for multi-component salts they usually do not (pyrvinium pamoate returns nothing at all via that path), which is why the pamoate rows stay at `investigate`.

---

## Verdict values

`id_cas_verdict` and `name_cas_verdict` each take:

- `match` — identical InChIKey
- `stereo_differs` — same skeleton, different stereo/isotope/charge layer
- `salt_differs` — different skeleton, but both sides share a desalted parent
- `skeleton_differs` — different molecule

Plus the reasons a comparison could not be made, which are **never findings about the compound**:

- `name_unresolved` — PubChem resolved no structure for the name
- `cas_unresolved` — PubChem resolved no structure for the CAS
- `chebi_unresolved` — malformed ID, or ChEBI has no such record
- `chebi_unreachable` — ChEBI could not be asked at all
- `name_ambiguous` — the name resolved to more structures than were compared, and none matched (see the cap below)
- `chebi_no_structure` — ChEBI has the entry but records no structure (class terms and R-group entries do not)
- `not_checked` — nothing was supplied for that check

`compare_structures()` additionally returns `not_comparable` when handed an empty structure on either side. That never reaches a row (`check_row` guards both callers), but the function is public, so it says "could not compare" rather than expressing it as a difference.

---

## Output columns

| Column | Contents |
|---|---|
| `review_rank` | 1 `investigate` · 2 `check` · 3 `unverified` · 4 `ok` — **sort on this** |
| `review` | The same answer in words |
| `id_cas_verdict`, `name_cas_verdict` | The two verdicts |
| `cas_inchikey` | Structure PubChem gives for the CAS |
| `cas_pubchem_name` | PubChem's preferred name for the CAS — useful for eyeballing a flagged row |
| `chebi_inchikey` | Structure ChEBI gives for the ID |
| `name_query` | **The exact string that resolved.** A difference must never be traceable to a query you cannot see. |
| `name_inchikey` | Structure(s) the name resolved to, pipe-separated |

---

## How names are resolved

Sheets write names loosely, so a small, auditable ladder is tried — whole-string forms first, first hit wins:

1. The cell as-is (`Scriptaid`, `PIK-75_HCl`)
2. Underscores as spaces (`Alexidine dihydrochloride`)
3. **Fallback only if neither resolved:** underscore-separated tokens, results unioned

**Salt information is never stripped.** `PIK-75 HCl` and `PIK-75` are different molecules with different InChIKeys, so dropping a `_HCl` would compare the wrong chemistry and report a difference of the tool's own manufacture. Any cell containing a counterion, salt, or hydration word (`SALT_TOKENS` in `client.py`) has the token fallback **disabled entirely** — reporting `name_unresolved` is better than comparing a molecule the row never claimed.

The token fallback exists for the `Vorinostat_SAHA` convention, where one cell holds two aliases for one compound. Results are unioned and matched against any candidate, so a stray token that resolves to something unrelated can only *mask* a finding, never invent one — the safe direction for a tool whose output drives edits.

---

## Package layout

```
bcp/structure_check/
├── __init__.py
├── __main__.py
├── cli.py       # argparse entrypoint
├── io.py        # CSV batch + single row
└── client.py    # resolvers, InChIKey comparison, verdicts
```

Composes the two sibling packages rather than re-implementing their HTTP: PubChem access comes from `chebi_lookup.client`, ChEBI access from `chebi_terms.client`. It inherits their retry and rate limiting.

**It does not inherit an answer-vs-outage distinction on the PubChem side, because there isn't one.** `chebi_lookup.get_with_retry` returns `None` for a 404 *and* for exhausted retries, so `cas_unresolved` and `name_unresolved` cannot tell "PubChem has no such compound" from "PubChem was unreachable". The ChEBI side can, because `fetch_compound` raises rather than returning `None` — which is why `chebi_unreachable` and `chebi_unresolved` are separate verdicts here rather than one.

### A check that never worked is not a clean sheet

Two ways a run can produce a complete-looking CSV that means nothing:

1. **Nothing was compared at all.** Point `--cas-column` at an existing but wrong column, or let PubChem go down, and every row lands on `unverified`.
2. **One side failed run-wide while the other passed.** This is the subtler one. If ChEBI is unreachable but the name column resolves fine, every row gets `chebi_unresolved`/`chebi_unreachable` on the ID side *and* a name match — so `review` says `ok`, `review_rank` is 4, and the sheet reads as clean with the ID question never answered on a single row.

So `RunSummary.degraded` fires on either, across at least 5 rows so a two-row spot check stays quiet, and the CLI exits **1**. `dead_checks` names the specific check that never succeeded, so the error says *which* question went unanswered rather than just that something went wrong. Findings themselves still exit 0 — they are the product.

A dead check is counted against the rows where that side was **attempted** — cells that were non-empty — not against every row. A partly-mapped sheet with blank ChEBI cells is the normal input for a tool whose job is checking ChEBI IDs, and measuring against total rows would let a single blank cell hide a total outage. A column that is entirely blank was never attempted, so it is not a dead check; those rows come out `unverified` and the first condition covers them.

**Programmatic use:**

```python
from structure_check import check_row

r = check_row(
    cas="22573-88-2", chebi_id="CHEBI:27391", name="Alexidine dihydrochloride"
)
print(r["review"], r["id_cas_verdict"], r["name_cas_verdict"])
# investigate match skeleton_differs
```

`compare_structures(reference_key, candidate_keys)` is pure — hand it InChIKeys you already have.

---

## Cost

Per distinct `(CAS, ChEBI ID, name)` triple: **2** PubChem calls for the CAS (CID, then a targeted `InChIKey,Title` property call — not `lookup_cas`'s four, whose xrefs and synonyms were being discarded), 1–4 for the name, and 1 ChEBI call. Plus 3 per structure for the desalted-parent check, which only runs where a difference was already found; successful parents are cached for the run, failures are not.

A row with neither a ChEBI ID nor a name costs nothing — the CAS is not resolved when there is nothing to compare it against. A name that resolves to more than `MAX_NAME_CANDIDATES` (10) structures is truncated, since each extra candidate costs three more requests during refinement. **Truncation can never produce a finding.** PubChem returns structures in CID order rather than relevance order, so the true match may sit past the cap; a truncated comparison that finds no match therefore reports `name_ambiguous`, not a difference. A match *inside* the compared set is still a match, since matching any candidate is sound. Either way the row records `(truncated: N structures)` in `name_query`, so a flagged row is auditable from the CSV alone. Identical triples are cached within a run, so repeated compounds cost nothing extra. Expect roughly 2–3 seconds per distinct row.

Findings exit **0** — they are the product, not a failure. Only a broken run (missing file, bad column, unwritable output) exits 1.

---

## Testing

**Default (offline, mocked HTTP):**

```bash
cd bcp
pytest tests/test_structure_check.py tests/test_structure_check_cli.py -v
```

**Live (opt-in, requires network):**

```bash
pytest -m "pubchem and chebi" tests/test_structure_check_live.py -v
```

The live tests pin the real findings this tool was built for, so a change in either upstream API that would silently flip a verdict fails here rather than in a curation sheet.
