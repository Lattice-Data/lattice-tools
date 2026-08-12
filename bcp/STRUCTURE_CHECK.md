## Structure check

Cross-check a curation sheet **by structure instead of by name**. Part of BCP tooling under `bcp/structure_check`.

Comparing chemical names as strings does not work. ChEBI calls `CHEBI:92401` *"6-(1,3-dioxo-2-benzo[de]isoquinolinyl)-N-hydroxyhexanamide"*; a sheet calls it *"Scriptaid"*. Both are correct, and no string rule separates that harmless case from a genuinely wrong row.

InChIKeys do work, because they are derived from the structure itself. So this tool resolves each identifier in a row to a structure **independently** — the CAS number and the compound name through PubChem, the ChEBI ID through ChEBI — and compares those.

The third of the three ChEBI tools here, and the only one that cross-checks identifiers against *each other*: [CHEBI_LOOKUP.md](CHEBI_LOOKUP.md) goes from a CAS number to a ChEBI ID, [CHEBI_TERMS.md](CHEBI_TERMS.md) goes from a ChEBI ID you already hold to its authoritative name and a correctness verdict.

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
| `ok` | At least one comparison succeeded and every comparison made agreed. A side that could not be compared — `name_ambiguous`, `name_unresolved`, `chebi_unreachable` — does not undo the side that did; check the verdict columns and `unasked` for what was left open. |
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
- `name_ambiguous` — the name resolved to more structures than were compared, and none matched (see the cap below)
- `chebi_no_structure` — ChEBI has the entry but records no structure (class terms and R-group entries do not)
- `not_checked` — nothing was supplied for that comparison. Either the row's own cell was blank, or the **CAS** was: the CAS is the pivot for both checks, so a blank one leaves nothing to compare a name or an ID *against*, and no request is made for either.

And separately, the three that mean **the upstream was never reached** — the only verdicts that say something about the network rather than about the compound:

- `chebi_unreachable` — ChEBI could not be asked at all
- `cas_unreachable` — PubChem could not be asked about the CAS
- `name_unreachable` — PubChem could not be asked about the name

That line matters more than it looks. "PubChem has no such CAS" is evidence about the compound; "PubChem never answered" is evidence about the afternoon. Only the second kind can be fixed by running again, and only the second kind is allowed to make the run itself untrustworthy.

`compare_structures()` additionally returns `not_comparable` when handed an empty structure on either side. That never reaches a row (`check_row` guards both callers), but the function is public, so it says "could not compare" rather than expressing it as a difference.

---

## Output columns

| Column | Contents |
|---|---|
| `review_rank` | 1 `investigate` · 2 `check` · 3 `unverified` · 4 `ok` — **sort on this** |
| `review` | The same answer in words |
| `unasked` | `id`, `name`, `both`, or empty — which requested check went unasked because an upstream was unreachable. **Not** a finding: these rows need a re-run, not a curator. |
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

Both sides can now tell an answer from a silence. `chebi_terms.fetch_compound` raises `ChebiUnavailableError`; on the PubChem side `chebi_lookup.get_with_retry_status` reports whether a request 404'd or ran out of retries, and the resolvers here raise `PubChemUnavailableError` when it was the latter. `get_with_retry` remains as a thin wrapper for callers that only want the payload — **same signature and same return semantics**, though its existing callers do inherit three behaviour changes: malformed-request 4xx are no longer retried, the pointless sleep after the final attempt is gone, and `cas_to_cid` now percent-encodes the CAS before it reaches the URL path (an injection fix, since these values come from untrusted sheets).

A 200 carrying something that is not this endpoint's JSON — a maintenance page, a moved endpoint — counts as unreachable too. Read as "no such compound" it would look exactly like a column full of junk.

### Is this run trustworthy?

A complete-looking CSV whose emptiness means nothing is the worst thing this tool could produce, so `RunSummary.degraded` answers that question and the CLI exits **1** when it is true. Findings themselves still exit 0 — they are the product, and the "N rows need attention" line is logged *before* the degraded exit so a network blip cannot swallow the pointer to real findings.

**`degraded` is three separate signals, not one predicate, and that is the whole design.** A single boolean here was wrong four times running, each fix correct for the input it targeted and blind to the next. The reason is that it was being asked two unrelated questions at once — *"was upstream up?"* and *"is this the right column?"* — which have different evidence, different natural units, and different remedies. A re-run fixes one and cannot fix the other.

| Signal | Fires when | Counted in |
|---|---|---|
| `outages` | on some side, **more rows went unanswered by an unreachable upstream than the upstream answered at all**, and at least 5 such rows | rows |
| `suspect_columns` | on some side, nothing resolved across **≥ 5 distinct values**, with no outage on *that* side | distinct values |
| `nothing_compared` | every row came out `unverified` — across at least 5 rows, or any number if an outage was involved | rows |

Why each is shaped the way it is:

- **A majority, not "any outage at all."** The sibling `chebi_terms` degrades on a single failed lookup, and can afford to: it makes one request per ID. This tool makes 3–44 per row against an API that rate-limits, and a throttle that outlasts its retries is indistinguishable from an outage. Firing on one occurrence would exit 1 on healthy large runs, and an exit code that cries wolf stops being read — which costs more than the case it was added to catch. A majority needs no tuned constant: once most of what was asked went unanswered, the output is not a verdict whatever the rest says.
- **Rows below that line are still not hidden.** They get `cas_unreachable`/`name_unreachable`/`chebi_unreachable` in their verdict column, a mark in **`unasked`**, and a WARNING naming the count. The `review` column deliberately does *not* change: `ok` means every comparison that was made agreed, which stays true when the other side was unreachable. Overloading the one documented sort key would bury a real `skeleton_differs` behind a transient 503.
- **Distinct values, not rows, for a wrong column.** A dose-response sheet repeats each compound thirty times; 500 rows of one junk string is a single mistake, not 500 pieces of evidence. Two retired ChEBI IDs on a 117-row sheet is a *finding* — exactly what this tool exists to surface — and must not be reported as a broken run.
- **The outage guard is per side.** A ChEBI column full of free text is rejected by a regex without a single request, so a PubChem outage is logically incapable of explaining it. Disqualifying globally would hide a real misconfiguration behind an unrelated blip.
- **A 404 is an answer, and counts as one.** The majority weighs unanswered rows against every row the upstream *replied* to, resolved or not. Counting only the rows that resolved let five throttled rows outweigh a hundred honest "no such compound" replies — and it misfired worst against a wrong column, where nothing resolves by construction, so the outage branch fired and suppressed the diagnosis that would have named the real problem.
- **The two floors count different things.** `MIN_OUTAGE_ROWS_FOR_DEGRADED` exists because the majority test degenerates against zero: with nothing answered, a *single* outage is a "majority", which is the normal shape of a sparse column (two ChEBI IDs on a 117-row sheet, both unlucky) and of a wrong one. What has to be substantial is the evidence of an outage, not the size of the sheet.
- **`nothing_compared` earns its place** because it is the only term that reasons about *comparisons performed* rather than identifiers resolved, and the two genuinely diverge. A ChEBI column of class terms resolves perfectly — ChEBI holds every record — and still compares nothing, because a class term carries no structure. Both other signals stay silent on that sheet. Its row floor lifts when an outage is in play, so a four-row batch against a dead ChEBI is not quieter than the same row would be in single mode — which exits 1 on one unanswered check. An outage is positive evidence that we failed to ask, so it needs no volume; absence of results *without* one is only evidence in bulk.

Blank cells never enter any of these counts. A partly-mapped sheet is the normal input for a tool whose job is checking ChEBI IDs, and a column that is entirely blank is not evidence about anything.

Two things are deliberately excluded from the outage signal. **Upstream answers** — `chebi_unresolved`, `chebi_no_structure`, `name_ambiguous` — are answers, not silences, and belong in the `Not compared` breakdown rather than in a predicate that exits 1 and tells the operator to check their network. And **the salt-form refinement**: it only ever downgrades a difference already found, so a failure there leaves a row flagged for a human, the safe direction. It is also by far the highest-volume caller, so counting it would exit 1 on a run whose only content is a genuine finding. It logs a warning instead.

### Giving up on one upstream without giving up on the run

A side that fails `MAX_CONSECUTIVE_OUTAGE_ROWS` (5) attempts in a row is treated as down and **stops being queried**: its rows are marked `unreachable` without a request, which is what they are, and the run carries on with the columns that still answer. Tracked per side, because a whole-row rule cannot see the case that costs most — with ChEBI dead and a live `--name-column` whose names match, every row reads `ok`, so a row-level streak resets on every row and the sheet walks all of itself at ~14s of backoff each to reach a verdict already settled by row 5.

It stops rather than *aborts* deliberately. Killing the run would throw away the name check on every remaining row because ChEBI is down; the exit code and the CSV say the same thing either way, so there is nothing to gain by discarding the work that still succeeds.

Being tripped is an *inference*, so it is re-tested every `OUTAGE_PROBE_INTERVAL` (20) attempts, and a successful probe clears it. Without that, five scattered blips on a sparsely mapped column would poison every remaining row with an `unreachable` no request was ever made for — an answer manufactured from a silence, which is the one thing this tool refuses to do.

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

Per distinct `(CAS, ChEBI ID, name)` triple: **2** PubChem calls for the CAS (CID, then a targeted `InChIKey,Title` property call — not `lookup_cas`'s four, whose xrefs and synonyms were being discarded), typically 1–2 for the name (more if neither whole-string form resolves and the token fallback tries each token in turn), and 1 ChEBI call. Plus 3 per structure for the desalted-parent check, which only runs where a difference was already found and the name was not truncated; successful parents are cached for the run, failures are not.

A row with neither a ChEBI ID nor a name costs nothing — the CAS is not resolved when there is nothing to compare it against — and so does a row whose **CAS** is blank or unresolvable, since the pivot is gone and neither comparison can be made whatever the other identifiers turn out to be. That second case is what a misdirected `--cas-column` looks like, and it is the difference between paying full price for every row of a broken run and paying for one lookup. A name that resolves to more than `MAX_NAME_CANDIDATES` (10) structures is truncated, since each extra candidate costs three more requests during refinement. **Truncation can never produce a finding.** PubChem returns structures in CID order rather than relevance order, so the true match may sit past the cap; a truncated comparison that finds no match therefore reports `name_ambiguous`, not a difference. A match *inside* the compared set is still a match, since matching any candidate is sound. Either way the row records `(truncated: N structures)` in `name_query`, so a flagged row is auditable from the CSV alone. Identical triples are cached within a run, so repeated compounds cost nothing extra — **except** a triple whose lookup hit an outage, which is retried rather than cached, so one unlucky moment does not decide every row that repeats it. Expect roughly 2–3 seconds per distinct row.

Findings exit **0** — they are the product, not a failure. Only a broken run (missing file, bad column, unwritable output, or a degraded run — see above) exits 1. Single mode holds the same contract from the other end: it exits 1 when *every* requested check went unasked because an upstream was unreachable, so `--cas X --name Y || alert` does not stay silent through an outage.

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
