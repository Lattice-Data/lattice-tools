## ChEBI identification

Establish a **verified ChEBI identifier** for every compound in a batch of experimental
perturbations, given only a name and a CAS number per row. Part of BCP tooling. There is
no CLI: this is a method run in phases against one spreadsheet at a time, and its
reusable parts live in `bcp/structure_check` (InChI layers, CAS validation, structure
comparison) and `bcp/chebi_terms` (ChEBI records, names, CAS accessions).

The constraint that forces the design: a name and a CAS number in adjacent cells look
like two facts about one compound, and often are not. Across two real batches, 25 rows
in one and 4 in the other had a name and a CAS that resolve to **unrelated molecules** —
`Haloprogin` carrying Asulam-sodium's number, `Calmidazolium chloride` carrying
chelerythrine chloride's. No amount of care with either identifier alone finds these.
Resolving both **independently** and comparing structures does, and it is the reason a
batch of this size can be done at all.

To go from a CAS number to a ChEBI ID for a single compound, see
[CHEBI_LOOKUP.md](CHEBI_LOOKUP.md). To check a ChEBI ID you already hold, see
[CHEBI_TERMS.md](CHEBI_TERMS.md). To cross-check a whole sheet's existing identifiers
against each other by structure, see [STRUCTURE_CHECK.md](STRUCTURE_CHECK.md). This
document is the layer above all three: it starts from a sheet with **no** trustworthy
identifier and produces one, or says precisely why it cannot.

Upstream services: [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest),
the ChEBI backend REST API, and ChEBI's Elasticsearch index. The EBI SOAP service and
`libChEBIpy` are both dead — do not reach for them.

---

## What "identified" has to mean

Two questions have to be settled in order, and the first gates the second.

| Question | Settled by |
|---|---|
| **Which molecule does this row denote?** | the name and the CAS resolving *separately* to the same structure |
| **Does ChEBI hold exactly that molecule?** | a full InChIKey match to the form the row names |

Keeping them apart is the point. A row that disagrees with itself about which compound
it is cannot have a ChEBI ID verified for it however good the ChEBI match looks, and
reporting "not in ChEBI" for such a row is a false claim about ChEBI. The verdict has to
name which of the two failed, because a curator does something completely different in
each case: one is a question for the experimental team, the other is a submission to
ChEBI.

**Identity needs two independent lines.** A row where only one side resolves is not
eligible, and this rule has a measured cost: it refuses `Thimerosal` → CHEBI:9546 when a
sheet gives no CAS, which is a false negative. Those rows are the cheapest for a human to
rescue — supplying the missing number is all it takes — so the rule is kept and the cost
is reported rather than absorbed.

---

## The phases

Exactly one data file persists between phases, gaining columns as it goes. Nothing is
forked, and no phase reads another's rendered output — a compound whose name contains a
separator character cannot corrupt a verdict.

| Phase | What it does | Network |
|---|---|---|
| 0 | load the sheet, drop non-perturbagens, deduplicate, match against finished batches | no |
| 2 | CAS → PubChem CID → InChIKey | yes |
| 3 | name → PubChem CID → InChIKey, **independently of the CAS** | yes |
| 4 | InChIKey → ChEBI by four routes | yes |
| 5 | adjudicate to a verdict per row | no |
| 6 | independent cross-checks over the verified set | yes |
| 7 | write the deliverables | no |

Phase 3 is the one that matters most. Phase 5 is deliberately separate from phase 4:
the policy is a judgement call that may be revisited, and re-deciding must not mean
re-querying.

### What carries between batches

Two things persist, and confusing them is expensive.

`chebi_answer_ledger.tsv` is **committed**, and is the reason each batch costs less
than the last: one row per distinct `(name, CAS)`, with the verdict, the identifier,
the structure it was matched on and the evidence. 1125 rows so far — 722 verified, 403
diagnosed dead ends. `chebi_ledger.py` reads, merges and writes it.

It is tracked while every other data file here is ignored, because it is the one that
**cannot be regenerated at any price**. A network cache is slow to rebuild;
these verdicts embody a round of 113 questions to the experimental team, a set of
hand-adjudicated promotions, and a stereochemistry correction that withdrew eight
published identifiers. Losing the file means asking the lab the same questions again.

The **network caches** are regenerable and live outside the repository, at
`~/.cache/lattice-bcp-chebi/{pubchem,chebi}.json` — 2362 and 4659 entries, about
14,100 HTTP requests and 45 minutes of wall clock inside PubChem's rate limits. Worth
keeping, not worth committing, and *not* worth trusting for SMILES: they were written
before the `SMILES`/`ConnectivitySMILES` rename was fixed, so every cached
`isomeric_smiles` and `canonical_smiles` is empty. InChI, InChIKey and formula are
unaffected, which is everything the method uses.

**A later batch may change an answer, but never silently.** `merge()` returns every
key whose verdict or identifier moved, so a run reports it. Those eight withdrawn
identifiers are precisely what that return value is for: an answer moving is a finding
about what has already been deposited, not bookkeeping.

### Deduplicate against finished batches before any network work

A second batch from the same screening library is mostly the first batch again. The
second real batch was 3597 rows, which collapsed to 1106 distinct `(name, CAS)` pairs, of
which **987 already had a verdict** — 89% answered for free, covering 2340 of 2612 source
rows. Only 119 needed the network. Across both batches there are 1125 distinct keys and
986 of them appear in both, so the overlap is the norm rather than a coincidence.

Match on the **full set** of CAS values a finished row carries, not just the sheet's
original: a completed row may hold the original, a mechanical repair, a lab correction
and the value actually queried. Comparing against the original alone reported 12
already-repaired numbers as fresh disagreements.

**A half match never inherits an answer.** Same name with a different CAS, or the same
CAS under a different name, is exactly the signal this method exists to interrogate. Send
it through the phases with the previous answer attached for comparison; do not resolve it
by eye. Nine such rows in the second batch produced two real findings — a doxorubicin row
that had swapped to the free base's CAS, and a row that had **reverted** to a wrong CAS
the lab had already corrected once.

---

## Locked policy decisions

- **Salts: strict.** Verified requires a full InChIKey match to the form the row names.
  The free base is not the hydrochloride. Lattice DB has no salt-form field, so the form
  has to be right in the identifier itself — which also means a salt absent from ChEBI is
  a submission candidate, never something to paper over with the parent's ID.
- **Stereochemistry: strict, with one exception that is not a loosening.** A *racemate*
  matched to a stereo-**unspecified** ChEBI entry is fine — that is how ChEBI represents a
  racemate. A *single isomer* matched to an unresolved parent is not. This is why
  `(±)-Nipecotic acid` → CHEBI:116931 was accepted and `AT 101`, which is
  (−)-gossypol, was **refused** CHEBI:28584.
- **A PubChem stereo duplicate is not a disagreement.** Where one PubChem record omits a
  stereo layer the other specifies, that is PubChem's cataloguing. Treat as agreement and
  anchor on the side that *defines* the stereochemistry.
- **Never auto-resolve a name/CAS conflict.** The name is usually right — 23 of 25 in the
  first batch — but `Compound W` resolves by name to salicylic acid, because that is a
  wart-remover brand, while its CAS correctly named the research compound. `Lylamine` was
  a corruption of *Leelamine* whose CAS belonged to *allylamine*. An earlier draft emitted
  `Compound W → 69-72-7` as a confirmed fix; shipping that would have been worse than
  shipping nothing.
- **A pre-filled ChEBI ID or ChEBI term column in the source is noise.** In the first
  batch these were wrong in every case checked, and for 11 rows the wrong term named the
  *same wrong compound* the wrong CAS resolved to — one error copied across several
  columns. Carry them as `*_untrusted`; never read them as input.

---

## Traps the check digit cannot catch

`structure_check.cas` implements this. **A valid check digit means the number is
well-formed, not that it is the right compound's number.** Two corruptions are
arithmetically invisible to it:

- **A spurious leading zero.** `0362-07-2` for `362-07-2`. Digits are weighted by position
  from the right, so a leading zero contributes `0 x weight` and shifts nothing: the two
  strings produce the *identical* check digit. Only PubChem catches it — 6 rows in the
  first batch.
- **A zero-padded check digit.** `42971-09-05` for `42971-09-5`. Worse than cosmetic,
  because it makes the string look like `YYYY-MM-DD` and a spreadsheet then converts it
  to a date. `2113-05-05 00:00:00` is m-chlorophenylbiguanide hydrochloride's
  `2113-05-5`, two corruptions deep. **41 numbers** in the second batch were zero-padded
  and 13 of those had become dates. The two are one causal chain, not two coincidences.

Two more shapes, both repairable, and each repair recorded so a curator can see that the
value looked up is not the value in the cell: **rotated segments** (`3-4-7689` for
`7689-03-4`, where 12 of 12 repairs validated afterwards — about one in 10^12 by chance)
and **mojibake separators** (`864461?31?4`).

**One shape that is deliberately not repaired.** A bare digit run such as `6857789` looks
like `6857-78-9`, and something upstream once made exactly that inference. The value was
a PubChem **CID**, and `6857-78-9` is a real registry number belonging to an unrelated
compound. Guessing does not recover a lost number, it invents a wrong one.

**A leading zero can survive every repair and still validate**, because the checksum
cannot see it. Fuzzing the repairs found three ways in: `00-00-0` sums to zero and its
check digit *is* zero, so a zero-filled placeholder cell passed; `007-00-1` came out of the
leading-zero repair as `07-00-1`, still corrupted but reported repaired; and a rotated
value hid one in its last segment, so `3-4-07689` rotated to `07689-03-4`, passed the
check digit, and was reported valid. A real registry number has no leading zero in its
first segment, and that is now the final guard rather than an assumption. Relatedly, a
pattern written with `\d` accepts *every* Unicode decimal and `int()` parses them, so
`٥٠-٠٠-٠` matched and satisfied the checksum — a string PubChem can never resolve.

**One ambiguity is irreducible, and is documented rather than guarded.** Because a padded
check digit makes a CAS look like a date, the converse holds too: a genuine date with a
single-digit day satisfies the checksum about one time in ten, so `2020-01-01 00:00:00`
becomes `2020-01-1` and validates. No string rule separates that from `2113-05-05`, whose
number really is year-shaped. The repair is applied because in a CAS column it is right far
more often than not, it is always reported, and PubChem settles it.

**A CAS recovered from PubChem synonyms is reliable only when the sheet's own number
fails its check digit.** Where the sheet's number is well-formed but merely unindexed, the
sheet is usually right and the synonym is a different form. Seven rows were demoted to
"unproven leads" on this basis, and vendor lookups later confirmed the sheet was right in
three of them (`Cromakalim`, `PSB 11`, `UB 165 fumarate`) — so the caution was correct.

---

## Resolving the name independently of the CAS

**Scope note, stated first because the rest of this section is easy to misread.** These
rules are **recorded here, not promoted.** `structure_check.client` supplies
`name_candidates()` and the `SALT_TOKENS` list; everything else below —
`_COUNTERION`, `fragments()`, `lost_apostrophe()`, `strip_racemic()`, the `GREEK`
transliteration map and `_CODE_NOISE` — lives only in a per-run script and does **not**
survive the deletion of those directories. Promoting the ladder is the obvious next piece
of this work, and until it happens this section is a specification to re-implement from,
not a description of code in the package.

A name ladder tries query forms most-faithful-first. What it must **never** do is widen
the query, because a widened query returns a confident wrong answer rather than no
answer.

- **Salt words are never stripped.** `PIK-75 HCl` and `PIK-75` are different structures;
  dropping the counterion manufactures a disagreement that is our own doing.
- **Stereodescriptors are never stripped as aliases.** An early draft turned
  `(R)-(-)-Apomorphine hydrochloride` into `Apomorphine hydrochloride`, which resolves to
  the racemate. Alias parentheticals like `(CAY10572)` come off; `(R)-` does not.
- **Racemic markers come off last and are flagged.** PubChem 404s on `±` outright:
  `(±)-Blebbistatin` is not found, `Blebbistatin` is CID 3476986. Dropping the marker
  widens "the racemate" to "unspecified", so those rungs are tried last, the row is marked
  `stereo_widened`, and the CAS still has to agree independently.
- **A Greek letter is as unsearchable as the `?` that replaced it.** PubChem 404s on
  `Pifithrin-μ` and resolves `Pifithrin-mu`. *Correcting* mojibake to the proper character
  therefore **breaks** the lookup unless the transliteration is offered as well — which
  cost a row that had previously resolved.
- **Fragment rungs are fenced hard.** Some cells hold two synonyms jammed together:
  `ABT 263 Navitoclax`, `Actinomycin D Dactinomycin`. Whitespace prefixes and suffixes are
  offered as last rungs — but **barred outright when the name mentions a counterion or a
  stereodescriptor**. Chopping `Doxorubicin hydrochloride` offers `hydrochloride`, which
  PubChem happily resolves to hydrochloric acid, and `Doxorubicin`, which is the free base.
- **A lost prime is repaired, a well-formed locant is not.** `4 -DEMETHYL...` is
  `4'-DEMETHYL...`: a digit followed by whitespace then a hyphen is not valid nomenclature
  in any system. Rewriting `2-Chloroadenosine` into `2'-Chloroadenosine` would invent a
  different compound, so only the whitespace form is touched.

Any rung that widens or guesses is recorded on the verdict, and both rows that won on a
fragment were confirmed independently by their CAS.

**PubChem 404s more than you expect,** and not always because the compound is unknown.
`Ivermectin` is a B1a/B1b mixture with no single-structure compound record at all. That is
a fact about the substance, not a lookup failure, and it stays unverified with that as the
stated reason.

---

## Four routes into ChEBI, and why absence needs all four

| Route | What it asks |
|---|---|
| A | PubChem's own ChEBI cross-reference |
| B | UniChem, InChIKey → ChEBI |
| C | **ChEBI's Elasticsearch index, queried with the raw InChIKey** |
| D | ChEBI by compound name — **for reason codes only** |

**Route A alone is not enough.** PubChem's cross-reference reached 59% of rows in the
first batch. Route C found an exact ChEBI match for **64 rows** route A missed, including
Auranofin, Cisplatin, Bestatin and Gemcitabine hydrochloride. Coverage varies a lot by
library: the same cross-reference reached **91%** of the second batch's resolved CAS
numbers (106 of 116), and 95% of its rows counting a hit on either the CAS or the name
side. That batch was mostly well-characterised clinical drugs.

**Route B contributed nothing and can be dropped.** It was asked about all 389 structures
routes A and C left empty and returned a ChEBI ID for **zero** of them. That is a clean
negative rather than a silent failure — 581 UniChem queries succeeded overall and 186
returned an ID, just none among the ones that mattered. Keep it only as independent
confirmation of an absence, and run it *only* where the cheap routes came up empty:
UniChem answers in about a second when it answers, but fails intermittently with a 30s
error whose retry costs roughly two minutes per structure.

**Route D never verifies anything.** A structure search cannot find a ChEBI entry that
carries **no structure**, and ChEBI holds plenty of those. Without asking by name, all of
them would be reported as "not in ChEBI", which is a different and false claim. It
reclassified 14 of 292 apparent absences into two more accurate buckets. A name match is
never a verdict.

### An outage must never read as an absence

A timeout silently becoming "no ChEBI entry exists" is the worst bug this method can
have, because it publishes a network hiccup as a curation finding. Per-query failures
return an explicit `UNREACHABLE` sentinel, are **never cached**, and are recorded as "not
asked" — but if more than about 10% of a batch fails, the run **aborts**, because that is
an outage and reporting a thousand unchecked absences is worse than stopping.

---

## Reading stereochemistry from the InChI layers

`structure_check.inchi` implements this, and it is the correction of a comparison that
shipped wrong.

Salt and stereo differences are discriminated from the **InChI layers**, not from
PubChem's desalted-parent lookup. `refine_skeleton_difference` costs three requests per
structure and misclassified **15 of 36** rows; comparing the principal component of the
formula layer got all 36. Note the plumbing limit: that function receives InChI**Keys**,
and a formula cannot be recovered from a 14-character connectivity hash, so the cheap path
is available only to callers holding the InChI string.

**Enumerate the layers you are willing to ignore, not the ones you remember to check.**
This is the single most expensive lesson here, and it had to be learned **twice**.
An earlier version compared `/t` and `/m`. It therefore could not see `/b` double-bond
geometry at all, and — because `/t`-equal-but-`/m`-opposite fell through to the wrong
branch — read a pair of **enantiomers** as a cataloguing duplicate. `/t` carries relative
parity and `/m` is the mirror flag; same `/t` with opposite `/m` *is* an enantiomer pair.

**Eight curated identifiers were wrong as a result**, and three were the exact
substitution the strict policy was written to prevent, with the policy right and the code
not: `Formoterol hemifumarate` — a racemate — had been given **arformoterol**'s ID, the
(R,R) enantiomer; `Nebivolol hydrochloride` had the **(R,S,S,S)** single-diastereomer
entry; `Dihydrosphingosine` had the *erythro* entry while its name resolves to
sphinganine. `CCMI` was a genuine E/Z pair verified as one compound. Two more were
improvements once anchoring was fixed — `PNU 74654` and `SU 5416` moved off
systematic-name duplicates onto ChEBI's named drug entries — and one row was rescued from
unverified.

The replacement then made the same mistake one level out. Having widened the stereo
comparison to all four of `/b`, `/t`, `/m` and `/s`, it still compared *only* the formula
layer and those four — so it never read `/c`, the connectivity table, or `/h`, the hydrogen
positions. Two molecules sharing a formula and differing only in constitution therefore
compared **identical**: ethanol and dimethyl ether, n-butane and isobutane, catechol and
hydroquinone. So did an InChI whose `/h` layer a spreadsheet cell had truncated, against
the intact string. That is the worst direction a structure comparison can fail in, and it
was worst precisely in this module's intended slot — the InChIKey path only asks for a
finer comparison once it already knows connectivity differs.

Three more layers were invisible for the same reason, all of which leave the formula layer
as the neutral parent's: `/i` (ethanol against ethanol-d6, a separate ChEBI entry with its
own CAS), and `/q` and `/p` (acetic acid against acetate). The comparison is now tiered —
principal component, then formula, then constitution, then ionisation, then isotope, then
stereo — so a difference can never be misfiled as the wrong *kind* of difference. Replayed
over all 2071 real name/CAS structure pairs from both batches, it reports **zero** pairs
with differing InChIKeys as identical.

What surfaced the first bug was **phase 3 and phase 5 disagreeing about the same pair**. Both compare
the same two structures with independently written code, and that redundancy is a
cross-check, not waste. It only pays if someone diffs the two answers.

**Silence must never render as agreement.** `classify_pair("", "")` returned `"identical"`
in the old version and escaped consequences only because its single caller guarded the
call site. Missing input now returns `not_comparable`.

**cisplatin and transplatin share an InChIKey** — standard InChI does not encode
square-planar geometry. Structure matching is blind there and only a name tiebreak can
decide. Do not remove that tiebreak as redundant with the structure comparison; it is the
one place where a name comparison does real work.

### One precisely bounded limit

`principal_component` ranks components by **heavy-atom count, then total atoms** — not by
formula-string length, and not by total atoms alone. Both weaker keys fail on real drugs.
String length ties on exactly the compound this is most needed for: pyrvinium pamoate is
`2C26H28N3.C23H16O6`, whose components are both eight characters once the multiplier is
dropped, so a string tiebreak returns the *pamoate counterion* — the same wrong answer
PubChem's desalting gives. Across 407 multi-component structures from a real run those two
rankings disagree on three. Total atoms alone then lets hydrogen outvote everything:
diatrizoate meglumine is `C11H9I3N2O4.C7H17NO5`, where the hydrogen-rich counterion wins
30 atoms to 29 while the iodinated parent wins 20 heavy atoms to 13 — so the salt compared
as a different compound from its own free acid.

For a **coordination complex or an organometallic** that InChI has split into fragments —
`6ClH.14H3N.2O.3Ru` for ruthenium red, `C6H5.C2H4O2.Hg` for phenylmercuric acetate — no
component is the parent and the answer is meaningless whichever key is used. Structure
*comparison* is unaffected, because the same ranking is applied to both sides and only the
results are compared.

---

## Choosing among several ChEBI entries for one structure

ChEBI legitimately holds more than one entry per structure, for three different reasons,
and one tiebreak will not do:

- **Zwitterion pairs** — `glycine` / `glycine zwitterion` share a standard InChIKey. The
  neutral entry is the compound; the zwitterion exists for reaction contexts.
- **Bulk-import duplicates beside a curated entry** — `CGS 15943` at 3 stars next to its
  systematic-name twin at 2.
- **Genuinely different molecules InChI cannot separate** — cisplatin and transplatin.

Resolve by: drop zwitterions → prefer the higher star rating → prefer an exact name match
→ prefer the lowest accession. If candidates still tie *and* their SMILES differ, InChI has
hidden a real difference and the row goes to a human rather than to a coin flip.

---

## Is a run trustworthy?

Two cross-checks, both run over the **whole** verified set rather than a sample, so the
result is a measured error rate and not an estimate.

| Signal | Why it is independent |
|---|---|
| ChEBI's own CAS Registry Numbers vs the row's | sourced from ChemIDplus and KEGG, not from PubChem — a different path entirely |
| the Lattice DB's existing ChEBI terms | chosen by a human curator, never read as input |

Measured: **371/371** consistent on the first batch, **450/450** on the second. And on
every compound where a curator had previously chosen an ID for the Lattice DB, the method
independently chose the same one — **48/48**, zero disagreements. That is the strongest
available evidence that the strict policy is calibrated rather than merely severe.

**A different registry number is not a different compound.** CAS issues more than one
number per substance and ChEBI may file the entry under the other one. Resolve ChEBI's
number and compare structures before calling it a finding: **17** apparent disagreements
in the first batch were all this, and 18 in the second. (The first batch originally
reported six; re-running the audit after the stereo corrections raised it, because three
rows moved out of the verified set and others moved in.)

### What the CAS cross-check cannot see

Stated against its own interest, because the 100% figure invites more confidence than it
earns: **a CAS cross-check validates which compound, not which stereoisomer.** CAS numbers
for stereo variants are frequently conflated in ChemIDplus and KEGG, the very sources
ChEBI imports from. Of the eight wrong identifiers above, six had `chebi_registers_no_cas`
— nothing to compare — and the two that "agreed" agreed on the skeleton. The audit was
sound and blind on precisely the axis where the defect lived. A second axis is not a
substitute for getting the comparison right.

---

## Cost

Phases 2 and 3 spend about four PubChem requests per row. PubChem's limits are 5 req/s and
400 req/min; three worker threads sit inside both, sustaining roughly 60 rows/min, so
budget about 40 minutes each for a cold 2500 rows. Phase 4 is faster. Run them backgrounded
and wait on the **process**, not on a log file — and never `until grep` a log you have not
truncated, or you will match the previous run's output and believe a job has finished when
it has not.

Caches are keyed by the literal query string and store negatives too, because "PubChem
does not know this CAS" is a finding worth not re-deriving. They are worth keeping between
batches: two real runs left about 7000 paid-for answers behind.

---

## What is deliberately out of scope

**Deciding a conflict.** Where a name and a CAS disagree, the method reports both
candidates and stops. The most valuable thing to get from the experimental team is the
**vendor catalogue number**, because it identifies a physical product where a name does
not. It corrected conclusions reached here twice, and caught a salt-form error nothing
else would have: `(R)-CR8` is the *trihydrochloride*, a third value different from both the
sheet's number and the one recovered from PubChem.

Their answers are authority about the vial, not about chemistry. One answer read
`Lylamine hydrochloride was used`, repeating the sheet's own misspelling — but the
catalogue number in the same reply, Tocris 2139, is *Leelamine* hydrochloride. The
catalogue number wins, because a name identifies only what someone typed.

**Submitting to ChEBI.** A compound whose identity is settled and which ChEBI does not
hold is a worklist entry, not a dead end — 352 and 361 across the two batches, of which
roughly 45% are incremental (ChEBI already has the parent; add the salt). Nothing is ever
sent outward automatically, and a row whose **identity** is unsettled is never a submission
candidate: submitting a structure we are unsure of would put our uncertainty into a public
database.

---

## Testing

Offline by default, from the `bcp` directory:

```bash
cd bcp
pytest tests/test_structure_check_inchi_layers.py tests/test_structure_check_cas_validation.py
```

Every lesson above that is enforceable is a test rather than a paragraph, and the fixtures
are the **real compounds that failed** — CCMI's opposing `/b` layers, dihydrosphingosine's
`/m` flip, formoterol's `/m` flip inside a multi-component salt, pyrvinium pamoate's tied
component lengths, ABT-702 as its mono- and di-hydrochloride, `2113-05-05 00:00:00`,
`6857789`. Replayed against the implementation they replaced, **12 of 14** structure
assertions fail — and the component-ranking assertions fail against *both* predecessors,
each of which got one real drug right and the other wrong: string length picked
diatrizoate correctly and pyrvinium's counterion, total atoms picked pyrvinium correctly
and meglumine.

This matters more than it sounds. The first batch's own record already contained the
sentence *"Stereochemistry: strict. Racemate-for-enantiomer is never acceptable."* It was
written, it was correct, and the code violated it for eight rows anyway. **Prose cannot
enforce an invariant.** Where a rule here can be expressed as an assertion, it belongs in
`tests/`, and what is left in this document is the part that genuinely cannot be.

Live API tests are opt-in and separate: `pytest -m pubchem` and `pytest -m chebi`.
