# Protospacer Set Signature: implementation spec

Implement this as described. The fixtures in `tests/fixtures/` and the golden digests
below are the contract. Do not modify the fixtures. Do not change the algorithm to make
a failing test pass; if a golden digest does not reproduce, the implementation is wrong.

## Purpose

Given a guide-template TSV, compute a deterministic string that is identical for any two
files containing the same **set** of protospacer sequences, and different otherwise. Row
order and all non-protospacer columns must not affect the result.

## Scope boundary (important)

This is a string-identity check on one column. It is **not** a biological equivalence
check.

Do NOT normalize sequence length. Do NOT strip a leading or trailing G. Do NOT
reverse-complement. Do NOT canonicalize strand or orientation. Do NOT map to genomic
coordinates. Two files whose protospacers differ by a trimmed base or a strand flip
contain different strings and MUST produce different signatures.

If you believe a normalization step is needed that is not in the algorithm below, stop
and raise it as a question rather than implementing it. Every such step is a way to
produce a false match, and a false match is the one failure mode this tool exists to
prevent.

## Algorithm

Exactly these steps, in this order. No additions.

1. Read the TSV with `encoding='utf-8-sig'` and `delimiter='\t'`. Files carry a UTF-8 BOM
   and use CRLF line endings. Use `newline=''` with the `csv` module so both are handled.
2. The first row is the header. Drop every subsequent row whose **first column** value,
   after `lstrip()`, starts with `#`. In this template the first column is named
   `property` and the rows `#comment` and `#example` must be excluded. Filter on
   first-column *position*, not on the literal header name `property`, and not on the
   literal strings `#comment` / `#example`.
   - Amendment (documented deviation from the original algorithm): a **fully blank
     row** -- one the `csv` reader yields as `[]` or as a list whose every cell is the
     empty string -- is skipped, not treated as an error. This covers the trailing empty
     line that Excel/Google Sheets exports routinely append, and short rows made entirely
     of dropped trailing tabs. A row that is short *and* has any non-empty cell is still a
     ragged row and must raise (see Errors). This is the only addition to the literal
     step list; it is captured here so the spec stays authoritative.
3. From each remaining row take the protospacer column. Default name
   `guide_protospacer`, overridable by the caller.
4. Apply `.strip()` then `.upper()` to the value.
5. Skip values that are empty after stripping.
6. Reject the file if any remaining value contains a character outside `{A, C, G, T}`.
   Raise. Do not coerce, do not warn and continue.
7. Collect into a `set` (deduplicate).
8. Sort lexicographically (plain Python `sorted()` on the uppercase ASCII strings).
9. `payload = '\n'.join(sorted_sequences)` -- LF only, no trailing newline.
10. `digest = sha256(payload.encode('utf-8')).hexdigest()[:32]`
11. Return `f'gsig1:set:n={len(sequences)}:{digest}'`

The `gsig1` prefix and the `n` field are part of the output contract. `gsig1` exists so a
future change to these rules cannot silently make old and new signatures comparable. `n`
exists so a mismatch is diagnosable without re-parsing the files.

## Errors

Raise, with the file path and the 1-based physical row number where applicable:

- protospacer column absent from the header
- non-ACGT character in a data value
- zero usable sequences after filtering (empty set) -- never return a signature for this
- file unreadable, not TSV-parseable, or not decodable as UTF-8 (a file saved as
  Windows-1252/Latin-1 raises rather than propagating a raw `UnicodeDecodeError`)
- a ragged row (shorter than the header) that is not a fully blank row, named by its
  1-based physical row number

## Interface

Library:

```python
signature(path: str | Path, column: str = 'guide_protospacer') -> str
protospacer_set(path: str | Path, column: str = 'guide_protospacer') -> set[str]
```

`protospacer_set` is exposed so callers can diff without reimplementing the parser.

CLI:

```
guidesig FILE [FILE ...] [--column NAME]
    one line per file: <signature>\t<path>

guidesig --compare FILE_A FILE_B [--column NAME]
    reports match or mismatch; on mismatch prints |A|, |B|, |A n B|, |A - B|, |B - A|
    and up to 10 example sequences from each difference
```

Exit non-zero on error, and on `--compare` mismatch.

Stdlib only: `csv`, `hashlib`, `argparse`, `pathlib`. No third-party runtime dependencies.
Follow the repo's existing package layout and test conventions if present; otherwise a
single module plus `tests/` is fine.

## Fixtures

Four synthetic files in `tests/fixtures/`. All sequences are randomly generated. They
carry no real library content, no real gene identifiers, and no real sample or order
identifiers. They reproduce the structural properties of real submissions:

| file | data rows | unique protospacers | properties exercised |
|---|---|---|---|
| `lib_one.tsv` | 63 | 60 | 19-mers; 3 duplicate protospacers under different `guide_id`s; one sequence and its exact reverse complement both present as separate guides; BOM; CRLF; `#comment` and `#example` rows |
| `lib_one_prepended_g.tsv` | 66 | 60 | same 60 sequences as `lib_one` with a `G` prepended (20-mers); different `guide_id`s; different row order; 6 duplicate rows |
| `lib_two_layout_a.tsv` | 60 | 60 | 8 columns; mixed 19-nt and 20-nt protospacers within one file |
| `lib_two_layout_b.tsv` | 61 | 60 | same protospacer set as `layout_a`; 15 columns in a different order; an extra `guide_rc_sequence` column; some values lowercase; some values padded with spaces; one row with an empty protospacer |

## Golden values

```
gsig1:set:n=60:62d08f4954bf50b7249c6277e2808b16  lib_one.tsv
gsig1:set:n=60:8a3af4c57186c0766609fd1841136e7c  lib_one_prepended_g.tsv
gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78  lib_two_layout_a.tsv
gsig1:set:n=60:0ae925412d48830d0e5bee973c24ea78  lib_two_layout_b.tsv
gsig1:set:n=60:429cecb0d8b27564cddc32dfde6972ff  lib_two_layout_b.tsv --column guide_rc_sequence
```

Three assertions carry most of the design intent:

- **`lib_two_layout_a` and `lib_two_layout_b` MUST match** (`0ae92541...`). 8 columns vs
  15, different column order, different row order, mixed case, padded whitespace, and an
  empty row, yet identical. This is the end-to-end proof that layout, order, case, and
  padding are irrelevant.
- **`lib_one` and `lib_one_prepended_g` MUST NOT match**, despite both having `n=60`.
  They differ only by a 5' G. If a change ever makes these match, length normalization
  was introduced in error. Note that equal `n` with unequal digest is the exact case a
  count-only check would miss.
- **`lib_two_layout_b` with `--column guide_rc_sequence` MUST differ** from the same file
  read with the default column. Confirms `--column` is honored and that a
  reverse-complement column is not interchangeable with the protospacer column.

## Tests

pytest. Build small synthetic inputs in `tmp_path` for the unit tests; use the fixtures
for the regression tests.

Invariance:
- row order permuted -> identical signature
- `guide_id`, target gene, PAM, coordinates all changed, protospacers held constant -> identical
- columns added, removed, or reordered in the header -> identical
- a protospacer present N times vs present once -> identical (set semantics)
- lowercase input -> identical to uppercase input
- leading/trailing spaces, and a stray `\r` on the value -> identical to clean input
- signature stable across two calls in the same process and across a re-read

Sensitivity:
- one sequence changed -> different signature
- one sequence added -> different signature, `n` increments by 1
- 19-mers vs the same sequences with a `G` prepended -> different signature
- a sequence and its reverse complement both present -> both retained, `n` counts both

Row filtering:
- `#comment` and `#example` sequences absent from the resulting set
- a file whose only data row is `#example` -> raises the empty-set error, and does not
  return a signature derived from the example sequence
- first column `# example` (space after the hash) -> still filtered
- a data row whose first column is empty -> retained, not filtered

Parsing:
- BOM present -> parses, and header lookup for the first column still works
- CRLF line endings -> no `\r` reaches the sequence set
- empty protospacer cell -> skipped, not an error
- ragged row shorter than the header -> raises, with the row number

Validation:
- value containing `N` -> raises, message names the row number
- value containing `U`, a digit, or a hyphen -> raises
- missing protospacer column -> raises

Serialization pinning:
- construct a two-sequence input (`AAAA`, `CCCC`), and assert the digest equals
  `sha256('AAAA\nCCCC')` computed inline in the test. This pins the join character and
  the no-trailing-newline rule rather than only checking self-consistency.

Regression:
- all five golden values above reproduce exactly
- the three named assertions above, written as their own explicit tests with comments
  stating why they matter

CLI:
- multi-file invocation prints one `<signature>\t<path>` line per file
- `--compare` on the `lib_two` pair reports a match and exits 0
- `--compare` on the `lib_one` pair reports a mismatch, prints the set-difference
  summary, and exits non-zero

## Deliverables

Module, CLI entry point, test suite, and a README section covering what the signature
does and does not guarantee, why length and orientation are deliberately not normalized,
and why the version prefix exists. Run the full suite and report results.
