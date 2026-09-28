# fastq_chunker

Split very large gzipped FASTQs that live in one S3 bucket into record-aligned
chunks below the SRA 100 GB per-file limit, writing the chunks straight to a
second bucket. Nothing is written to local disk. Mates of a set (read1, read2,
read3, index1, index2) are split independently with the same reads-per-chunk
value, so chunk *k* of every mate holds the same reads.

The data path per file is one shell pipeline:

    s3io get SRC | pigz -dc | split -l LINES --filter='pigz -c | s3io put DST/STEM.$FILE.fastq.gz' - part

Python only sits at the two S3 endpoints (through `s3fs`); counting, decompression
and compression are GNU `split` and `pigz`. That is what makes a 500 GB file take
about an hour rather than a day.

## Input

Data portal JSON, as exported from the Lattice portal:

- `sequence_file_set` objects supply the grouping: which files are mates. Only
  the read slots (`read1`, `read2`, `read3`, `index1`, `index2`) are used; CRAM
  slots are ignored, and a set with no read slots is skipped.
- `sequence_file` objects supply `s3_uri`, `file_size` and `read_count`.

Each may be a single object, a JSON list, or a search result carrying `@graph`.
A set only embeds the `@id` of its reads, so both must be supplied; a file that
belongs to more than one set is an error, since the two sets would want it split
at different reads-per-chunk values.

## Commands

Run from `bcp/`.

```
python -m fastq_chunker plan --file-sets sets.json --files files.json \
    --dst s3://dst-bucket/run42_chunks/ --out plan.json [--target-gb 80] [--names-only]
python -m fastq_chunker run --plan plan.json --workers 4 --pigz-threads 8 --log-dir logs/
python -m fastq_chunker verify --plan plan.json --level quick
python -m fastq_chunker verify --plan plan.json --level full
python -m fastq_chunker batch --plan plan.json --run-manifest run_manifest.tsv --out batches.tsv
```

`plan` touches no network. A file below the 100 GB limit is legal as-is and is
never split; otherwise the set is cut into `ceil(largest / 80 GB)` chunks and
reads per chunk is rounded up to a multiple of one million. The plan refuses any
set whose estimated largest chunk is not below the limit.

`run` preflights (tools present, every source exists with the planned size,
destination writable), skips files whose chunks are already complete, deletes
leftovers before re-splitting a file, and runs `--workers` pipelines at once.
Each chunk gets a `<chunk>.md5` sidecar written during upload. A failed file has
its partial chunks removed and the run continues; the exit code is non-zero and
the `--only ...` command to rerun is printed. `--copy-singletons` server-side
copies files that need no splitting into `--dst` as well (no sidecar, since
nothing reads their bytes). Run under `tmux`, and export
`AWS_MAX_ATTEMPTS=10 AWS_RETRY_MODE=adaptive` for long jobs.

`verify --level quick` checks every chunk exists below the limit, flags size
drift beyond 15 percent of the estimate, checks the sidecar, decompresses the
first record of each chunk for sanity, and confirms the first read ID of chunk
*k* is identical across mates. `--level full` costs about as much as the run:
it decompresses every chunk, checks `lines % 4 == 0`, that the read count equals
the plan, that the MD5 matches the sidecar, and that first and last read IDs
agree across mates and do not repeat across chunk boundaries. Any `FAIL` means
re-running that file with `--force`.

`batch` is organisational only: it packs whole sets, in plan order, into
submission batches under `--batch-limit-gb` (default 5000), using actual chunk
sizes from the run manifest when given.

## Outputs

- `<stem>.partNNN.fastq.gz` and `<stem>.partNNN.fastq.gz.md5` under `--dst`
- `plan.json`, `run_manifest.tsv`, `verify_report.tsv`, `batches.tsv`, `logs/<stem>.log`

## Requirements

`pigz` and GNU `split` with `--filter` (coreutils 8.13+). On macOS:
`brew install pigz coreutils` (the latter provides `gsplit`, which `run` finds
on its own). Python deps are in `bcp/requirements.txt`; S3 credentials come
from the environment (IRSA on the JupyterHub pods).

## Tests

`tests/test_fastq_chunker_*.py`. Unit tests need nothing. The S3 endpoint tests
run against an in-process `moto` server. The integration test drives the real
pipeline over `file://` URLs on a synthetic quadruple (`test-data`), then checks
the ground truth independently of the tool: concatenating a file's chunks
decompresses to the original byte for byte, and every sidecar matches its chunk.
It skips, with a message, when `pigz` or GNU `split` is missing.

## Known limitations

- Assumes four-line FASTQ records; `verify` catches a violation after the fact.
- Assumes mates hold identical read counts in identical order; `verify` checks
  first and last IDs per chunk, not every read.
- Resume granularity is a whole file.
