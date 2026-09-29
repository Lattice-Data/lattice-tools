# fastq_chunker

Split very large gzipped FASTQs that live in one S3 bucket into record-aligned
chunks below the SRA 100 GB per-file limit, writing the chunks straight to a
second bucket. Nothing is written to local disk. Mates of a set (read1, read2,
read3, index1, index2) are split independently with the same reads-per-chunk
value, so chunk *k* of every mate holds the same reads.

The data path per file is one shell pipeline (values abbreviated; `run
--dry-run` prints the exact commands):

    s3io get --concurrency 8 SRC \
    | pigz -dc                                  # or: rapidgzip -d -c -P T
    | split -l LINES --numeric-suffixes=1 -a 3 \
        --filter='s3io put "DST/STEM.$FILE.fastq.gz" --upload-concurrency 4 --compress pigz -c -p T -6' \
        - part

Python sits only at the two S3 endpoints: `get` fetches ranges of the source
object with several requests in flight, and `put` runs the compressor and
uploads its output as a multipart upload with several parts in flight. Record
counting is GNU `split -l` (a chunk boundary is always a line boundary, and
FASTQ records are four lines), inflation is `pigz` or `rapidgzip`, compression
is `pigz`. No Python code touches the FASTQ bytes in between. Measured on a
96-core pod, a 500 GB file takes roughly 2 to 3 hours; see Performance below
for what was measured and how the first version's 11 hours were brought down.

## Input

Data portal JSON, as exported from the Lattice portal:

- `sequence_file_set` objects supply the grouping: which files are mates. Only
  the read slots (`read1`, `read2`, `read3`, `index1`, `index2`) are used; CRAM
  slots are ignored, and a set with no read slots is skipped.
- `sequence_file` objects supply `s3_uri`, `file_size` and `read_count`.

Each may be a single object, a JSON list, or a search result carrying `@graph`.
A set only embeds the `@id` of its reads, so both must be supplied; a file that
belongs to more than one set is an error, since the two sets would want it split
at different reads-per-chunk values. Two files with the same basename (a sample
name reused across runs, say) are also an error: every chunk lands under one
destination prefix, so their chunk names would collide.

## Commands

Run from `bcp/`.

```
python -m fastq_chunker plan --file-sets sets.json --files files.json \
    --dst s3://dst-bucket/run42_chunks/ --out plan.json [--target-gb 80] [--names-only]
python -m fastq_chunker run --plan plan.json --workers 4 --decompressor rapidgzip --log-dir logs/
python -m fastq_chunker verify --plan plan.json --level quick
python -m fastq_chunker verify --plan plan.json --level full
python -m fastq_chunker batch --plan plan.json --run-manifest run_manifest.tsv --out batches.tsv
```

`plan` touches no network. A file below the 100 GB limit is legal as-is and is
never split; otherwise the set is cut into `ceil(largest / 80 GB)` chunks and
reads per chunk is rounded up to a multiple of one million. The plan refuses any
set whose estimated largest chunk is not below the limit. The estimate assumes
the chunks compress like the original; on real vendor data (Psomagen, pilot of
2026-09-28) `pigz -6` output was 15 to 22 percent smaller than that, so an
80 GB target lands near 65 GB and `--target-gb 90` still leaves the full
margin. Check the pilot's `size_drift` before relying on that for a new vendor.

`run` preflights (tools present, every source exists with the planned size,
destination writable), skips files whose chunks are already complete, deletes
leftovers before re-splitting a file, and runs `--workers` pipelines at once.
Each pipeline reads its source with `--read-concurrency` (default 8) range
requests in flight, inflates with `pigz -dc` or, with `--decompressor rapidgzip`,
in parallel, and uploads each chunk as a multipart upload with
`--upload-concurrency` (default 4) parts in flight. `put` runs `pigz -c` itself
and aborts the upload if the compressor fails, so a chunk object is never a
truncated success. Each chunk's `<chunk>.md5` sidecar (MD5 computed during
upload) is written only after every chunk of that file has landed, so on resume
"chunk plus sidecar" means a chunk from a completed pipeline. A failed file has
its partial chunks removed and the run continues; the exit code is non-zero and
the `--only ...` command to rerun is printed. If `split` produces a different
number of chunks than the plan names, which is what a wrong portal `read_count`
causes, the file is treated the same way. `--copy-singletons` server-side
copies files that need no splitting into `--dst` as well (no sidecar, since
nothing reads their bytes). Run under `tmux` or `nohup`, and export
`AWS_MAX_ATTEMPTS=10 AWS_RETRY_MODE=adaptive` for long jobs. A pipeline that is
killed outright (not one that fails) can leave an incomplete multipart upload
behind, which S3 bills for and does not list; give the destination bucket a
lifecycle rule that aborts incomplete multipart uploads after a day.

`verify --level quick` checks every chunk exists below the limit, flags size
drift (more than 15 percent over the estimate, or more than 30 percent under),
checks the sidecar, decompresses the first record of each chunk for sanity, and
confirms the first read ID of chunk *k* is identical across mates. `--level
full` costs about as much as the run: it decompresses every chunk, checks
`lines % 4 == 0`, that the read count equals the plan, that the MD5 matches the
sidecar (a missing sidecar is a `FAIL`, not a crash), and that first and last
read IDs agree across mates and do not repeat across chunk boundaries. Any
`FAIL` means re-running that file with `--force`. `verify` needs `pigz` only.

`batch` is organisational only: it packs whole sets, in plan order, into
submission batches under `--batch-limit-gb` (default 5000), using actual chunk
sizes from the run manifest when given.

## Performance

Measured in the JupyterHub pod (96 cores) on one real Psomagen set: 20 GB
across R1/R2/I1/I2, 274 M reads, largest file (R2) 9.9 GB, four pipelines in
parallel, 2026-09-28/29. Each row adds one change to the previous:

| change | set wall | R2 | R2, compressed | 500 GB file |
|---|---|---|---|---|
| first pilot (s3fs streaming, `pigz -dc`, 8 threads) | 13 m 23 s | 791 s | 12.5 MB/s | ~11 h |
| parallel ranged reads (`--read-concurrency 8`) | 8 m 43 s | 514 s | 19 MB/s | ~7 h |
| `--decompressor rapidgzip` | 8 m 32 s | 498 s | 20 MB/s | ~7 h |
| concurrent multipart uploads (`--upload-concurrency 4`) | 7 m 17 s | 427 s | 23 MB/s | ~6 h |
| `--pigz-threads 16` | 3 m 49 s | 223 s | 44 MB/s | ~3 h |

What each step removed:

- s3fs streams `open("rb")` one range at a time on one connection: 32 MB/s
  alone, under 10 MB/s with three other streams. `get` fetches 64 MiB ranges
  with several in flight; peak memory per `get` is about `--read-concurrency`
  blocks (512 MiB at the default).
- Single-threaded `pigz -dc` inflates at ~340 MB/s uncompressed; `rapidgzip`
  did 640 MB/s from disk and works from a pipe. It matters most for the small,
  highly compressible index reads (236 s to 81 s). `verify --level full`
  always inflates with `pigz`, so the check runs a different decoder than the
  run.
- s3fs uploads each part inside `write()`, stalling the whole pipeline behind
  it for the duration (62 MB/s per stream, measured with a 2 GB probe), so
  compute and upload ran serially. `put` uploads parts through boto3 with
  `--upload-concurrency` in flight; compression continues while parts upload.
- After those three, the pipeline is compression-bound: 150 bp reads with
  quality strings compress at ~19 MB/s uncompressed per pigz thread at level 6,
  and doubling the threads halved R2's time. `run`'s default is
  `cores / workers - 1` threads (23 on that pod with 4 workers), so do not pass
  `--pigz-threads` there unless you want fewer; expect roughly 2 to 2.5 hours
  per 500 GB file, and a quadruple with two 500 GB mates in about the same wall
  time since all four run at once. `--gzip-level 5` buys another ~25 percent
  for ~3 percent larger chunks.

Every run above produced byte-identical chunks (same sizes and MD5s as the
first pilot), and quick and full verify passed on each.

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
