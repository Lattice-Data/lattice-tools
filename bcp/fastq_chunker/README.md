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
copies files that need no splitting into `--dst` as well, `--copy-workers`
files (default 4) at a time, and preflight then checks their sources too. On
S3 each file of 64 MiB or more is copied as 128 MiB parts, `--copy-parts`
(default 32) of them at once. One copy request moves only a few MB/s, so the
parts are what make it fast: in the TOOLS-382 pilot a 13.4 GB file took 48
minutes one part at a time and 2 minutes with 32 at once. Each copy
is checked after it lands: its size must match the plan, and when DB2
recorded a CRC64NVME for it (`crc64nvme_base64`), S3 must hold a full-object
CRC64NVME equal to it. A copy S3 has no such checksum for (one uploaded some
other way, or with a COMPOSITE checksum) cannot be compared, so it counts as
wrong: a rerun copies it again and the gate fails it. Only without a DB2 CRC,
or off S3, is the size alone compared. The manifest's `checked` column says
which was compared (`size` or `size+crc64nvme`). A failed copy is reported and
the other copies go on; what is at its target is removed only once a check
found it wrong, never because a request failed, so a transient error leaves a
good copy in place. A rerun copies only what is missing or wrong; the hint a
failed run prints is the flags to add, `--copy-singletons` included. `--only`
naming a copied file without `--copy-singletons` is refused. A copy
gets no sidecar, since nothing reads its bytes. A plan with nothing to split
only copies, and needs neither `pigz` nor `split`. Run under `tmux` or `nohup`, and export
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

## From db2_flattener's SRA sheets

db2_flattener writes a submittable SRA file sheet per batch
(`sra_file_set{N}.csv`) whose `filename…` cells hold DB2's predicted names
(`sra_file_names` on SequenceFile: a file's own name, or its chunk names once
its set reaches the limit), plus `sra_sequence_fileset_set{N}.csv`
(`library_ID,sequence_file_set_uuid`, one row per SequenceFileSet). The
destination then holds exactly what the sheet names: chunks for a file that is
split (never the original), the original for a file that is not.

```
# where DB2 credentials are (DB2_PROD_KEY, DB2_PROD_SECRET, DB2_PROD_SERVER)
python -m fastq_chunker resolve --sra-sheet sra_file_set1.csv \
    --set-sheet sra_sequence_fileset_set1.csv --out-dir set1/
python -m fastq_chunker plan --file-sets set1/sets.json --files set1/files.json \
    --dst s3://dst-bucket/batch1/ --out set1/plan.json

# where S3 access is (the pod); plan.json carries every s3_uri
python -m fastq_chunker run --plan plan.json --copy-singletons --decompressor rapidgzip
python -m fastq_chunker verify --plan plan.json --level full
python -m fastq_chunker check-submission --plan plan.json \
    --sra-sheet sra_file_set1.csv --verify-report verify_report.tsv
```

`resolve` is read-only. It fetches the listed sets by uuid, and the files in
their read slots, and writes the `sets.json`/`files.json` that `plan` already
reads. It refuses, listing every problem at once, when:

- a uuid does not resolve, or a set or file is deleted;
- a set is listed twice, or a file belongs to two listed sets;
- DB2 gives a file no `sra_file_names` (which it does when it would refuse the
  set);
- a library's `sra_file_names`, across all its sets, differ from that row's
  `filename…` cells (order within a row is not compared);
- the batch is over `--batch-limit-gb` (default 5000).

Credentials follow the repo's `--mode` convention (default `db2_prod`).

`plan` enforces the name contract whenever the file objects carry
`sra_file_names`. Its own names, the chunk names or an unsplit file's
name, must equal DB2's, in order, or nothing is written. So the chunker
produces exactly the predicted names or refuses. A different `--target-gb` is
accepted only where it yields the same names.

`check-submission` is the gate before upload. It is read-only and reads the
destination itself, not the manifest. It requires that:

- the sheet's names are exactly the plan's;
- the destination holds exactly those objects, plus `.md5` sidecars of chunks
  and nothing else (no pre-split original, no probe, no sub-prefix);
- every copy passes the copy check above;
- every chunk is below the limit, has its sidecar, and passed
  `verify --level full` (read count and md5) with no FAIL for it, for its
  file's part (`boundary_distinct`), or for its part across mates
  (`first_id_aligned`, `last_id_aligned`; verify writes these by part label,
  `part001`, not by chunk name);
- the verify report, when given, has no FAIL anywhere and ends in an overall
  PASS, so a row the mapping misses, or a report cut short, still fails;
- the batch in actual bytes is within the limit.

It writes `submission_report.tsv` and exits non-zero on any FAIL.

## Copying a set into the SRA bucket, step by step

This puts one set's files into its own folder in the SRA bucket, so that the
folder holds exactly what the set's SRA sheet names. It does not submit
anything to SRA. Every step runs in a terminal on a Biohub JupyterHub pod. The
pod's role is what may read the vendor buckets and write the SRA bucket, so no
AWS login is needed. Values written `<like-this>` (the bucket, and where DB2
credentials come from) are not in this public repo; get them from the person
who handed you this.

You need:

- **The flattener's two CSVs for the set:** `sra_file_set<N>.csv` and
  `sra_sequence_fileset_set<N>.csv`.
- **DB2 credentials:** a key and secret for production DB2. They are only
  read; nothing is written to DB2.
- **lattice-tools on the pod**, at `~/lattice-tools`, up to date:

  ```bash
  cd ~/lattice-tools && git checkout main && git pull --ff-only && git log --oneline -1
  ```

Keep each set's working files in `~/fastq_batches/`, never inside
`~/lattice-tools`: the sheets and plans carry collaborators' file and library
names, and the repo is public.
Commands below are run one at a time. If you open a new terminal, set the
variables from steps 1 and 3 again (`SET`, `B`, `NAME`, `SRA`) before going on.

**1. Name the set's folder.** Upload the two CSVs with the JupyterHub file
browser into `~/fastq_batches/incoming/`, then set `SET` to the sheets'
`set<N>` and make the folder. Its name is the date, the set, and the first
8 characters of the sheet's sha256, so two versions of a sheet never share a
folder:

```bash
mkdir -p ~/fastq_batches/incoming
```
```bash
SET=set1
NAME="$(date -u +%Y%m%d)-$SET-$(sha256sum ~/fastq_batches/incoming/sra_file_$SET.csv | cut -c1-8)"; B=~/fastq_batches/$NAME; mkdir -p $B && mv ~/fastq_batches/incoming/sra_file_$SET.csv ~/fastq_batches/incoming/sra_sequence_fileset_$SET.csv $B/ && echo $NAME
```

Write the printed name down. Run this step once per set: on a later day it
would make a new name, and a rerun must use the old folder.

**2. Give the shell your DB2 credentials.** Typed at the prompts, they stay
out of the shell history and are gone when the terminal closes. The server is
`<db2-server-url>`:

```bash
read -rp "DB2 server: " DB2_PROD_SERVER; read -rp "DB2 key: " DB2_PROD_KEY; read -rsp "DB2 secret: " DB2_PROD_SECRET; echo; export DB2_PROD_SERVER DB2_PROD_KEY DB2_PROD_SECRET
```

**3. Resolve the set against DB2 and plan it.** Set `SRA` to the bucket, then
run `resolve`:

```bash
SRA=s3://<sra-bucket>
```
```bash
cd ~/lattice-tools/bcp && python -m fastq_chunker resolve --sra-sheet $B/sra_file_$SET.csv --set-sheet $B/sra_sequence_fileset_$SET.csv --out-dir $B
```

It prints `… libraries, … sets, … files, … GB; the SRA sheet matches DB2's
sra_file_names`. If it refuses instead, it lists every problem and writes
nothing. Stop there and send the output to whoever made the sheets: the sheets
and DB2 disagree, and copying would not fix that. Then plan:

```bash
cd ~/lattice-tools/bcp && python -m fastq_chunker plan --file-sets $B/sets.json --files $B/files.json --dst $SRA/$NAME/ --out $B/plan.json > $B/plan_summary.txt && head -1 $B/plan_summary.txt
```

It prints `G groups, S to split; …`. With `0 to split`, every file is copied
as it is. With more, the files of 100 GB or more are split into chunks, and
step 6 is needed.

**4. Check where it will write.** The first line must be `$SRA/$NAME/`. On a
first run the listing after it must be empty; if it shows files, stop, since
something else is already in that folder. On a rerun of the same set it shows
the files copied so far.

```bash
python -c "import json,sys; print(json.load(open(sys.argv[1]))['dst'])" $B/plan.json; aws s3 ls $SRA/$NAME/
```

**5. Start the copy.** It runs in the background, so closing the browser does
not stop it. Copies have run at 13 to 27 GB a minute (a 1.3 TB set took 100
minutes); splitting a 500 GB file takes 2 to 3 hours.

```bash
cd ~/lattice-tools/bcp && AWS_MAX_ATTEMPTS=10 AWS_RETRY_MODE=adaptive nohup sh -c "date; python -u -m fastq_chunker run --plan $B/plan.json --copy-singletons --log-dir $B/logs --manifest $B/run_manifest.tsv; date" > $B/run.log 2>&1 &
```

Check on it whenever you like. This prints the files copied so far, any
failures, and the last line of the log:

```bash
grep -c '^copied' $B/run.log; grep FAILED $B/run.log; tail -1 $B/run.log
```

It is finished when the last line of `run.log` is a date. If any line says
`FAILED`, or the run stopped partway, start step 5 again exactly as written:
files already copied and checked are skipped, and only what is missing or
wrong is redone.

**6. Only if step 3 said anything to split: verify the chunks.** This
decompresses every chunk and checks its read count, its md5, and that mates
line up. It also runs in the background:

```bash
cd ~/lattice-tools/bcp && nohup sh -c "date; python -u -m fastq_chunker verify --plan $B/plan.json --level full --out $B/verify_report.tsv; date" > $B/verify.log 2>&1 &
```

**7. Check the folder.** This is the gate. It reads the folder itself and
changes nothing:

```bash
cd ~/lattice-tools/bcp && python -m fastq_chunker check-submission --plan $B/plan.json --sra-sheet $B/sra_file_$SET.csv --out $B/submission_report.tsv
```

If step 6 ran, add `--verify-report $B/verify_report.tsv`. The set is done
when this prints `PASS: …, 0 failed`, with every copy checked by
`size+crc64nvme`. A `FAIL` line names the file and the check. Most are fixed
by running step 5 again, then this step; anything else, send it on.

**8. Record the summary**: start and end time, files, GB, and how each copy
was checked. Send it with step 7's output to whoever handed you the set:

```bash
grep -v '^nohup' $B/run.log | head -1; tail -1 $B/run.log; awk -F'\t' 'NR>1 {sub(/\r$/, ""); n++; gb+=$3/1e9; c[$7]++; s[$5]++} END {printf "%d files, %.1f GB; status:", n, gb; for (k in s) printf " %s=%d", k, s[k]; printf "; checked:"; for (k in c) printf " %s=%d", k, c[k]; print ""}' $B/run_manifest.tsv
```

Do not create folders in the bucket from the S3 console: the marker it leaves
fails step 7.

**Removing a set's folder**, when it was only a test or has to be redone from
scratch. Deletion is permanent. Check the count and the names first; every
line must start with `(dryrun) delete: $SRA/$NAME/`:

```bash
aws s3 rm $SRA/$NAME/ --recursive --dryrun | wc -l
```
```bash
aws s3 rm $SRA/$NAME/ --recursive --dryrun | head -3
```

Then delete it, and confirm the folder is gone. The listing should not show
`$NAME/`:

```bash
aws s3 rm $SRA/$NAME/ --recursive
```
```bash
aws s3 ls $SRA/
```

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
- copies of unsplit files under their own names, with `--copy-singletons`
- `sets.json`, `files.json` (from `resolve`), `plan.json`, `run_manifest.tsv`,
  `verify_report.tsv`, `submission_report.tsv`, `batches.tsv`, `logs/<stem>.log`

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
