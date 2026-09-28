"""End to end on synthetic FASTQs over file:// URLs: test-data, plan, run, verify, batch.

Needs pigz and GNU split (``brew install pigz coreutils`` on macOS); skipped
with a message when either is missing. The ground truth is independent of the
tool: concatenating a file's chunks must decompress to the original, byte for
byte, and every sidecar MD5 must match the chunk it names.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import stat
from pathlib import Path

import pytest

from fastq_chunker.cli import main
from fastq_chunker.plan import load_plan
from fastq_chunker.run import RunError, find_tools

try:
    TOOLS = find_tools()
except RunError as e:
    TOOLS = None
    SKIP_REASON = str(e)

pytestmark = pytest.mark.skipif(
    TOOLS is None, reason=(SKIP_REASON if TOOLS is None else "")
)

READS = 6_000
TARGET = 60_000
LIMIT = 75_000


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    src = (tmp_path / "src").as_uri() + "/"
    dst = (tmp_path / "dst").as_uri() + "/"
    assert (
        main(
            [
                "test-data",
                "--reads",
                str(READS),
                "--out-prefix",
                src,
                "--meta-dir",
                "meta",
            ]
        )
        == 0
    )
    assert (
        main(
            [
                "plan",
                "--file-sets",
                "meta/sets.json",
                "--files",
                "meta/files.json",
                "--dst",
                dst,
                "--target-bytes",
                str(TARGET),
                "--limit-bytes",
                str(LIMIT),
                "--round-to",
                "1",
                "--out",
                "plan.json",
            ]
        )
        == 0
    )
    return tmp_path


def run_args(*extra):
    return [
        "run",
        "--plan",
        "plan.json",
        "--workers",
        "4",
        "--pigz-threads",
        "2",
        *extra,
    ]


def manifest_statuses(path="run_manifest.tsv"):
    rows = [line.split("\t") for line in Path(path).read_text().splitlines()[1:]]
    return {(r[0], r[1]): r[4] for r in rows}


def test_end_to_end(workspace, capsys):
    plan = load_plan(Path("plan.json"))
    (g,) = plan.groups
    assert g.action == "split" and g.n_chunks > 1
    src_dir, dst_dir = workspace / "src", workspace / "dst"

    assert main(run_args()) == 0
    statuses = manifest_statuses()
    assert set(statuses.values()) == {"ok"}
    assert len(statuses) == 4 * g.n_chunks

    # independent ground truth: chunks concatenate back to the original
    for fp in g.files:
        original = gzip.decompress((src_dir / fp.file.filename).read_bytes())
        joined = b"".join((dst_dir / c.name).read_bytes() for c in fp.chunks)
        assert gzip.decompress(joined) == original
        for c in fp.chunks:
            data = (dst_dir / c.name).read_bytes()
            assert len(data) < LIMIT
            assert (
                dst_dir / (c.name + ".md5")
            ).read_text() == f"{hashlib.md5(data).hexdigest()}  {c.name}\n"
            assert gzip.decompress(data).count(b"\n") == 4 * c.reads
    # no stray objects: exactly the planned chunks, their sidecars, nothing else
    assert sorted(p.name for p in dst_dir.iterdir()) == sorted(
        n for fp in g.files for c in fp.chunks for n in (c.name, c.name + ".md5")
    )

    assert (
        main(
            ["verify", "--plan", "plan.json", "--level", "quick", "--out", "quick.tsv"]
        )
        == 0
    )
    assert (
        main(
            [
                "verify",
                "--plan",
                "plan.json",
                "--level",
                "full",
                "--workers",
                "4",
                "--out",
                "full.tsv",
            ]
        )
        == 0
    )
    full = [line.split("\t") for line in Path("full.tsv").read_text().splitlines()]
    checks = {row[3] for row in full[1:-1]}
    assert {
        "read_count",
        "md5",
        "first_id_aligned",
        "last_id_aligned",
        "boundary_distinct",
        "mate_alignment",
    } <= checks
    assert full[-1][4] == "PASS"
    assert all(row[4] != "FAIL" for row in full[1:])

    assert (
        main(
            [
                "batch",
                "--plan",
                "plan.json",
                "--run-manifest",
                "run_manifest.tsv",
                "--batch-limit-gb",
                "1",
                "--out",
                "batches.tsv",
            ]
        )
        == 0
    )
    rows = Path("batches.tsv").read_text().splitlines()
    assert len(rows) == 1 + 4 * g.n_chunks
    assert {r.split("\t")[0] for r in rows[1:]} == {"1"}

    # resume: a second run skips every file
    before = {p.name: p.stat().st_mtime_ns for p in dst_dir.iterdir()}
    assert main(run_args()) == 0
    assert set(manifest_statuses().values()) == {"skipped"}
    assert {p.name: p.stat().st_mtime_ns for p in dst_dir.iterdir()} == before

    # --only with --force touches just that file's chunks
    r1 = next(fp for fp in g.files if fp.file.slot == "read1")
    assert main(run_args("--only", r1.file.filename, "--force")) == 0
    statuses = manifest_statuses()
    assert set(statuses) == {(r1.file.filename, c.name) for c in r1.chunks}
    assert set(statuses.values()) == {"ok"}
    after = {p.name: p.stat().st_mtime_ns for p in dst_dir.iterdir()}
    changed = {n for n in after if after[n] != before[n]}
    assert changed == {n for c in r1.chunks for n in (c.name, c.name + ".md5")}
    for c in r1.chunks:
        assert (dst_dir / c.name).stat().st_size < LIMIT

    # a chunk without its sidecar is not trusted: that file is re-split
    i1 = next(fp for fp in g.files if fp.file.slot == "index1")
    (dst_dir / (i1.chunks[2].name + ".md5")).unlink()
    assert main(run_args()) == 0
    statuses = manifest_statuses()
    assert {s for (f, _), s in statuses.items() if f == i1.file.filename} == {"ok"}
    assert {s for (f, _), s in statuses.items() if f != i1.file.filename} == {"skipped"}
    assert (dst_dir / (i1.chunks[2].name + ".md5")).exists()
    capsys.readouterr()


def test_verify_detects_a_corrupted_chunk(workspace, capsys):
    plan = load_plan(Path("plan.json"))
    (g,) = plan.groups
    assert main(run_args()) == 0
    r2 = next(fp for fp in g.files if fp.file.slot == "read2")
    victim = workspace / "dst" / r2.chunks[1].name
    r1 = next(fp for fp in g.files if fp.file.slot == "read1")
    # swap in R1's chunk 3 under R2's chunk 2 name: valid gzip, wrong reads
    victim.write_bytes((workspace / "dst" / r1.chunks[2].name).read_bytes())
    # and drop the last two records of I1's chunk 4: valid gzip, too few reads
    i1 = next(fp for fp in g.files if fp.file.slot == "index1")
    short = workspace / "dst" / i1.chunks[3].name
    lines = gzip.decompress(short.read_bytes()).splitlines(keepends=True)
    short.write_bytes(gzip.compress(b"".join(lines[:-8])))
    assert (
        main(["verify", "--plan", "plan.json", "--level", "full", "--out", "v.tsv"])
        == 1
    )
    rows = [r.split("\t") for r in Path("v.tsv").read_text().splitlines()[1:]]
    failed = {(r[2], r[3]) for r in rows if r[4] == "FAIL"}
    assert (r2.chunks[1].name, "md5") in failed
    assert ("part002", "mate_alignment") in failed
    assert ("part002", "first_id_aligned") in failed
    assert (r2.chunks[1].name, "header") not in failed
    assert (i1.chunks[3].name, "read_count") in failed
    assert (i1.chunks[3].name, "lines_mod_4") not in failed
    assert ("part004", "last_id_aligned") in failed
    out = capsys.readouterr().out
    assert "FAIL:" in out


def test_failure_leaves_no_partial_chunks(workspace, capsys, monkeypatch):
    # a pigz that works for its first three calls (decompress, chunk 1, chunk 2)
    # and then dies, so two chunks are uploaded before the pipeline fails
    counter = workspace / "pigz_calls"
    monkeypatch.setenv("FAKE_PIGZ_COUNTER", str(counter))
    fake = workspace / "fake_pigz"
    fake.write_text(
        "#!/bin/sh\n"
        'n=$(cat "$FAKE_PIGZ_COUNTER" 2>/dev/null || echo 0)\n'
        'n=$((n + 1)); echo "$n" > "$FAKE_PIGZ_COUNTER"\n'
        'if [ "$n" -gt 3 ]; then exit 3; fi\n'
        f'exec {TOOLS.pigz} "$@"\n'
    )
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    rc = main(run_args("--workers", "1", "--pigz", str(fake), "--log-dir", "logs"))
    assert rc == 1
    assert set(manifest_statuses().values()) == {"failed"}
    dst_dir = workspace / "dst"
    assert not dst_dir.exists() or not any(dst_dir.iterdir())
    out = capsys.readouterr().out
    assert "rerun with: --only" in out
    logs = {log.name: log.read_text() for log in (workspace / "logs").glob("*.log")}
    assert len(logs) == 4
    # the first file got two chunks out before its compressor died
    uploaded = [log for log in logs.values() if log.count(".fastq.gz\t") == 2]
    assert len(uploaded) == 1
    assert "exited 3" in uploaded[0]


def test_dry_run_prints_commands_without_touching_dst(workspace, capsys):
    assert main(run_args("--dry-run")) == 0
    out = capsys.readouterr().out
    assert out.count("set -o pipefail") == 4
    assert "--filter=" in out and "$FILE" in out and "--compress" in out
    assert not (workspace / "dst").exists()


def test_preflight_rejects_stale_plan(workspace, capsys):
    plan = json.loads(Path("plan.json").read_text())
    plan["groups"][0]["files"][0]["file"]["size_bytes"] += 1
    Path("plan.json").write_text(json.dumps(plan))
    assert main(run_args()) == 2
    assert "size mismatch" in capsys.readouterr().err
    assert not (workspace / "dst").exists()


def test_only_unknown_name_is_an_error(workspace, capsys):
    assert main(run_args("--only", "nope.fastq.gz")) == 2
    assert "not in plan" in capsys.readouterr().err


@pytest.mark.parametrize("factor, status", [(2, "missing"), (0.5, "unplanned")])
def test_wrong_read_count_is_detected_and_cleaned_up(
    tmp_path, monkeypatch, capsys, factor, status
):
    """A portal read_count that is wrong makes split emit the wrong number of
    chunks with exit 0; the run must notice and leave nothing behind."""
    monkeypatch.chdir(tmp_path)
    src = (tmp_path / "src").as_uri() + "/"
    dst = (tmp_path / "dst").as_uri() + "/"
    assert (
        main(
            [
                "test-data",
                "--reads",
                str(READS),
                "--out-prefix",
                src,
                "--meta-dir",
                "meta",
            ]
        )
        == 0
    )
    files = json.loads(Path("meta/files.json").read_text())
    for f in files["@graph"]:
        f["read_count"] = int(READS * factor)
    Path("meta/files.json").write_text(json.dumps(files))
    assert (
        main(
            [
                "plan",
                "--file-sets",
                "meta/sets.json",
                "--files",
                "meta/files.json",
                "--dst",
                dst,
                "--target-bytes",
                str(TARGET),
                "--limit-bytes",
                str(LIMIT),
                "--round-to",
                "1",
                "--out",
                "plan.json",
            ]
        )
        == 0
    )
    assert main(run_args()) == 1
    statuses = manifest_statuses()
    assert status in set(statuses.values())
    assert "ok" in set(statuses.values())
    dst_dir = tmp_path / "dst"
    assert not any(dst_dir.iterdir())
    assert "rerun with: --only" in capsys.readouterr().out
