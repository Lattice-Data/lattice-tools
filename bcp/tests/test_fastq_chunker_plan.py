"""Tests for fastq_chunker.plan and fastq_chunker.portal: no network, no pigz."""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

import pytest

from fastq_chunker import portal
from fastq_chunker.cli import main
from fastq_chunker.plan import (
    GB,
    LIMIT_BYTES,
    TARGET_BYTES,
    FastqFile,
    PlanError,
    chunk_count,
    chunk_name,
    load_plan,
    make_plan,
    plan_group,
    reads_per_chunk,
    stem_of,
    write_plan,
)

FIXTURES = Path(__file__).parent / "fixtures" / "fastq_chunker"
SET1 = "/sequence_file_sets/11111111-1111-1111-1111-111111111111/"
SET2 = "/sequence_file_sets/22222222-2222-2222-2222-222222222222/"


def fq(filename, size_bytes, reads, group="G", slot="read1"):
    return FastqFile(
        group=group,
        slot=slot,
        filename=filename,
        s3_uri=f"s3://src/{filename}",
        size_bytes=size_bytes,
        reads=reads,
    )


def lib01():
    r = 2_000_000_000
    return [
        fq("LIB01_S1_L001_I1_001.fastq.gz", 31 * GB, r, slot="index1"),
        fq("LIB01_S1_L001_I2_001.fastq.gz", 29 * GB, r, slot="index2"),
        fq("LIB01_S1_L001_R1_001.fastq.gz", 480 * GB, r, slot="read1"),
        fq("LIB01_S1_L001_R2_001.fastq.gz", 510 * GB, r, slot="read2"),
    ]


@pytest.mark.parametrize(
    "filename, stem",
    [
        ("a_R1_001.fastq.gz", "a_R1_001"),
        ("a.fq.gz", "a"),
        ("dots.in.name.fastq.gz", "dots.in.name"),
    ],
)
def test_stem_of(filename, stem):
    assert stem_of(filename) == stem


@pytest.mark.parametrize("bad", ["a.fastq", "a.fq", "a.gz", ".fastq.gz", "a.bam"])
def test_stem_of_rejects(bad):
    with pytest.raises(PlanError):
        stem_of(bad)


def test_chunk_name_pads():
    assert chunk_name("s", 3, 3) == "s.part003.fastq.gz"
    assert chunk_name("s", 1234, 4) == "s.part1234.fastq.gz"


@pytest.mark.parametrize(
    "size, expected",
    [
        (99 * GB, 1),
        (LIMIT_BYTES - 1, 1),
        (LIMIT_BYTES, 2),
        (159 * GB, 2),
        (161 * GB, 3),
        (510 * GB, 7),
    ],
)
def test_chunk_count_never_splits_a_legal_file(size, expected):
    assert chunk_count(size, TARGET_BYTES, LIMIT_BYTES) == expected


def test_plan_group_worked_example():
    g = plan_group(lib01())
    assert g.action == "split"
    assert g.n_chunks == 7
    assert g.reads_per_chunk == 286_000_000
    assert g.lines_per_chunk == 1_144_000_000
    assert g.suffix_width == 3
    r2 = next(fp for fp in g.files if fp.file.slot == "read2")
    assert [c.reads for c in r2.chunks] == [286_000_000] * 6 + [284_000_000]
    assert r2.chunks[0].name == "LIB01_S1_L001_R2_001.part001.fastq.gz"
    assert r2.chunks[-1].name == "LIB01_S1_L001_R2_001.part007.fastq.gz"
    assert r2.chunks[0].est_bytes == round(510 * GB * 286 / 2000)
    assert all(c.est_bytes < LIMIT_BYTES for fp in g.files for c in fp.chunks)
    i1 = next(fp for fp in g.files if fp.file.slot == "index1")
    assert i1.chunks[0].est_bytes == round(31 * GB * 286 / 2000)
    assert i1.chunks[0].name == "LIB01_S1_L001_I1_001.part001.fastq.gz"
    for fp in g.files:
        assert sum(c.reads for c in fp.chunks) == 2_000_000_000


def test_plan_group_below_limit_is_skipped():
    g = plan_group([fq("LIB02_R1.fastq.gz", 95 * GB, 400_000_000)])
    assert g.action == "skip"
    assert g.n_chunks == 1
    assert g.reads_per_chunk == 400_000_000
    assert g.files[0].chunks == []
    assert g.largest_est_bytes() == 95 * GB


def test_rounding_falls_back_when_it_would_empty_the_last_chunk():
    # 3 chunks of 10 reads; rounding to 1,000,000 would put all 10 in chunk 1.
    assert reads_per_chunk(10, 3, 1_000_000) == 4
    assert reads_per_chunk(10, 3, 1) == 4
    assert reads_per_chunk(2_000_000_000, 7, 1_000_000) == 286_000_000
    assert reads_per_chunk(7, 1, 1_000_000) == 7


def test_read_count_mismatch_raises():
    files = lib01()
    files[0] = fq(
        "LIB01_S1_L001_I1_001.fastq.gz", 31 * GB, 1_999_999_999, slot="index1"
    )
    with pytest.raises(PlanError, match="read counts differ"):
        plan_group(files)


def test_target_at_or_above_limit_raises():
    with pytest.raises(PlanError, match="below limit"):
        plan_group(lib01(), target_bytes=LIMIT_BYTES, limit_bytes=LIMIT_BYTES)


def test_estimate_at_limit_raises():
    # 150 GB in 3 reads: 2 chunks of 2+1 reads, and 2/3 of 150 GB is a 100 GB
    # chunk, not below the 100 GB limit.
    with pytest.raises(PlanError, match="not below"):
        plan_group([fq("x.fastq.gz", 150 * GB, 3)], round_to=1)


def test_more_chunks_than_reads_raises():
    with pytest.raises(PlanError, match="cannot fill 3 chunks"):
        plan_group([fq("x.fastq.gz", 200 * GB, 2)], round_to=1)


def test_mixed_groups_raise():
    with pytest.raises(PlanError, match="more than one group"):
        plan_group(
            [fq("a.fastq.gz", 1, 1, group="A"), fq("b.fastq.gz", 1, 1, group="B")]
        )


def test_plan_invariants_hold_over_random_inputs():
    rng = random.Random(353)
    for _ in range(500):
        size = rng.randint(1, 2_000 * GB)
        reads = rng.randint(1, 5_000_000_000)
        round_to = rng.choice([1, 1_000, 1_000_000])
        try:
            g = plan_group([fq("f.fastq.gz", size, reads)], round_to=round_to)
        except PlanError:
            # only reachable when the reads are too few for the chunk count
            assert size >= LIMIT_BYTES
            assert reads < 2 * math.ceil(size / TARGET_BYTES) * 1_000_000
            continue
        n, k = g.reads_per_chunk, g.n_chunks
        assert (k - 1) * n < reads <= k * n
        if g.action == "split":
            chunks = g.files[0].chunks
            assert len(chunks) == k
            assert sum(c.reads for c in chunks) == reads
            assert all(c.reads > 0 for c in chunks)
            assert all(c.est_bytes < LIMIT_BYTES for c in chunks)
            assert len({c.name for c in chunks}) == k


def test_plan_round_trips_through_json(tmp_path):
    plan = make_plan(
        [plan_group(lib01())], "s3://dst/run/", TARGET_BYTES, LIMIT_BYTES, 1_000_000
    )
    write_plan(plan, tmp_path / "plan.json")
    back = load_plan(tmp_path / "plan.json")
    assert back == plan
    assert back.dst == "s3://dst/run/"
    assert back.groups[0].files[3].chunk_uri(
        back.dst, back.groups[0].files[3].chunks[0]
    ) == ("s3://dst/run/LIB01_S1_L001_R2_001.part001.fastq.gz")


# ---- portal JSON -----------------------------------------------------------


def load_fixture_objects():
    return (
        portal.load_objects(FIXTURES / "sets.json"),
        portal.load_objects(FIXTURES / "files.json"),
    )


def test_load_objects_accepts_object_list_and_graph(tmp_path):
    obj = {"@id": "/sequence_files/x/", "uuid": "x"}
    for payload in (obj, [obj], {"@graph": [obj]}):
        p = tmp_path / "o.json"
        p.write_text(json.dumps(payload))
        assert portal.load_objects(p) == [obj]
    (tmp_path / "bad.json").write_text("42")
    with pytest.raises(PlanError):
        portal.load_objects(tmp_path / "bad.json")


def test_build_groups_from_fixture():
    sets, files = load_fixture_objects()
    groups = portal.build_groups(sets, files)
    # the cram-only set contributes nothing; the other two are kept in set order
    assert [g[0].group for g in groups] == [SET1, SET2]
    quad = {f.slot: f for f in groups[0]}
    assert set(quad) == {"index1", "index2", "read1", "read2"}
    assert quad["read2"].filename == "LIB01_S1_L001_R2_001.fastq.gz"
    assert quad["read2"].size_bytes == 510 * GB
    assert quad["read2"].reads == 2_000_000_000
    assert quad["read2"].s3_uri.endswith("/raw/LIB01_S1_L001_R2_001.fastq.gz")
    assert groups[1][0].slot == "read1"
    assert portal.labels(sets)[SET1] == "labalpha:LIB01_S1_L001"


@pytest.mark.parametrize(
    "value",
    [
        "/sequence_files/aaaa0003-0000-0000-0000-000000000003/",
        "aaaa0003-0000-0000-0000-000000000003",
        {"@id": "/sequence_files/aaaa0003-0000-0000-0000-000000000003/"},
        {"uuid": "aaaa0003-0000-0000-0000-000000000003"},
    ],
)
def test_ref_id_forms(value):
    assert (
        portal.ref_id(value) == "/sequence_files/aaaa0003-0000-0000-0000-000000000003/"
    )


def test_build_groups_reports_every_problem_at_once():
    sets, files = load_fixture_objects()
    by_id = {f["@id"]: f for f in files}
    r1 = by_id["/sequence_files/aaaa0003-0000-0000-0000-000000000003/"]
    del r1["read_count"]
    r2 = by_id["/sequence_files/aaaa0004-0000-0000-0000-000000000004/"]
    r2["file_format"] = "cram"
    i1 = by_id["/sequence_files/aaaa0001-0000-0000-0000-000000000001/"]
    i1["s3_uri"] = "s3://example-bucket/x/LIB01_I1.bam"
    files = [
        f
        for f in files
        if not f["@id"].endswith("aaaa0002-0000-0000-0000-000000000002/")
    ]
    with pytest.raises(PlanError) as e:
        portal.build_groups(sets, files)
    msg = str(e.value)
    assert (
        "read1 (/sequence_files/aaaa0003-0000-0000-0000-000000000003/): missing read_count"
        in msg
    )
    assert (
        "read2 (/sequence_files/aaaa0004-0000-0000-0000-000000000004/): file_format is 'cram'"
        in msg
    )
    assert (
        "index1 (/sequence_files/aaaa0001-0000-0000-0000-000000000001/): unsupported extension"
        in msg
    )
    assert (
        "index2: file /sequence_files/aaaa0002-0000-0000-0000-000000000002/ not supplied"
        in msg
    )
    assert msg.count("\n") == 3


def test_file_in_two_sets_is_an_error():
    sets, files = load_fixture_objects()
    sets[1]["read2"] = sets[0]["read1"]
    with pytest.raises(PlanError, match="belongs to more than one set"):
        portal.build_groups(sets, files)


def test_zero_read_count_is_reported_as_missing():
    sets, files = load_fixture_objects()
    files[0]["read_count"] = 0
    with pytest.raises(PlanError, match="missing read_count"):
        portal.build_groups(sets, files)


# ---- CLI --------------------------------------------------------------------


def test_cli_plan_writes_plan_and_summary(tmp_path, capsys):
    out = tmp_path / "plan.json"
    rc = main(
        [
            "plan",
            "--file-sets",
            str(FIXTURES / "sets.json"),
            "--files",
            str(FIXTURES / "files.json"),
            "--dst",
            "s3://dst-bucket/run_chunks",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    plan = load_plan(out)
    assert plan.dst == "s3://dst-bucket/run_chunks/"
    assert [g.action for g in plan.groups] == ["split", "skip"]
    assert plan.groups[0].label == "labalpha:LIB01_S1_L001"
    text = capsys.readouterr().out
    assert "2 groups, 1 to split" in text
    assert f"plan written to {out}" in text
    assert "7 chunks x 286,000,000 reads" in text
    assert "LIB01_S1_L001_R2_001.part007.fastq.gz" in text
    assert "labalpha:LIB02: skip" in text


def test_cli_names_only(tmp_path, capsys):
    rc = main(
        [
            "plan",
            "--file-sets",
            str(FIXTURES / "sets.json"),
            "--files",
            str(FIXTURES / "files.json"),
            "--dst",
            "s3://dst/",
            "--out",
            str(tmp_path / "p.json"),
            "--names-only",
        ]
    )
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 28
    # slot order follows the schema: read1, read2, read3, index1, index2
    assert lines[0] == "LIB01_S1_L001_R1_001.part001.fastq.gz"
    assert lines[7] == "LIB01_S1_L001_R2_001.part001.fastq.gz"
    assert lines[-1] == "LIB01_S1_L001_I2_001.part007.fastq.gz"
    assert all(line.endswith(".fastq.gz") for line in lines)


def test_cli_tiny_targets_force_many_chunks(tmp_path):
    rc = main(
        [
            "plan",
            "--file-sets",
            str(FIXTURES / "sets.json"),
            "--files",
            str(FIXTURES / "files.json"),
            "--dst",
            "s3://dst/",
            "--out",
            str(tmp_path / "p.json"),
            "--target-bytes",
            str(50 * GB),
            "--limit-bytes",
            str(60 * GB),
            "--round-to",
            "1",
        ]
    )
    assert rc == 0
    plan = load_plan(tmp_path / "p.json")
    assert plan.groups[0].n_chunks == 11
    assert plan.groups[1].action == "split"
    assert plan.groups[1].n_chunks == 2
    assert plan.groups[1].reads_per_chunk == 200_000_000


def test_cli_reports_plan_errors_on_stderr(tmp_path, capsys):
    sets, files = load_fixture_objects()
    (tmp_path / "sets.json").write_text(json.dumps(sets))
    (tmp_path / "files.json").write_text(json.dumps(files[1:]))
    rc = main(
        [
            "plan",
            "--file-sets",
            str(tmp_path / "sets.json"),
            "--files",
            str(tmp_path / "files.json"),
            "--dst",
            "s3://dst/",
            "--out",
            str(tmp_path / "p.json"),
        ]
    )
    assert rc == 2
    assert "not supplied" in capsys.readouterr().err
    assert not (tmp_path / "p.json").exists()


def test_too_few_reads_for_the_chunk_count_is_a_plan_error():
    # 5 reads over 4 chunks: 2 per chunk fills 3 chunks, never 4
    with pytest.raises(PlanError, match="cannot be spread"):
        reads_per_chunk(5, 4, 1)
    with pytest.raises(PlanError):
        plan_group([fq("x.fastq.gz", 300 * GB, 5)], round_to=1)


def test_same_basename_in_two_sets_is_rejected():
    a = plan_group([fq("S1_L001_R1_001.fastq.gz", 200 * GB, 10**9, group="A")])
    b = plan_group([fq("S1_L001_R1_001.fastq.gz", 150 * GB, 10**9, group="B")])
    b.files[0].file = FastqFile(
        **{
            **b.files[0].file.__dict__,
            "s3_uri": "s3://src/runB/S1_L001_R1_001.fastq.gz",
        }
    )
    with pytest.raises(PlanError, match="share a stem") as e:
        make_plan([a, b], "s3://dst/", TARGET_BYTES, LIMIT_BYTES, 1)
    assert (
        "s3://src/S1_L001_R1_001.fastq.gz, s3://src/runB/S1_L001_R1_001.fastq.gz"
        in str(e.value)
    )
    # a skipped singleton collides with a split file's stem just the same
    small = plan_group([fq("S1_L001_R1_001.fq.gz", 10 * GB, 10**8, group="C")])
    with pytest.raises(PlanError, match="share a stem"):
        make_plan([a, small], "s3://dst/", TARGET_BYTES, LIMIT_BYTES, 1)


def test_stem_shaped_like_another_files_chunk_is_rejected():
    a = plan_group([fq("S1.fastq.gz", 200 * GB, 10**9, group="A")])
    b = plan_group([fq("S1.part2.fastq.gz", 10 * GB, 10**8, group="B")])
    with pytest.raises(PlanError, match="looks like a chunk of 'S1'"):
        make_plan([a, b], "s3://dst/", TARGET_BYTES, LIMIT_BYTES, 1)
    # the other way round is fine: 'S1' is not shaped like a chunk of anything
    c = plan_group([fq("S1.partial.fastq.gz", 10 * GB, 10**8, group="C")])
    make_plan([a, c], "s3://dst/", TARGET_BYTES, LIMIT_BYTES, 1)
