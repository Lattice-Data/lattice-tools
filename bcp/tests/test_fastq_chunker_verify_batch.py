"""Unit tests for verify helpers and batch packing; no pigz, no network."""

from __future__ import annotations

import gzip
import io
import shlex

import pytest

from fastq_chunker import batch, s3io, verify
from fastq_chunker.plan import (
    GB,
    FastqFile,
    PlanError,
    load_plan,
    make_plan,
    plan_group,
    write_plan,
)
from fastq_chunker.run import (
    MANIFEST_COLUMNS,
    RunError,
    Tools,
    decompress_command,
    find_tools,
    parse_put_lines,
    pipeline_command,
)


@pytest.mark.parametrize(
    "header, rid",
    [
        (
            "@A00123:45:HXYZ:1:1101:1000:2000 1:N:0:ACGT",
            "@A00123:45:HXYZ:1:1101:1000:2000",
        ),
        (
            "@A00123:45:HXYZ:1:1101:1000:2000 2:N:0:ACGT",
            "@A00123:45:HXYZ:1:1101:1000:2000",
        ),
        ("@SRR000001.1/1", "@SRR000001.1"),
        ("@SRR000001.1/2\n", "@SRR000001.1"),
        ("@UG-1-2-3-4", "@UG-1-2-3-4"),
    ],
)
def test_read_id_default_rule(header, rid):
    assert verify.read_id(header) == rid


def test_read_id_regex():
    assert verify.read_id("@x:1:2 mate=2", r"^(@[^ ]+) mate=\d") == "@x:1:2"
    with pytest.raises(ValueError):
        verify.read_id("@nomatch", r"^(zzz)")


def test_header_sanity():
    assert verify.header_sanity(["@r", "ACGT", "+", "FFFF"]) is None
    assert "four lines" in verify.header_sanity(["@r", "ACGT"])
    assert "line 1" in verify.header_sanity(["r", "ACGT", "+", "FFFF"])
    assert "line 3" in verify.header_sanity(["@r", "ACGT", "-", "FFFF"])
    assert "lengths" in verify.header_sanity(["@r", "ACGT", "+", "FFF"])


def test_first_record_reads_only_the_head(tmp_path):
    body = b"".join(f"@r{i} x\nACGT\n+\nFFFF\n".encode() for i in range(50_000))
    p = tmp_path / "x.fastq.gz"
    p.write_bytes(gzip.compress(body))
    fs = s3io.make_fs(p.as_uri())
    assert verify.first_record(fs, p.as_uri()) == ["@r0 x", "ACGT", "+", "FFFF"]


@pytest.mark.parametrize("block", [1, 3, 7, 64, 4096])
def test_count_stream_first_line_and_last_record_header(block):
    data = b"@r1 x\nACGT\n+\nFFFF\n@r2 y\nACGA\n+\nFF:F\n@r3 z\nTTTT\n+\n,,,,\n"
    lines, first, last = verify.count_stream(io.BytesIO(data), block)
    assert (lines, first, last) == (12, "@r1 x", "@r3 z")


def test_count_stream_edge_cases():
    # a partial trailing line is not counted and not a candidate header
    lines, first, last = verify.count_stream(io.BytesIO(b"@a\nA\n+\nF\npartial"), 2)
    assert (lines, first, last) == (4, "@a", "@a")
    assert verify.count_stream(io.BytesIO(b"a\nb\n"), 4) == (2, "a", "")
    assert verify.count_stream(io.BytesIO(b""), 4) == (0, "", "")


def test_parse_put_lines_ignores_noise():
    out = "junk\nx.part001.fastq.gz\t123\tabc\nnot\tnumber\tz\n"
    assert parse_put_lines(out) == [("x.part001.fastq.gz", 123, "abc")]


def lib():
    r = 2_000_000_000
    return [
        FastqFile(
            "S", "read1", "L_R1_001.fastq.gz", "s3://src/L_R1_001.fastq.gz", 480 * GB, r
        ),
        FastqFile(
            "S", "read2", "L_R2_001.fastq.gz", "s3://src/L_R2_001.fastq.gz", 510 * GB, r
        ),
    ]


def test_pipeline_command_keeps_FILE_for_split():
    g = plan_group(lib())
    tools = Tools(pigz="/usr/bin/pigz", split="/usr/bin/split", python="/py")
    cmd = pipeline_command(tools, g, g.files[1], "s3://dst/run/", threads=8, level=6)
    assert cmd.startswith("set -o pipefail\n")
    assert (
        "/py -m fastq_chunker.s3io get --concurrency 8 s3://src/L_R2_001.fastq.gz | /usr/bin/pigz -dc | "
        in cmd
    )
    assert "-l 1144000000 --numeric-suffixes=1 -a 3 --filter=" in cmd
    assert cmd.rstrip().endswith("- part")
    # the filter is single-quoted so bash leaves $FILE for split's shell, and
    # contains no pipe: the compressor runs inside put
    filt = shlex.split(cmd.split("--filter=")[1].split(" - part")[0])[0]
    assert filt.count("$FILE") == 1
    words = shlex.split(filt.replace("$FILE", "part001"))
    assert words[:4] == ["/py", "-m", "fastq_chunker.s3io", "put"]
    assert words[4] == "s3://dst/run/L_R2_001.part001.fastq.gz"
    assert words[5:] == [
        "--upload-concurrency",
        "4",
        "--compress",
        "/usr/bin/pigz",
        "-c",
        "-p",
        "8",
        "-6",
    ]


def two_group_plan(dst="s3://dst/run/"):
    small = [
        FastqFile(
            "T",
            "read1",
            "T_R1.fastq.gz",
            "s3://src/T_R1.fastq.gz",
            95 * GB,
            400_000_000,
        )
    ]
    return make_plan(
        [plan_group(lib()), plan_group(small)], dst, 80 * GB, 100 * GB, 1_000_000
    )


def test_batches_keep_a_set_together_and_use_actual_sizes(tmp_path):
    plan = two_group_plan()
    manifest = tmp_path / "run_manifest.tsv"
    rows = ["\t".join(MANIFEST_COLUMNS)]
    rows.append(
        "L_R1_001.fastq.gz\tL_R1_001.part001.fastq.gz\t70000000000\tabc\tok\t1.0"
    )
    manifest.write_text("\n".join(rows) + "\n")
    batches = batch.build_batches(plan, manifest, limit_bytes=5_000 * GB)
    assert len(batches) == 1
    items = {i.name: i for i in batches[0].items}
    assert items["L_R1_001.part001.fastq.gz"].bytes == 70 * GB
    assert (
        items["L_R1_001.part002.fastq.gz"].bytes
        == plan.groups[0].files[0].chunks[1].est_bytes
    )
    assert items["T_R1.fastq.gz"].uri == "s3://src/T_R1.fastq.gz"
    assert (
        items["L_R1_001.part001.fastq.gz"].uri
        == "s3://dst/run/L_R1_001.part001.fastq.gz"
    )
    assert len(batches[0].items) == 15

    # the split set is ~990 GB; a 1 TB limit puts the 95 GB singleton in batch 2
    batches = batch.build_batches(plan, None, limit_bytes=1_000 * GB)
    assert [sorted({i.group for i in b.items}) for b in batches] == [["S"], ["T"]]
    out = tmp_path / "batches.tsv"
    batch.write_batches(out, batches)
    lines = out.read_text().splitlines()
    assert lines[0] == "batch\tgroup\tlabel\tname\tbytes\turi"
    assert lines[-1].startswith("2\tT\tT\tT_R1.fastq.gz\t95000000000\t")

    with pytest.raises(PlanError, match="alone is"):
        batch.build_batches(plan, None, limit_bytes=500 * GB)


def test_write_report_marks_overall(tmp_path):
    checks = [verify.Check("g", "f", "c", "exists", "PASS", "1")]
    assert verify.write_report(tmp_path / "r.tsv", checks) is True
    assert (tmp_path / "r.tsv").read_text().splitlines()[-1].split("\t")[4] == "PASS"
    checks.append(verify.Check("g", "f", "c", "md5", "FAIL", "x"))
    assert verify.write_report(tmp_path / "r.tsv", checks) is False
    assert (tmp_path / "r.tsv").read_text().splitlines()[-1].split("\t")[4] == "FAIL"


def test_plan_json_survives_cli_round_trip(tmp_path):
    p = tmp_path / "plan.json"
    write_plan(two_group_plan(), p)
    assert load_plan(p).groups[1].action == "skip"


@pytest.mark.parametrize(
    "size, est, warns",
    [
        (100, 100, False),
        (114, 100, False),
        (116, 100, True),
        (80, 100, False),
        (71, 100, False),
        (69, 100, True),
        (50, 0, False),
    ],
)
def test_drift_is_asymmetric(size, est, warns):
    detail = verify.drift_detail(size, est)
    assert bool(detail) is warns
    if warns:
        assert f"vs estimate {est}" in detail and "%" in detail


def test_pipeline_command_with_rapidgzip_and_concurrency():
    g = plan_group(lib())
    tools = Tools(
        pigz="/usr/bin/pigz", split="/usr/bin/split", python="/py", rapidgzip="/rg"
    )
    cmd = pipeline_command(tools, g, g.files[0], "s3://dst/", 8, 6, "rapidgzip", 12)
    assert (
        "s3io get --concurrency 12 s3://src/L_R1_001.fastq.gz | /rg -d -c -P 8 | "
        in cmd
    )
    assert "-dc" not in cmd
    default = pipeline_command(tools, g, g.files[0], "s3://dst/", 8, 6)
    assert "--concurrency 8 " in default and "| /usr/bin/pigz -dc |" in default


def test_decompress_command_requires_located_rapidgzip():
    tools = Tools(pigz="/p", split="/s")
    assert decompress_command(tools, "pigz", 4) == "/p -dc"
    with pytest.raises(RunError, match="rapidgzip"):
        decompress_command(tools, "rapidgzip", 4)


def test_find_tools_rejects_unknown_decompressor(monkeypatch):
    with pytest.raises(RunError, match="unknown decompressor"):
        find_tools(decompressor="zstd")
    monkeypatch.setattr(
        "fastq_chunker.run.shutil.which",
        lambda name: None if name == "rapidgzip" else f"/fake/{name}",
    )
    monkeypatch.setattr(
        "fastq_chunker.run.find_split", lambda explicit=None: "/fake/split"
    )
    with pytest.raises(RunError, match="pip install rapidgzip"):
        find_tools(decompressor="rapidgzip")


def test_pipeline_quotes_hostile_stem_and_dst():
    stem = 'we$ird `x` "q" name_R1_001'
    files = [
        FastqFile(
            "S",
            "read1",
            stem + ".fastq.gz",
            "s3://src/" + stem + ".fastq.gz",
            480 * GB,
            10**9,
        )
    ]
    g = plan_group(files)
    tools = Tools(pigz="/p", split="/s", python="/py")
    cmd = pipeline_command(tools, g, g.files[0], "s3://d$t/run it/", 8, 6)
    filt_quoted = cmd.split("--filter=")[1].split(" - part")[0]
    # unwrap bash's quoting of the filter argument and check what split's shell will see
    filt = shlex.split(filt_quoted)[0]
    words = shlex.split(filt.replace("$FILE", "part001"))
    assert words[4] == "s3://d$t/run it/" + stem + ".part001.fastq.gz"
    assert filt.count("$FILE") == 1
    assert (
        shlex.split("echo " + cmd.split(" get ")[1].split(" | ")[0])[-1]
        == "s3://src/" + stem + ".fastq.gz"
    )


def test_run_plan_default_tools_honour_the_decompressor(monkeypatch, tmp_path):
    from fastq_chunker import run as run_mod
    from fastq_chunker.run import RunOptions, run_plan

    seen = {}

    def fake_find_tools(pigz=None, split=None, decompressor="pigz", rapidgzip=None):
        seen["decompressor"] = decompressor
        return Tools(pigz="/p", split="/s", python="/py", rapidgzip="/rg")

    monkeypatch.setattr(run_mod, "find_tools", fake_find_tools)
    plan = two_group_plan()
    rc = run_plan(
        plan, RunOptions(decompressor="rapidgzip", dry_run=True), out=io.StringIO()
    )
    assert rc == 0
    assert seen == {"decompressor": "rapidgzip"}
