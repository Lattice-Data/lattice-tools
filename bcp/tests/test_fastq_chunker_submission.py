"""check-submission over file://: the destination must be exactly the sheet."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from fastq_chunker import s3io, submission, verify
from fastq_chunker.cli import main
from fastq_chunker.plan import GB, FastqFile, make_plan, plan_group, write_plan
from fastq_chunker.run import PROBE_NAME
from fastq_chunker.s3io import make_fs


def fq(src: Path, name, slot, size, reads=10_000_000, group="G"):
    return FastqFile(group, slot, name, (src / name).as_uri(), size, reads)


def make_world(tmp_path: Path):
    """One unsplit set (copied) and one split set (2 chunks per mate)."""
    src, dst = tmp_path / "src", tmp_path / "dst"
    src.mkdir()
    dst.mkdir()
    small = [
        fq(src, "s_R1.fastq.gz", "read1", 11, group="S"),
        fq(src, "s_R2.fastq.gz", "read2", 12, group="S"),
    ]
    big = [
        fq(src, "b_R1.fastq.gz", "read1", 150 * GB, group="B"),
        fq(src, "b_R2.fastq.gz", "read2", 150 * GB, group="B"),
    ]
    plan = make_plan(
        [plan_group(small), plan_group(big)],
        dst.as_uri() + "/",
        80 * GB,
        100 * GB,
        1_000_000,
    )
    sheet = {
        "LIBS": ["s_R1.fastq.gz", "s_R2.fastq.gz"],
        "LIBB": [
            "b_R1.part001.fastq.gz",
            "b_R2.part001.fastq.gz",
            "b_R1.part002.fastq.gz",
            "b_R2.part002.fastq.gz",
        ],
    }
    return plan, sheet, dst


def stage(dst: Path, plan, sidecars=True):
    """Put everything the plan produces into dst, as a good run would."""
    for g in plan.groups:
        for fp in g.files:
            if g.action == "split":
                for c in fp.chunks:
                    (dst / c.name).write_bytes(b"chunk")
                    if sidecars:
                        (dst / (c.name + ".md5")).write_text(f"x  {c.name}\n")
            else:
                (dst / fp.file.filename).write_bytes(b"x" * fp.file.size_bytes)


def verify_report(
    path: Path, plan, fail: str | None = None, checks=("read_count", "md5")
):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group", "filename", "chunk", "check", "status", "detail"])
        for g in plan.split_files():
            grp, fp = g
            for c in fp.chunks:
                for chk in checks:
                    w.writerow([grp.group, fp.file.filename, c.name, chk, "PASS", ""])
                if c.name == fail:
                    w.writerow(
                        [
                            grp.group,
                            fp.file.filename,
                            c.name,
                            "boundary_distinct",
                            "FAIL",
                            "x",
                        ]
                    )
        w.writerow(["*", "*", "*", "overall", "FAIL" if fail else "PASS", ""])
    return path


def full_report(path: Path, plan, monkeypatch, misalign=None, repeat=()):
    """A report from the real ``verify.full``, over chunks whose reads are faked.

    ``misalign``: a chunk whose first read differs from its mate's;
    ``repeat``: chunks whose last read starts their file's next chunk too.
    """
    fs = make_fs(plan.dst)
    chunks = {c.name: c for _, fp in plan.split_files() for c in fp.chunks}

    def fake_full_one(pigz, url):
        c = chunks[s3io.basename(url)]
        md5, _ = s3io.read_sidecar(fs, url)
        first = "@MATE_X" if c.name == misalign else f"@READ_{c.index}"
        last = f"@READ_{c.index + 1}" if c.name in repeat else f"@LAST_{c.index}"
        return verify.FullResult(c.reads * 4, first, last, md5)

    monkeypatch.setattr(verify, "full_one", fake_full_one)
    verify.write_report(path, verify.full(plan, fs, "pigz", workers=1))
    return path


def fails(rows):
    return {(r.name, r.check) for r in rows if r.status == "FAIL"}


def run_check(plan, sheet, dst, report=None, **kw):
    return submission.check(plan, sheet, make_fs(plan.dst), report, **kw)


def test_a_complete_destination_passes(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    rows = run_check(plan, sheet, dst, verify_report(tmp_path / "v.tsv", plan))
    assert fails(rows) == set()
    assert {r.check for r in rows} >= {"copy", "verified_full", "batch_limit"}


def test_missing_and_extra_objects_fail(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    (dst / "s_R2.fastq.gz").unlink()
    (dst / "b_R1.fastq.gz").write_bytes(b"pre-split original")  # never uploaded
    (dst / PROBE_NAME).write_bytes(b"ok")
    (dst / "stray.fastq.gz.md5").write_text("x")
    (dst / "logs").mkdir()
    got = fails(run_check(plan, sheet, dst, verify_report(tmp_path / "v.tsv", plan)))
    assert got == {
        ("s_R2.fastq.gz", "exists"),
        ("b_R1.fastq.gz", "nothing_extra"),
        (PROBE_NAME, "nothing_extra"),
        ("stray.fastq.gz.md5", "nothing_extra"),
        ("logs/", "nothing_extra"),
    }


def test_a_wrong_size_copy_fails(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    (dst / "s_R1.fastq.gz").write_bytes(b"short")
    rows = run_check(plan, sheet, dst, verify_report(tmp_path / "v.tsv", plan))
    assert fails(rows) == {("s_R1.fastq.gz", "copy")}


def test_sheet_and_plan_must_agree(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    sheet["LIBB"] = sheet["LIBB"][:3] + ["b_R2.part003.fastq.gz"]
    sheet["LIBS"].append("s_R1.fastq.gz")
    got = fails(run_check(plan, sheet, dst, verify_report(tmp_path / "v.tsv", plan)))
    assert ("b_R2.part003.fastq.gz", "plan_matches_sheet") in got
    assert ("b_R2.part002.fastq.gz", "plan_matches_sheet") in got
    assert ("s_R1.fastq.gz", "sheet_unique") in got


def test_chunks_need_sidecars_and_a_full_verify(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan, sidecars=False)
    got = fails(run_check(plan, sheet, dst, None))
    chunks = sheet["LIBB"]
    assert got == {(c, "sidecar") for c in chunks} | {
        (c, "verified_full") for c in chunks
    }


def test_a_quick_only_or_failing_verify_report_is_not_enough(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    quick = verify_report(tmp_path / "q.tsv", plan, checks=("exists", "header"))
    assert fails(run_check(plan, sheet, dst, quick)) == {
        (c, "verified_full") for c in sheet["LIBB"]
    }
    # read count alone is not a full verify: the md5 must have been compared too
    half = verify_report(tmp_path / "h.tsv", plan, checks=("read_count",))
    assert fails(run_check(plan, sheet, dst, half)) == {
        (c, "verified_full") for c in sheet["LIBB"]
    }
    bad = verify_report(tmp_path / "f.tsv", plan, fail="b_R2.part002.fastq.gz")
    assert fails(run_check(plan, sheet, dst, bad)) == {
        ("b_R2.part002.fastq.gz", "verified_full"),
        ("(verify report)", "verify_clean"),
    }


def test_a_real_full_verify_report_passes(tmp_path, monkeypatch):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    rows = run_check(
        plan, sheet, dst, full_report(tmp_path / "v.tsv", plan, monkeypatch)
    )
    assert fails(rows) == set()
    assert ("(verify report)", "verify_clean") in {(r.name, r.check) for r in rows}


@pytest.mark.parametrize(
    "kw, verify_fails, bad_chunks",
    [
        # verify writes these with filename * and chunk part001: the FAIL is
        # about the part across mates, so both mates' chunks fail
        (
            {"misalign": "b_R2.part001.fastq.gz"},
            {("*", "part001", "first_id_aligned")},
            {"b_R1.part001.fastq.gz", "b_R2.part001.fastq.gz"},
        ),
        # in both mates, so the mates still agree and only boundary_distinct
        # fails: written per file, with chunk part001
        (
            {"repeat": ("b_R1.part001.fastq.gz", "b_R2.part001.fastq.gz")},
            {
                ("b_R1.fastq.gz", "part001", "boundary_distinct"),
                ("b_R2.fastq.gz", "part001", "boundary_distinct"),
            },
            {"b_R1.part001.fastq.gz", "b_R2.part001.fastq.gz"},
        ),
    ],
)
def test_part_level_verify_fails_reach_the_chunks(
    tmp_path, monkeypatch, kw, verify_fails, bad_chunks
):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    report = full_report(tmp_path / "v.tsv", plan, monkeypatch, **kw)
    with open(report, newline="") as f:
        got = {
            (r["filename"], r["chunk"], r["check"])
            for r in csv.DictReader(f, delimiter="\t")
            if r["status"] == "FAIL" and r["check"] != "overall"
        }
    assert got == verify_fails  # the rows verify really writes for this
    assert fails(run_check(plan, sheet, dst, report)) == {
        (c, "verified_full") for c in bad_chunks
    } | {("(verify report)", "verify_clean")}


def test_any_verify_fail_or_no_overall_pass_fails_the_batch(tmp_path, monkeypatch):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    good = full_report(tmp_path / "v.tsv", plan, monkeypatch)
    lines = good.read_text().splitlines()
    # a FAIL naming nothing the plan has still counts
    odd = tmp_path / "odd.tsv"
    odd.write_text("\n".join(lines[:-1] + ["Z\t*\tpart009\tx\tFAIL\t", lines[-1]]))
    # a report cut short has no overall line
    cut = tmp_path / "cut.tsv"
    cut.write_text("\n".join(lines[:-1]) + "\n")
    for report in (odd, cut):
        assert fails(run_check(plan, sheet, dst, report)) == {
            ("(verify report)", "verify_clean")
        }


def test_a_chunk_at_the_limit_fails(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)  # every chunk is 5 bytes
    report = verify_report(tmp_path / "v.tsv", plan)
    plan.limit_bytes = 6
    assert fails(run_check(plan, sheet, dst, report)) == set()
    plan.limit_bytes = 5  # a chunk must be below the limit, not at it
    assert fails(run_check(plan, sheet, dst, report)) == {
        (c, "size_below_limit") for c in sheet["LIBB"]
    }


def test_batch_limit_uses_actual_bytes(tmp_path):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    report = verify_report(tmp_path / "v.tsv", plan)
    total = 11 + 12 + 4 * len(b"chunk")
    assert fails(run_check(plan, sheet, dst, report, batch_limit_bytes=total)) == set()
    assert fails(run_check(plan, sheet, dst, report, batch_limit_bytes=total - 1)) == {
        ("(batch)", "batch_limit")
    }


def test_command_exit_codes_and_report(tmp_path, capsys):
    plan, sheet, dst = make_world(tmp_path)
    stage(dst, plan)
    write_plan(plan, tmp_path / "plan.json")
    sheet_path = tmp_path / "sra_file_set1.csv"
    rows = ["library_ID,filename,filename2,filename3,filename4"]
    for lib, names in sheet.items():
        rows.append(lib + "," + ",".join(names + [""] * (4 - len(names))))
    sheet_path.write_text("\n".join(rows) + "\n")
    args = [
        "check-submission",
        "--plan",
        str(tmp_path / "plan.json"),
        "--sra-sheet",
        str(sheet_path),
        "--out",
        str(tmp_path / "r.tsv"),
    ]
    assert (
        main(args + ["--verify-report", str(verify_report(tmp_path / "v.tsv", plan))])
        == 0
    )
    assert "PASS:" in capsys.readouterr().out
    (dst / "s_R1.fastq.gz").unlink()
    assert main(args + ["--verify-report", str(tmp_path / "v.tsv")]) == 1
    out = capsys.readouterr().out
    assert "FAIL s_R1.fastq.gz exists" in out
    assert "s_R1.fastq.gz\texists\tFAIL" in (tmp_path / "r.tsv").read_text()
