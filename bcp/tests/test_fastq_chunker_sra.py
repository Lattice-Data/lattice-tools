"""Tests for fastq_chunker.sra and fastq_chunker.db2: SRA sheets as input, the
name contract, and the resolve command. No network: DB2 is a fake search."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastq_chunker import db2, sra
from fastq_chunker.cli import cmd_resolve, main
from fastq_chunker.plan import GB, PlanError, make_plan, plan_group
from fastq_chunker.portal import build_groups

SET_A = "11111111-0000-0000-0000-00000000000a"
SET_B = "11111111-0000-0000-0000-00000000000b"
SET_C = "11111111-0000-0000-0000-00000000000c"


def fuuid(n: int) -> str:
    return f"aaaaaaaa-0000-0000-0000-{n:012d}"


def fobj(n, name, size=2 * GB, reads=1_000_000, names=None, **extra):
    o = {
        "@id": f"/sequence_files/{fuuid(n)}/",
        "uuid": fuuid(n),
        "file_format": "fastq",
        "file_size": size,
        "read_count": reads,
        "s3_uri": f"s3://example-bucket/raw/{name}",
        "status": "current",
        "sra_file_names": [name] if names is None else names,
    }
    o.update(extra)
    return o


def sobj(set_uuid, *file_numbers, **extra):
    slots = ("read1", "read2", "read3", "index1", "index2")
    o = {
        "@id": f"/sequence_file_sets/{set_uuid}/",
        "uuid": set_uuid,
        "status": "current",
        "summary": f"labalpha:{set_uuid[-1]}",
    }
    for slot, n in zip(slots, file_numbers):
        o[slot] = f"/sequence_files/{fuuid(n)}/"
    o.update(extra)
    return o


def world():
    """Two libraries: LIB1 with one set (R1/R2), LIB2 with two sets (two runs)."""
    files = {
        fuuid(1): fobj(1, "run1-LIB1_R1_001.fastq.gz"),
        fuuid(2): fobj(2, "run1-LIB1_R2_001.fastq.gz"),
        fuuid(3): fobj(3, "run1-LIB2_R1_001.fastq.gz"),
        fuuid(4): fobj(4, "run1-LIB2_R2_001.fastq.gz"),
        fuuid(5): fobj(5, "run2-LIB2_R1_001.fastq.gz"),
        fuuid(6): fobj(6, "run2-LIB2_R2_001.fastq.gz"),
    }
    sets = {
        SET_A: sobj(SET_A, 1, 2),
        SET_B: sobj(SET_B, 3, 4),
        SET_C: sobj(SET_C, 5, 6),
    }
    sheet = {
        "LIB1": ["run1-LIB1_R1_001.fastq.gz", "run1-LIB1_R2_001.fastq.gz"],
        "LIB2": [
            "run1-LIB2_R1_001.fastq.gz",
            "run1-LIB2_R2_001.fastq.gz",
            "run2-LIB2_R1_001.fastq.gz",
            "run2-LIB2_R2_001.fastq.gz",
        ],
    }
    rows = [("LIB1", SET_A), ("LIB2", SET_B), ("LIB2", SET_C)]
    return sheet, rows, sets, files


# --- sheets -----------------------------------------------------------------


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_read_sra_sheet_keeps_filename_columns_in_order_and_ignores_filetype(tmp_path):
    p = write(
        tmp_path / "s.csv",
        "sample_name,library_ID,filetype,filename,filename2,filename3\n"
        "s1,LIB1,fastq,a_R1.fastq.gz,a_R2.fastq.gz,\n"
        's1,LIB2,fastq,"b_R1.fastq.gz",b_R2.fastq.gz,c.fastq.gz\n',
    )
    assert sra.read_sra_sheet(p) == {
        "LIB1": ["a_R1.fastq.gz", "a_R2.fastq.gz"],
        "LIB2": ["b_R1.fastq.gz", "b_R2.fastq.gz", "c.fastq.gz"],
    }


def test_read_sra_sheet_rejects_repeated_library(tmp_path):
    p = write(
        tmp_path / "s.csv", "library_ID,filename\nLIB1,a.fastq.gz\nLIB1,b.fastq.gz\n"
    )
    with pytest.raises(PlanError, match="LIB1 repeats"):
        sra.read_sra_sheet(p)


def test_read_sra_sheet_needs_library_and_filename_columns(tmp_path):
    with pytest.raises(PlanError, match="no library_ID"):
        sra.read_sra_sheet(write(tmp_path / "a.csv", "lib,filename\nx,y\n"))
    with pytest.raises(PlanError, match="no filename columns"):
        sra.read_sra_sheet(write(tmp_path / "b.csv", "library_ID,filetype\nx,y\n"))


def test_read_fileset_sheet(tmp_path):
    p = write(
        tmp_path / "f.csv",
        f"library_ID,sequence_file_set_uuid\nLIB1,{SET_A.upper()}\nLIB2,{SET_B}\n\n",
    )
    assert sra.read_fileset_sheet(p) == [("LIB1", SET_A), ("LIB2", SET_B)]


@pytest.mark.parametrize(
    "body, match",
    [
        ("library_ID,uuid\nLIB1,x\n", "header must be"),
        (
            f"library_ID,sequence_file_set_uuid\nLIB1,{SET_A}\nLIB2,{SET_A}\n",
            "already listed",
        ),
        ("library_ID,sequence_file_set_uuid\nLIB1,not-a-uuid\n", "not a uuid"),
        ("library_ID,sequence_file_set_uuid\n", "no rows"),
        (
            f"library_ID,sequence_file_set_uuid\nLIB1,{SET_A},extra\n",
            "expected 2 columns",
        ),
    ],
)
def test_read_fileset_sheet_rejects(tmp_path, body, match):
    with pytest.raises(PlanError, match=match):
        sra.read_fileset_sheet(write(tmp_path / "f.csv", body))


# --- resolve ----------------------------------------------------------------


def test_resolve_matches_and_returns_objects_plan_reads():
    sheet, rows, sets, files = world()
    r = sra.resolve(sheet, rows, sets, files)
    assert [s["uuid"] for s in r.sets] == [SET_A, SET_B, SET_C]
    assert len(r.files) == 6
    assert r.n_libraries == 2
    assert r.total_bytes == 12 * GB
    # the output feeds the existing portal path unchanged
    assert len(build_groups(r.sets, r.files)) == 3


def test_resolve_ignores_sheet_order_within_a_library():
    sheet, rows, sets, files = world()
    sheet["LIB2"] = list(reversed(sheet["LIB2"]))
    sra.resolve(sheet, rows, sets, files)


def errors(sheet, rows, sets, files, **kw) -> str:
    with pytest.raises(PlanError) as e:
        sra.resolve(sheet, rows, sets, files, **kw)
    return str(e.value)


def test_resolve_reports_a_name_missing_from_the_sheet():
    sheet, rows, sets, files = world()
    sheet["LIB2"].remove("run2-LIB2_R2_001.fastq.gz")
    msg = errors(sheet, rows, sets, files)
    assert "LIB2: names differ" in msg
    assert "predicted but not in the sheet: run2-LIB2_R2_001.fastq.gz" in msg


def test_resolve_reports_a_name_db2_does_not_predict():
    sheet, rows, sets, files = world()
    sheet["LIB1"].append("stray.fastq.gz")
    assert "in the sheet but not predicted: stray.fastq.gz" in errors(
        sheet, rows, sets, files
    )


def test_resolve_compares_chunk_names_when_db2_predicts_a_split():
    sheet, rows, sets, files = world()
    chunks = ["run1-LIB1_R1_001.part001.fastq.gz", "run1-LIB1_R1_001.part002.fastq.gz"]
    files[fuuid(1)]["sra_file_names"] = chunks
    # the sheet still has the unsplit name: a stale sheet is caught
    assert "LIB1: names differ" in errors(sheet, rows, sets, files)
    sheet["LIB1"] = chunks + ["run1-LIB1_R2_001.fastq.gz"]
    sra.resolve(sheet, rows, sets, files)


def test_resolve_reports_a_duplicate_name_inside_one_row():
    sheet, rows, sets, files = world()
    sheet["LIB1"].append("run1-LIB1_R1_001.fastq.gz")
    assert "in the sheet but not predicted: run1-LIB1_R1_001.fastq.gz" in errors(
        sheet, rows, sets, files
    )


def test_resolve_reports_libraries_missing_from_either_sheet():
    sheet, rows, sets, files = world()
    sheet["LIB3"] = ["x.fastq.gz"]
    rows.append(("LIB9", "11111111-0000-0000-0000-0000000000ff"))
    msg = errors(sheet, rows, sets, files)
    assert "LIB3: in the SRA sheet but has no set" in msg
    assert "LIB9: listed in the set sheet but not in the SRA sheet" in msg


def test_resolve_reports_missing_and_deleted_objects_together():
    sheet, rows, sets, files = world()
    del sets[SET_C]
    files[fuuid(1)]["status"] = "deleted"
    del files[fuuid(2)]
    msg = errors(sheet, rows, sets, files)
    assert f"set {SET_C} not found in DB2" in msg
    assert f"file {fuuid(1)} is deleted" in msg
    assert f"file {fuuid(2)} not found in DB2" in msg


def test_resolve_rejects_a_deleted_set_and_a_set_without_reads():
    sheet, rows, sets, files = world()
    sets[SET_A]["status"] = "deleted"
    for slot in ("read1", "read2"):
        del sets[SET_B][slot]
    msg = errors(sheet, rows, sets, files)
    assert f"set {SET_A} is deleted" in msg
    assert f"set {SET_B} has no read slots" in msg


def test_resolve_rejects_a_file_in_two_listed_sets():
    sheet, rows, sets, files = world()
    sets[SET_C]["read2"] = f"/sequence_files/{fuuid(4)}/"
    assert f"file {fuuid(4)} is also in set {SET_B}" in errors(sheet, rows, sets, files)


def test_resolve_rejects_a_file_db2_gives_no_names():
    sheet, rows, sets, files = world()
    files[fuuid(3)]["sra_file_names"] = []
    del files[fuuid(4)]["sra_file_names"]
    msg = errors(sheet, rows, sets, files)
    assert f"file {fuuid(3)} has no sra_file_names" in msg
    assert f"file {fuuid(4)} has no sra_file_names" in msg


def test_resolve_enforces_the_batch_limit():
    sheet, rows, sets, files = world()
    sra.resolve(sheet, rows, sets, files, batch_limit_bytes=12 * GB)
    assert "above the 11 GB limit" in errors(
        sheet, rows, sets, files, batch_limit_bytes=11 * GB
    )


def test_slot_refs_accepts_embedded_objects():
    s = {"read1": {"@id": f"/sequence_files/{fuuid(1)}/"}, "index1": fuuid(2)}
    assert sra.slot_refs(s) == [("read1", fuuid(1)), ("index1", fuuid(2))]


# --- name contract ----------------------------------------------------------


def contract_plan(files: list[dict]):
    sets = [sobj(SET_A, 1, 2)]
    groups = [plan_group(m) for m in build_groups(sets, files)]
    return make_plan(groups, "s3://dst/p/", 80 * GB, 100 * GB, 1_000_000)


def test_contract_passes_unsplit_and_split_files():
    files = [fobj(1, "a_R1.fastq.gz"), fobj(2, "a_R2.fastq.gz")]
    sra.name_contract(contract_plan(files), files)
    big = [
        fobj(
            1,
            "a_R1.fastq.gz",
            size=150 * GB,
            reads=10_000_000,
            names=["a_R1.part001.fastq.gz", "a_R1.part002.fastq.gz"],
        ),
        fobj(
            2,
            "a_R2.fastq.gz",
            size=150 * GB,
            reads=10_000_000,
            names=["a_R2.part001.fastq.gz", "a_R2.part002.fastq.gz"],
        ),
    ]
    sra.name_contract(contract_plan(big), big)


def test_contract_refuses_a_prediction_the_plan_does_not_make():
    files = [
        fobj(
            1, "a_R1.fastq.gz", names=["a_R1.part001.fastq.gz", "a_R1.part002.fastq.gz"]
        ),
        fobj(2, "a_R2.fastq.gz"),
    ]
    with pytest.raises(PlanError, match="a_R1.fastq.gz: plan names"):
        sra.name_contract(contract_plan(files), files)


def test_contract_refuses_chunk_order_or_count_drift():
    big = [
        fobj(
            1,
            "a_R1.fastq.gz",
            size=250 * GB,
            reads=10_000_000,
            names=["a_R1.part001.fastq.gz", "a_R1.part002.fastq.gz"],
        ),
        fobj(
            2,
            "a_R2.fastq.gz",
            size=250 * GB,
            reads=10_000_000,
            names=["a_R2.part001.fastq.gz", "a_R2.part002.fastq.gz"],
        ),
    ]
    # 250 GB at an 80 GB target is 4 chunks; DB2 says 2
    with pytest.raises(PlanError, match="differ from DB2"):
        sra.name_contract(contract_plan(big), big)


def test_contract_refuses_the_right_chunks_in_the_wrong_order():
    def big(n, stem):
        names = [f"{stem}.part002.fastq.gz", f"{stem}.part001.fastq.gz"]
        return fobj(n, f"{stem}.fastq.gz", size=150 * GB, reads=10_000_000, names=names)

    files = [big(1, "a_R1"), big(2, "a_R2")]
    # 150 GB is 2 chunks, as DB2 says, but DB2 lists them part002 first
    with pytest.raises(PlanError, match="differ from DB2"):
        sra.name_contract(contract_plan(files), files)


def test_contract_needs_names_on_every_file_once_any_has_them():
    files = [fobj(1, "a_R1.fastq.gz"), fobj(2, "a_R2.fastq.gz")]
    del files[1]["sra_file_names"]
    with pytest.raises(PlanError, match="a_R2.fastq.gz: no sra_file_names"):
        sra.name_contract(contract_plan(files), files)


def test_contract_is_off_for_exports_without_names():
    files = [fobj(1, "a_R1.fastq.gz"), fobj(2, "a_R2.fastq.gz")]
    for f in files:
        del f["sra_file_names"]
    sra.name_contract(contract_plan(files), files)


def test_plan_command_applies_the_contract(tmp_path):
    files = [
        fobj(1, "a_R1.fastq.gz", names=["wrong.fastq.gz"]),
        fobj(2, "a_R2.fastq.gz"),
    ]
    (tmp_path / "sets.json").write_text(json.dumps([sobj(SET_A, 1, 2)]))
    (tmp_path / "files.json").write_text(json.dumps(files))
    out = tmp_path / "plan.json"
    rc = main(
        [
            "plan",
            "--file-sets",
            str(tmp_path / "sets.json"),
            "--files",
            str(tmp_path / "files.json"),
            "--dst",
            "s3://dst/p/",
            "--out",
            str(out),
        ]
    )
    assert rc == 2
    assert not out.exists()


# --- db2 client and the resolve command -------------------------------------


class FakeSearch:
    def __init__(self, objects: list[dict]):
        self.objects = objects
        self.calls: list[dict] = []

    def __call__(self, params: dict) -> list[dict]:
        self.calls.append(params)
        want = set(params["uuid"])
        kind = {
            "SequenceFileSet": "/sequence_file_sets/",
            "SequenceFile": "/sequence_files/",
        }
        prefix = kind[params["type"]]
        return [
            o for o in self.objects if o["uuid"] in want and o["@id"].startswith(prefix)
        ]


def test_fetch_by_uuid_batches_and_drops_duplicates(monkeypatch):
    monkeypatch.setattr(db2, "BATCH", 2)
    _, _, sets, _ = world()
    search = FakeSearch(list(sets.values()))
    got = db2.fetch_by_uuid(search, "SequenceFileSet", [SET_A, SET_B, SET_A, SET_C])
    assert set(got) == {SET_A, SET_B, SET_C}
    assert [len(c["uuid"]) for c in search.calls] == [2, 1]
    assert all(c["frame"] == "object" and c["limit"] == "all" for c in search.calls)


def test_client_needs_all_three_variables(monkeypatch):
    monkeypatch.setenv("DB2_TEST_KEY", "k")
    monkeypatch.delenv("DB2_TEST_SECRET", raising=False)
    monkeypatch.delenv("DB2_TEST_SERVER", raising=False)
    with pytest.raises(PlanError, match="DB2_TEST_SECRET, DB2_TEST_SERVER"):
        db2.Client("db2_test", session=object())


class FakeResponse:
    def __init__(self, status_code, body=None):
        self.status_code, self.body = status_code, body or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self.body


class FakeSession:
    def __init__(self, response):
        self.response, self.headers = response, {}

    def get(self, url, params, timeout):
        return self.response


@pytest.mark.parametrize(
    "response, want",
    [
        # snovault answers a search that finds nothing with 404
        (FakeResponse(404, {"@graph": []}), []),
        (FakeResponse(200, {"@graph": [{"uuid": "u"}]}), [{"uuid": "u"}]),
    ],
)
def test_client_search_reads_the_graph_and_treats_404_as_empty(
    monkeypatch, response, want
):
    for k in ("KEY", "SECRET", "SERVER"):
        monkeypatch.setenv(f"DB2_TEST_{k}", "x")
    client = db2.Client("db2_test", session=FakeSession(response))
    assert client.search({"type": "SequenceFile"}) == want


def test_client_search_raises_on_other_errors(monkeypatch):
    for k in ("KEY", "SECRET", "SERVER"):
        monkeypatch.setenv(f"DB2_TEST_{k}", "x")
    client = db2.Client("db2_test", session=FakeSession(FakeResponse(500)))
    with pytest.raises(RuntimeError, match="HTTP 500"):
        client.search({"type": "SequenceFile"})


class FakeClient:
    def __init__(self, objects):
        self.search = FakeSearch(objects)


def sheets(tmp_path, sheet, rows):
    sra_path = tmp_path / "sra_file_set1.csv"
    width = max(len(v) for v in sheet.values())
    cols = ["filename"] + [f"filename{i}" for i in range(2, width + 1)]
    lines = ["sample_name,library_ID," + ",".join(cols)]
    for lib, names in sheet.items():
        lines.append(f"s,{lib}," + ",".join(names + [""] * (width - len(names))))
    sra_path.write_text("\n".join(lines) + "\n")
    set_path = tmp_path / "sra_sequence_fileset_set1.csv"
    set_path.write_text(
        "library_ID,sequence_file_set_uuid\n"
        + "".join(f"{lib},{u}\n" for lib, u in rows)
    )
    return sra_path, set_path


def resolve_args(tmp_path, sra_path, set_path, batch_limit_gb=5000.0):
    import argparse

    return argparse.Namespace(
        sra_sheet=sra_path,
        set_sheet=set_path,
        mode="db2_test",
        out_dir=tmp_path / "out",
        batch_limit_gb=batch_limit_gb,
    )


def test_resolve_command_writes_inputs_that_plan_accepts(tmp_path, capsys):
    sheet, rows, sets, files = world()
    sra_path, set_path = sheets(tmp_path, sheet, rows)
    client = FakeClient(list(sets.values()) + list(files.values()))
    assert cmd_resolve(resolve_args(tmp_path, sra_path, set_path), client=client) == 0
    assert "2 libraries, 3 sets, 6 files, 12 GB" in capsys.readouterr().out
    out = tmp_path / "out"
    rc = main(
        [
            "plan",
            "--file-sets",
            str(out / "sets.json"),
            "--files",
            str(out / "files.json"),
            "--dst",
            "s3://dst/p/",
            "--out",
            str(out / "plan.json"),
        ]
    )
    assert rc == 0
    plan = json.loads((out / "plan.json").read_text())
    assert [g["action"] for g in plan["groups"]] == ["skip"] * 3


def test_resolve_command_writes_nothing_on_a_mismatch(tmp_path):
    sheet, rows, sets, files = world()
    sheet["LIB1"] = sheet["LIB1"][:1]
    sra_path, set_path = sheets(tmp_path, sheet, rows)
    client = FakeClient(list(sets.values()) + list(files.values()))
    with pytest.raises(PlanError, match="LIB1: names differ"):
        cmd_resolve(resolve_args(tmp_path, sra_path, set_path), client=client)
    assert not (tmp_path / "out").exists()


def test_resolve_command_applies_the_batch_limit(tmp_path):
    sheet, rows, sets, files = world()
    sra_path, set_path = sheets(tmp_path, sheet, rows)
    client = FakeClient(list(sets.values()) + list(files.values()))
    args = resolve_args(tmp_path, sra_path, set_path, batch_limit_gb=11.0)
    with pytest.raises(PlanError, match="above the 11 GB limit"):
        cmd_resolve(args, client=client)
    assert not (tmp_path / "out").exists()
