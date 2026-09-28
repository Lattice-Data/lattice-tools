"""Synthetic quadruple (I1, I2, R1, R2) with shared, deterministic read names.

Writes the gzipped FASTQs to any fsspec prefix and emits ``sets.json`` and
``files.json`` shaped like the data portal, so ``plan`` consumes them unchanged.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import random
from pathlib import Path

from . import s3io

SLOTS = (
    ("index1", "I1", 8),
    ("index2", "I2", 8),
    ("read1", "R1", 150),
    ("read2", "R2", 150),
)
SET_UUID = "00000000-0000-4000-8000-00000000f5e7"


def read_name(i: int, mate: int) -> str:
    # whitespace-separated comment, so the read-ID rule in verify is exercised
    return f"@TEST:1:FC:1:1:{i}:{i} {mate}:N:0:ACGTACGT"


def fastq_bytes(reads: int, length: int, mate: int, seed: int, level: int = 6) -> bytes:
    rng = random.Random(seed)
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", compresslevel=level, mtime=0) as gz:
        for i in range(1, reads + 1):
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            qual = "".join(rng.choice("FFFF:,") for _ in range(length))
            gz.write(f"{read_name(i, mate)}\n{seq}\n+\n{qual}\n".encode())
    return buf.getvalue()


def file_uuid(stem: str) -> str:
    h = hashlib.md5(stem.encode()).hexdigest()
    return f"{h[:8]}-{h[8:12]}-4{h[13:16]}-8{h[17:20]}-{h[20:32]}"


def generate(
    reads: int,
    out_prefix: str,
    meta_dir: Path,
    stem: str = "TEST_S1_L001",
    seed: int = 353,
) -> tuple[Path, Path]:
    out_prefix = out_prefix.rstrip("/") + "/"
    fs = s3io.make_fs(out_prefix)
    set_id = f"/sequence_file_sets/{SET_UUID}/"
    files = []
    fset = {
        "@id": set_id,
        "uuid": SET_UUID,
        "aliases": [f"synthetic:{stem}"],
        "summary": f"synthetic:{stem}",
        "run_cardinality": "paired-end-with-dual-index",
    }
    for k, (slot, tag, length) in enumerate(SLOTS):
        name = f"{stem}_{tag}_001.fastq.gz"
        data = fastq_bytes(
            reads, length, mate=1 if tag in ("R1", "I1") else 2, seed=seed + k
        )
        url = out_prefix + name
        with fs.open(url, "wb") as f:
            f.write(data)
        uuid = file_uuid(name)
        fset[slot] = {"@id": f"/sequence_files/{uuid}/"}
        files.append(
            {
                "@id": f"/sequence_files/{uuid}/",
                "uuid": uuid,
                "file_format": "fastq",
                "file_size": len(data),
                "read_count": reads,
                "s3_uri": url,
                "sequence_file_sets": [set_id],
            }
        )
    meta_dir.mkdir(parents=True, exist_ok=True)
    sets_path = meta_dir / "sets.json"
    files_path = meta_dir / "files.json"
    sets_path.write_text(json.dumps({"@graph": [fset]}, indent=2) + "\n")
    files_path.write_text(json.dumps({"@graph": files}, indent=2) + "\n")
    return sets_path, files_path
