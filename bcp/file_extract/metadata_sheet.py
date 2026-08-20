"""Fetch the Lattice sample-template tab for Scale extract."""

from __future__ import annotations

import csv
import io
import json
import re
from typing import Callable
from urllib.request import Request, urlopen

from bs4 import BeautifulSoup

from .scale_wells import ScaleExtractError, parse_rt_index_cell

SAMPLE_TEMPLATE_TAB = "sample template"

_BOOTSTRAP_RE = re.compile(r"var bootstrapData = (.*?)};")


def get_tab_gid(
    sheet_id: str,
    tab_name: str = SAMPLE_TEMPLATE_TAB,
    *,
    opener: Callable = urlopen,
) -> str:
    """Resolve a tab name to its numeric gid via the public spreadsheet HTML."""
    sheet_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}"
    req = Request(sheet_url, headers={"User-Agent": "Magic Browser"})
    with opener(req) as resp:
        html = resp.read().decode("utf-8", errors="replace")
    soup = BeautifulSoup(html, "html.parser")
    tab_ids: dict[str, str] = {}
    for script in soup.find_all("script"):
        match = _BOOTSTRAP_RE.search(str(script))
        if not match:
            continue
        data = json.loads(match.group()[20:-1])
        for entry in data.get("changes", {}).get("topsnapshot", []):
            parts = str(entry[1]).split('"')
            if len(parts) > 5:
                tab_ids[parts[5]] = parts[1]
    if tab_name not in tab_ids:
        raise ScaleExtractError(
            f"tab {tab_name!r} not found in Google Sheet {sheet_id}"
        )
    return tab_ids[tab_name]


def fetch_sheet_csv(
    sheet_id: str,
    gid: str,
    *,
    opener: Callable = urlopen,
) -> str:
    url = (
        f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid={gid}"
    )
    req = Request(url, headers={"User-Agent": "Magic Browser"})
    with opener(req) as resp:
        return resp.read().decode("utf-8")


def fetch_sample_template(
    sheet_id: str,
    *,
    tab_name: str = SAMPLE_TEMPLATE_TAB,
    opener: Callable = urlopen,
) -> str:
    """Download the sample template tab as CSV text."""
    gid = get_tab_gid(sheet_id, tab_name, opener=opener)
    return fetch_sheet_csv(sheet_id, gid, opener=opener)


def parse_sample_template(
    csv_text: str, *, experiment: str | None = None
) -> list[dict[str, str]]:
    """Return one row per well with sample_name, RT_index, and normalized well.

    Rows whose first column starts with ``#`` are Lattice comment rows and
    are skipped entirely. When ``experiment`` is set, only rows whose
    ``experiment_name`` equals that value (after strip) are kept.
    """
    reader = csv.DictReader(io.StringIO(csv_text))
    if not reader.fieldnames or "RT_index" not in reader.fieldnames:
        raise ScaleExtractError("sample template must include an RT_index column")
    wanted = None
    if experiment is not None:
        wanted = experiment.strip()
        if not wanted:
            raise ScaleExtractError("experiment must not be empty")
        if "experiment_name" not in reader.fieldnames:
            raise ScaleExtractError(
                "sample template must include an experiment_name column"
            )
    first_col = reader.fieldnames[0]
    rows: list[dict[str, str]] = []
    for raw in reader:
        first_val = (raw.get(first_col) or "").strip()
        if first_val.startswith("#"):
            continue
        if wanted is not None:
            experiment_name = (raw.get("experiment_name") or "").strip()
            if experiment_name != wanted:
                continue
        rt_index = (raw.get("RT_index") or "").strip()
        if not rt_index:
            continue
        sample_name = (raw.get("sample_name") or "").strip()
        for well in parse_rt_index_cell(rt_index):
            rows.append(
                {
                    "sample_name": sample_name,
                    "RT_index": rt_index,
                    "well": well,
                }
            )
    if wanted is not None and not rows:
        raise ScaleExtractError(
            f"no sample template rows with experiment_name {wanted!r}"
        )
    return rows


def sheet_wells_from_csv(csv_text: str, *, experiment: str | None = None) -> set[str]:
    return {
        row["well"] for row in parse_sample_template(csv_text, experiment=experiment)
    }
