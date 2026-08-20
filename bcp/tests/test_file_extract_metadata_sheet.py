"""Tests for file_extract.metadata_sheet."""

from __future__ import annotations

import json

import pytest

from file_extract.metadata_sheet import fetch_sample_template, get_tab_gid
from file_extract.scale_wells import ScaleExtractError


def _bootstrap_html(tab_name: str, gid: str) -> bytes:
    payload = {
        "changes": {
            "topsnapshot": [
                [0, f'xxxx "{gid}" yyyy "ignored" zzzz "{tab_name}" extra'],
            ]
        }
    }
    # get_tab_gid slices the full match: "var bootstrapData = " is 20 chars,
    # then json, then "};" — [20:-1] drops the trailing semicolon.
    script = f"var bootstrapData = {json.dumps(payload)};"
    return f"<html><script>{script}</script></html>".encode("utf-8")


class _FakeResponse:
    def __init__(self, body: bytes) -> None:
        self._body = body

    def read(self) -> bytes:
        return self._body

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, *args: object) -> None:
        return None


def test_get_tab_gid_reads_sample_template() -> None:
    html = _bootstrap_html("sample template", "12345")

    def opener(_req: object) -> _FakeResponse:
        return _FakeResponse(html)

    assert get_tab_gid("sheet-id", opener=opener) == "12345"


def test_get_tab_gid_missing_tab() -> None:
    html = _bootstrap_html("other tab", "1")

    def opener(_req: object) -> _FakeResponse:
        return _FakeResponse(html)

    with pytest.raises(ScaleExtractError, match="sample template"):
        get_tab_gid("sheet-id", opener=opener)


def test_fetch_sample_template_uses_gid_then_csv() -> None:
    html = _bootstrap_html("sample template", "99")
    csv_body = b"sample_name,RT_index\nfoo,SCALEQUANT-A1\n"
    calls: list[str] = []

    def opener(req: object) -> _FakeResponse:
        url = getattr(req, "full_url", None) or getattr(req, "get_full_url")()
        calls.append(url)
        if "export" in url:
            return _FakeResponse(csv_body)
        return _FakeResponse(html)

    text = fetch_sample_template("sheet-id", opener=opener)
    assert "SCALEQUANT-A1" in text
    assert any("gid=99" in url for url in calls)
