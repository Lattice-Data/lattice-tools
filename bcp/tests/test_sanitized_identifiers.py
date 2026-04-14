from __future__ import annotations

from pathlib import Path
import re


TEST_ROOT = Path(__file__).parent


BLOCKED_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b(?:weissman|trapnell|hamazaki|marson|lange|ucsf|califano)\b", re.IGNORECASE
    ),
    re.compile(r"\bNVUS(?!0000000000-)\d{10}-\d+\b"),
    re.compile(r"\bAN(?!0000000[0-9])\d{8}\b"),
    re.compile(r"/ORPROJ1/"),
    re.compile(r"\b(?:liguo|pennyyang)\b", re.IGNORECASE),
)


def _iter_test_text_files() -> list[Path]:
    return sorted(
        p
        for p in TEST_ROOT.rglob("*")
        if p.is_file()
        and p.suffix.lower() in {".py", ".csv", ".tsv", ".json", ".html"}
        and p.name != "test_sanitized_identifiers.py"
    )


def test_test_assets_are_sanitized() -> None:
    violations: list[str] = []

    for path in _iter_test_text_files():
        content = path.read_text(encoding="utf-8")
        for pattern in BLOCKED_PATTERNS:
            for match in pattern.finditer(content):
                violations.append(
                    f"{path.relative_to(TEST_ROOT)} :: {pattern.pattern} :: {match.group(0)}"
                )

    assert not violations, "Unsanitized test identifiers found:\n" + "\n".join(
        violations
    )
