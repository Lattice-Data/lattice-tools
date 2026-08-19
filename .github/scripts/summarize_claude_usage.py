#!/usr/bin/env python3
"""Summarize Claude Code OTLP telemetry captured during a GitHub Actions run.

Reads the newline-delimited OTLP JSON written by the collector's file exporter
and emits a markdown report of token usage per model, plus Claude Code's own
list-price USD estimate.

Usage:
    summarize_claude_usage.py <telemetry.jsonl> <report.md>

Exits 0 even when nothing was captured, so a telemetry failure never fails the
review job. A run that exported nothing produces a report saying so rather than
a zero total, because zero tokens and no data mean very different things.
"""

import json
import sys
from collections import defaultdict

TOKEN_METRIC = "claude_code.token.usage"
COST_METRIC = "claude_code.cost.usage"

# OTLP aggregationTemporality. The wire enum is an int, but protojson - which
# the collector's file exporter uses - serializes enums by name, so the same
# field arrives as "AGGREGATION_TEMPORALITY_DELTA" from a real collector and as
# 1 from a hand-built fixture. Accept both spellings of each.
DELTA = frozenset([1, "1", "AGGREGATION_TEMPORALITY_DELTA"])
CUMULATIVE = frozenset([2, "2", "AGGREGATION_TEMPORALITY_CUMULATIVE"])

# Claude Code prefixes its event names, but the attribute has been seen both
# with and without the prefix. Compare against the bare suffix.
EVENT_PREFIX = "claude_code."
API_REQUEST_EVENT = "api_request"

# Order and display names for the `type` attribute on the token counter.
TOKEN_TYPES = [
    ("input", "Input"),
    ("output", "Output"),
    ("cacheRead", "Cache read"),
    ("cacheCreation", "Cache write"),
]


def attr_value(value):
    """Unwrap an OTLP AnyValue into a plain Python scalar."""
    if not isinstance(value, dict):
        return None
    if "stringValue" in value:
        return value["stringValue"]
    # int64 crosses the wire as a JSON string under protojson, so coerce.
    if "intValue" in value:
        return int(value["intValue"])
    if "doubleValue" in value:
        return float(value["doubleValue"])
    if "boolValue" in value:
        return bool(value["boolValue"])
    return None


def attrs_to_dict(attrs):
    out = {}
    for a in attrs or []:
        key = a.get("key")
        if key is not None:
            out[key] = attr_value(a.get("value", {}))
    return out


def as_number(value):
    """Coerce an attribute to float, or None if it is not numeric."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def datapoint_value(dp):
    if "asInt" in dp:
        return float(int(dp["asInt"]))
    if "asDouble" in dp:
        return float(dp["asDouble"])
    return 0.0


def iter_metric_points(record):
    """Yield (metric_name, attributes, value, temporality) for every datapoint."""
    for rm in record.get("resourceMetrics") or []:
        for sm in rm.get("scopeMetrics") or []:
            for metric in sm.get("metrics") or []:
                name = metric.get("name")
                body = metric.get("sum")
                if not isinstance(body, dict):
                    body = metric.get("gauge")
                if not isinstance(body, dict):
                    continue
                temporality = body.get("aggregationTemporality")
                for dp in body.get("dataPoints") or []:
                    yield (
                        name,
                        attrs_to_dict(dp.get("attributes")),
                        datapoint_value(dp),
                        temporality,
                    )


def event_name(log_record, attrs):
    """The event's name, from either the LogRecord field or its attributes."""
    name = log_record.get("eventName") or attrs.get("event.name") or ""
    if not isinstance(name, str):
        return ""
    if name.startswith(EVENT_PREFIX):
        name = name[len(EVENT_PREFIX) :]
    return name


def iter_log_events(record):
    """Yield (event_name, attributes) for every log record."""
    for rl in record.get("resourceLogs") or []:
        for sl in rl.get("scopeLogs") or []:
            for lr in sl.get("logRecords") or []:
                attrs = attrs_to_dict(lr.get("attributes"))
                yield event_name(lr, attrs), attrs


def series_key(name, attrs):
    """Identity of one metric time series: name plus its attribute set."""
    return (name, tuple(sorted((k, str(v)) for k, v in attrs.items())))


def collect(path):
    """Return (token_totals, cost_totals, api_requests, parse_errors)."""
    # Claude Code is pinned to delta temporality, where summing every datapoint
    # of a series gives the run total. If a series arrives cumulative anyway,
    # the last value is already the total, so take the max instead of summing.
    series_points = defaultdict(list)
    series_temporality = {}
    series_attrs = {}
    api_requests = []
    parse_errors = 0

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                parse_errors += 1
                continue
            if not isinstance(record, dict):
                parse_errors += 1
                continue
            for name, attrs, value, temporality in iter_metric_points(record):
                if name not in (TOKEN_METRIC, COST_METRIC):
                    continue
                key = series_key(name, attrs)
                series_points[key].append(value)
                series_temporality[key] = temporality
                series_attrs[key] = attrs
            for name, attrs in iter_log_events(record):
                if name == API_REQUEST_EVENT:
                    api_requests.append(attrs)

    token_totals = defaultdict(float)  # (model, type) -> tokens
    cost_totals = defaultdict(float)  # model -> usd
    for key, values in series_points.items():
        name, _ = key
        attrs = series_attrs[key]
        if series_temporality[key] in CUMULATIVE:
            total = max(values)
        else:
            total = sum(values)
        model = attrs.get("model") or "unknown"
        if name == TOKEN_METRIC:
            token_totals[(model, attrs.get("type") or "unknown")] += total
        else:
            cost_totals[model] += total

    return token_totals, cost_totals, api_requests, parse_errors


def human(n):
    return f"{int(round(n)):,}"


def render(token_totals, cost_totals, api_requests, parse_errors):
    lines = ["<!-- claude-cost-report -->", "### Claude review usage", ""]

    if not token_totals and not api_requests:
        lines += [
            "No telemetry was captured for this run.",
            "",
            "Either the review step never called the model, or the OTLP exporter "
            "failed to reach the collector. Check the collector logs in the job "
            "output. This is not a zero-token review.",
        ]
        return "\n".join(lines) + "\n"

    models = sorted({model for model, _ in token_totals})
    header = "| Model | " + " | ".join(label for _, label in TOKEN_TYPES) + " | Total |"
    divider = "|---" * (len(TOKEN_TYPES) + 2) + "|"
    lines += [header, divider]

    grand_total = 0.0
    for model in models:
        cells = []
        row_total = 0.0
        for type_key, _ in TOKEN_TYPES:
            value = token_totals.get((model, type_key), 0.0)
            row_total += value
            cells.append(human(value))
        grand_total += row_total
        lines.append(f"| `{model}` | " + " | ".join(cells) + f" | {human(row_total)} |")

    if len(models) > 1:
        # Blank cells for every token-type column, so the row lines up with the
        # header however many types are listed above.
        blanks = " | ".join([""] * len(TOKEN_TYPES))
        lines.append(f"| **All models** | {blanks} | **{human(grand_total)}** |")

    lines.append("")

    total_cost = sum(cost_totals.values())
    if total_cost:
        lines.append(
            f"Claude Code's own list-price estimate for this run: "
            f"**${total_cost:,.4f}**."
        )

    if api_requests:
        durations = [as_number(r.get("duration_ms")) for r in api_requests]
        wall = sum(d for d in durations if d is not None) / 1000.0
        lines.append(f"{len(api_requests)} API requests, {wall:,.1f}s of model time.")

    lines += [
        "",
        "<sub>Cost metrics are approximations. This workflow authenticates with a "
        "subscription OAuth token, so the figure above is a list-rate valuation of "
        "usage drawn against plan limits, not an amount billed. It also excludes "
        "GitHub Actions minutes.</sub>",
    ]

    if parse_errors:
        lines.append(
            f"<sub>{parse_errors} telemetry line(s) could not be parsed.</sub>"
        )

    return "\n".join(lines) + "\n"


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    telemetry_path, report_path = sys.argv[1], sys.argv[2]
    try:
        token_totals, cost_totals, api_requests, parse_errors = collect(telemetry_path)
    except OSError as exc:
        # A missing or unreadable file means the collector never wrote anything.
        # Report that, do not fail the job.
        print(f"could not read {telemetry_path}: {exc}", file=sys.stderr)
        token_totals, cost_totals, api_requests, parse_errors = {}, {}, [], 0

    report = render(token_totals, cost_totals, api_requests, parse_errors)
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
