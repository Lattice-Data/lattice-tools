#!/usr/bin/env python3
"""
Extract Lattice portal profile schemas and classify properties.

By default, lists every profile type from GET /profiles?format=json&frame=object,
then fetches and summarizes each schema. Use --profile to extract one or more
types only.

Mirrors logic in:
  - functions/ProfileRegistry.js + ProfileSlugParse.js (profile list)
  - functions/Profile.js (property classification)

Usage:
  export LATTICE_ACCESS_KEY=...
  export LATTICE_SECRET_KEY=...
  python scripts/extract_lattice_profiles.py -o lattice_profiles.json

  # All profiles (default):
  python scripts/extract_lattice_profiles.py --endpoint https://api.sandbox.lattice-data.org

  # Single profile:
  python scripts/extract_lattice_profiles.py --profile document -o document.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import requests
from requests.auth import HTTPBasicAuth

PROFILE_ID_RE = re.compile(r"/profiles/([^/?#]+)\.json$")

# Same list as Profile.js ADMIN_OR_SYSTEM_PERMISSIONS
ADMIN_OR_SYSTEM_PERMISSIONS = frozenset({
    "add_unvalidated",
    "edit_unvalidated",
    "expand",
    "impersonate",
    "import_items",
    "index",
    "submit_for_any",
    "view_raw",
})

DEFAULT_ENDPOINT = "https://lattice-api-dev.demo.lattice-data.org/"


def fix_profile_json_text(raw_text: str) -> str:
    """Fix known bad regex escapes returned by some portal profile responses."""
    return (
        raw_text.replace("\\\\S\\\\:", "\\\\S:")
        .replace("\\\\-", "-")
    )


def _collect_slug_from_id_value(id_val: str, slugs: set[str]) -> None:
    match = PROFILE_ID_RE.search(id_val)
    if match:
        slugs.add(match.group(1))


def parse_profile_slugs(payload: Any) -> list[str]:
    """Collect profile slugs from /profiles listing payloads.

    Supports both response shapes used by Lattice APIs:
    - JSON-LD: nested objects with ``@id`` (e.g. ``.../profiles/foo.json``)
    - Framed object: profile schemas with ``$id`` (e.g. ``/profiles/foo.json``)
    """
    slugs: set[str] = set()

    def visit(node: Any) -> None:
        if node is None:
            return
        if isinstance(node, list):
            for item in node:
                visit(item)
            return
        if not isinstance(node, dict):
            return

        for id_key in ("@id", "$id"):
            id_val = node.get(id_key)
            if isinstance(id_val, str):
                _collect_slug_from_id_value(id_val, slugs)

        for value in node.values():
            visit(value)

    visit(payload)
    return sorted(slugs)


class LatticeProfileClient:
    def __init__(
        self,
        endpoint: str,
        access_key: str | None = None,
        secret_key: str | None = None,
        timeout: int = 60,
    ) -> None:
        self.endpoint = endpoint.rstrip("/")
        self.session = requests.Session()
        self.session.headers["Accept"] = "application/json"
        if access_key and secret_key:
            self.session.auth = HTTPBasicAuth(access_key, secret_key)
        self.timeout = timeout

    def fetch_profile_slug_list(self) -> list[str]:
        url = f"{self.endpoint}/profiles?format=json&frame=object&limit=all"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return parse_profile_slugs(response.json())

    def fetch_profile_hierarchy(self) -> dict[str, Any]:
        url = f"{self.endpoint}/profiles?format=json&frame=object&limit=all"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        return response.json()["_hierarchy"]["Item"]

    def fetch_profile(self, slug: str) -> dict[str, Any]:
        url = f"{self.endpoint}/profiles/{slug}?format=json"
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        fixed = fix_profile_json_text(response.text)
        return json.loads(fixed)


def is_required_prop(profile: dict[str, Any], prop: str) -> bool:
    if "required" in profile:
        required = profile["required"]
        return isinstance(required, list) and prop in required
    if "anyOf" in profile:
        for sub in profile["anyOf"]:
            if isinstance(sub, dict) and is_required_prop(sub, prop):
                return True
    return False


def get_prop_def(profile: dict[str, Any], prop: str) -> dict[str, Any] | None:
    props = profile.get("properties")
    if not isinstance(props, dict):
        return None
    val = props.get(prop)
    return val if isinstance(val, dict) else None


def get_prop_type(profile: dict[str, Any], prop: str) -> str | None:
    prop_def = get_prop_def(profile, prop)
    if not prop_def:
        return None
    typ = prop_def.get("type")
    return typ if isinstance(typ, str) else None


def get_link_to(profile: dict[str, Any], prop: str) -> str | None:
    prop_def = get_prop_def(profile, prop)
    if not prop_def:
        return None
    if "linkTo" in prop_def and isinstance(prop_def["linkTo"], str):
        return prop_def["linkTo"]
    items = prop_def.get("items")
    if isinstance(items, dict) and isinstance(items.get("linkTo"), str):
        return items["linkTo"]
    return None


def is_searchable_prop(profile: dict[str, Any], prop: str) -> bool:
    return get_link_to(profile, prop) is not None


def is_identifying_prop(profile: dict[str, Any], prop: str) -> bool:
    identifying = profile.get("identifyingProperties", [])
    return isinstance(identifying, list) and prop in identifying


def is_readonly_prop(profile: dict[str, Any], prop: str) -> bool:
    prop_def = get_prop_def(profile, prop)
    return bool(prop_def and prop_def.get("readonly") is True)


def is_not_submittable_prop(profile: dict[str, Any], prop: str) -> bool:
    prop_def = get_prop_def(profile, prop)
    return bool(prop_def and prop_def.get("notSubmittable") is True)


def has_element_type(profile: dict[str, Any], prop: str) -> bool:
    prop_def = get_prop_def(profile, prop)
    return bool(prop_def and prop_def.get("items"))


def get_elements_type(profile: dict[str, Any], prop: str) -> str | None:
    prop_def = get_prop_def(profile, prop)
    if not prop_def:
        return None
    if has_element_type(profile, prop):
        items = prop_def.get("items")
        if isinstance(items, dict) and isinstance(items.get("type"), str):
            return items["type"]
    return None


def has_do_not_submit_comment(profile: dict[str, Any], prop: str) -> bool:
    prop_def = get_prop_def(profile, prop)
    if not prop_def:
        return False
    comment = prop_def.get("comment")
    return isinstance(comment, str) and comment.lower().startswith("do not submit")


def is_admin_or_system_prop(profile: dict[str, Any], prop: str) -> bool:
    prop_def = get_prop_def(profile, prop)
    if not prop_def:
        return False
    permission = prop_def.get("permission")
    return isinstance(permission, str) and permission in ADMIN_OR_SYSTEM_PERMISSIONS


def is_postable_prop(profile: dict[str, Any], prop: str, for_admin: bool = False) -> bool:
    if get_prop_def(profile, prop) is None:
        return False
    if for_admin:
        return not is_not_submittable_prop(profile, prop)
    return (
        is_required_prop(profile, prop)
        or is_identifying_prop(profile, prop)
        or (
            not is_readonly_prop(profile, prop)
            and not is_not_submittable_prop(profile, prop)
            and not has_do_not_submit_comment(profile, prop)
            and not is_admin_or_system_prop(profile, prop)
        )
    )


def get_profile_schema_version(profile: dict[str, Any]) -> Any:
    schema_version = get_prop_def(profile, "schema_version")
    if schema_version and "default" in schema_version:
        return schema_version["default"]
    return None


def get_dependency_comment(profile: dict[str, Any], prop: str) -> str | None:
    dependent = profile.get("dependentSchemas")
    if not isinstance(dependent, dict):
        return None
    dep = dependent.get(prop)
    if isinstance(dep, dict) and isinstance(dep.get("comment"), str):
        return dep["comment"]
    return None


@dataclass
class PropertyInfo:
    name: str
    type: str | None = None
    element_type: str | None = None
    required: bool = False
    optional: bool = True
    identifying: bool = False
    link_to: str | None = None
    searchable: bool = False
    readonly: bool = False
    not_submittable: bool = False
    admin_or_system: bool = False
    do_not_submit_comment: bool = False
    postable: bool = False
    postable_admin: bool = False
    enum: list[Any] | None = None
    default: Any = None
    title: str | None = None
    description: str | None = None
    comment: str | None = None
    permission: str | None = None
    attachment: bool = False
    dependency_comment: str | None = None


@dataclass
class ProfileSummary:
    slug: str
    title: str | None = None
    schema_version: Any = None
    json_schema_uri: str | None = None
    identifying_properties: list[str] = field(default_factory=list)
    required_properties: list[str] = field(default_factory=list)
    properties: list[PropertyInfo] = field(default_factory=list)


def summarize_profile(slug: str, profile: dict[str, Any]) -> ProfileSummary:
    props = profile.get("properties", {})
    if not isinstance(props, dict):
        props = {}

    property_infos: list[PropertyInfo] = []
    required_props: list[str] = []

    for prop_name in sorted(props.keys()):
        prop_def = props[prop_name]
        if not isinstance(prop_def, dict):
            continue

        required = is_required_prop(profile, prop_name)
        if required:
            required_props.append(prop_name)

        link_to = get_link_to(profile, prop_name)
        enum_vals = prop_def.get("enum")
        if enum_vals is not None and not isinstance(enum_vals, list):
            enum_vals = None

        property_infos.append(
            PropertyInfo(
                name=prop_name,
                type=get_prop_type(profile, prop_name),
                element_type=get_elements_type(profile, prop_name),
                required=required,
                optional=not required,
                identifying=is_identifying_prop(profile, prop_name),
                link_to=link_to,
                searchable=is_searchable_prop(profile, prop_name),
                readonly=is_readonly_prop(profile, prop_name),
                not_submittable=is_not_submittable_prop(profile, prop_name),
                admin_or_system=is_admin_or_system_prop(profile, prop_name),
                do_not_submit_comment=has_do_not_submit_comment(profile, prop_name),
                postable=is_postable_prop(profile, prop_name, for_admin=False),
                postable_admin=is_postable_prop(profile, prop_name, for_admin=True),
                enum=enum_vals,
                default=prop_def.get("default"),
                title=prop_def.get("title") if isinstance(prop_def.get("title"), str) else None,
                description=(
                    prop_def.get("description")
                    if isinstance(prop_def.get("description"), str)
                    else None
                ),
                comment=prop_def.get("comment") if isinstance(prop_def.get("comment"), str) else None,
                permission=(
                    prop_def.get("permission")
                    if isinstance(prop_def.get("permission"), str)
                    else None
                ),
                attachment=bool(prop_def.get("attachment") is True),
                dependency_comment=get_dependency_comment(profile, prop_name),
            )
        )

    identifying = profile.get("identifyingProperties", [])
    if not isinstance(identifying, list):
        identifying = []

    return ProfileSummary(
        slug=slug,
        title=profile.get("title") if isinstance(profile.get("title"), str) else None,
        schema_version=get_profile_schema_version(profile),
        json_schema_uri=(
            profile.get("$schema") if isinstance(profile.get("$schema"), str) else None
        ),
        identifying_properties=[str(x) for x in identifying],
        required_properties=required_props,
        properties=property_infos,
    )


def summaries_to_json(summaries: list[ProfileSummary], hierarchy: dict[str, Any], endpoint: str) -> str:
    payload = {
        "endpoint": endpoint,
        "profiles": [asdict(summary) for summary in summaries],
        "profile_count": len(summaries),
        "hierarchy": hierarchy
    }
    return json.dumps(payload, indent=2, sort_keys=False, default=str)


def digest(file_name: os.PathLike | str) -> str:
    with open(file_name, "rb") as f:
        # Read the contents of the file in chunks
        chunk_size = 1024
        hasher = hashlib.sha256()
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Lattice profile schemas as JSON. "
            "Default: discover all profile types from /profiles, then fetch each schema."
        ),
    )
    parser.add_argument(
        "--endpoint",
        default=os.environ.get("LATTICE_ENDPOINT", DEFAULT_ENDPOINT),
        help=f"Lattice API base URL (default: {DEFAULT_ENDPOINT})",
    )
    parser.add_argument(
        "--access-key",
        default=os.environ.get("DEMO_KEY"),
        help="Portal access key (or set LATTICE_ACCESS_KEY)",
    )
    parser.add_argument(
        "--secret-key",
        default=os.environ.get("DEMO_SECRET"),
        help="Portal secret key (or set LATTICE_SECRET_KEY)",
    )
    parser.add_argument(
        "--profile",
        action="append",
        dest="profiles",
        metavar="SLUG",
        help=(
            "Fetch only these profile slugs (repeatable). "
            "Omit to discover and fetch all types from /profiles."
        ),
    )
    parser.add_argument(
        "--output", "-o",
        help="Write JSON output to this file instead of stdout",
    )
    parser.add_argument(
        "--include-raw",
        action="store_true",
        help="Include full raw profile JSON under each summary",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress messages on stderr",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.access_key or not args.secret_key:
        print(
            "Error: credentials required. Set LATTICE_ACCESS_KEY and LATTICE_SECRET_KEY "
            "or pass --access-key / --secret-key.",
            file=sys.stderr,
        )
        return 2

    client = LatticeProfileClient(
        endpoint=args.endpoint,
        access_key=args.access_key,
        secret_key=args.secret_key,
    )

    if args.profiles:
        slugs = args.profiles
        if not args.quiet:
            print(
                f"Fetching {len(slugs)} profile(s): {', '.join(slugs)}",
                file=sys.stderr,
            )
    else:
        profiles_url = f"{client.endpoint}/profiles?format=json&frame=object&limit=all"
        if not args.quiet:
            print(f"Discovering profile types from {profiles_url}", file=sys.stderr)
        slugs = client.fetch_profile_slug_list()
        if not args.quiet:
            print(f"Found {len(slugs)} profile type(s)", file=sys.stderr)

    if not slugs:
        print("No profile slugs found.", file=sys.stderr)
        return 1

    summaries: list[ProfileSummary] = []
    raw_by_slug: dict[str, dict[str, Any]] = {}

    for i, slug in enumerate(slugs, start=1):
        if not args.quiet:
            print(f"[{i}/{len(slugs)}] Fetching {slug}", file=sys.stderr)
        profile = client.fetch_profile(slug)
        summaries.append(summarize_profile(slug, profile))
        if args.include_raw:
            raw_by_slug[slug] = profile

    hierarchy = client.fetch_profile_hierarchy()
    output = summaries_to_json(summaries, hierarchy, args.endpoint)
    if args.include_raw:
        payload = json.loads(output)
        for item in payload["profiles"]:
            item["raw_profile"] = raw_by_slug[item["slug"]]
        output = json.dumps(payload, indent=2, sort_keys=False, default=str)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(output)
        print(f"Wrote {len(summaries)} profile(s) to {args.output}", file=sys.stderr)
        old_hash = digest(Path("json_profiles_current.json"))
        new_hash = digest(Path(args.output))
        print(f"Hash of old JSON: {old_hash}")
        print(f"Hash of generated JSON: {new_hash}")
        if old_hash == new_hash:
            Path(args.output).unlink()
            print(f"No changes to new JSON, removing {args.output}...")
    else:
        print(output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
