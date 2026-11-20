"""Tests for re-sequencing order parsing, the 10x layout check, and routing."""

from __future__ import annotations

import pytest

from file_extract.models import ListedObject
from file_extract.orders import (
    PROCESSED_LAYOUT,
    RAW_LAYOUT,
    OrderError,
    OrderPrefix,
    check_layout,
    check_single_processed,
    group_id_for,
    order_id_from_prefix,
    parse_order_prefixes,
    route,
    scan_order_folders,
    split_order_values,
)
from file_extract.sheets import validate_cro_order

BUCKET = "czi-novogene"
MAIN = f"s3://{BUCKET}/smith-proj/ORDER_A"


def _orders(*reseq: str, main_cro_order: str | None = "ORDER_A"):
    return parse_order_prefixes(MAIN, list(reseq), main_cro_order=main_cro_order)


def _objects(prefix: str, *rel_keys: str) -> list[ListedObject]:
    return [ListedObject(key=f"{prefix}{rel}", size_bytes=1) for rel in rel_keys]


# --- parsing -------------------------------------------------------------


def test_parse_main_first_then_each_reseq_order() -> None:
    orders = _orders(
        f"s3://{BUCKET}/smith-proj/ORDER_B,s3://{BUCKET}/smith-proj/ORDER_C"
    )
    assert [o.cro_order for o in orders] == ["ORDER_A", "ORDER_B", "ORDER_C"]
    assert [o.is_primary for o in orders] == [True, False, False]
    assert orders[1].prefix == "smith-proj/ORDER_B/"
    assert orders[1].uri == f"s3://{BUCKET}/smith-proj/ORDER_B/"


def test_parse_accepts_space_and_comma_separated_values() -> None:
    orders = _orders(
        f"s3://{BUCKET}/smith-proj/ORDER_B/",
        f"s3://{BUCKET}/smith-proj/ORDER_C, s3://{BUCKET}/smith-proj/ORDER_D",
    )
    assert [o.cro_order for o in orders[1:]] == ["ORDER_B", "ORDER_C", "ORDER_D"]


def test_parse_main_id_defaults_to_last_folder() -> None:
    orders = _orders(f"s3://{BUCKET}/smith-proj/ORDER_B", main_cro_order=None)
    assert orders[0].cro_order == "ORDER_A"


def test_parse_runs_order_ids_through_validator() -> None:
    seen: list[str] = []

    def validate(value: str) -> str:
        seen.append(value)
        return validate_cro_order(value)

    parse_order_prefixes(
        MAIN,
        [f"s3://{BUCKET}/smith-proj/ORDER_B"],
        main_cro_order="ORDER_A",
        validate_order_id=validate,
    )
    assert seen == ["ORDER_B"]


@pytest.mark.parametrize(
    ("values", "message"),
    [
        ([], "at least one order"),
        ([""], "empty values"),
        ([f"s3://{BUCKET}/smith-proj/ORDER_B,"], "empty values"),
        (["smith-proj/ORDER_B"], "Invalid --resequencing-orders"),
        ([f"s3://{BUCKET}"], "a bucket alone"),
        (["s3://other-bucket/smith-proj/ORDER_B"], "same bucket"),
        ([MAIN], "more than once, or is the main order"),
        ([f"{MAIN}/"], "more than once, or is the main order"),
        ([f"{MAIN}/GRP_001"], "overlaps"),
        ([f"s3://{BUCKET}/smith-proj"], "overlaps"),
        (
            [
                f"s3://{BUCKET}/smith-proj/ORDER_B",
                f"s3://{BUCKET}/smith-proj/ORDER_B",
            ],
            "more than once",
        ),
        (
            [
                f"s3://{BUCKET}/smith-proj/ORDER_B",
                f"s3://{BUCKET}/smith-proj/ORDER_B/sub",
            ],
            "overlaps",
        ),
        (
            [
                f"s3://{BUCKET}/proj-1/ORDER_B",
                f"s3://{BUCKET}/proj-2/ORDER_B",
            ],
            "order ID 'ORDER_B'",
        ),
        ([f"s3://{BUCKET}/other-proj/ORDER_A"], "order ID 'ORDER_A'"),
    ],
)
def test_parse_rejects(values: list[str], message: str) -> None:
    with pytest.raises(OrderError, match=message):
        _orders(*values)


def test_parse_rejects_reseq_named_like_main_cro_order() -> None:
    # The main folder and --cro-order can differ; both names are reserved.
    with pytest.raises(OrderError, match="order ID 'NVUS-15'"):
        parse_order_prefixes(
            MAIN,
            [f"s3://{BUCKET}/other/NVUS-15"],
            main_cro_order="NVUS-15",
        )


def test_parse_rejects_main_order_at_bucket_root() -> None:
    with pytest.raises(OrderError, match="whole bucket"):
        parse_order_prefixes(f"s3://{BUCKET}", [f"s3://{BUCKET}/x/ORDER_B"])


def test_split_order_values_keeps_empty_tokens_for_the_caller() -> None:
    assert split_order_values(["a,b", "c"]) == ["a", "b", "c"]
    assert split_order_values(["a,"]) == ["a", ""]
    assert split_order_values(None) == []


def test_order_id_from_prefix() -> None:
    assert order_id_from_prefix("smith-proj/ORDER_B/") == "ORDER_B"
    assert order_id_from_prefix("ORDER_B") == "ORDER_B"


# --- GroupID -------------------------------------------------------------


@pytest.mark.parametrize(
    ("key", "prefix", "expected"),
    [
        ("p/ORDER_A/GRP_001/raw/x.fastq.gz", "p/ORDER_A/", "GRP_001"),
        ("p/ORDER_A/GRP_001/processed/a/b.h5", "p/ORDER_A/", "GRP_001"),
        ("p/ORDER_A/x.fastq.gz", "p/ORDER_A/", ""),
        ("GRP_001/raw/x.fastq.gz", "", "GRP_001"),
    ],
)
def test_group_id_for(key: str, prefix: str, expected: str) -> None:
    assert group_id_for(key, prefix) == expected


def test_group_id_for_rejects_key_outside_prefix() -> None:
    with pytest.raises(ValueError):
        group_id_for("p/ORDER_B/GRP/raw/x", "p/ORDER_A/")


# --- 10x layout check ----------------------------------------------------


def _order(name: str, primary: bool = False) -> OrderPrefix:
    return OrderPrefix(BUCKET, f"smith-proj/{name}/", name, is_primary=primary)


def test_raw_layout_accepts_10x_shape() -> None:
    a = _order("ORDER_A", primary=True)
    check_layout(
        [(a, _objects(a.prefix, "GRP_001/raw/RUN1-GRP_001_R1_001.fastq.gz"))],
        RAW_LAYOUT,
    )


def test_layout_errors_on_file_directly_in_order_folder() -> None:
    a = _order("ORDER_A", primary=True)
    with pytest.raises(OrderError, match="no GroupID folder") as info:
        check_layout([(a, _objects(a.prefix, "RUN1_R1_001.fastq.gz"))], RAW_LAYOUT)
    assert "smith-proj/ORDER_A/RUN1_R1_001.fastq.gz" in str(info.value)


@pytest.mark.parametrize("folder", ["raw", "processed", "per_sample_outs"])
def test_layout_errors_on_structural_group_id(folder: str) -> None:
    a = _order("ORDER_A", primary=True)
    with pytest.raises(OrderError, match="must be the order folder"):
        check_layout([(a, _objects(a.prefix, f"{folder}/x/y.fastq.gz"))], RAW_LAYOUT)


def test_layout_errors_on_scale_sci_shape_in_any_order() -> None:
    a = _order("ORDER_A", primary=True)
    b = _order("ORDER_B")
    listings = [
        (a, _objects(a.prefix, "GRP_001/raw/ok.cram")),
        (b, _objects(b.prefix, "RNA3_098/raw/426971/426971-GRP01_QSR-1.cram")),
    ]
    with pytest.raises(OrderError, match="only the 10x layout") as info:
        check_layout(listings, RAW_LAYOUT)
    assert "ORDER_B/RNA3_098/raw/426971/" in str(info.value)
    assert "ok.cram" not in str(info.value)


def test_layout_error_lists_first_paths_and_counts_the_rest() -> None:
    a = _order("ORDER_A", primary=True)
    rel = [f"G{i}/raw/{i}/f.cram" for i in range(8)]
    with pytest.raises(OrderError) as info:
        check_layout([(a, _objects(a.prefix, *rel))], RAW_LAYOUT)
    message = str(info.value)
    assert message.startswith("8 file(s) are not in the 10x layout")
    assert "... and 3 more" in message


def test_layout_error_reports_each_kind_of_problem() -> None:
    a = _order("ORDER_A", primary=True)
    objects = _objects(a.prefix, "loose.cram", "raw/x.cram", "G/raw/1/x.cram")
    with pytest.raises(OrderError) as info:
        check_layout([(a, objects)], RAW_LAYOUT)
    message = str(info.value)
    assert "no GroupID folder" in message
    assert "must be the order folder" in message
    assert "only the 10x layout" in message


def test_processed_layout() -> None:
    a = _order("ORDER_A", primary=True)
    ok = "GRP_001/processed/cellranger/run/outs/per_sample_outs/S1/m.h5"
    check_layout([(a, _objects(a.prefix, ok))], PROCESSED_LAYOUT)
    with pytest.raises(OrderError, match="only the 10x layout"):
        check_layout([(a, _objects(a.prefix, "GRP_001/other/m.h5"))], PROCESSED_LAYOUT)


# --- routing -------------------------------------------------------------


def _route(main_rel, *reseq):
    a = _order("ORDER_A", primary=True)
    listings = [(a, _objects(a.prefix, *main_rel))]
    for name, rel in reseq:
        order = _order(name)
        listings.append((order, _objects(order.prefix, *rel)))
    return route(listings)


def _keys(output_set) -> list[str]:
    return [routed.obj.key for routed in output_set.files]


def test_route_known_group_is_added_to_main_and_main_files_stay() -> None:
    routing = _route(
        ["GRP_001/raw/RUN1_R1.fastq.gz", "GRP_001/raw/RUN1_R2.fastq.gz"],
        ("ORDER_B", ["GRP_001/raw/RUN2_R1.fastq.gz", "GRP_001/raw/RUN2_R2.fastq.gz"]),
    )
    assert len(routing.output_sets) == 1
    main = routing.output_sets[0]
    assert main.owner.cro_order == "ORDER_A"
    assert _keys(main) == [
        "smith-proj/ORDER_A/GRP_001/raw/RUN1_R1.fastq.gz",
        "smith-proj/ORDER_A/GRP_001/raw/RUN1_R2.fastq.gz",
        "smith-proj/ORDER_B/GRP_001/raw/RUN2_R1.fastq.gz",
        "smith-proj/ORDER_B/GRP_001/raw/RUN2_R2.fastq.gz",
    ]
    assert [r.order.cro_order for r in main.files] == [
        "ORDER_A",
        "ORDER_A",
        "ORDER_B",
        "ORDER_B",
    ]
    assert routing.reports[1].added_to_main == {"GRP_001": 2}
    assert routing.reports[1].new_groups == {}
    assert routing.warnings == []


def test_route_new_group_gets_its_own_set() -> None:
    routing = _route(
        ["GRP_001/raw/a.fastq.gz"],
        ("ORDER_B", ["GRP_003/raw/b.fastq.gz"]),
    )
    assert [s.owner.cro_order for s in routing.output_sets] == ["ORDER_A", "ORDER_B"]
    assert _keys(routing.output_sets[0]) == [
        "smith-proj/ORDER_A/GRP_001/raw/a.fastq.gz"
    ]
    assert _keys(routing.output_sets[1]) == [
        "smith-proj/ORDER_B/GRP_003/raw/b.fastq.gz"
    ]
    assert routing.reports[1].new_groups == {"GRP_003": 1}


def test_route_splits_a_mixed_order_and_skips_orders_with_nothing_new() -> None:
    routing = _route(
        ["GRP_001/raw/a.fastq.gz", "GRP_002/raw/a2.fastq.gz"],
        ("ORDER_B", ["GRP_001/raw/b.fastq.gz", "GRP_003/raw/b3.fastq.gz"]),
        ("ORDER_C", ["GRP_002/raw/c.fastq.gz"]),
    )
    assert [s.owner.cro_order for s in routing.output_sets] == ["ORDER_A", "ORDER_B"]
    assert _keys(routing.output_sets[0])[-2:] == [
        "smith-proj/ORDER_B/GRP_001/raw/b.fastq.gz",
        "smith-proj/ORDER_C/GRP_002/raw/c.fastq.gz",
    ]
    assert _keys(routing.output_sets[1]) == [
        "smith-proj/ORDER_B/GRP_003/raw/b3.fastq.gz"
    ]
    assert routing.reports[0].group_ids == ["GRP_001", "GRP_002"]
    assert [r.order.cro_order for r in routing.reports] == [
        "ORDER_A",
        "ORDER_B",
        "ORDER_C",
    ]


def test_route_two_orders_with_new_groups_make_two_extra_sets() -> None:
    routing = _route(
        ["GRP_001/raw/a.fastq.gz"],
        ("ORDER_B", ["GRP_001/raw/b.fastq.gz", "GRP_003/raw/b3.fastq.gz"]),
        ("ORDER_C", ["GRP_001/raw/c.fastq.gz", "GRP_004/raw/c4.fastq.gz"]),
    )
    assert [s.owner.cro_order for s in routing.output_sets] == [
        "ORDER_A",
        "ORDER_B",
        "ORDER_C",
    ]


def test_route_warns_and_records_an_empty_reseq_order() -> None:
    routing = _route(["GRP_001/raw/a.fastq.gz"], ("ORDER_B", []))
    assert [o.cro_order for o in routing.empty_orders] == ["ORDER_B"]
    assert any("ORDER_B" in w and "matched no files" in w for w in routing.warnings)
    assert len(routing.output_sets) == 1


def test_route_warns_when_every_group_is_new() -> None:
    routing = _route(["GRP_001/raw/a.fastq.gz"], ("ORDER_B", ["GRP_009/raw/b.cram"]))
    assert any("check that this is the right order" in w for w in routing.warnings)


def test_route_warns_on_near_miss_group_ids_but_treats_them_as_new() -> None:
    routing = _route(
        ["GRP-001/raw/a.fastq.gz", "GRP_002/raw/a2.fastq.gz"],
        ("ORDER_B", ["grp_001/raw/b.fastq.gz", "GRP_002/raw/b2.fastq.gz"]),
    )
    assert routing.reports[1].new_groups == {"grp_001": 1}
    assert any(
        "'grp_001'" in w and "'GRP-001'" in w and "treated as a new GroupID" in w
        for w in routing.warnings
    )


def test_route_warns_when_main_order_is_empty() -> None:
    routing = _route([], ("ORDER_B", ["GRP_001/raw/b.fastq.gz"]))
    assert any("main order ORDER_A" in w for w in routing.warnings)
    assert [s.owner.cro_order for s in routing.output_sets] == ["ORDER_A", "ORDER_B"]
    assert routing.output_sets[0].files == []


def test_route_requires_main_order_first() -> None:
    with pytest.raises(ValueError):
        route([(_order("ORDER_B"), [])])


# --- GroupID folders and processed/ ------------------------------------------


def _folders(order, *rel):
    return scan_order_folders(order, _objects(order.prefix, *rel))


def test_scan_order_folders_counts_raw_only_groups() -> None:
    a = OrderPrefix(BUCKET, "proj/ORDER_A/", "ORDER_A", is_primary=True)
    folders = _folders(
        a,
        "GRP_001/raw/x_R1_001.fastq.gz",
        "GRP_002/processed/outs/m.h5",
        "md5.txt",
        "raw/stray.fastq.gz",
    )
    assert folders.group_ids == {"GRP_001", "GRP_002"}
    assert folders.processed == {"GRP_002"}


def test_check_single_processed_allows_one_processed_per_shared_group() -> None:
    a = OrderPrefix(BUCKET, "proj/ORDER_A/", "ORDER_A", is_primary=True)
    b = OrderPrefix(BUCKET, "proj/ORDER_B/", "ORDER_B")
    c = OrderPrefix(BUCKET, "proj/ORDER_C/", "ORDER_C")
    check_single_processed(
        [
            _folders(a, "GRP_001/raw/x.fastq.gz", "GRP_002/processed/m.h5"),
            _folders(b, "GRP_001/processed/m.h5", "GRP_003/processed/m.h5"),
            # GRP_003 is not the main order's, so two processed/ is not checked.
            _folders(c, "GRP_002/raw/y.fastq.gz", "GRP_003/processed/m.h5"),
        ]
    )


def test_check_single_processed_errors_and_names_every_folder() -> None:
    a = OrderPrefix(BUCKET, "proj/ORDER_A/", "ORDER_A", is_primary=True)
    b = OrderPrefix(BUCKET, "proj/ORDER_B/", "ORDER_B")
    c = OrderPrefix(BUCKET, "proj/ORDER_C/", "ORDER_C")
    with pytest.raises(OrderError) as exc:
        check_single_processed(
            [
                _folders(a, "GRP_001/raw/x.fastq.gz", "GRP_002/processed/m.h5"),
                _folders(b, "GRP_001/processed/m.h5"),
                _folders(c, "GRP_001/processed/m.h5", "GRP_002/processed/m.h5"),
            ]
        )
    message = str(exc.value)
    assert message.startswith("2 GroupID(s)")
    for uri in (
        f"s3://{BUCKET}/proj/ORDER_B/GRP_001/processed/",
        f"s3://{BUCKET}/proj/ORDER_C/GRP_001/processed/",
        f"s3://{BUCKET}/proj/ORDER_A/GRP_002/processed/",
        f"s3://{BUCKET}/proj/ORDER_C/GRP_002/processed/",
    ):
        assert uri in message
    assert f"s3://{BUCKET}/proj/ORDER_A/GRP_001/processed/" not in message


def test_route_main_group_ids_route_raw_only_groups_to_main() -> None:
    a = OrderPrefix(BUCKET, "proj/ORDER_A/", "ORDER_A", is_primary=True)
    b = OrderPrefix(BUCKET, "proj/ORDER_B/", "ORDER_B")
    routing = route(
        [
            (a, _objects(a.prefix, "GRP_002/processed/m.h5")),
            (b, _objects(b.prefix, "GRP_001/processed/m.h5")),
        ],
        main_group_ids={"GRP_001", "GRP_002"},
    )
    assert len(routing.output_sets) == 1
    assert [r.obj.key for r in routing.output_sets[0].files] == [
        "proj/ORDER_A/GRP_002/processed/m.h5",
        "proj/ORDER_B/GRP_001/processed/m.h5",
    ]


# --- include_new (--resequencing-contains-new) -------------------------------


def test_route_include_new_false_leaves_new_groups_out() -> None:
    a = _order("ORDER_A", primary=True)
    b = _order("ORDER_B")
    routing = route(
        [
            (a, _objects(a.prefix, "GRP_001/raw/a.fastq.gz")),
            (
                b,
                _objects(b.prefix, "GRP_001/raw/b.fastq.gz", "GRP_003/raw/b3.fastq.gz"),
            ),
        ],
        include_new=False,
    )
    assert [s.owner.cro_order for s in routing.output_sets] == ["ORDER_A"]
    assert _keys(routing.output_sets[0]) == [
        "smith-proj/ORDER_A/GRP_001/raw/a.fastq.gz",
        "smith-proj/ORDER_B/GRP_001/raw/b.fastq.gz",
    ]
    assert [r.obj.key for r in routing.left_out] == [
        "smith-proj/ORDER_B/GRP_003/raw/b3.fastq.gz"
    ]
    # The report still counts what was left out.
    assert routing.reports[1].new_groups == {"GRP_003": 1}
    assert routing.include_new is False


def test_route_include_new_false_warnings_say_left_out() -> None:
    a = _order("ORDER_A", primary=True)
    b = _order("ORDER_B")
    routing = route(
        [
            (a, _objects(a.prefix, "GRP-001/raw/a.fastq.gz")),
            (b, _objects(b.prefix, "GRP_001/raw/b.fastq.gz")),
        ],
        include_new=False,
    )
    assert any("are left out" in w for w in routing.warnings)
    assert any("treated as a new GroupID and left out" in w for w in routing.warnings)
