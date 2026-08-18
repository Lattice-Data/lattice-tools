"""Tests for order_qa.spec: resolving a spec to exactly one order prefix."""

import pytest

from order_qa.spec import (
    OrderSpec,
    SpecError,
    classify_order_shape,
    is_non_order_name,
    resolve_spec,
)


class TestResolveSpec:
    def test_bare_three_segment_form(self):
        spec = resolve_spec("czi-psomagen/orglab-alpha/AN00000001")
        assert spec.bucket == "czi-psomagen"
        assert spec.project == "orglab-alpha"
        assert spec.order == "AN00000001"
        assert spec.prefix == "orglab-alpha/AN00000001/"
        assert spec.s3_path == "s3://czi-psomagen/orglab-alpha/AN00000001"

    def test_s3_uri_form(self):
        assert (
            resolve_spec("s3://czi-novogene/somelab-atlas/NVUS0000000000-36").key
            == "czi-novogene/somelab-atlas/NVUS0000000000-36"
        )

    def test_trailing_slash_and_whitespace_tolerated(self):
        assert (
            resolve_spec("  czi-psomagen/orglab-alpha/AN00000001/  ").order
            == "AN00000001"
        )

    def test_deeper_path_resolves_to_its_order_not_the_tail(self):
        """A pasted object key must resolve to the order that contains it.

        Reading project/order backwards from the tail turned this into
        'czi-psomagen/raw/x.fastq.gz' and QA'd a prefix nobody named.
        """
        spec = resolve_spec(
            "czi-psomagen/orglab-alpha/AN00000001/raw/GEX/x_R1_001.fastq.gz"
        )
        assert spec.project == "orglab-alpha"
        assert spec.order == "AN00000001"
        assert spec.prefix == "orglab-alpha/AN00000001/"

    def test_prefix_has_trailing_slash_so_it_cannot_match_a_sibling(self):
        """AN00000002 must not list AN000000020."""
        assert resolve_spec("b/p/AN00000002").prefix.endswith("/")

    @pytest.mark.parametrize("bad", ["", "   ", "http://x/y/z", "czi-psomagen/proj"])
    def test_rejects_unusable_specs(self, bad):
        with pytest.raises(SpecError):
            resolve_spec(bad)

    @pytest.mark.parametrize(
        "name", ["jsonlds", "metadata-tsv", "test", "test_orglab-alpha"]
    )
    def test_rejects_known_non_orders(self, name):
        with pytest.raises(SpecError, match="not an order"):
            resolve_spec(f"czi-psomagen/orglab-beta/{name}")

    def test_bare_order_id_is_ambiguous_rather_than_guessed(self):
        """Three watched projects means a bare ID names three possible orders."""
        with pytest.raises(SpecError, match="ambiguous"):
            resolve_spec("AN00000001")


class TestNonOrderNames:
    @pytest.mark.parametrize(
        "name",
        ["jsonlds", "metadata-tsv", "test", "test.txt", "test_orglab-alpha", "TEST"],
    )
    def test_known_non_orders(self, name):
        assert is_non_order_name(name)

    @pytest.mark.parametrize(
        "name", ["AN00000001", "NVUS0000000000-36", "REF3", "testis-atlas"]
    )
    def test_real_orders_are_not_filtered(self, name):
        """'testis-atlas' starts with 'test' but is not a smoke test."""
        assert not is_non_order_name(name)


class TestOrderShape:
    @pytest.mark.parametrize(
        "order,shape",
        [
            ("AN00000001", "psomagen_AN"),
            ("AN00000002", "psomagen_AN"),
            ("NVUS0000000000-36", "novogene_NVUS"),
            ("REF3", "unrecognized"),
            ("AN0002835", "unrecognized"),
        ],
    )
    def test_classification(self, order, shape):
        assert classify_order_shape(order) == shape

    def test_unrecognized_shape_still_resolves(self):
        """A new vendor format must not be refused; it is reported, not enforced."""
        spec = resolve_spec("czi-ultima/somelab/WHATEVER-2026")
        assert spec.order_shape == "unrecognized"
        assert spec.order == "WHATEVER-2026"


class TestOrderSpecIdentity:
    def test_identity_includes_bucket_and_project(self):
        """Order ID alone does not name a thing to QA."""
        a = OrderSpec("czi-psomagen", "proj-a", "AN1", "psomagen_AN")
        b = OrderSpec("czi-novogene", "proj-a", "AN1", "psomagen_AN")
        assert a.key != b.key
