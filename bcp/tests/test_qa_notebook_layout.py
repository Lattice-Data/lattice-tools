"""
Structural invariants of qa.ipynb.

Why this module exists
---------------------
The SeaHub cells were originally inserted in the middle of the notebook, between
the "extra raw files" cell and the CellRanger cell. That moved CellRanger
validation from cell index 13 to 17, and a partially-executed run therefore
produced a ``*_errors.txt`` missing every CellRanger finding -- which was
reported as a 10x regression and cost a day to trace back to layout rather than
code.

The fix is that all SeaHub cells form one contiguous block at the very end,
tagged ``seahub``. These tests pin that, so inserting an assay-specific cell in
the middle of the shared pipeline fails loudly instead of silently renumbering
everything after it.

They check cell ORDER and dataflow only. Whether a reordered notebook still runs
end to end is a separate question that only executing it can answer.
"""

from __future__ import annotations

import ast
import json
import os

import pytest

NOTEBOOK_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "qa.ipynb")
SEAHUB_TAG = "seahub"


@pytest.fixture(scope="module")
def notebook():
    with open(NOTEBOOK_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def cells(notebook):
    return notebook["cells"]


def _source(cell) -> str:
    return "".join(cell["source"])


def _tags(cell) -> list[str]:
    return cell.get("metadata", {}).get("tags", [])


def _seahub_indices(cells) -> list[int]:
    return [i for i, c in enumerate(cells) if SEAHUB_TAG in _tags(c)]


def _bound_names(source: str) -> set[str]:
    """Every name a cell binds: assignment, loop target, comprehension, import."""
    tree = ast.parse(source)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            names.add(node.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names.update((a.asname or a.name).split(".")[0] for a in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
    return names


def _loaded_names(source: str) -> set[str]:
    tree = ast.parse(source)
    return {
        node.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
    }


def _index_of(cells, needle: str) -> int:
    matches = [i for i, c in enumerate(cells) if needle in _source(c)]
    assert len(matches) == 1, (
        f"expected exactly one cell containing {needle!r}, got {matches}"
    )
    return matches[0]


class TestSeahubBlockIsTrailing:
    def test_the_seahub_block_is_non_empty_and_tagged(self, cells):
        assert _seahub_indices(cells), (
            "no cell carries the 'seahub' tag; SeaHub cells must be tagged so "
            "this invariant is checkable"
        )

    def test_the_seahub_block_is_contiguous_and_at_the_end(self, cells):
        indices = _seahub_indices(cells)
        assert indices == list(range(indices[0], len(cells))), (
            f"seahub-tagged cells are at {indices} but the block must run to the "
            f"end of the notebook ({len(cells)} cells). An assay-specific cell in "
            "the middle of the shared pipeline renumbers every cell after it."
        )

    def test_cellranger_and_scale_validation_precede_the_seahub_block(self, cells):
        first_seahub = _seahub_indices(cells)[0]
        cellranger = _index_of(cells, "### Go through cellranger logs")
        scale_first = _index_of(cells, "### Scale cb_tag validation")

        assert cellranger < first_seahub
        assert scale_first < first_seahub

    def test_no_untagged_cell_is_seahub_specific(self, cells):
        """A SeaHub cell that forgot its tag would defeat every test above.

        Five cells legitimately mention SeaHub outside the block: cell 0 imports
        the SeaHub modules for the whole notebook, and the parameter block, the
        processed-validation override, the Q30 skip set and the wafer table's
        discovered_wafers comment each reference the assay. Anything else naming
        SeaHub belongs in the tagged block.
        """
        allowed_untagged = {
            0,
            _index_of(cells, "# Choose data source:"),
            _index_of(cells, "### Gather all raw + processed file listings"),
            _index_of(cells, "### Validate sublibrary Q30"),
            _index_of(cells, "### Wafer-level summary:"),
        }
        offenders = [
            i
            for i, c in enumerate(cells)
            if SEAHUB_TAG not in _tags(c)
            and i not in allowed_untagged
            and "seahub" in _source(c).lower()
        ]
        assert offenders == [], f"untagged SeaHub-specific cells at {offenders}"


class TestSharedPipelineIsIndependentOfTheSeahubBlock:
    def test_no_shared_cell_reads_a_name_bound_only_in_the_seahub_block(self, cells):
        """The dataflow condition that makes the trailing placement legal.

        If a shared cell consumed something only the SeaHub block defines, moving
        the block to the end would break every non-SeaHub run.

        Names the shared cells bind for themselves don't count, however common
        they are in the SeaHub block too -- a loop variable called ``group`` in
        both places is not a dependency. Hence the AST rather than a regex, which
        also stops comments and string literals registering as uses.
        """
        seahub = set(_seahub_indices(cells))
        shared = [
            c
            for i, c in enumerate(cells)
            if i not in seahub and c["cell_type"] == "code"
        ]

        seahub_binds = set()
        for i in seahub:
            if cells[i]["cell_type"] == "code":
                seahub_binds |= _bound_names(_source(cells[i]))
        assert seahub_binds, "the SeaHub block binds no names; the check is vacuous"

        shared_binds: set[str] = set()
        for cell in shared:
            shared_binds |= _bound_names(_source(cell))

        seahub_only = seahub_binds - shared_binds
        offenders = {
            i: sorted(_loaded_names(_source(cells[i])) & seahub_only)
            for i, c in enumerate(cells)
            if i not in seahub
            and c["cell_type"] == "code"
            and _loaded_names(_source(c)) & seahub_only
        }
        assert offenders == {}, f"shared cells depend on SeaHub-only names: {offenders}"

    def test_the_processed_validation_override_stays_in_the_gather_cell(self, cells):
        """``validate_processed = False`` for SeaHub must precede the branches.

        Cells 13-16 branch on validate_processed, so this override cannot move
        into the trailing SeaHub block even though it is SeaHub-specific.
        """
        gather = _index_of(cells, "### Gather all raw + processed file listings")
        source = _source(cells[gather])

        # Quote-agnostic: notebooks are excluded from ruff (pyproject.toml), so
        # their quote style is whatever the author used and must not be asserted.
        assert "raw_assay ==" in source
        assert "seahub_sci" in source
        assert "validate_processed = False" in source
        assert gather < _seahub_indices(cells)[0]


class TestParameterBlock:
    def test_the_untrimmed_source_parameter_is_a_list(self, cells):
        """One experiment spans several vendor orders, so it cannot be a scalar."""
        params = _source(cells[_index_of(cells, "# Choose data source:")])

        assert "untrimmed_s3_paths = []" in params
        assert "untrimmed_s3_path =" not in params

    def test_the_seahub_block_gates_on_seahub_active(self, cells):
        """Every SeaHub code cell must self-skip for other assays."""
        for i in _seahub_indices(cells):
            if cells[i]["cell_type"] != "code":
                continue
            assert "seahub_active" in _source(cells[i]), f"cell {i} has no gate"


class TestCellZero:
    def test_every_import_in_the_first_cell_resolves(self, cells):
        """Cell 0 is shared by every assay: an ImportError here kills the notebook.

        Only the import statements are executed -- constructing the boto3 client
        would require AWS configuration that a unit test should not need.
        """
        source = _source(cells[0])
        tree = ast.parse(source)
        imports = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]

        assert imports, "cell 0 has no imports; the notebook layout changed"

        module = ast.Module(body=imports, type_ignores=[])
        exec(compile(module, filename="qa.ipynb#cell0", mode="exec"), {})

    def test_the_first_cell_is_imports_only(self, cells):
        """Keeps cell 0 cheap to validate and free of side effects beyond the client."""
        tree = ast.parse(_source(cells[0]))
        non_imports = [
            n for n in tree.body if not isinstance(n, (ast.Import, ast.ImportFrom))
        ]
        assert [type(n).__name__ for n in non_imports] == ["Assign"], (
            "cell 0 should contain imports plus the single s3client assignment"
        )
