## graph_db2

Visualize the DB2 reference graph. Any object in the database is a node; any
`/object_type/uuid/` reference in its JSON is an edge.

There are two ways to draw it:

| | Entry point | Shape |
|---|---|---|
| **Interactive explorer** | `python -m graph_db2` | Dash + Cytoscape app. Starts at one object and grows only where you click. |
| **Static snapshot** | `graphing.py` in a notebook | Walks the whole reachable graph up front and writes a self-contained pyvis HTML file. |

The explorer is the one to reach for when you don't know how big the graph is.
It resolves neighbors on demand, so the seed can sit anywhere in the database
without materializing everything behind it. The static path is better when you
want a file to hand to someone.

---

## Requirements

- **Python 3.11+**
- **`db2_flattener`**, installed from the sibling
  [db2-flattener](https://github.com/Lattice-Data/db2-flattener) repo. Supplies
  `DB2Gatherer` for batched reports, `SchemaIDs` for parsing object paths, and
  `constants.yaml` for the per-server profile config.
- **DB2 API credentials** as environment variables.
- **`dash`** and **`dash-cytoscape`** for the explorer; **`pyvis`** for the
  static snapshot.

### Install

From the `db2-flattener` checkout:

```bash
pip install -e ".[dev]"
```

Then the visualization dependencies:

```bash
pip install dash dash-cytoscape pyvis
```

Verified against `dash` 4.4.1, `dash-cytoscape` 1.0.2, `pyvis` 0.3.2. Dash 3+
is required — the app calls `app.run()`, not the removed `run_server()`.

### Credentials

Variable names follow the `--mode` value, which must start with `DB2_`:

```bash
export DB2_PROD_KEY=...
export DB2_PROD_SECRET=...
export DB2_PROD_SERVER=https://api.data.lattice-data.org/

export DB2_DEMO_KEY=...
export DB2_DEMO_SECRET=...
export DB2_DEMO_SERVER=https://lattice-api-dev.demo.lattice-data.org/
```

If you keep these in a conda env, `conda env config vars list` shows what is
set. A missing triple exits with a message naming the variables it wanted.

Both servers above have an entry in the flattener's `constants.yaml`. Pointing
`--mode` at anything else raises a `KeyError` listing the known endpoints; use
`--fetch-new` to pull profile schemas from the instance directly instead.

---

## Running the explorer

Run from the **`bcp`** directory — the package uses relative imports, so it has
to be invoked as a module, not as a loose script:

```bash
cd bcp
python -m graph_db2 --help
```

```bash
python -m graph_db2
python -m graph_db2 --mode db2_demo --seed /matrix_file_sets/<uuid>/
```

Then open <http://localhost:8050>.

| Flag | Description |
|------|-------------|
| `--seed` | Object path to start from. Defaults to a sample on `db2_prod`; on any other mode, defaults to an empty canvas. |
| `--mode` | DB2 instance (default `db2_prod`). Must start with `db2_`. |
| `--fetch-new` | Fetch profile schemas from the instance instead of reading `constants.yaml`. |
| `--port` | Default `8050`. |
| `--debug` | Dash debug mode with hot reload. |

`--seed` accepts any of these spellings, and normalizes to the canonical form:

```
/object_type/uuid/
object_type/uuid
object_type/uuid/
/object_type/uuid
```

A seed only means something on the server it came from, so `--mode` **without**
`--seed` deliberately starts empty rather than carrying a path across
deployments. Anything malformed is rejected up front with
`Expected 'object_type/uuid', got '...'`.

### Using it

- **Click a node** to resolve its neighbors and fold them into the graph. A
  dashed ring means "not expanded yet"; a solid ring means already resolved.
- **Node colours** come from the `NodeColor` enum in `models.py`, keyed by the
  abstract class (so `Tissue`, `CellLine`, and `Organoid` all read as
  `Biosample`). Unmapped types fall back to grey.
- **Wide fans collapse.** A fan larger than the "draw up to N neighbors"
  budget (default `DRAW_BUDGET = 25`) puts each oversized type behind one
  placeholder like `RawMatrixFile × 512`. Clicking a placeholder opens a
  searchable multi-select in the side panel — tick as many members as you want
  and only those land on the canvas. "Fan out all N" draws the lot.

  This is not cosmetic: a 512-wide fan lays out roughly 23,000px tall and is
  unreadable at any zoom, and drawing it costs 512 label fetches.
- **Layout** is picked from the graph's shape on load — `concentric` when one
  node touches ≥80% of the others (a star), `dagre` otherwise (lineage
  chains). Change it from the dropdown at any time; expansions never override
  your choice.
- **"Show types"** at the bottom of the panel toggles whole node types.
- The status line reports real counts — `64 drawn`, `512 RawMatrixFile
  grouped` — and turns red on failure, since an empty canvas otherwise looks
  identical to a silent 404.

### Cost per click

One authenticated GET for the clicked object, plus one batched report per
neighbor type actually drawn. Grouped types cost nothing until you open them,
because an object's type is readable from its path without fetching it.

Objects referenced by types in `EXCLUDED_SCHEMAS` (`labs`, `users`, `terms`,
`controlled_terms`, …) are skipped, which is why a MatrixFileSet whose JSON
names a lab and a submitter can still show only its files.

---

## Static pyvis snapshot

`graphing.py` keeps the original whole-graph approach. Call
`fetch_all_references_from_cache()` first — it batch-loads every reachable
object into the `LatticeNode` cache, and skipping it makes the adjacency walk
request one object at a time (~40× slower).

Beyond visualization, this can be used to generate the full graph for other
potential QA/inspection tools that make use of graph nodes and edges.

```python
from graph_db2.connection import get_connection
from graph_db2.graphing import (
    build_adjacency_dict,
    create_nodes_for_pyvis,
    find_edges,
    graph_dict,
)
from graph_db2.models import LatticeNode, fetch_all_references_from_cache
from graph_db2.schema import create_configs

from db2_flattener.gather.gatherer import DB2Gatherer
from pyvis import network as net

mode = "db2_prod"
gatherer = DB2Gatherer(get_connection(mode), create_configs(mode))

# LatticeNode.mode is a ClassVar, not a constructor argument: set it once and
# every node in the process talks to that server
LatticeNode.mode = mode

start = LatticeNode("/matrix_file_sets/<uuid>/")
start.object_json  # seeds the cache
fetch_all_references_from_cache(start, gatherer)  # batch-load the rest
build_adjacency_dict(start)

graph = net.Network(notebook=True, cdn_resources="remote", select_menu=True)
for node in create_nodes_for_pyvis(graph_dict):
    graph.add_node(**node)
graph.add_edges(find_edges(graph_dict))
graph.show("graph.html")
```

`graph_dict` is a module-level global that accumulates across calls — reset it
between unrelated graphs.

---

## Module layout

| File | Role |
|------|------|
| `cli.py` / `__main__.py` | Argument parsing and entry point for `python -m graph_db2`. |
| `explorer.py` | Dash app: layout, stylesheet, callbacks, layout heuristic. |
| `cyto_elements.py` | Graph logic with no Dash dependency — path normalization, one-hop expansion, grouping, element construction. |
| `models.py` | `LatticeNode` (lazy, API-backed, class-level cache), `NodeColor`, batch-request helpers. |
| `graphing.py` | Whole-graph walk and pyvis element construction. |
| `connection.py` | Cached `Connection` per mode. |
| `schema.py` | Builds `Configs` from `constants.yaml` or live profiles. |
| `constants.py` | Defaults, `EXCLUDED_SCHEMAS`, `ABSTRACT_MAPPING`. |
| `graphing_playground.ipynb` | Scratch notebook for the pyvis path. |

`cyto_elements.py` holds no Dash imports on purpose, so expansion and grouping
can be driven from a notebook or a test without starting a server.

---

## Known limitations

- **Edges are undirected.** `LatticeNode.get_ids()` walks the JSON on both
  sides of a reference, so a Library and its FileSet each name the other.
  Choosing a direction from discovery order would draw arrows that don't mean
  lineage. True lineage arrows would need specific fields (`derived_from`,
  `libraries`) rather than the generic reference walk.
- **Single process only.** `LatticeNode._cache` and the explorer's
  `_fully_fetched` set are class- and module-level globals, and they are not
  keyed by mode. One Dash worker, one mode per process. Adding an in-app mode
  switcher would need a mode dimension on both, or prod and demo objects with
  colliding paths will cross-contaminate.
- **Re-expanding a node whose group you already fanned out** re-creates the
  placeholder.
- **Full URLs are not accepted** as a seed — only `object_type/uuid` paths, not
  `https://api.data.lattice-data.org/matrix_file_sets/<uuid>/`.
- The dev server is Flask's. Fine for a local tool; put it behind a real WSGI
  server if it ever gets shared.
