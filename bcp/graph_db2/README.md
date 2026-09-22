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
| `--seed` | Object to start from, as a path, an alias or a bare uuid. Defaults to a sample on `db2_prod`; on any other mode, defaults to an empty canvas. |
| `--mode` | DB2 instance (default `db2_prod`). Must start with `db2_`. |
| `--fetch-new` | Fetch profile schemas from the instance instead of reading `constants.yaml`. |
| `Hold View` | Toolbar checkbox (not a flag). Stops the canvas refitting on every draw — see [Using it](#using-it). |
| `Hold Layout` | Toolbar checkbox (not a flag). Stops the layout re-running, so drawn nodes stay put — see [Using it](#using-it). |
| `columns by type` | Layout dropdown entry (not a flag). Vertical columns in pipeline order — see [Using it](#using-it). |
| `--port` | Default `8050`. |
| `--debug` | Dash debug mode with hot reload. |

`--seed` takes any of the three things a DB2 instance will resolve, with or
without surrounding slashes and whitespace:

```
/object_type/uuid/          matrix_file_sets/<uuid>
/alias/                     test-lab:my_matrixfileset
/uuid/                      <uuid>
```

An alias is usually what a curator has to hand, and a bare uuid is what a
report or a Jira ticket carries — but neither can be a node id, because the
graph is keyed by `@id`. So the seed is requested once, the profile that comes
back is cached, and its `@id` becomes the canonical form; the seed box and the
status line both show that rather than what was typed.

Because only the server can tell an alias from a uuid, almost nothing is
rejected on shape — a bare word could be an alias. The exceptions are an empty
seed, a full URL (`https://api.data.lattice-data.org/...`), and an empty path
segment (`matrix_file_sets//<uuid>`), each of which would otherwise become a
request for something the user never typed. Everything else that fails — an
unknown alias, an object the key cannot see, a collection endpoint like
`matrix_file_sets` that answers `200` without naming one object — comes back as
a single `GraphDB2Error` naming the seed and the server it was tried on.

A seed only means something on the server it came from, so `--mode` **without**
`--seed` deliberately starts empty rather than carrying a path across
deployments — and resolution always goes to the `--mode` server, never the
default one.

### Using it

- **Click a node** to resolve its neighbors and fold them into the graph. A
  dashed ring means "not expanded yet"; a solid ring means already resolved.
  This is a property of the current canvas, not of the fetch cache, so pressing
  **Load** again makes every node clickable from scratch — cheaply, since the
  profiles are still cached.
- **`Hold View`** stops the canvas re-centering and zooming to fit every time
  nodes are drawn. Off by default. On a large graph the refit is what loses your
  place: you zoom into one corner, click a node, and the canvas snaps back out
  to show everything. Tick the box and the viewport stays exactly where you put
  it; untick it to fit the whole graph again.

  Turning it on does not strand you. dash-cytoscape fits the viewport itself
  when new elements land *entirely* off screen, so loading a seed somewhere else
  on the canvas still snaps to it.

  It holds your window, not the positions inside it — for those, tick
  `Hold Layout` as well.
- **`Hold Layout`** stops the graph moving under you. Off by default, and
  independent of `Hold View`: that one holds the viewport, this one holds the
  arrangement. Normally every add or remove re-runs the chosen layout, and
  dagre re-flows the *whole* graph when one node arrives — so expanding a leaf
  rearranges everything you were reading. Tick the box and nodes and edges
  already drawn stay exactly where they are, including any you dragged there
  by hand.

  New nodes have to go somewhere, and cytoscape drops a node with no position
  at the origin, so every expansion places its own: a ring around the node you
  clicked, spilling onto further rings for a wide fan. (This happens whether
  or not the box is ticked — a running layout is free to overrule the ring,
  but it does not run until 100ms after the nodes land, which is long enough
  to watch them pile into the corner.) They are not laid out — nothing stops
  one landing on top of something else — so a held session eventually wants
  tidying.

  Three things still re-arrange the graph, all of them things you asked for:
  picking a layout from the dropdown, unticking the box, and pressing **Load**
  (a new seed has no positions to hold, so it gets one layout run). Holding
  also outranks `columns by type`, which otherwise re-columns on every draw.

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

  **The ticks are the canvas.** A member that is drawn shows as ticked, and
  unticking it takes it back off — so the same picker prunes a fan as well as
  builds one, including after "Fan out all N", which comes back with every box
  ticked. The picker applies the *difference* between its selection and the
  canvas rather than the selection itself, so only newly ticked members cost a
  fetch and re-firing with an unchanged selection does nothing.

  Removal takes the member and its edges, and nothing else. Anything you
  expanded *from* that member stays put as an isolated node rather than being
  cascaded away — losing a subtree to an untick would be far worse than seeing
  one float.
- **Layout** is picked from the graph's shape on load — `concentric` when one
  node touches ≥80% of the others (a star), `dagre` otherwise (lineage
  chains). Change it from the dropdown at any time; expansions never override
  your choice. Picking one always re-runs it, `Hold Layout` or not.
- **`columns by type`** draws the graph as vertical columns in pipeline
  order, left to right:

  ```
  Donor → Biosample → GeneticModification/Treatment/ExperimentalCondition →
  Library → SequenceFile → SequenceFileSet →
  RawMatrixFile/ProcessedMatrixFile/TabularFile → MatrixFileSet
  ```

  A column is a *legend bucket*, not an `api_name`, so `Tissue`, `CellLine`
  and `Organoid` line up together the same way they share a colour, and a new
  Biosample subtype needs no code change. Anything the legend does not map
  gets a trailing column of its own rather than being hidden in someone
  else's. Empty columns are packed out, so a graph with no libraries has no
  gutter where they would have been.

  Rows are not alphabetical: after an initial sort by label, two barycentre
  passes pull each node level with its neighbours in the columns either side.
  Without them a column of 30 files faces its donors in an unrelated order and
  every edge crosses every other.

  Unlike the others this one is computed in Python (`column_positions()`) and
  handed to cytoscape as a `preset` — none of the bundled layouts can group by
  an attribute, and dagre ranks by distance from a root, which puts a Tissue
  and a SequenceFile in the same column whenever the path lengths happen to
  match. Being a positions map, it has to be rebuilt whenever the canvas
  changes, so it re-runs on every draw — which also means a node you drag
  snaps back on the next one. Tick `Hold Layout` if you want your own
  arrangement to stick.
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
- **An expansion places new nodes, it does not lay them out.** They ring the
  node you expanded without consulting the rest of the canvas, so they can
  land on top of something already drawn. Under `Hold Layout` that is where
  they stay: a layout that arranged only the new nodes would need the whole
  graph as fixed constraints, which none of the bundled cytoscape layouts
  take.
- **`columns by type` re-columns on every draw**, which costs a round trip
  and a positions map for the whole graph each time, and undoes anything you
  dragged. It is a `preset`, so there is nothing for cytoscape to re-run
  incrementally. Tick `Hold Layout` to stop it.
- **Unticking a member removes it even if another expansion drew it too.** The
  picker's ticks mean "on the canvas", and a node is on the canvas once
  regardless of how many paths led to it. Unticking it also removes the edges
  those other expansions contributed; re-expanding the neighbor puts them back.
- **Full URLs are not accepted** as a seed — a path, alias or uuid, not
  `https://api.data.lattice-data.org/matrix_file_sets/<uuid>/`.
- **A seed costs one extra request the first time.** Resolving it is a GET, and
  the full profile fetch is another. Only a path already in the cache — which
  is every node on the canvas — resolves for free, so clicking around does not
  pay for it. An alias is not a cache key, so re-Loading the same alias
  re-resolves it.
- The dev server is Flask's. Fine for a local tool; put it behind a real WSGI
  server if it ever gets shared.
