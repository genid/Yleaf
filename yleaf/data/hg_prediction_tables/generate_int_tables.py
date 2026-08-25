"""Regenerate the per-backbone expected-state tables (``<node>_int.txt``).

Each tree directory ships one ``<node>_int.txt`` per backbone node listed in its
``Intermediates.txt``. ``predict_haplogroup.get_qc1_score()`` reads the table of
the candidate's most specific backbone ancestor and compares the observed state
of every backbone node against the expected one.

Generation rule
---------------
For backbone node ``X``, the expected state of backbone node ``Y`` is ``D`` when
``Y`` is an ancestor of ``X`` or ``Y == X``; otherwise ``A``. Rows are written in
the order the nodes appear in ``Intermediates.txt``, one ``name<TAB>state`` pair
per line.

Because every table enumerates *all* backbone nodes, adding a node to
``Intermediates.txt`` invalidates every existing table — they all need the new
row. Run this script after any change to a tree's ``Intermediates.txt``.

Usage
-----
    python generate_int_tables.py --tree isogg           # rewrite the tables
    python generate_int_tables.py --tree isogg --check   # report, write nothing

A node listed in ``Intermediates.txt`` that is absent from the tree JSON (the
ISOGG tree has one such phantom, ``A``, kept for backwards compatibility) is
only ever derived with respect to itself.

The rule above reproduces the shipped ISOGG tables byte-identically. The YFull
tables record the same states but list all ``D`` rows before all ``A`` rows, so
regenerating them would reorder every row without changing any state; the FTDNA
tables encode different expectations altogether (that tree is flat: all 20
backbone nodes are direct children of the root). Only ``--tree isogg`` is
maintained by this script today.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

HERE = Path(__file__).resolve().parent

# tree name -> (tables directory, tree JSON file)
TREES = {
    "yfull": ("major_tables", "tree.json"),
    "yfull_v10": ("yfull_v10_major_tables", "yfull_v10_tree.json"),
    "ftdna": ("ftdna_major_tables", "ftdna_tree.json"),
    "isogg": ("isogg_major_tables", "isogg_tree.json"),
}


def read_backbone(tables_dir: Path) -> List[str]:
    """Backbone node names in file order. Mirrors read_backbone_groups(): lines
    containing '~' (approximate placement) are not backbone nodes."""
    names = []
    with open(tables_dir / "Intermediates.txt") as f:
        for line in f:
            if "~" in line:
                continue
            name = line.strip()
            if name and name not in names:
                names.append(name)
    return names


def read_parents(tree_path: Path) -> Tuple[Dict[str, str], str]:
    """Return (child -> parent, root) for an adjacency-list tree JSON."""
    adjacency = json.loads(tree_path.read_text())
    parents = {}
    for parent, children in adjacency.items():
        for child in children or []:
            parents[child] = parent
    roots = [node for node in adjacency if node not in parents]
    if len(roots) != 1:
        sys.exit(f"expected exactly one root in {tree_path.name}, found {roots}")
    return parents, roots[0]


def render_table(node: str, backbone: List[str], parents: Dict[str, str], root: str) -> str:
    """The <node>_int.txt content: every backbone node with its expected state."""
    derived = set()
    current = node
    while current is not None and current != root:
        derived.add(current)
        current = parents.get(current)
    return "".join(f"{name}\t{'D' if name in derived else 'A'}\n" for name in backbone)


def read_table(path: Path) -> Dict[str, str]:
    states = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        name, state = line.split("\t")
        states[name] = state
    return states


def check_existing_rows_unchanged(tables_dir: Path, tables: Dict[str, str]) -> None:
    """Regression guard: for every table that already exists, the states it
    currently records must be preserved by the regenerated version. Only the
    rows that were already present are compared — new backbone nodes add rows
    but must never change an existing one."""
    for node, new_text in tables.items():
        path = tables_dir / f"{node}_int.txt"
        if not path.exists():
            continue
        old_states = read_table(path)
        new_states = {}
        for line in new_text.splitlines():
            name, state = line.split("\t")
            new_states[name] = state
        for name, state in old_states.items():
            if new_states.get(name) != state:
                sys.exit(f"REGRESSION: {path.name} row '{name}' was '{state}' but "
                         f"would become '{new_states.get(name)}' — refusing to write")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", required=True, choices=sorted(TREES),
                        help="which tree's tables to regenerate")
    parser.add_argument("--check", action="store_true",
                        help="report what would change without writing")
    args = parser.parse_args()

    subdir, tree_file = TREES[args.tree]
    tables_dir = HERE / subdir
    backbone = read_backbone(tables_dir)
    parents, root = read_parents(HERE / tree_file)

    missing = [node for node in backbone if node not in parents]
    if missing:
        print(f"note: backbone nodes absent from {tree_file}: {missing}")

    tables = {node: render_table(node, backbone, parents, root) for node in backbone}
    check_existing_rows_unchanged(tables_dir, tables)

    created, changed = [], []
    for node, text in tables.items():
        path = tables_dir / f"{node}_int.txt"
        if not path.exists():
            created.append(path.name)
        elif path.read_text() != text:
            changed.append(path.name)

    print(f"{args.tree}: {len(backbone)} backbone nodes")
    print(f"  to create : {len(created)} {created}")
    print(f"  to rewrite: {len(changed)}")

    if args.check:
        print("  --check given: nothing written")
        return 0

    for node, text in tables.items():
        (tables_dir / f"{node}_int.txt").write_text(text)
    print(f"  wrote {len(tables)} tables in {tables_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
