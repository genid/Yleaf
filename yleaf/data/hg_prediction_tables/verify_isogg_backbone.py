"""Verify the ISOGG backbone covers every lineage in the ISOGG tree.

``predict_haplogroup.get_qc1_score()`` scores a candidate by walking its path to
the root until it finds a backbone node; if it finds none, QC1 is forced to 0 and
the candidate is discarded, so the sample is predicted ``NA``.

The ISOGG ``Intermediates.txt`` used to list only the 20 single-letter clades
A-T. 49 nodes — the deep-African A0/A1/A00 region and the megahaplogroup
internal nodes BT, CT, CF and DE — have no single-letter ancestor and could
therefore never be predicted. Adding ``A0-T`` (the ancestor of every non-root
node) closes the gap.

This script asserts that state and the table invariants that go with it. Run it
after regenerating the tables with ``generate_int_tables.py``:

    PYTHONPATH=<repo> python verify_isogg_backbone.py

Exit code 0 means every assertion held.
"""

import json
import sys
from pathlib import Path

from yleaf import predict_haplogroup as ph
from yleaf import yleaf_constants as C
from yleaf.tree import Tree

TABLES = C.HG_PREDICTION_FOLDER
ISOGG_TABLES = TABLES / "isogg_major_tables"
ROOT = 'ROOT (Y-Chromosome "Adam")'
LETTERS = set("ABCDEFGHIJKLMNOPQRST")
FAILED = []


def check(label, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    if not condition:
        FAILED.append(label)


def read_tree():
    adjacency = json.load(open(TABLES / C.ISOGG_TREE_FILE))
    parents = {child: parent for parent, kids in adjacency.items() for child in (kids or [])}
    return set(adjacency) | set(parents), parents


def path_to_root(node, parents):
    path = [node]
    while path[-1] in parents:
        path.append(parents[path[-1]])
    return path


def main():
    nodes, parents = read_tree()
    ph.read_backbone_groups(TABLES / C.ISOGG_TREE_FILE)
    backbone = set(ph.MAIN_HAPLO_GROUPS)

    print(f"\n== 1. Backbone set as loaded ({len(backbone)} entries)")
    print(f"  {sorted(backbone)}")
    check("backbone is the 20 single-letter clades plus A0-T",
          backbone == LETTERS | {"A0-T"})

    print("\n== 2. QC1 finds a backbone ancestor for every lineage")
    tree = Tree(TABLES / C.ISOGG_TREE_FILE)

    def qc1(name):
        node = tree.get(name)
        path = [node.name]
        while node.parent is not None:
            node = node.parent
            path.append(node.name)
        return ph.get_qc1_score(path, {}), path

    # previously hard-zero: deep-African lineages and the megahaplogroup nodes
    for name in ("A0b", "A1a", "A0", "A1b1a1", "BT", "CT", "CF", "DE"):
        score, path = qc1(name)
        hit = next((v for v in path if v in backbone), None)
        check(f"{name:9s} QC1={score!r:5s} backbone_hit={hit!r}",
              score >= 0.95 and hit is not None)
    # unaffected lineages must be unchanged
    for name in ("B2b1", "R1b1a1b", "E1b1a1"):
        score, path = qc1(name)
        hit = next((v for v in path if v in backbone), None)
        check(f"{name:9s} QC1={score!r:5s} backbone_hit={hit!r} (unchanged)",
              score == 1.0 and hit in LETTERS)

    print("\n== 3. Every tree node has a backbone ancestor")
    orphans = [n for n in nodes
               if n != ROOT and not any(v in backbone for v in path_to_root(n, parents))]
    check("no node is unreachable from the backbone", not orphans,
          f"n={len(orphans)} of {len(nodes)}" + (f" {sorted(orphans)[:12]}" if orphans else ""))

    print("\n== 4. 'A' is a phantom: in the backbone set but not in the tree")
    missing = sorted(b for b in backbone if b not in nodes)
    check("'A' is the only backbone entry absent from isogg_tree.json",
          missing == ["A"], f"missing={missing}")
    positions = C.DATA_FOLDER / "hg38" / "isogg_positions_hg38.txt"
    if positions.exists():
        labelled_a = sum(1 for line in open(positions)
                         if len(line.split("\t")) > 2 and line.split("\t")[2] == "A")
        check("no marker in isogg_positions_hg38.txt is labelled exactly 'A'",
              labelled_a == 0, f"count={labelled_a}")

    print("\n== 5. Expected-state tables: one per backbone node")
    have = sorted(f.name[:-len("_int.txt")] for f in ISOGG_TABLES.glob("*_int.txt"))
    check(f"every backbone node has an _int.txt ({len(have)} files)",
          set(have) == backbone, f"missing={sorted(backbone - set(have))}")

    print("\n== 6. Tables follow the generation rule")
    print("  For backbone X: expected state of backbone Y is 'D' iff Y is an")
    print("  ancestor of X or Y == X; otherwise 'A'.")
    ok = 0
    for name in sorted(backbone):
        table = {}
        for line in open(ISOGG_TABLES / f"{name}_int.txt"):
            if line.strip():
                key, state = line.strip().split("\t")
                table[key] = state
        ancestors = set(path_to_root(name, parents)) if name in nodes else {name}
        expected = {y: ("D" if y in ancestors else "A") for y in backbone}
        if expected == table:
            ok += 1
        else:
            diff = {k: (table.get(k), expected.get(k))
                    for k in set(table) | set(expected) if table.get(k) != expected.get(k)}
            print(f"    mismatch {name}: {diff}")
    check(f"rule reproduces all tables ({ok}/{len(backbone)})", ok == len(backbone))

    print("\n== 7. Other trees are untouched")
    for subdir, expected_n in (("major_tables", 56), ("yfull_v10_major_tables", 56),
                               ("ftdna_major_tables", 20), ("isogg_major_tables", 21)):
        directory = TABLES / subdir
        n_names = len([l for l in open(directory / "Intermediates.txt")
                       if l.strip() and "~" not in l])
        n_int = len(list(directory.glob("*_int.txt")))
        check(f"{subdir:24s} Intermediates={n_names:3d}  _int.txt={n_int:3d}",
              n_names == expected_n)

    print("\n" + "=" * 68)
    if FAILED:
        print(f"{len(FAILED)} CHECK(S) FAILED: {FAILED}")
        return 1
    print("All checks passed: ISOGG backbone covers the whole tree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
