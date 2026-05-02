#!/usr/bin/env python

"""
Mixture analysis for Yleaf: detect multiple Y-chromosome contributors in a sample.

Algorithm
---------
1. Load "Below base majority" positions from the .fmf file — these are positions
   where both alleles have meaningful read support (heterozygous under mixture).
2. Compute the derived-allele fraction (der_frac) per position:
     state='D'  → der_frac = called_perc / 100   (derived is majority)
     state='A'  → der_frac = 1 - called_perc / 100  (derived is minority)
3. Find the Lowest Common Ancestor (LCA) of all het-position-bearing nodes — this
   is the deepest shared ancestor of all contributors.
4. Recursively decompose the tree below the LCA: at each branching point, child
   subtrees that contain het positions represent distinct contributor paths.
5. For each uniquely-assigned leaf branch:
   - Haplogroup = deepest node in the tree reachable via het positions
   - Ratio = median(der_frac) across the branch's het positions
   - QC = fraction of positions whose der_frac is within 2×MAD of the ratio

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass

from yleaf.tree import Tree

LOG = logging.getLogger("yleaf_logger")


@dataclass
class MixtureContributor:
    rank: int
    haplogroup: str
    ratio: float
    markers: int
    qc: float


@dataclass
class MixtureResult:
    tree: str
    common_ancestor: str
    contributors: List[MixtureContributor]


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------

def _get_ancestors(node_name: str, tree: Tree) -> set:
    """All ancestor names of node_name, not including node_name itself."""
    ancestors = set()
    node = tree.get(node_name)
    while node.parent is not None:
        ancestors.add(node.parent.name)
        node = node.parent
    return ancestors


def _get_subtree_nodes(node_name: str, tree: Tree) -> set:
    """All descendant names of node_name (recursive), not including node_name."""
    result = set()
    queue = list(tree.get(node_name).children)
    while queue:
        name = queue.pop()
        if name not in tree.node_mapping:
            continue
        result.add(name)
        queue.extend(tree.get(name).children)
    return result


def _find_lca(node_names: List[str], tree: Tree) -> str:
    """Deepest node that is an ancestor-or-self of every node in node_names."""
    if not node_names:
        raise ValueError("node_names is empty")
    ancestor_sets = [_get_ancestors(n, tree) | {n} for n in node_names]
    common = ancestor_sets[0]
    for s in ancestor_sets[1:]:
        common &= s
    if not common:
        # Fall back to the root's first child
        root_key = 'ROOT (Y-Chromosome "Adam")'
        children = tree.get(root_key).children
        return children[0] if children else root_key
    return max(common, key=lambda n: tree.get(n).depth)


# ---------------------------------------------------------------------------
# Core decomposition
# ---------------------------------------------------------------------------

def _decompose_subtree(
    node_name: str,
    het_df: pd.DataFrame,
    tree: Tree,
    max_n: int,
) -> List[Tuple[str, pd.DataFrame]]:
    """
    Recursively decompose the subtree rooted at node_name into leaf contributors.

    Returns a list of (haplogroup, unique_het_positions) pairs — one per contributor.
    het_df is pre-filtered to positions within this subtree (plus node_name itself).
    """
    if het_df.empty or max_n == 0:
        return []

    node = tree.get(node_name)
    branches: List[Tuple[str, set, pd.DataFrame]] = []
    for child_name in node.children:
        if child_name not in tree.node_mapping:
            continue
        child_nodes = {child_name} | _get_subtree_nodes(child_name, tree)
        child_het = het_df[het_df["haplogroup"].isin(child_nodes)]
        if not child_het.empty:
            branches.append((child_name, child_nodes, child_het))

    if not branches:
        # No children have het positions: current node is the deepest reachable.
        return [(node_name, het_df)]

    if len(branches) == 1:
        # Single branch: keep descending (this node is still on the shared path).
        return _decompose_subtree(branches[0][0], branches[0][2], tree, max_n)

    # Multiple branches: each is a separate contributor (or sub-group of contributors).
    result: List[Tuple[str, pd.DataFrame]] = []
    for child_name, _, child_het in branches:
        if len(result) >= max_n:
            break
        sub = _decompose_subtree(child_name, child_het, tree, max_n - len(result))
        result.extend(sub)
    return result


# ---------------------------------------------------------------------------
# Ratio and QC
# ---------------------------------------------------------------------------

def _der_fracs(het_df: pd.DataFrame) -> np.ndarray:
    """Derived-allele fraction per position."""
    fracs = np.where(
        het_df["state"].values == "D",
        het_df["called_perc"].values / 100.0,
        1.0 - het_df["called_perc"].values / 100.0,
    )
    return fracs.astype(float)


def _estimate_ratio(het_df: pd.DataFrame) -> float:
    return float(np.median(_der_fracs(het_df)))


def _compute_qc(het_df: pd.DataFrame, ratio: float) -> float:
    if het_df.empty:
        return 0.0
    fracs = _der_fracs(het_df)
    deviations = np.abs(fracs - ratio)
    mad = float(np.median(deviations))
    threshold = 2 * mad if mad > 0 else 0.1
    return float(np.sum(deviations <= threshold) / len(fracs))


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def analyze_mixture(
    out_file: Path,
    fmf_file: Path,
    tree_name: str,
    tree_file: Path,
    reads_threshold: int,
    max_contributors: int = 5,
) -> Optional[MixtureResult]:
    """
    Analyse .out + .fmf for a single sample/tree pair and return a MixtureResult,
    or None if no mixture signal is detected.
    """
    # --- Load het positions (below base majority, non-NA state, sufficient reads) ---
    try:
        fmf_df = pd.read_csv(fmf_file, sep="\t")
    except Exception as e:
        LOG.warning(f"Mixture: could not read {fmf_file}: {e}")
        return None

    het_df = fmf_df[
        (fmf_df["Description"] == "Below base majority") &
        (fmf_df["state"].isin(["A", "D"]))
    ].copy()

    if het_df.empty:
        LOG.debug("Mixture: no heterozygous positions found — no mixture signal.")
        return None

    # For BAM, reads = total; for VCF, reads = called (majority) reads.
    # Either way, filter on reads >= threshold (consistent with normal pipeline).
    het_df = het_df[het_df["reads"] >= reads_threshold]
    if het_df.empty:
        LOG.debug("Mixture: all het positions below reads threshold.")
        return None

    # --- Load tree ---
    try:
        tree = Tree(tree_file)
    except Exception as e:
        LOG.warning(f"Mixture: could not load tree {tree_file}: {e}")
        return None

    # Filter to nodes that exist in the tree
    valid = het_df["haplogroup"].isin(tree.node_mapping)
    het_df = het_df[valid]
    if het_df.empty:
        LOG.debug("Mixture: no het positions with recognised tree nodes.")
        return None

    # --- Remove isolated noise nodes ---
    # A node is isolated if no other het node is its ancestor or descendant.
    # Scattered single-position noise spreads across unrelated branches and
    # pulls the LCA to ROOT, drowning real contributor signal.
    het_node_set = set(het_df["haplogroup"].unique())
    connected = set()
    for n in het_node_set:
        relatives = _get_ancestors(n, tree) | _get_subtree_nodes(n, tree)
        if relatives & het_node_set:
            connected.add(n)
    if connected != het_node_set:
        LOG.debug(
            f"Mixture: dropped {len(het_node_set) - len(connected)} isolated noise node(s); "
            f"{len(connected)} connected node(s) remain."
        )
    het_df = het_df[het_df["haplogroup"].isin(connected)]
    if het_df.empty:
        LOG.debug("Mixture: no connected het positions after noise filtering.")
        return None

    het_nodes = list(het_df["haplogroup"].unique())
    LOG.debug(f"Mixture: {len(het_df)} het positions across {len(het_nodes)} nodes after filtering.")

    # --- Find LCA of all het nodes ---
    try:
        lca = _find_lca(het_nodes, tree)
    except Exception as e:
        LOG.warning(f"Mixture: LCA computation failed: {e}")
        return None

    LOG.debug(f"Mixture: LCA = {lca}")

    # --- Decompose subtree into contributor branches ---
    lca_nodes = {lca} | _get_subtree_nodes(lca, tree)
    lca_het = het_df[het_df["haplogroup"].isin(lca_nodes)]
    decomposed = _decompose_subtree(lca, lca_het, tree, max_contributors)

    if not decomposed:
        return None

    if len(decomposed) < 2:
        LOG.debug("Mixture: only one branch found — insufficient evidence for mixture.")
        return None

    # --- Build contributor objects, sorted by descending ratio ---
    raw: List[Tuple[str, float, int, float]] = []
    for hg, unique_het in decomposed:
        ratio = _estimate_ratio(unique_het)
        qc = _compute_qc(unique_het, ratio)
        raw.append((hg, ratio, len(unique_het), qc))

    # Normalise ratios to sum to 1
    total_ratio = sum(r[1] for r in raw)
    if total_ratio > 0:
        raw = [(hg, r / total_ratio, m, q) for hg, r, m, q in raw]
    raw.sort(key=lambda x: x[1], reverse=True)

    contributors = [
        MixtureContributor(rank=i + 1, haplogroup=hg, ratio=round(ratio, 4),
                           markers=markers, qc=round(qc, 4))
        for i, (hg, ratio, markers, qc) in enumerate(raw)
    ]

    LOG.info(f"Mixture: detected {len(contributors)} contributors (tree={tree_name}, LCA={lca})")
    for c in contributors:
        LOG.info(f"  Contributor {c.rank}: {c.haplogroup}  ratio={c.ratio:.2%}  "
                 f"markers={c.markers}  QC={c.qc:.3f}")

    return MixtureResult(tree=tree_name, common_ancestor=lca, contributors=contributors)


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def write_mix_file(mix_file: Path, result: MixtureResult) -> None:
    with open(mix_file, "w") as f:
        f.write("#mixture_version\t1.0\n")
        f.write(f"#tree\t{result.tree}\n")
        f.write(f"#common_ancestor\t{result.common_ancestor}\n")
        f.write(f"#n_contributors\t{len(result.contributors)}\n")
        f.write("contributor\thaplogroup\tratio\tmarkers\tqc\n")
        for c in result.contributors:
            f.write(f"{c.rank}\t{c.haplogroup}\t{c.ratio}\t{c.markers}\t{c.qc}\n")


def read_mix_file(mix_file: Path) -> Optional[MixtureResult]:
    """Parse a .mix file and return a MixtureResult, or None on failure."""
    try:
        meta: Dict[str, str] = {}
        contributors: List[MixtureContributor] = []
        with open(mix_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    parts = line.lstrip("#").split("\t", 1)
                    if len(parts) == 2:
                        meta[parts[0]] = parts[1]
                elif line.startswith("contributor\t"):
                    continue  # header row
                else:
                    parts = line.split("\t")
                    if len(parts) >= 5:
                        contributors.append(MixtureContributor(
                            rank=int(parts[0]),
                            haplogroup=parts[1],
                            ratio=float(parts[2]),
                            markers=int(parts[3]),
                            qc=float(parts[4]),
                        ))
        return MixtureResult(
            tree=meta.get("tree", ""),
            common_ancestor=meta.get("common_ancestor", ""),
            contributors=contributors,
        )
    except Exception as e:
        LOG.warning(f"Could not read mix file {mix_file}: {e}")
        return None
