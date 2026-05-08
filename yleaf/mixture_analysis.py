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
    branch_root: Optional[str] = None,
) -> List[Tuple[str, pd.DataFrame, str]]:
    """
    Recursively decompose the subtree rooted at node_name into leaf contributors.

    Returns a list of (haplogroup, unique_het_positions, branch_root) triples.
    branch_root is the first child of the split node on this branch — used by the
    caller to limit the path-coherence check to below the split.
    het_df is pre-filtered to positions within this subtree (plus node_name itself).
    """
    if het_df.empty:
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
        return [(node_name, het_df, branch_root or node_name)]

    if len(branches) == 1:
        # Single branch: keep descending (this node is still on the shared path).
        return _decompose_subtree(branches[0][0], branches[0][2], tree, branch_root)

    # Multiple branches: split found here. Process largest branches first so that
    # real contributors (many het positions) fill slots before single-position noise.
    branches.sort(key=lambda x: len(x[2]), reverse=True)
    result: List[Tuple[str, pd.DataFrame, str]] = []
    for child_name, _, child_het in branches:
        sub = _decompose_subtree(child_name, child_het, tree,
                                 branch_root=child_name)
        result.extend(sub)
    return result


def _path_has_called_intermediate(
    candidate: str,
    branch_root: str,
    tree: Tree,
    called_nodes: set,
) -> bool:
    """
    Walk ancestors from candidate up to (but not including) branch_root.
    Return True if any intermediate node appears in called_nodes (has a non-het
    majority call), which contradicts the candidate being a real contributor.
    """
    node = tree.get(candidate)
    while node.parent is not None and node.name != branch_root:
        if node.name != candidate and node.name in called_nodes:
            return True
        node = node.parent
    # Also check branch_root itself
    if branch_root in called_nodes:
        return True
    return False


def _het_fraction(node_name: str, het_per_node: "pd.Series", out_per_node: "pd.Series") -> float:
    """Fraction of observed positions for node_name that are heterozygous."""
    n_het = het_per_node.get(node_name, 0)
    n_tot = n_het + out_per_node.get(node_name, 0)
    return n_het / n_tot if n_tot > 0 else 0.0


def _path_coherent_to_lca(
    node_name: str,
    lca: str,
    tree: Tree,
    het_per_node: "pd.Series",
    out_per_node: "pd.Series",
) -> bool:
    """
    Return True if every observed intermediate node between node_name and lca has
    het fraction > 0.5.  Nodes with zero total observations are unobserved and
    treated as neutral (they do not block the path).
    """
    node = tree.get(node_name)
    if node.parent is None or node.name == lca:
        return True
    node = node.parent
    while node is not None and node.name != lca:
        n_het = het_per_node.get(node.name, 0)
        n_tot = n_het + out_per_node.get(node.name, 0)
        if n_tot > 0 and n_het / n_tot <= 0.5:
            return False
        if node.parent is None:
            break
        node = node.parent
    return True


def _find_deepest_coherent_ancestor(
    candidate: str,
    lca: str,
    tree: Tree,
    het_per_node: "pd.Series",
    out_per_node: "pd.Series",
) -> Optional[str]:
    """
    Starting from candidate, climb toward lca.
    Return the deepest node (inclusive of candidate) whose path to lca is coherent
    (every intermediate has het fraction > 0.5).  Returns None if no such node exists.
    """
    node = tree.get(candidate)
    while node is not None and node.name != lca:
        if _path_coherent_to_lca(node.name, lca, tree, het_per_node, out_per_node):
            return node.name
        if node.parent is None:
            break
        node = node.parent
    return None


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
) -> Optional[MixtureResult]:
    """
    Analyse .out + .fmf for a single sample/tree pair and return a MixtureResult,
    or None if no mixture signal is detected.
    """
    # --- Load called (non-het majority) positions from .out file ---
    out_df = pd.DataFrame()
    called_nodes: set = set()
    try:
        out_df = pd.read_csv(out_file, sep="\t")
        called_nodes = set(out_df["haplogroup"].dropna().unique())
    except Exception as e:
        LOG.warning(f"Mixture: could not read {out_file}: {e}")

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

    # --- Filter het nodes for LCA computation: majority of observed positions must be het ---
    # A node qualifies only if more than half its total observed positions (called + het)
    # are heterozygous.  Universal ancestor nodes (A0-T, A1, …) have mostly clear derived
    # calls and will be excluded, preventing them from pulling the LCA above the real split.
    out_per_node = out_df["haplogroup"].value_counts() if not out_df.empty else pd.Series(dtype=int)
    het_per_node = het_df["haplogroup"].value_counts()
    lca_het_nodes = [
        n for n in het_nodes
        if (het_per_node.get(n, 0) + out_per_node.get(n, 0)) > 0
        and het_per_node.get(n, 0) / (het_per_node.get(n, 0) + out_per_node.get(n, 0)) > 0.5
    ]
    if not lca_het_nodes:
        lca_het_nodes = het_nodes  # fallback if filter removes everything

    dropped = set(het_nodes) - set(lca_het_nodes)
    if dropped:
        LOG.debug(f"Mixture: excluded {len(dropped)} node(s) from LCA computation (het fraction ≤ 0.5): {dropped}")

    # --- Find LCA from het-majority nodes only ---
    try:
        lca = _find_lca(lca_het_nodes, tree)
    except Exception as e:
        LOG.warning(f"Mixture: LCA computation failed: {e}")
        return None

    # --- Validate LCA: must have ≥2 direct children with het-majority positions ---
    # If the computed LCA has fewer than 2 qualifying children, descend into the
    # single qualifying child until we reach a proper split point.
    # If 0 children qualify directly, the LCA is too high (pulled up by noise);
    # descend into the child whose subtree contains the most het nodes and try again.
    lca_het_nodes_set = set(lca_het_nodes)
    while True:
        qualifying_children = [
            c for c in tree.get(lca).children
            if c in tree.node_mapping and _het_fraction(c, het_per_node, out_per_node) > 0.5
        ]
        if len(qualifying_children) >= 2:
            break
        if len(qualifying_children) == 1:
            lca = qualifying_children[0]
        else:
            best_child = max(
                (c for c in tree.get(lca).children if c in tree.node_mapping),
                key=lambda c: len((_get_subtree_nodes(c, tree) | {c}) & lca_het_nodes_set),
                default=None,
            )
            if best_child is None or not ((_get_subtree_nodes(best_child, tree) | {best_child}) & lca_het_nodes_set):
                LOG.debug("Mixture: no het signal below LCA — no mixture detected.")
                return None
            lca = best_child

    LOG.debug(f"Mixture: LCA = {lca}")

    # --- Decompose subtree into contributor branches ---
    lca_nodes = {lca} | _get_subtree_nodes(lca, tree)
    lca_het = het_df[het_df["haplogroup"].isin(lca_nodes)]
    decomposed = _decompose_subtree(lca, lca_het, tree)

    if not decomposed:
        return None

    # --- Keep only deepest node per lineage (drop ancestors of other results) ---
    hg_names = {hg for hg, _, _ in decomposed}
    decomposed = [
        (hg, het, br) for hg, het, br in decomposed
        if not (_get_subtree_nodes(hg, tree) & hg_names)
    ]

    # --- Resolve each branch to its deepest coherent ancestor ---
    # A branch is coherent when every node on the path from the candidate back to
    # the LCA has het fraction > 0.5 (majority of observed positions are heterozygous).
    # Incoherent branches are climbed toward the LCA until a coherent ancestor is found.
    # Multiple branches that collapse to the same ancestor are merged automatically.
    resolved: Dict[str, pd.DataFrame] = {}
    for hg, _unique_het, _branch_root in decomposed:
        if _path_coherent_to_lca(hg, lca, tree, het_per_node, out_per_node):
            coherent_hg = hg
        else:
            coherent_hg = _find_deepest_coherent_ancestor(hg, lca, tree, het_per_node, out_per_node)
            if coherent_hg is None:
                LOG.debug(f"Mixture: discarding {hg} — no coherent ancestor found below LCA.")
                continue
            if coherent_hg != hg:
                LOG.debug(f"Mixture: {hg} → promoted to coherent ancestor {coherent_hg}")

        if coherent_hg in called_nodes:
            LOG.debug(f"Mixture: discarding {coherent_hg} — has called positions (incoherent with mixture)")
            continue

        if coherent_hg in resolved:
            continue  # het positions for this ancestor already collected

        coherent_nodes = {coherent_hg} | _get_subtree_nodes(coherent_hg, tree)
        branch_het = lca_het[lca_het["haplogroup"].isin(coherent_nodes)]
        if not branch_het.empty:
            resolved[coherent_hg] = branch_het

    decomposed_pairs = list(resolved.items())

    if len(decomposed_pairs) < 2:
        LOG.debug("Mixture: fewer than 2 coherent branches — insufficient evidence for mixture.")
        return None

    # --- Build contributor objects, sorted by descending ratio ---
    raw: List[Tuple[str, float, int, float]] = []
    for hg, unique_het in decomposed_pairs:
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
