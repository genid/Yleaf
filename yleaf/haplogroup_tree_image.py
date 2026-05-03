#!/usr/bin/env python

"""
Script for drawing haplogroups in the haplogroup tree. Is designed to be used with Yleaf but can also be used with
Yfull haplogroups

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import re
import json
from pathlib import Path
from typing import List, Dict, Tuple, Set
import graphviz
from collections import defaultdict
import argparse
import logging

from yleaf.tree import Tree
from yleaf import yleaf_constants

LOG = logging.getLogger("yleaf_logger")

# One color per major haplogroup (A–T), plus ROOT and SAMPLE.
# Medium-saturation fills — readable with black graphviz text.
_HG_COLORS: Dict[str, str] = {
    'ROOT':   '#d0d0d0',
    'A':      '#e05050',
    'A0':     '#f07030',
    'B':      '#d4a017',
    'C':      '#7fbf1f',
    'D':      '#3cb44b',
    'E':      '#38b8c8',
    'F':      '#4363d8',
    'G':      '#7b2eb4',
    'H':      '#c032b6',
    'I':      '#909090',
    'J':      '#9b1a1a',
    'K':      '#8b6914',
    'L':      '#6b8000',
    'M':      '#2d8070',
    'N':      '#2a2a9b',
    'O':      '#9b4dca',
    'P':      '#e07840',
    'Q':      '#2da04a',
    'R':      '#8060d0',
    'S':      '#b09020',
    'T':      '#508050',
    'SAMPLE': '#6bcb77',
    'OTHER':  '#e0e0e0',
}


def _get_major_hg(name: str) -> str:
    name = name.strip().split('\n')[0]
    if not name or name.startswith('ROOT'):
        return 'ROOT'
    if name.upper().startswith('A0'):
        return 'A0'
    c = name[0].upper()
    return c if c.isalpha() else 'OTHER'


def _html_decode(s: str) -> str:
    return (s.replace('&amp;', '&').replace('&lt;', '<')
             .replace('&gt;', '>').replace('&quot;', '"')
             .replace('&#10;', '\n').replace('&#9;', '\t'))


def _inject_hg_attrs(svg: str, node_to_hg: Dict[str, str]) -> str:
    """Inject data-hg attributes onto <g class="node"> and <g class="edge"> elements."""

    def _patch_node(m: re.Match) -> str:
        title = _html_decode(m.group(1))
        hg = node_to_hg.get(title, _get_major_hg(title))
        return m.group(0).replace('class="node"', f'class="node" data-hg="{hg}"', 1)

    def _patch_edge(m: re.Match) -> str:
        raw = m.group(1)
        child_raw = raw.split('-&gt;')[-1] if '-&gt;' in raw else raw
        child = _html_decode(child_raw)
        hg = node_to_hg.get(child, _get_major_hg(child))
        return m.group(0).replace('class="edge"', f'class="edge" data-hg="{hg}"', 1)

    svg = re.sub(
        r'(<g\s[^>]*class="node"[^>]*>)\s*<title>([\s\S]*?)</title>',
        lambda m: m.group(0).replace('class="node"', 'class="node" data-hg="{}"'.format(
            node_to_hg.get(_html_decode(m.group(2)), _get_major_hg(_html_decode(m.group(2))))), 1),
        svg,
    )
    svg = re.sub(
        r'(<g\s[^>]*class="edge"[^>]*>)\s*<title>([\s\S]*?)</title>',
        lambda m: m.group(0).replace('class="edge"', 'class="edge" data-hg="{}"'.format(
            node_to_hg.get(_html_decode(m.group(2).split('-&gt;')[-1]),
                           _get_major_hg(_html_decode(m.group(2).split('-&gt;')[-1])))), 1),
        svg,
    )
    return svg


def _write_html(svg: str, output_path: Path, hg_colors: Dict[str, str]) -> None:
    colors_js = json.dumps(hg_colors)
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Yleaf Haplogroup Tree</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: sans-serif; display: flex; height: 100vh;
          background: #1a1a2e; color: #e0e0e0; overflow: hidden; }}
  #sidebar {{ width: 195px; flex-shrink: 0; display: flex; flex-direction: column;
              background: #16213e; border-right: 1px solid #0f3460; }}
  #sidebar-header {{ padding: 10px 14px 8px;
                     border-bottom: 1px solid #0f3460; flex-shrink: 0; }}
  #sidebar-header h2 {{ font-size: 0.78rem; text-transform: uppercase;
                        letter-spacing: 0.12em; color: #6a8ac0; margin-bottom: 8px; }}
  .select-all-row {{ display: flex; align-items: center; gap: 7px;
                     font-size: 0.80rem; cursor: pointer; }}
  .select-all-row input {{ cursor: pointer; accent-color: #4a9eff; flex-shrink: 0; }}
  #sidebar-scroll {{ flex: 1; overflow-y: auto; padding: 10px 14px; }}
  #sidebar-scroll h3 {{ font-size: 0.65rem; text-transform: uppercase;
                        letter-spacing: 0.1em; color: #506070; margin-bottom: 8px; }}
  .filter-item {{ display: flex; align-items: center; gap: 7px;
                  padding: 3px 0; cursor: pointer; font-size: 0.80rem; }}
  .filter-item input {{ cursor: pointer; accent-color: #4a9eff; flex-shrink: 0; }}
  .color-dot {{ width: 10px; height: 10px; border-radius: 50%;
               border: 1px solid rgba(255,255,255,0.25); flex-shrink: 0; }}
  .btn-row {{ display: flex; gap: 6px; padding: 10px 14px 12px;
              border-top: 1px solid #0f3460; flex-shrink: 0; }}
  .btn {{ flex: 1; padding: 5px 4px; border: 1px solid #1a4070;
          background: #0f2a50; color: #80b0e0; border-radius: 5px;
          cursor: pointer; font-size: 0.72rem; }}
  .btn:hover {{ background: #1a4a8a; }}
  #svg-wrap {{ flex: 1; overflow: hidden; }}
  #svg-wrap svg {{ width: 100%; height: 100%; }}
  .svgpanzoom-control-icons {{ filter: invert(0.8); }}
  @media print {{
    #sidebar {{ display: none; }}
    body {{ display: block; height: auto; background: white; }}
    #svg-wrap {{ width: 100vw; height: 100vh; }}
    #svg-wrap svg {{ width: 100%; height: 100%; }}
    .svgpanzoom-control-icons {{ display: none; }}
  }}
</style>
</head>
<body>
<div id="sidebar">
  <div id="sidebar-header">
    <h2>Yleaf Haplogroup Tree</h2>
    <label class="select-all-row">
      <input type="checkbox" id="select-all-cb" checked onchange="toggleAll(this.checked)">
      <span>Select all</span>
    </label>
  </div>
  <div id="sidebar-scroll">
    <h3>Filter by major haplogroup</h3>
    <div id="filters"></div>
  </div>
  <div class="btn-row">
    <button class="btn" onclick="fitVisible()">Fit view</button>
    <button class="btn" onclick="savePdf()">Save PDF</button>
  </div>
</div>
<div id="svg-wrap">
{svg}
</div>
<script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom@3.6.1/dist/svg-pan-zoom.min.js"></script>
<script>
const COLOR_MAP = {colors_js};
let panZoom = null;

window.addEventListener('load', function () {{
  const svgEl = document.querySelector('#svg-wrap svg');
  if (svgEl) {{
    svgEl.setAttribute('width', '100%');
    svgEl.setAttribute('height', '100%');
    panZoom = svgPanZoom(svgEl, {{
      zoomEnabled: true,
      controlIconsEnabled: true,
      fit: true,
      center: true,
      minZoom: 0.05,
      maxZoom: 100,
    }});
  }}

  // Build filter sidebar
  const nodeEls = document.querySelectorAll('g.node[data-hg]');
  const groupsPresent = [...new Set([...nodeEls].map(n => n.dataset.hg))].sort();
  const filtersDiv = document.getElementById('filters');
  for (const hg of groupsPresent) {{
    const color = COLOR_MAP[hg] || '#888';
    const lbl = document.createElement('label');
    lbl.className = 'filter-item';
    lbl.innerHTML =
      '<input type="checkbox" checked onchange="toggleHg(\\''+hg+'\\', this.checked)">' +
      '<span class="color-dot" style="background:' + color + '"></span>' +
      '<span>' + hg + '</span>';
    filtersDiv.appendChild(lbl);
  }}
}});

function toggleHg(hg, visible) {{
  document.querySelectorAll('g.node[data-hg="' + hg + '"], g.edge[data-hg="' + hg + '"]')
    .forEach(el => el.style.display = visible ? '' : 'none');
  updateSelectAllState();
  fitVisible();
}}

function toggleAll(checked) {{
  document.querySelectorAll('g.node[data-hg], g.edge[data-hg]')
    .forEach(el => el.style.display = checked ? '' : 'none');
  document.querySelectorAll('#filters input').forEach(cb => cb.checked = checked);
  fitVisible();
}}

function updateSelectAllState() {{
  const all = document.querySelectorAll('#filters input');
  const checked = [...all].filter(cb => cb.checked).length;
  const cb = document.getElementById('select-all-cb');
  cb.indeterminate = checked > 0 && checked < all.length;
  cb.checked = checked === all.length;
}}

function fitVisible() {{
  if (!panZoom) return;
  const svgEl = document.querySelector('#svg-wrap svg');
  const visible = [...document.querySelectorAll('g.node[data-hg]')]
    .filter(el => el.style.display !== 'none');
  if (visible.length === 0) return;

  // Compute union bounding box in SVG coordinates
  const ctm = svgEl.getScreenCTM();
  const inv = ctm.inverse();
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (const el of visible) {{
    const r = el.getBoundingClientRect();
    // Convert screen corners to SVG coordinate space
    const tl = new DOMPoint(r.left, r.top).matrixTransform(inv);
    const br = new DOMPoint(r.right, r.bottom).matrixTransform(inv);
    minX = Math.min(minX, tl.x, br.x);
    minY = Math.min(minY, tl.y, br.y);
    maxX = Math.max(maxX, tl.x, br.x);
    maxY = Math.max(maxY, tl.y, br.y);
  }}

  const pad = 30;
  const bw = maxX - minX + pad * 2;
  const bh = maxY - minY + pad * 2;
  const wrap = document.getElementById('svg-wrap');
  const scaleX = wrap.clientWidth / bw;
  const scaleY = wrap.clientHeight / bh;
  const newZoom = Math.min(scaleX, scaleY, panZoom.getMaxZoom());

  panZoom.zoom(newZoom);
  // Center on the midpoint of the bounding box
  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;
  panZoom.pan({{
    x: wrap.clientWidth / 2 - cx * newZoom,
    y: wrap.clientHeight / 2 - cy * newZoom,
  }});
}}

function savePdf() {{
  window.print();
}}
</script>
</body>
</html>"""
    output_path.write_text(html, encoding='utf-8')


def main(namespace: argparse.Namespace = None):
    if namespace is None:
        namespace = get_arguments()
    LOG.info("Starting with drawing haplogroups")
    haplogroups, sample_mapping = read_input_file(namespace.input)
    if len(haplogroups) == 0:
        LOG.warning("No haplogroups found in provided input file.")
        return
    add_main_haplogroups(haplogroups)
    tree_file = getattr(namespace, 'tree_file', None)
    partial_haplogroup_dict = haplogroup_tree_dict(haplogroups, tree_file)
    if namespace.collapse_mode:
        edge_mapping = collapse_tree_dict(partial_haplogroup_dict, sample_mapping)
    else:
        edge_mapping = {}
    make_dendrogram(partial_haplogroup_dict, sample_mapping, edge_mapping, namespace.outfile)
    LOG.info("Finished drawing haplogroups")


def get_arguments() -> argparse.Namespace:
    """Get the arguments provided by the user to this script"""
    parser = argparse.ArgumentParser(description="Erasmus MC: Genetic Identification\n Haplogroup tree creation")

    parser.add_argument("-i", "--input", required=True,
                        help="prediction file (.hg file) containing the predicted haplogroups for different samples. "
                             "Alternalively you can provide a file containing samples in the first column and "
                             "haplogroups in the second.",
                        metavar="FILE")
    parser.add_argument("-c", "--collapse_mode", help="Add this flag to compress the haplogroup tree image and remove"
                                                      " all uninformative haplogroups from it.",
                        action="store_true")
    parser.add_argument("-o", "--outfile", required=True, help="Output file name", metavar="FILE")

    args = parser.parse_args()
    return args


def read_input_file(
    file: str
) -> Tuple[List[str], Dict[str, List[str]]]:
    haplogroups = []
    sample_mapping = defaultdict(list)
    with open(file) as f:
        f.readline()
        for line in f:
            values = line.split("\t")
            sample, haplogroup = values[:2]
            if haplogroup == 'NA':
                continue
            haplogroup = haplogroup.split("*")[0]
            haplogroups.append(haplogroup)
            sample_mapping[haplogroup].append(sample)
    return haplogroups, dict(sample_mapping)


def add_main_haplogroups(haplogroups: List[str]):
    haplogroups.extend([chr(nr) for nr in range(66, 85)])
    haplogroups.append("A0-T")
    haplogroups.append("A00")


def haplogroup_tree_dict(
    haplogroups: List[str],
    tree_file=None
) -> Dict[str, Set[str]]:
    if tree_file is None:
        tree_file = yleaf_constants.DATA_FOLDER / yleaf_constants.HG_PREDICTION_FOLDER / yleaf_constants.TREE_FILE
    tree = Tree(tree_file)
    partial_tree_dict = defaultdict(set)
    for name in haplogroups:
        path = []
        try:
            node = tree.get(name)
        except KeyError:
            continue
        while node is not None:
            path.append(node.name)
            node = node.parent
        parent = path.pop()
        while len(path) > 0:
            child = path.pop()
            partial_tree_dict[parent].add(child)
            parent = child
    return dict(partial_tree_dict)


def collapse_tree_dict(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]]
) -> Dict[Tuple[str, str], str]:
    edge_mapping = {}
    remove_keys = set()
    for parent, children in partial_tree_dict.items():
        if parent in remove_keys:
            continue
        orig_parent = parent
        for child in children.copy():
            total_collapsed = 0
            partial_tree_dict[orig_parent].remove(child)
            while True:
                if not can_collapse(child, partial_tree_dict, sample_mapping):
                    break
                remove_keys.add(child)
                child = next(iter(partial_tree_dict[child]))
                total_collapsed += 1
            partial_tree_dict[orig_parent].add(child)
            if total_collapsed > 0:
                edge_mapping[(orig_parent, child)] = str(total_collapsed)
    for name in remove_keys:
        del partial_tree_dict[name]
    return edge_mapping


def can_collapse(
    name: str,
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]]
) -> bool:
    if name not in partial_tree_dict:
        return False
    if len(partial_tree_dict[name]) != 1:
        return False
    if name == 'A0-T' or '-' not in name:
        return False
    if name in sample_mapping:
        return False
    return True


def make_dendrogram(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]],
    edge_mapping: Dict[Tuple[str, str], str],
    output_file,
):
    dot = graphviz.Digraph()
    dot.attr(ratio="compress", bgcolor="white")
    covered_nodes = set()
    covered_edges = set()

    # Map every node label (incl. sample nodes) → major haplogroup for JS filtering
    node_to_hg: Dict[str, str] = {}

    for parent, children in partial_tree_dict.items():
        parent_hg = _get_major_hg(parent)
        parent_color = _HG_COLORS.get(parent_hg, _HG_COLORS['OTHER'])
        if parent not in covered_nodes:
            dot.node(parent, parent, shape="box", style="filled",
                     fillcolor=parent_color, fontcolor="black")
            node_to_hg[parent] = parent_hg
            covered_nodes.add(parent)

        for child in children:
            child_hg = _get_major_hg(child)
            child_color = _HG_COLORS.get(child_hg, _HG_COLORS['OTHER'])
            if child not in covered_nodes:
                dot.node(child, child, shape="box", style="filled",
                         fillcolor=child_color, fontcolor="black")
                node_to_hg[child] = child_hg
                covered_nodes.add(child)

            edge_name = (parent, child)
            if edge_name not in covered_edges:
                label = edge_mapping.get(edge_name, '')
                dot.edge(parent, child, weight="1", label=label)
                covered_edges.add(edge_name)

                if child in sample_mapping:
                    sample_label = '\n'.join(sample_mapping[child])
                    dot.node(sample_label, sample_label, shape="box",
                             style="filled", fillcolor=_HG_COLORS['SAMPLE'],
                             fontcolor="black")
                    node_to_hg[sample_label] = child_hg
                    dot.edge(child, sample_label, weight="10")

    # Generate SVG in memory (no temp files)
    svg_bytes = dot.pipe(format='svg')
    svg = svg_bytes.decode('utf-8')

    # Strip the XML declaration so the SVG embeds cleanly in HTML
    svg = re.sub(r'<\?xml[^?]*\?>', '', svg).strip()

    # Inject data-hg attributes for JS filtering
    svg = _inject_hg_attrs(svg, node_to_hg)

    # Write self-contained HTML
    html_path = Path(str(output_file) + '.html')
    _write_html(svg, html_path, _HG_COLORS)
    LOG.info(f"Haplogroup tree written to {html_path}")


if __name__ == '__main__':
    main()
