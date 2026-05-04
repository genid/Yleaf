#!/usr/bin/env python

"""
Script for drawing haplogroups in the haplogroup tree. Is designed to be used with Yleaf but can also be used with
Yfull haplogroups

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

import json
from pathlib import Path
from typing import List, Dict, Tuple, Set
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


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Yleaf Haplogroup Tree</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: sans-serif; display: flex; height: 100vh;
         background: #1a1a2e; color: #e0e0e0; overflow: hidden; }
  #sidebar { width: 195px; flex-shrink: 0; display: flex; flex-direction: column;
             background: #16213e; border-right: 1px solid #0f3460; }
  #sidebar-header { padding: 10px 14px 8px; border-bottom: 1px solid #0f3460; flex-shrink: 0; }
  #sidebar-header h2 { font-size: 0.78rem; text-transform: uppercase;
                       letter-spacing: 0.12em; color: #6a8ac0; margin-bottom: 8px; }
  .select-all-row { display: flex; align-items: center; gap: 7px;
                    font-size: 0.80rem; cursor: pointer; }
  .select-all-row input { cursor: pointer; accent-color: #4a9eff; flex-shrink: 0; }
  #sidebar-scroll { flex: 1; overflow-y: auto; padding: 10px 14px; }
  #sidebar-scroll h3 { font-size: 0.65rem; text-transform: uppercase;
                       letter-spacing: 0.1em; color: #506070; margin-bottom: 8px; }
  .filter-item { display: flex; align-items: center; gap: 7px;
                 padding: 3px 0; cursor: pointer; font-size: 0.80rem; }
  .filter-item input { cursor: pointer; accent-color: #4a9eff; flex-shrink: 0; }
  .color-dot { width: 10px; height: 10px; border-radius: 50%;
               border: 1px solid rgba(255,255,255,0.25); flex-shrink: 0; }
  .btn-row { display: flex; gap: 6px; padding: 10px 14px 12px;
             border-top: 1px solid #0f3460; flex-shrink: 0; }
  .btn { flex: 1; padding: 5px 4px; border: 1px solid #1a4070;
         background: #0f2a50; color: #80b0e0; border-radius: 5px;
         cursor: pointer; font-size: 0.72rem; }
  .btn:hover { background: #1a4a8a; }
  .btn:disabled { opacity: 0.5; cursor: default; }
  #svg-wrap { flex: 1; overflow: hidden; position: relative; }
  #svg-wrap svg { width: 100%; height: 100%; }
  #loading { position: absolute; inset: 0; display: flex; align-items: center;
             justify-content: center; font-size: 0.9rem; color: #6a8ac0;
             background: #1a1a2e; }
  .svgpanzoom-control-icons { filter: invert(0.8); }
  @media print {
    #sidebar { display: none; }
    body { display: block; height: auto; background: white; }
    #svg-wrap { width: 100vw; height: 100vh; }
    #svg-wrap svg { width: 100%; height: 100%; }
    .svgpanzoom-control-icons { display: none; }
  }
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
    <button class="btn" id="fit-btn" onclick="fitView()">Fit view</button>
    <button class="btn" id="pdf-btn" onclick="savePdf()">Save PDF</button>
  </div>
</div>
<div id="svg-wrap">
  <div id="loading">Rendering tree…</div>
</div>

<script src="https://cdn.jsdelivr.net/npm/viz.js@2.1.2/viz.js"></script>
<script src="https://cdn.jsdelivr.net/npm/viz.js@2.1.2/full.render.js"></script>
<script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom@3.6.1/dist/svg-pan-zoom.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js"></script>
<script>
const TREE = TREE_DATA_PLACEHOLDER;

// Build adjacency maps
const childrenOf = {};
const parentOf   = {};
for (const [from, to, lbl] of TREE.edges) {
  if (!childrenOf[from]) childrenOf[from] = [];
  childrenOf[from].push([to, lbl]);
  parentOf[to] = from;
}

// Bottom-up DFS: a node is visible if its hg is checked OR any descendant is visible
function computeVisible(checkedHgs) {
  const visited = new Set();
  const visible = new Set();
  function visit(id) {
    if (visited.has(id)) return visible.has(id);
    visited.add(id);
    const node = TREE.nodes[id];
    const kids = childrenOf[id] || [];
    // Visit ALL children — never short-circuit
    let anyChildVis = false;
    for (const [c] of kids) {
      if (visit(c)) anyChildVis = true;
    }
    if (checkedHgs.has(node.hg) || anyChildVis) {
      visible.add(id);
      return true;
    }
    return false;
  }
  const roots = Object.keys(TREE.nodes).filter(id => !(id in parentOf));
  for (const r of roots) visit(r);
  return visible;
}

function dotEsc(s) {
  return s.replace(/\\\\/g, '\\\\\\\\')
          .replace(/"/g, '\\\\"')
          .replace(/\\n/g, '\\\\n');
}

function buildDot(visible) {
  let dot = 'digraph {\\n  ratio=compress bgcolor=white\\n';
  for (const id of visible) {
    const n = TREE.nodes[id];
    const lbl = dotEsc(n.label);
    dot += `  "${id}" [label="${lbl}" shape=box style=filled fillcolor="${n.color}" fontcolor=black]\\n`;
  }
  for (const [from, to, lbl] of TREE.edges) {
    if (!visible.has(from) || !visible.has(to)) continue;
    const lattr = lbl ? ` [label="${dotEsc(lbl)}"]` : '';
    dot += `  "${from}" -> "${to}"${lattr}\\n`;
  }
  dot += '}';
  return dot;
}

let viz = null;
let panZoom = null;
let rendering = false;
let pendingChecked = null;

async function render(checkedHgs) {
  if (rendering) { pendingChecked = checkedHgs; return; }
  rendering = true;
  document.getElementById('loading').style.display = 'flex';
  document.getElementById('fit-btn').disabled = true;

  const visible = computeVisible(checkedHgs);
  const dot = buildDot(visible);

  try {
    const svgStr = await viz.renderString(dot, { format: 'svg' });
    const wrap = document.getElementById('svg-wrap');
    if (panZoom) { panZoom.destroy(); panZoom = null; }
    wrap.innerHTML = '<div id="loading" style="display:none"></div>' + svgStr;
    const svgEl = wrap.querySelector('svg');
    svgEl.setAttribute('width', '100%');
    svgEl.setAttribute('height', '100%');
    panZoom = svgPanZoom(svgEl, {
      zoomEnabled: true, controlIconsEnabled: true, fit: true, center: true,
      minZoom: 0.01, maxZoom: 100,
    });
  } catch(e) {
    console.error('viz.js render error:', e);
  }

  document.getElementById('loading').style.display = 'none';
  document.getElementById('fit-btn').disabled = false;
  rendering = false;

  if (pendingChecked) {
    const next = pendingChecked;
    pendingChecked = null;
    render(next);
  }
}

function buildSidebar() {
  const filtersDiv = document.getElementById('filters');
  const hgGroups = Object.keys(TREE.hg_colors)
    .filter(h => h !== 'ROOT' && h !== 'SAMPLE' && h !== 'OTHER').sort();
  for (const hg of hgGroups) {
    const color = TREE.hg_colors[hg];
    const lbl = document.createElement('label');
    lbl.className = 'filter-item';
    lbl.innerHTML =
      '<input type="checkbox" checked onchange="onFilterChange()">' +
      '<span class="color-dot" style="background:' + color + '"></span>' +
      '<span>' + hg + '</span>';
    filtersDiv.appendChild(lbl);
  }
}

function getCheckedHgs() {
  const checked = new Set(['ROOT', 'SAMPLE', 'OTHER']);
  document.querySelectorAll('#filters input').forEach(cb => {
    const label = cb.parentElement.querySelector('span:last-child').textContent;
    if (cb.checked) checked.add(label);
  });
  return checked;
}

function onFilterChange() {
  updateSelectAllState();
  render(getCheckedHgs());
}

function toggleAll(checked) {
  document.querySelectorAll('#filters input').forEach(cb => cb.checked = checked);
  render(getCheckedHgs());
}

function updateSelectAllState() {
  const all = document.querySelectorAll('#filters input');
  const n = [...all].filter(cb => cb.checked).length;
  const cb = document.getElementById('select-all-cb');
  cb.indeterminate = n > 0 && n < all.length;
  cb.checked = n === all.length;
}

function fitView() {
  if (panZoom) { panZoom.fit(); panZoom.center(); }
}

function savePdf() {
  const svgEl = document.querySelector('#svg-wrap svg');
  if (!svgEl) return;
  const vb = svgEl.viewBox.baseVal;
  const w = vb.width || svgEl.clientWidth || 1200;
  const h = vb.height || svgEl.clientHeight || 800;

  const clone = svgEl.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width', w);
  clone.setAttribute('height', h);

  const svgStr = new XMLSerializer().serializeToString(clone);
  const svgData = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgStr)));

  const img = new Image();
  img.onload = function() {
    const scale = 2;
    const canvas = document.createElement('canvas');
    canvas.width = w * scale; canvas.height = h * scale;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#1a1a2e';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0, w, h);
    const { jsPDF } = window.jspdf;
    const doc = new jsPDF({ orientation: w > h ? 'l' : 'p', unit: 'px', format: [w, h] });
    doc.addImage(canvas.toDataURL('image/jpeg', 0.92), 'JPEG', 0, 0, w, h);
    doc.save('yleaf_haplogroup_tree.pdf');
  };
  img.onerror = () => alert('PDF export failed. Use File → Print → Save as PDF instead.');
  img.src = svgData;
}

window.addEventListener('load', async function() {
  viz = new Viz();
  buildSidebar();
  await render(getCheckedHgs());
});
</script>
</body>
</html>
"""


def _write_html(tree_data: dict, output_path: Path) -> None:
    tree_json = json.dumps(tree_data, ensure_ascii=False)
    html = _HTML_TEMPLATE.replace('TREE_DATA_PLACEHOLDER', tree_json)
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
    node_ids: Dict[str, str] = {}
    nodes: Dict[str, dict] = {}
    edges = []
    covered_nodes: Set[str] = set()
    covered_edges: Set[Tuple[str, str]] = set()

    def get_id(name: str) -> str:
        if name not in node_ids:
            node_ids[name] = str(len(node_ids))
        return node_ids[name]

    for parent, children in partial_tree_dict.items():
        parent_hg = _get_major_hg(parent)
        if parent not in covered_nodes:
            nodes[get_id(parent)] = {"label": parent, "hg": parent_hg,
                                     "color": _HG_COLORS.get(parent_hg, _HG_COLORS['OTHER']),
                                     "is_sample": False}
            covered_nodes.add(parent)

        for child in children:
            child_hg = _get_major_hg(child)
            if child not in covered_nodes:
                nodes[get_id(child)] = {"label": child, "hg": child_hg,
                                        "color": _HG_COLORS.get(child_hg, _HG_COLORS['OTHER']),
                                        "is_sample": False}
                covered_nodes.add(child)

            if (parent, child) not in covered_edges:
                label = edge_mapping.get((parent, child), '')
                edges.append([get_id(parent), get_id(child), label])
                covered_edges.add((parent, child))

                if child in sample_mapping:
                    sample_label = '\n'.join(sorted(sample_mapping[child]))
                    if sample_label not in covered_nodes:
                        nodes[get_id(sample_label)] = {"label": sample_label, "hg": child_hg,
                                                       "color": _HG_COLORS['SAMPLE'],
                                                       "is_sample": True}
                        covered_nodes.add(sample_label)
                    edges.append([get_id(child), get_id(sample_label), ""])

    hg_colors = {hg: _HG_COLORS.get(hg, _HG_COLORS['OTHER'])
                 for hg in sorted(set(n["hg"] for n in nodes.values()))}
    tree_data = {"nodes": nodes, "edges": edges, "hg_colors": hg_colors}

    html_path = Path(str(output_file) + '.html')
    _write_html(tree_data, html_path)
    LOG.info(f"Haplogroup tree written to {html_path}")


if __name__ == '__main__':
    main()
