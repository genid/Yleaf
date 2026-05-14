#!/usr/bin/env python

"""
Script for drawing haplogroups in the haplogroup tree. Is designed to be used with Yleaf but can also be used with
Yfull haplogroups

Developed at Erasmus Medical Center Department of Genetic Identification

License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

Autor: Bram van Wersch
"""

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
  body { font-family: sans-serif; display: flex; flex-direction: column; height: 100vh;
         background: #1a1a2e; color: #e0e0e0; overflow: hidden; }
  #toolbar { display: flex; align-items: center; gap: 6px; padding: 6px 10px;
             background: #16213e; border-bottom: 1px solid #0f3460; flex-shrink: 0; flex-wrap: wrap; }
  #toolbar h2 { font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.12em;
                color: #6a8ac0; margin-right: 6px; white-space: nowrap; }
  .tab-btn { padding: 4px 10px; border: 1px solid #1a4070; background: #0f2a50;
             color: #80b0e0; border-radius: 4px; cursor: pointer; font-size: 0.75rem; }
  .tab-btn:hover { background: #1a4a8a; }
  .tab-btn.active { background: #1a4a8a; border-color: #4a9eff; color: #fff; }
  .spacer { flex: 1; }
  .btn { padding: 4px 10px; border: 1px solid #1a4070; background: #0f2a50;
         color: #80b0e0; border-radius: 4px; cursor: pointer; font-size: 0.75rem; white-space: nowrap; }
  .btn:hover { background: #1a4a8a; }
  #svg-wrap { flex: 1; overflow: hidden; position: relative; }
  .svg-container { display: none; width: 100%; height: 100%; position: absolute; inset: 0; }
  .svg-container.active { display: block; }
  .svg-container > svg { width: 100%; height: 100%; }
  .svgpanzoom-control-icons { filter: invert(0.8); }
  @media print {
    #toolbar { display: none; }
    body { display: block; height: auto; background: white; }
    #svg-wrap { width: 100vw; height: 100vh; }
  }
</style>
</head>
<body>
<div id="toolbar">
  <h2>Yleaf Haplogroup Tree</h2>
  TAB_BUTTONS_PLACEHOLDER
  <span class="spacer"></span>
  <button class="btn" onclick="fitView()">Fit view</button>
  <button class="btn" onclick="savePdf()">Save PDF</button>
</div>
<div id="svg-wrap">
SVG_CONTAINERS_PLACEHOLDER
</div>
PDF_LINKS_PLACEHOLDER
<script>SVG_PAN_ZOOM_JS_PLACEHOLDER</script>
<script>
const panZooms = {};
let activeTab = null;

function initPanZoom(id) {
  if (panZooms[id]) return;
  const svgEl = document.querySelector('#tab-' + id + ' > svg');
  if (svgEl) {
    panZooms[id] = svgPanZoom(svgEl, {
      zoomEnabled: true, controlIconsEnabled: true, fit: true, center: true,
      minZoom: 0.01, maxZoom: 100,
    });
  }
}

function switchTab(id) {
  document.querySelectorAll('.svg-container').forEach(c => c.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById('tab-' + id).classList.add('active');
  document.querySelector('[data-tab="' + id + '"]').classList.add('active');
  activeTab = id;
  initPanZoom(id);
}

function fitView() {
  if (activeTab && panZooms[activeTab]) { panZooms[activeTab].fit(); panZooms[activeTab].center(); }
}

function savePdf() {
  if (!activeTab) return;
  const link = document.getElementById('pdf-' + activeTab);
  if (link) link.click();
}

window.addEventListener('load', function() {
  const first = document.querySelector('.tab-btn');
  if (first) switchTab(first.dataset.tab);
});
</script>
</body>
</html>
"""


import re as _re

_DATA_DIR = Path(__file__).parent / 'data'


def _dot_esc(s: str) -> str:
    return s.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n')


def _build_dot_string(nodes: dict, edges: list) -> str:
    lines = ['digraph {', '  ratio=compress bgcolor=white']
    for nid, n in nodes.items():
        lines.append(
            f'  "{nid}" [label="{_dot_esc(n["label"])}" shape=box style=filled '
            f'fillcolor="{n["color"]}" fontcolor=black]'
        )
    for from_id, to_id, lbl in edges:
        lattr = f' [label="{_dot_esc(lbl)}"]' if lbl else ''
        lines.append(f'  "{from_id}" -> "{to_id}"{lattr}')
    lines.append('}')
    return '\n'.join(lines)


def _render_dot(dot_str: str):
    """Return (svg_string, pdf_bytes) for a DOT graph string."""
    import graphviz
    src = graphviz.Source(dot_str)
    svg = src.pipe(format='svg').decode('utf-8')
    svg = _re.sub(r'<\?xml[^?]*\?>', '', svg)
    svg = _re.sub(r'<!DOCTYPE[^>]*>', '', svg)
    pdf_bytes = src.pipe(format='pdf')
    return svg.strip(), pdf_bytes


def _write_html(svgs: list, output_path: Path) -> None:
    """svgs: list of (label, svg_string, pdf_bytes) tuples; first entry shown by default."""
    import base64
    tab_buttons = '\n  '.join(
        f'<button class="tab-btn" data-tab="{lbl.lower()}" '
        f'onclick="switchTab(\'{lbl.lower()}\')">{lbl}</button>'
        for lbl, _, _pdf in svgs
    )
    containers = '\n'.join(
        f'<div id="tab-{lbl.lower()}" class="svg-container">{svg}</div>'
        for lbl, svg, _pdf in svgs
    )
    pdf_links = '\n'.join(
        f'<a id="pdf-{lbl.lower()}" href="data:application/pdf;base64,{base64.b64encode(pdf_bytes).decode()}" '
        f'download="yleaf_tree_{lbl}.pdf" style="display:none"></a>'
        for lbl, _svg, pdf_bytes in svgs
    )
    spz = (_DATA_DIR / 'svg-pan-zoom.min.js').read_text(encoding='utf-8')
    html = _HTML_TEMPLATE
    html = html.replace('TAB_BUTTONS_PLACEHOLDER', tab_buttons)
    html = html.replace('SVG_CONTAINERS_PLACEHOLDER', containers)
    html = html.replace('PDF_LINKS_PLACEHOLDER', pdf_links)
    html = html.replace('SVG_PAN_ZOOM_JS_PLACEHOLDER', spz)
    output_path.write_text(html, encoding='utf-8')


def main(namespace: argparse.Namespace = None):
    if namespace is None:
        namespace = get_arguments()
    LOG.info("Starting with drawing haplogroups")
    haplogroups, sample_mapping = read_input_file(namespace.input)
    if len(haplogroups) == 0:
        LOG.warning("No haplogroups found in provided input file.")
        return
    original_haplogroups = list(haplogroups)
    add_main_haplogroups(haplogroups)
    tree_file = getattr(namespace, 'tree_file', None)
    partial_haplogroup_dict = haplogroup_tree_dict(haplogroups, tree_file)
    if namespace.collapse_mode:
        edge_mapping = collapse_tree_dict(partial_haplogroup_dict, sample_mapping)
    else:
        edge_mapping = {}
    make_dendrogram(partial_haplogroup_dict, sample_mapping, edge_mapping, namespace.outfile,
                    original_haplogroups=original_haplogroups, tree_file=tree_file,
                    collapse_mode=namespace.collapse_mode)
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


def _build_nodes_edges(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]],
    edge_mapping: Dict[Tuple[str, str], str],
) -> Tuple[dict, list]:
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
    return nodes, edges


def make_dendrogram(
    partial_tree_dict: Dict[str, Set[str]],
    sample_mapping: Dict[str, List[str]],
    edge_mapping: Dict[Tuple[str, str], str],
    output_file,
    original_haplogroups: List[str] = None,
    tree_file=None,
    collapse_mode: bool = False,
):
    svgs = []

    # Full tree
    LOG.info("Rendering full tree SVG (this may take a while for large datasets)…")
    all_nodes, all_edges = _build_nodes_edges(partial_tree_dict, sample_mapping, edge_mapping)
    svgs.append(('All', *_render_dot(_build_dot_string(all_nodes, all_edges))))

    # Per-major-haplogroup trees
    if original_haplogroups:
        by_major: Dict[str, List[str]] = {}
        for hg in original_haplogroups:
            major = _get_major_hg(hg)
            if major not in ('ROOT', 'OTHER'):
                by_major.setdefault(major, []).append(hg)

        for major in sorted(by_major):
            LOG.info(f"Rendering haplogroup {major} SVG…")
            sub_hgs = list(by_major[major])
            add_main_haplogroups(sub_hgs)
            sub_partial = haplogroup_tree_dict(sub_hgs, tree_file)
            sub_edge_map = collapse_tree_dict(sub_partial, sample_mapping) if collapse_mode else {}
            sub_nodes, sub_edges = _build_nodes_edges(sub_partial, sample_mapping, sub_edge_map)
            svgs.append((major, *_render_dot(_build_dot_string(sub_nodes, sub_edges))))

    html_path = Path(str(output_file) + '.html')
    _write_html(svgs, html_path)
    LOG.info(f"Haplogroup tree written to {html_path}")


if __name__ == '__main__':
    main()
