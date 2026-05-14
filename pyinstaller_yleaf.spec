# PyInstaller spec for standalone Yleaf releases (Linux, macOS, Windows).
# Bundles the yleaf Python package and platform-specific bioinformatics tool
# binaries (samtools, minimap2, bcftools) so no external dependencies are needed.
#
# Build with:
#   python -m PyInstaller pyinstaller_yleaf.spec
#
# On Linux/macOS: samtools/minimap2/bcftools must be on PATH (install via conda bioconda).
# On Windows:     set YLEAF_BUNDLED_TOOLS to a directory containing the .exe + .dll files
#                 (populated by the GitHub Actions CI workflow via MSYS2 MINGW64).
#
# Large full-genome reference files (full_reference.fa) are intentionally excluded;
# users configure those paths via yleaf/config.txt.

import os
import shutil
import glob
import sys
from pathlib import Path

YLEAF_PKG = Path(SPECPATH) / 'yleaf'

# ── Data files ──────────────────────────────────────────────────────────────

datas = []

# Per-build reference data (positions, BED files, chrY.fa, snp_data.csv)
for build in ['hg19', 'hg38', 't2t']:
    build_dir = YLEAF_PKG / 'data' / build
    for f in sorted(build_dir.iterdir()):
        if f.name != 'full_reference.fa' and f.is_file():
            datas.append((str(f), f'yleaf/data/{build}'))

# Haplogroup prediction tables and tree JSON files (all subdirs)
datas.append((
    str(YLEAF_PKG / 'data' / 'hg_prediction_tables'),
    'yleaf/data/hg_prediction_tables',
))

# Config and chromosome naming convention
datas.append((str(YLEAF_PKG / 'config.txt'), 'yleaf'))
datas.append((str(YLEAF_PKG / 'chr_naming_convention.txt'), 'yleaf'))

# JS files for haplogroup tree HTML output
for js in ['svg-pan-zoom.min.js', 'viz.js', 'full.render.js']:
    datas.append((str(YLEAF_PKG / 'data' / js), 'yleaf/data'))

# ── External tool binaries ───────────────────────────────────────────────────

binaries = []
tools_dir = os.environ.get('YLEAF_BUNDLED_TOOLS', '')

if tools_dir and os.path.isdir(tools_dir):
    # Windows CI path: staged directory of .exe + .dll files from MSYS2 MINGW64
    for f in glob.glob(os.path.join(tools_dir, '*.exe')):
        binaries.append((f, '.'))
    for f in glob.glob(os.path.join(tools_dir, '*.dll')):
        binaries.append((f, '.'))
else:
    # Linux / macOS: locate tools on PATH (installed via conda bioconda)
    for tool in ['samtools', 'minimap2', 'bcftools']:
        path = shutil.which(tool)
        if path:
            binaries.append((path, '.'))
        else:
            print(f"WARNING: '{tool}' not found on PATH — it will not be bundled.")

# ── Analysis ─────────────────────────────────────────────────────────────────

a = Analysis(
    [str(YLEAF_PKG / 'Yleaf.py')],
    pathex=[str(YLEAF_PKG.parent)],
    binaries=binaries,
    datas=datas,
    hiddenimports=[
        'yleaf.predict_haplogroup',
        'yleaf.tree',
        'yleaf.haplogroup_tree_image',
        'yleaf.download_reference',
        'yleaf.yleaf_constants',
        'multiprocessing.pool',
        'multiprocessing.managers',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='yleaf',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
    argv_emulation=False,
    target_arch=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='yleaf',
)
