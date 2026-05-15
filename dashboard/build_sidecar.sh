#!/usr/bin/env bash
# Builds the Yleaf PyInstaller sidecar and places it in src-tauri/binaries/
# with the Tauri sidecar naming convention: yleaf-<rust-target-triple>
#
# Usage: ./build_sidecar.sh
# Run from anywhere; script resolves its own directory.

set -euo pipefail

YLEAF_PYTHON=/home/arwin2/miniconda3/envs/yleaf/bin/python
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINARIES_DIR="$SCRIPT_DIR/src-tauri/binaries"
DIST_DIR=/tmp/yleaf-sidecar-dist
BUILD_DIR=/tmp/yleaf-sidecar-build

TARGET=$(rustc -vV | sed -n 's/host: //p')
echo "Target triple: $TARGET"

if ! "$YLEAF_PYTHON" -c "import PyInstaller" 2>/dev/null; then
    echo "Installing PyInstaller into yleaf conda env..."
    "$YLEAF_PYTHON" -m pip install pyinstaller
fi

echo "Running PyInstaller..."
"$YLEAF_PYTHON" -m PyInstaller \
    --distpath "$DIST_DIR" \
    --workpath "$BUILD_DIR" \
    --noconfirm \
    "$SCRIPT_DIR/pyinstaller_yleaf.spec"

mkdir -p "$BINARIES_DIR"
cp "$DIST_DIR/yleaf" "$BINARIES_DIR/yleaf-$TARGET"
chmod +x "$BINARIES_DIR/yleaf-$TARGET"

echo "Sidecar ready: $BINARIES_DIR/yleaf-$TARGET"
