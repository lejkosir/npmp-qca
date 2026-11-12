#!/usr/bin/env bash
# Builds the headless simulator engine/bin/qca_sim (run inside WSL Ubuntu).
# Needs: git, gcc, pkg-config, libglib2.0-dev
#
# The QCADesigner-E source is not stored in this repository. On the first run it is
# cloned from GitHub at a fixed commit into engine/QCADesigner-E/ (ignored by git).
# Our only addition to it is engine/main_qca_sim.c; no upstream file is modified.
set -euo pipefail
# -O0 on purpose: at -O2 the legacy 32-bit-era code misbehaves (design loads with 0 inputs).

UPSTREAM_URL="https://github.com/FSillT/QCADesigner-E.git"
UPSTREAM_COMMIT="3bc90302bddacebfff838c74fb91cc56500387ea"   # 2025-08-17, master

HERE="$(cd "$(dirname "$0")" && pwd)"
QDE="$HERE/QCADesigner-E"
SRC="$QDE/QCADesignerE/src"
OUT="$HERE/bin"

if [ ! -d "$QDE/.git" ]; then
  git clone --quiet "$UPSTREAM_URL" "$QDE"
fi
git -C "$QDE" checkout --quiet "$UPSTREAM_COMMIT"
mkdir -p "$OUT"

SOURCES=(
  fileio.c fileio_helpers.c custom_widgets.c
  simulation.c simulation_data.c bistable_simulation.c
  coherence_vector.c coherence_vector_energy.c
  design.c exp_array.c hsl2rgb.c generic_utils.c graph_dialog_widget_data.c
  objects/object_helpers.c
  objects/QCADCell.c objects/QCADCompoundDO.c objects/QCADDesignObject.c
  objects/QCADDOContainer.c objects/QCADLabel.c objects/QCADLayer.c
  objects/QCADRuler.c objects/QCADStretchyObject.c objects/QCADSubstrate.c
)

cd "$SRC"
gcc -O0 -g -w -I"$SRC" -I"$SRC/objects" -DSTDIO_FILEIO -DPACKAGE=\"QCADesigner\" -DVERSION=\"2.2\" \
  -DGETTEXT_PACKAGE=\"QCADesigner\" -DPACKAGE_DATA_DIR=\"/usr/local/share\" \
  -DPACKAGE_LOCALE_DIR=\"/usr/local/share/locale\" -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 \
  $(pkg-config --cflags glib-2.0 gobject-2.0) \
  "${SOURCES[@]}" "$HERE/main_qca_sim.c" \
  -o "$OUT/qca_sim" \
  $(pkg-config --libs glib-2.0 gobject-2.0 gthread-2.0) -lm

echo "built $OUT/qca_sim"
