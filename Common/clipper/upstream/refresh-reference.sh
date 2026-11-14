#!/bin/sh
# Regenerate native reference geometry from a local checkout, never from a moving branch.
set -eu
if [ "$#" -ne 1 ]; then
  echo "Usage: sh $0 /path/to/PolyClipping-Clipper2" >&2
  exit 1
fi
upstream=$(cd "$1" && pwd)
expected=884ba821abc002d895ce72dcc21fe139c207229f
actual=$(git -C "$upstream" rev-parse HEAD)
if [ "$actual" != "$expected" ]; then
  echo "Expected upstream $expected; got $actual. Audit changes before updating the pin." >&2
  exit 1
fi
if ! git -C "$upstream" diff --quiet HEAD -- CPP Tests; then
  echo "Upstream source/test files have uncommitted changes; use the clean pinned revision." >&2
  exit 1
fi
repo=$(cd "$(dirname "$0")/../../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
compiler=${CXX:-c++}
"$compiler" -std=c++17 -O2 -ffp-contract=off -DUSINGZ \
  -I"$upstream/CPP/Clipper2Lib/include" -I"$upstream/CPP/Utils" \
  "$repo/Common/clipper/upstream/GenerateReference.cpp" \
  "$upstream/CPP/Utils/ClipFileLoad.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.engine.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.offset.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.rectclip.cpp" -o "$work/reference"
"$compiler" -std=c++17 -O2 -ffp-contract=off -DUSINGZ \
  -I"$upstream/CPP/Clipper2Lib/include" \
  "$repo/Common/clipper/upstream/GenerateOffsetReference.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.engine.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.offset.cpp" \
  "$upstream/CPP/Clipper2Lib/src/clipper.rectclip.cpp" -o "$work/offsets"
"$work/reference" "$upstream/Tests/Polygons.txt" > "$work/polygons.json"
"$work/reference" "$upstream/Tests/Lines.txt" > "$work/lines.json"
"$work/offsets" > "$work/offsets.json"
python3 - "$work" "$repo/UnitTests/TestData/Clipper2" <<'PY'
import json, pathlib, sys
for name in ('polygons', 'lines', 'offsets'):
    data = json.loads((pathlib.Path(sys.argv[1]) / (name + '.json')).read_text())
    # One case per line keeps fixture diffs readable without duplicating test inputs.
    output = '[\n' + ',\n'.join(json.dumps(case, separators=(',', ':')) for case in data) + '\n]\n'
    (pathlib.Path(sys.argv[2]) / (name + '.json')).write_text(output)
PY
