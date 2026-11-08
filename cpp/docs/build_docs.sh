#!/bin/bash
# Build the C++ user guide: compile and run the example programs against the installed library,
# then run Doxygen over the pages and the headers. Logs stay next to this script.
#
#   bash cpp/docs/build_docs.sh                     # library installed under cpp/install
#   TREX_INSTALL_PREFIX=/opt/trex bash cpp/docs/build_docs.sh
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
PREFIX=${TREX_INSTALL_PREFIX:-$HERE/../install}
STATUS=0

# The version shown in the page header, read from the library's CMake project.
TREX_VERSION=$(awk '/^project\(TRexSelector/ {f=1} f && /VERSION/ {print $2; exit}' "$HERE/../CMakeLists.txt")
export TREX_VERSION

echo "=== configure and build the examples ($(date +%H:%M:%S))"
cmake -S "$HERE/examples" -B "$HERE/examples/build" \
      -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_BUILD_TYPE=Release \
      > "$HERE/_build_examples.log" 2>&1 \
  && cmake --build "$HERE/examples/build" --parallel 4 >> "$HERE/_build_examples.log" 2>&1
rc=$?
echo "build exit $rc"
if [ "$rc" -ne 0 ]; then
    grep -n -E "error|Error" "$HERE/_build_examples.log" | head -20
    echo "=== done, exit 1 ($(date +%H:%M:%S))"
    exit 1
fi

echo "=== run the examples ($(date +%H:%M:%S))"
rm -rf "$HERE/examples/output"
cmake --build "$HERE/examples/build" --target run_examples > "$HERE/_run_examples.log" 2>&1
rc=$?
echo "run exit $rc, $(ls "$HERE/examples/output" 2>/dev/null | wc -l | tr -d ' ') output sections"
[ "$rc" -ne 0 ] && { tail -20 "$HERE/_run_examples.log"; STATUS=1; }

# The figures are drawn from the examples' data with matplotlib. Without a Python that has it, the
# committed PNG files under pages/figures stay as they are.
echo "=== figures ($(date +%H:%M:%S))"
PY=${TREX_DOCS_PYTHON:-python3}
if "$PY" -c "import matplotlib, numpy" > /dev/null 2>&1; then
    "$PY" "$HERE/figures/make_figures.py" > "$HERE/_build_figures.log" 2>&1
    rc=$?
    echo "figures exit $rc ($PY)"
    [ "$rc" -ne 0 ] && { tail -20 "$HERE/_build_figures.log"; STATUS=1; }
else
    echo "no matplotlib for $PY, the committed figures stay (set TREX_DOCS_PYTHON)"
fi

echo "=== doxygen $TREX_VERSION ($(date +%H:%M:%S))"
mkdir -p "$HERE/_build"
# Doxygen 1.18.0 as released dies at random with a bus error, see the README for the patch.
(cd "$HERE" && doxygen Doxyfile > "$HERE/_build_doxygen.log" 2>&1)
rc=$?
echo "doxygen exit $rc, $(wc -l < "$HERE/_build/doxygen_warnings.log" | tr -d ' ') warning lines"
[ "$rc" -ne 0 ] && STATUS=1

# Doxygen's bibliography style prints no DOI, so the links are added to the citation list afterwards.
if [ "$rc" -eq 0 ] && [ -f "$HERE/_build/html/citelist.html" ]; then
    python3 "$HERE/tools/add_doi_links.py" "$HERE/_build/html/citelist.html" "$HERE/../../references.bib" \
        || STATUS=1
fi

echo "=== done, exit $STATUS ($(date +%H:%M:%S))"
exit $STATUS
