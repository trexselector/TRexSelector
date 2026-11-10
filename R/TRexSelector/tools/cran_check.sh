#!/usr/bin/env bash
# Vendor check, source tarball and R CMD check --as-cran of the R package, the checks before a CRAN upload.
# The remote incoming check stays off (it fetches the package's web links), the timing thresholds are CRAN's.
#   bash R/TRexSelector/tools/cran_check.sh [out dir]
# Default out dir: notes/cran_upload/build_<date> (gitignored), so earlier tarballs stay.
set -u
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT="$(cd "$PKG_DIR/../.." && pwd)"
OUT=${1:-$ROOT/notes/cran_upload/build_$(date +%Y-%m-%d)}
mkdir -p "$OUT"
cd "$OUT" || exit 1

finish() { echo "=== done, exit $1 ($(date +%H:%M:%S))"; exit "$1"; }

echo "=== vendor check ($(date +%H:%M:%S))"
bash "$SCRIPT_DIR/bundle_for_cran.sh" --check > "$OUT/vendor_check.log" 2>&1
rc=$?
echo "vendor check exit $rc"
tail -2 "$OUT/vendor_check.log"
[ "$rc" -ne 0 ] && finish 1

echo "=== R CMD build ($(date +%H:%M:%S))"
rm -f "$OUT"/TRexSelector_*.tar.gz
R CMD build "$PKG_DIR" > "$OUT/build.log" 2>&1
rc=$?
echo "build exit $rc"
tail -2 "$OUT/build.log"
[ "$rc" -ne 0 ] && { grep -m5 -n -i 'error' "$OUT/build.log"; finish 1; }
TARBALLS=("$OUT"/TRexSelector_*.tar.gz)
TARBALL=${TARBALLS[0]}
echo "tarball: $(basename "$TARBALL"), $(wc -c < "$TARBALL" | tr -d ' ') bytes"
shasum -a 256 "$TARBALL" | awk '{print "sha256:", $1}'

echo "=== R CMD check --as-cran ($(date +%H:%M:%S))"
export _R_CHECK_CRAN_INCOMING_REMOTE_=false
export _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
export _R_CHECK_TESTS_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
export _R_CHECK_VIGNETTE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
rm -rf "$OUT/TRexSelector.Rcheck"
R CMD check --as-cran "$TARBALL" > "$OUT/check.log" 2>&1
rc=$?
echo "check exit $rc"
grep -E '^Status:|WARNING|NOTE|ERROR' "$OUT/TRexSelector.Rcheck/00check.log" | head -20
[ "$rc" -ne 0 ] && finish 1
finish 0
