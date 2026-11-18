# TRexSelector 2.0.2

## Summary of the update

A patch release of 2.0.1, which is on CRAN since 2026-11-12.
List of corrections:

* `TRexDASelector` with the methods `"AR1"` and `"EQUI"` recomputed its penalty at every step of the
  path. The deflated occurrences were then not monotone, and the estimate of the false discovery
  proportion was too small. The two methods take the penalty of the last step now, as the other
  methods of the selector always did.

* `TRexDASelector` with a greedy solver exceeded a target FDR of 0.2 on strongly correlated designs
  (0.21 to 0.27 at a correlation of 0.9 in a Monte Carlo study). It runs these solvers with the
  exchangeable-tie band now, which was an option that the user had to set, and holds the target there
  (0.13 to 0.17).

* `TRexTikhonovSelector` ran by default on a coordinate-descent solver that can stop at its sweep
  budget before it has converged. On the example of the vignette its default run selected 33 variables
  with a false discovery proportion of 0.70 at a target of 0.1. The selector runs on the homotopy
  solver of the same problem now.

* A second `select()` on the same selector object ran on the design at its original scale and returned
  a wrong selection without an error. It is an error now.

* `TRexBiobankScreeningSelector` drew one set of dummies for all phenotypes of a seeded run, and
  `TRexScreeningSelector` reported an FDR estimate of 0 for an empty selection.

* A coordinate-descent solve that had not converged gave one warning, and the selection was built on
  the iterate it had. It is an error now.

Beside the corrections the coordinate-descent solver family is deprecated ahead of its removal in
2.1.0, the screening selector gains a nearest-neighbour variant, a read-only memory-mapped design is
accepted wherever nothing is written into it, and the help pages and vignettes are corrected where
they disagreed with the code.
The Rcpp bindings compile as one unit, so the package builds 37 objects
where 2.0.1 built 48 and installs in less time. NEWS.md lists every change. The package has no
reverse dependencies on CRAN.

## Test environments

* macOS 27.2 (Apple silicon), R 4.6.1: `R CMD check --as-cran`, also with the CPU-to-elapsed
  thresholds CRAN applies (`_R_CHECK_*_TIMING_CPU_TO_ELAPSED_THRESHOLD_` at 2.5).
* R-devel (2026-11-08 r90859, GCC 16.2) in a rocker/r-devel container (Debian, x86_64):
  `R CMD check --as-cran`.
* Ubuntu 24.04 (arm64), R 4.3.3, GCC 13.3: install, the testthat suite and `R CMD check`.
* clang 14 with libc++ 14 (Ubuntu 24.04 container), the closest available match to CRAN's macOS
  builders for R-release on x86_64 and for R-oldrel: the 37 sources of the package compile without a
  warning.
* win-builder, R-release (R 4.6.2) and R-devel (2026-11-16 r90914, GCC 14.3): `R CMD check --as-cran`,
  the vignettes built there as well.
* The package's GitHub Actions jobs: install, the testthat suite and `R CMD check` with the vignettes
  built and re-built on Ubuntu, macOS and Windows (R release), the sequential build without OpenMP on
  Ubuntu, and Windows with the daily R-devel snapshot.

## R CMD check results

0 errors | 0 warnings | 0 notes on macOS. win-builder reports the incoming-feasibility NOTE alone on
R-release and on R-devel, "Days since last update: 6", which the corrected results above explain.

Two NOTEs in the R-devel container, both from the container itself:

* checking compilation flags used ... NOTE: non-portable flags `-Wdate-time`,
  `-Werror=format-security` and `-Wformat`. They come from the Makeconf of Debian's R-devel build,
  not from the package's configure or Makevars.
* checking HTML version of manual ... NOTE: math rendering skipped, as V8 is unavailable there.

Two NOTEs on Ubuntu:

* checking installed package size ... NOTE: installed size is 7.8 MB, sub-directories of 1 MB or more:
  libs 5.1 MB, R 2.0 MB. The compiled backend (the Eigen template instantiations of the solvers)
  accounts for the first. The shared object is stripped on Linux (src/Makevars.in).
* checking DESCRIPTION meta-information ... NOTE: 'LinkingTo' for 'BH' is unused as it has no
  'include' directory. An artifact of Ubuntu's r-cran-bh binary, which ships no headers. The CRAN BH
  package has them, and the sources include the Boost headers through it.

Compilation emits warnings from RcppEigen's Eigen headers and none from the package's own sources:
`-Wunused-but-set-variable` under clang (SparseLU, "variable 'count' set but not used" and similar)
and `-Wignored-attributes` under GCC ("ignoring attributes on template argument"). They originate in
the upstream headers.

## System requirements and OpenMP

A C++20 compiler (GCC 10, Clang 12, Apple Clang 14 or newer, Rtools on Windows) and R 4.3 or newer.
OpenMP is optional: `configure` uses R's `SHLIB_OPENMP_CXXFLAGS` where R provides them (Linux,
Windows). On macOS without them it looks for an `omp.h` (Homebrew libomp or the mac.r-project.org
OpenMP package) and links R's bundled libomp, and without one it compiles sequentially. Under R-devel
every compile line of the install log carries `-fopenmp`. The examples run in under 0.3 s each, and
the vignettes limit themselves to two threads.

## Downstream dependencies

There are no reverse dependencies on CRAN.
