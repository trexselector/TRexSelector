# TRexSelector 2.0.3

## Summary of the update

A patch release of 2.0.2, which is on CRAN since 2026-11-18. It fixes the ERROR that the check flavor
r-devel-linux-x86_64-fedora-clang reports for 2.0.2.

* The selectors lost the class and the message of an error that was raised inside their random
  experiments, on systems where the package is built with libc++ while the R process also holds the GNU
  C++ library. R then reported "c++ exception (unknown reason)", which failed the test
  `test-trex_gvs.R:769` on that flavor. The library carried such an error out of its OpenMP region by
  `std::exception_ptr`. It records the class and the message now and throws a fresh exception after the
  region. The defect was present since 2.0.0. Results do not change.

* The tests ask for two threads at most. One test of 2.0.2 ran a solver step with four threads, which
  the same flavor reported as "OMP: Warning #96: Cannot form a team with 4 threads, using 2 instead".

NEWS.md lists the change. The package has no reverse dependencies on CRAN.

## Test environments

* macOS 27.2 (Apple silicon), R 4.6.1: `R CMD check --as-cran`, also with the CPU-to-elapsed
  thresholds CRAN applies (`_R_CHECK_*_TIMING_CPU_TO_ELAPSED_THRESHOLD_` at 2.5).
* R-devel (2026-11-08 r90859, GCC 16.2) in a rocker/r-devel container (Debian, x86_64):
  `R CMD check --as-cran`.
* clang 18 with libc++ 18 in an Ubuntu 24.04 container (arm64), R 4.3.3, with the flags of the
  fedora-clang flavor (`-O3 -Wall -pedantic -frtti -Wp,-D_FORTIFY_SOURCE=3`), Rcpp built the same way and
  an R process that also holds the GNU C++ library. This is the configuration of the reported ERROR.
  With 2.0.2 the container reproduces the message of the flavor. With 2.0.3 the failing call carries its
  message, and the testthat suite and the examples pass under `OMP_THREAD_LIMIT=2` without a request for
  more threads.
* The same container with clang and the GNU C++ library, and with the sequential build without OpenMP:
  the testthat suite and the examples pass.
* Ubuntu 24.04 (arm64), R 4.3.3, GCC 13.3: install, the testthat suite and `R CMD check`.
* clang 14 with libc++ 14 (Ubuntu 24.04 container), the closest available match to CRAN's macOS
  builders for R-release on x86_64 and for R-oldrel: the 37 sources of the package compile without a
  warning.
* win-builder, R-release (R 4.6.2) and R-devel (2026-11-17 r90926, GCC 14.3): `R CMD check --as-cran`,
  the vignettes built there as well.
* The 20 containers that the R-hub project provides for CRAN's Linux flavors and additional checks
  (x86_64, R-devel 2026-11-18 r90935 or R 4.6.2), each with the check settings of its image and
  `OMP_THREAD_LIMIT=2`:
  * `clang23` (clang 23 with libc++ and the flags of the fedora-clang flavor): `R CMD check` shows the
    reported ERROR for 2.0.2, word for word, and Status OK for 2.0.3.
  * `ubuntu-clang`, `ubuntu-gcc16`, `ubuntu-next` and `ubuntu-release` (the counterparts of debian-clang,
    debian-gcc, r-patched-linux and r-release-linux), `debian-clang`, `gcc16`, `atlas`, `blis`, `mkl`,
    `nold`, `lto`, `donttest` and `nosuggests`: no ERROR and no WARNING.
  * `clang-asan`, `clang-ubsan` and `gcc-asan`: Status OK, and no sanitizer report in any output.
  * `valgrind`: Status OK, no invalid access and no use of an uninitialised value. The one record,
    368 bytes "possibly lost" below `pthread_create` in libgomp, is the thread-local storage of the
    OpenMP worker thread. A ten-line OpenMP program without code of the package shows it there as well.
  * `rchk`: no report for a source file of the package. The reported lines are in Rcpp's headers.
  * `vnu`: the HTML of the help pages validates (errors 0, warnings 0).
* The package's GitHub Actions jobs: install, the testthat suite and `R CMD check` with the vignettes
  built and re-built on Ubuntu, macOS and Windows (R release), the sequential build without OpenMP on
  Ubuntu, and Windows with the daily R-devel snapshot.

## R CMD check results

0 errors | 0 warnings | 0 notes on macOS. win-builder reports the incoming-feasibility NOTE alone on
R-release and on R-devel, "Days since last update: 2".

Two NOTEs in the R-devel container, both from the container itself:

* checking compilation flags used ... NOTE: non-portable flags `-Wdate-time`,
  `-Werror=format-security` and `-Wformat`. They come from the Makeconf of Debian's R-devel build,
  not from the package's configure or Makevars.
* checking HTML version of manual ... NOTE: math rendering skipped, as V8 is unavailable there.

The R-hub containers report one NOTE at most, from the container itself: the math rendering NOTE above
in `ubuntu-clang`, `ubuntu-gcc16`, `ubuntu-next`, `ubuntu-release` and `lto`, and in `debian-clang` the
non-portable flag `-Wp,-D_FORTIFY_SOURCE=3` of the container's Makevars.

Two NOTEs on Ubuntu:

* checking installed package size ... NOTE: installed size is 7.8 MB, sub-directories of 1 MB or more:
  libs 5.1 MB, R 2.1 MB. The compiled backend (the Eigen template instantiations of the solvers)
  accounts for the first. The shared object is stripped on Linux (src/Makevars.in).
* checking DESCRIPTION meta-information ... NOTE: 'LinkingTo' for 'BH' is unused as it has no
  'include' directory. An artifact of Ubuntu's r-cran-bh binary, which ships no headers. The CRAN BH
  package has them, and the sources include the Boost headers through it.

Compilation emits warnings from RcppEigen's Eigen headers and from the Boost headers of BH, and none
from the package's own sources, also under `-Wall -pedantic` with clang 18: `-Wunused-but-set-variable`
under clang (SparseLU, "variable 'count' set but not used" and similar) and `-Wignored-attributes`
under GCC ("ignoring attributes on template argument"). They originate in the upstream headers.

The days since the last update are few because 2.0.2 is from 2026-11-18. This update answers the mail of
the CRAN Team, which asks for a correction of the ERROR on the check page of 2.0.2 before 2026-12-11.

## System requirements and OpenMP

A C++20 compiler (GCC 10, Clang 12, Apple Clang 14 or newer, Rtools on Windows) and R 4.3 or newer.
OpenMP is optional: `configure` uses R's `SHLIB_OPENMP_CXXFLAGS` where R provides them (Linux,
Windows). On macOS without them it looks for an `omp.h` (Homebrew libomp or the mac.r-project.org
OpenMP package) and links R's bundled libomp, and without one it compiles sequentially. Under R-devel
every compile line of the install log carries `-fopenmp`. The examples run in under 0.3 s each, and
the vignettes limit themselves to two threads.

## Downstream dependencies

There are no reverse dependencies on CRAN.
