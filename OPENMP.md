# OpenMP on your platform

The TRexSelector library and its R and Python packages run their parallel regions with OpenMP (Open
Multi-Processing): the random experiments of a selector, the passes of the T-solvers over the design, and
the scalers, the clustering and the machine learning (ML) methods. A build without OpenMP runs the same
work in sequence and selects the same variables. This page says how each package gets OpenMP on each
platform and how to check it.

## The check

| Package | Call | Without OpenMP |
| --- | --- | --- |
| R | `trex_get_max_threads()` | returns 1 |
| Python | `trexselector.utils.get_max_threads()` | returns 1 |
| C++ | `find_package(TRexSelector)` sets `TRexSelector_WITH_OPENMP` | OFF |

## R

- **Windows.** The binary on the Comprehensive R Archive Network (CRAN) is built with Rtools and carries
  OpenMP.
- **Linux.** A source install compiles with the OpenMP flags of R. The configure script of the package
  reads them from R's Makeconf, since R 4.6 no longer serves `SHLIB_OPENMP_CXXFLAGS` through
  `R CMD config`.
- **macOS.** CRAN builds R for macOS without OpenMP, so the CRAN binary of the package runs on one
  thread. A build from source uses the `omp.h` header of Homebrew's libomp together with the OpenMP
  runtime that R itself ships.

  ```sh
  brew install libomp
  ```

  ```r
  install.packages("TRexSelector", type = "source")
  trex_get_max_threads()
  ```

  Without the header the package compiles single-threaded and `trex_get_max_threads()` returns 1.
- `TREX_OPENMP=no` in the environment of the installation forces the single-threaded build on every
  platform.

The sources behind these rules are the OpenMP sections of
[R Installation and Administration](https://cran.r-project.org/doc/manuals/r-release/R-admin.html) and
[Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html), and CRAN's page on
[OpenMP for R on macOS](https://mac.r-project.org/openmp/).

## Python

Every wheel on the Python Package Index (PyPI) bundles the OpenMP runtime of its platform, so the package
runs on several threads out of the box. A build from source takes OpenMP from the toolchain. On macOS
install Homebrew's LLVM and libomp (`brew install llvm libomp`), which the build picks up, because Apple
clang has no OpenMP. The switch below forces the single-threaded build.

```sh
pip install -C cmake.define.TREX_WITH_OPENMP=OFF Python/trexselector
```

## C++

CMake finds OpenMP through its FindOpenMP module and reports "OpenMP not found: building
single-threaded" otherwise. On macOS `brew install llvm libomp` provides both, and the build prefers
Homebrew's LLVM. `-DTREX_WITH_OPENMP=OFF` forces the single-threaded build. The
[C++ README](cpp/README.md) lists the dependencies per platform and describes the thread budgets.

## The CRAN rule

CRAN checks packages on at most two cores. The tests and the vignettes of the R package cap the thread
count with `trex_set_num_threads(min(2L, trex_get_max_threads()))` for their duration and restore it at
the end, the idiom for any package that calls TRexSelector from its own tests.
