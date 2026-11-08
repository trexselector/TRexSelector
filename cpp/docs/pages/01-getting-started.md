# Getting started {#getting_started}

[TOC]

## What the library does

The Terminating-Random Experiments (T-Rex) selector selects variables while controlling the false
discovery rate (FDR) at a target level chosen by the user, and it does so with a high true positive rate
(TPR) \cite machkour2025trex. It applies to low- and high-dimensional designs and is tailored to the
high-dimensional problem. It takes a design matrix \f$\boldsymbol{X}\f$ with \f$n\f$ observations of
\f$p\f$ variables, a response vector \f$\boldsymbol{y}\f$ and the target FDR, and runs \f$K\f$ random
experiments in which the real variables compete with computer-generated dummy variables, random columns
with finite mean and variance that act as flagged null variables. In every experiment the dummies are
appended to the design, a forward selection solver walks its path and stops as soon as \f$T\f$ dummies
have entered, and the share of experiments in which a variable was selected is its vote. Because a dummy
enters as easily as a null variable, the included dummies give an estimate of the false discoveries, and
a calibration over the number of dummies, the stopping time and the voting level selects as many
variables as the FDR target allows, with a finite-sample guarantee that holds for a broad class of
forward selection solvers. The scheme is not tied to the linear model and scales to millions of
variables. The library currently covers the supervised linear model with a response vector and a design
matrix.

The C++ library is the core of the TRexSelector project. The R package and the Python package are
bindings to it, so the three surfaces select the same variables on the same data. The library is
written in C++20 on top of Eigen and installs as static libraries with CMake package files.

The library currently has five groups of tools.

- **T-Rex selectors.** The core TRexSelector and its variants for dependent variables
  (TRexDASelector), grouped variables (TRexGVSSelector), fast screening (ScreenTRexSelector,
  BiobankScreenTRex), sparse principal component analysis (PCA) (TRexSPCA) and a general Tikhonov
  penalty (TRexTikhonovSelector), under `src/trex_selector_methods/ad/`.
- **T-solvers.** The twenty terminating forward selection solvers the selectors run inside their
  random experiments, in four families (least angle regression, orthogonal matching pursuit,
  adaptive forward stepwise and coordinate descent based) under `src/tsolvers/ad/`.
- **Machine learning (ML) methods.** Scalers, ridge and elastic-net regression with cross-validation,
  singular value decomposition (SVD), PCA and hierarchical clustering under `src/ml_methods/`.
- **Memory-mapped matrices.** Disk-backed designs are accepted everywhere an Eigen map is accepted.
  The algorithms are agnostic to the data representation, so the library processes designs larger
  than the random-access memory (RAM) of the machine out of core, under `src/utils/memmap/`.
- **Utilities.** The selection metrics (false discovery proportion, true positive proportion,
  precision and recall), the synthetic data generators and the threading functions under `src/utils/`.

## Installation

The library is built from source with CMake. It needs a C++20 compiler, Eigen3, Boost and cereal,
and OpenMP for its parallel regions. The README of the library lists the dependencies per platform
and the CMake presets. The local install puts the headers, the static libraries and the CMake package
files under `cpp/install/`:

```bash
cmake --workflow --preset local-install
cmake --install build/release-portable
```

A downstream project finds the library with `find_package()` and links one target, which brings
Eigen, cereal and the other dependencies with it:

```cmake
find_package(TRexSelector REQUIRED)
target_link_libraries(my_app PRIVATE TRexSelector::trex_selector_methods)
```

The example programs of this guide are built exactly like that against the local install and run
when the documentation is built. Their output is shown under the code.

The library builds on Eigen \cite guennebaud2010eigen for its linear algebra, on Boost
\cite boost2026libraries with Boost.Interprocess \cite gaztanaga2008interprocess for the memory mapping,
on cereal \cite grant2017cereal for the serialization, and on OpenMP \cite openmp2021spec for its
parallel regions.

## A first run

The headers of the selector and the metrics, and the names used on this page:

\snippet{trimleft} getting_started.cpp includes

A small helper prints index vectors:

\snippet{trimleft} getting_started.cpp helpers

A Gaussian data set with \f$n = 50\f$ observations, \f$p = 100\f$ variables and three active
variables, drawn with a fixed seed. Indices are 0-based.

\snippet{trimleft} getting_started.cpp data
\output{getting_started.data.txt}

The selector is a class. It takes Eigen maps of the design and the response, so the data stay where
they are and the design is used in place. Construct it, run it with `select()` and read the results
from its getters. A fixed `seed` makes the dummy draws reproducible. The default `seed = -1` draws
them from the hardware entropy source, which is the right choice for Monte Carlo studies, see the T-Rex
selector page.

\snippet{trimleft} getting_started.cpp first-run
\output{getting_started.first-run.txt}

The two metrics of variable selection compare the selection with the true support. The false
discovery proportion (FDP) is the share of false positives among the selected variables and the true
positive proportion (TPP) is the share of the true support that was found. Their expectations are
the FDR and the TPR.

\snippet{trimleft} getting_started.cpp metrics
\output{getting_started.metrics.txt}

## The result fields

Every selector exposes the calibration it arrived at.

\snippet{trimleft} getting_started.cpp result-fields
\output{getting_started.result-fields.txt}

`getFDPHatMat()` returns the FDP estimate per stopping time (rows) and voting level (columns),
`getPhiMat()` the relative occurrence of every variable over the \f$K\f$ experiments after each stopping
time, and `getVotingGrid()` the voting levels. The T-Rex selector page walks through them.

\snippet{trimleft} getting_started.cpp matrices
\output{getting_started.matrices.txt}

## T-Rex variants of the library

| Situation | Selector | Page |
| --- | --- | --- |
| Independent or modestly dependent variables | TRexSelector | T-Rex selector |
| Strongly dependent variables (autoregressive, equi-correlated, clustered) | TRexDASelector | T-Rex+DA |
| Groups of correlated variables that should enter together | TRexGVSSelector | T-Rex GVS |
| Very large \f$p\f$, a fast first pass before the full selector | ScreenTRexSelector | Screen-T-Rex |
| Many phenotypes against one genotype matrix | BiobankScreenTRex | Screen-T-Rex |
| Sparse principal components | TRexSPCA | T-Rex SPCA |
| A structured penalty (smoothness, graph Laplacian, groups) | TRexTikhonovSelector | T-Rex Tikhonov |

The solver inside the random experiments is set in `TRexControlParameter::solver_type`. T-LARS
(terminating least angle regression) is the default and the solver of the T-Rex paper. T-LASSO adds
the lasso modification, T-ENET the elastic net, the orthogonal matching pursuit (OMP) family and the
adaptive forward stepwise family are greedy and faster per step, and the coordinate descent (CD)
family computes exact penalized minimizers at every dummy crossing. The T-solvers page describes the
families, and every solver can be run on its own.

\snippet{trimleft} getting_started.cpp solver
\output{getting_started.solver.txt}

## Citing the library

The citation metadata of the repository is in its `CITATION.cff`. The methods are published in
\cite machkour2025trex (T-Rex), \cite machkour2025datrex (T-Rex+DA), \cite machkour2022gvs and
\cite machkour2023ien (T-Rex GVS), \cite machkour2023screen (Screen-T-Rex), \cite machkour2024spca
(T-Rex SPCA), \cite scheidt2023bigtrex (memory mapping) and \cite scheidt2026greedy (the terminating
greedy pursuits).
