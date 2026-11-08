# Working without data copies {#no_copies}

[TOC]

## The rule

The library works on the caller's data wherever possible. An Eigen matrix is bound through an
`Eigen::Map` of its storage, a memory-mapped file through its mapping, and the preprocessing a
method needs is applied in place and undone afterwards. A design of several gigabytes is therefore
held once, by the caller or by the file, and not a second time by the library. This page states the
contract that makes it safe. The map is Eigen's view of foreign storage \cite guennebaud2010eigen,
and a memory-mapped file comes from Boost.Interprocess \cite gaztanaga2008interprocess.

## Pass double matrices

The library binds double storage in column-major order, the layout Eigen reads without a copy. The
constructors take `Eigen::Map<Eigen::MatrixXd>` and `Eigen::Map<Eigen::VectorXd>`, so the type is
settled at compile time. An integer or single-precision matrix is cast to doubles by the caller,
once, a row-major matrix is copied into a column-major one, once, and the map of that storage passes
through untouched. To keep a large design copy-free, create it as a column-major matrix of doubles
from the start. The headers and names of this page:

\snippet{trimleft} no_copies.cpp includes

\snippet{trimleft} no_copies.cpp double
\output{no_copies.double.txt}

## Selectors: in place at construction, back at the end

A selector centers and scales `X` in place when it is constructed and writes the original values
back when `select()` returns. In between, `X` holds the working design. `restore()` writes back
early for a selector that will not run, and a selector that is destroyed without running writes back
in its destructor.

\snippet{trimleft} no_copies.cpp selector
\output{no_copies.selector.txt}

## T-solvers: writers and readers

A T-solver with `normalize = true` or `intercept = true` in its SolverOptions is a writer. It
transforms `X` and `D` in place, keeps the statistics, and holds a claim on both matrices until
`restore()` or its destructor writes them back. A second writer on a claimed matrix is refused,
because a second transformation would move the data by rounding. The `restore()` of a T-solver
takes the three maps it was constructed on.

\snippet{trimleft} no_copies.cpp writers
\output{no_copies.writers.txt}

A solver with `normalize = false` and `intercept = false` is a reader. It runs on the data as it
finds them and never writes, so any number of readers can share one design. A pipeline scales the
design once with a scaler and constructs readers.

\snippet{trimleft} no_copies.cpp readers
\output{no_copies.readers.txt}

`y` is never written. The solver centers an internal copy of the response.

## Loading a saved solver

`load()` restores a solver into one constructed on the same data in the same state. Load into a
reader on pre-scaled data, or into a writer constructed with the same options on the same raw data,
which then keeps the duty to write back.

## Mapped files are never read into memory

A memory-mapped matrix is bound through the same interface as an Eigen matrix. The library touches
the pages it needs and the file keeps the data, so the working set of a run is the library's own
buffers, never a copy of the design. The accessor `getMap()` of the handle is a view of the mapping
as well, and the element access reads from the file, so no accessor copies the design either.

## The low-memory engines

The ridge and elastic-net classes have a `low_memory` engine that avoids every allocation of the
size of `X`.

| Class | `low_memory = false` | `low_memory = true` |
| --- | --- | --- |
| `ridge_cv_svd` | each fold gathers and standardizes its training rows and takes a thin singular value decomposition (SVD), about three to four copies of `X` in flight | the dual (kernel) form, \f$n \times n\f$ kernels accumulated from column blocks read in place, memory \f$O(n^2)\f$ |
| `enet_gaussian` | coordinate descent on a centered and standardized working copy, one \f$n \times p\f$ matrix | coordinate descent on `X` itself, no allocation of the size of `X` |
| `enet_cv_ccd` | one working copy per fit | every fold on `X` itself with a training mask |

The two engines compute the same estimator on the same grid and agree to rounding. The selectors
that tune \f$\lambda_2\f$ by cross-validation (GVS, SPCA, Tikhonov) use these engines, and the elastic-net
tuner of the sparse principal component selector keeps its peak memory under 0.3 times the size of `X`.

\snippet{trimleft} no_copies.cpp low-memory
\output{no_copies.low-memory.txt}

## Measuring the footprint

`memoryBytes()` reports what a T-solver holds on the heap. `X` and `D` are borrowed and not
counted.

\snippet{trimleft} no_copies.cpp bytes
\output{no_copies.bytes.txt}
