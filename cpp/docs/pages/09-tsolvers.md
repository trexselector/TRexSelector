# The T-solvers {#tsolvers}

[TOC]

## Terminating forward selection

A T-solver is a forward selection algorithm on the enlarged design \f$[\boldsymbol X, \boldsymbol
D]\f$ of predictors and dummies that stops as soon as `T_stop` dummies are active. The T-Rex
selectors run one inside every random experiment, and every T-solver is available as a class of its
own, for a single terminated path, for a study of the solvers, or as a building block. The
terminated form of each algorithm follows the T-Rex paper \cite machkour2025trex, the greedy family
the terminating greedy pursuits \cite scheidt2026greedy. The classes live in four namespaces under
`trex::tsolvers::ad::linear_model`, one per family, and each family has an umbrella header that
includes its solvers.

| Family | Header | Classes | Algorithm |
| --- | --- | --- | --- |
| LARS | `lars_based.hpp` | TLARS_Solver, TLASSO_Solver, TSTEPWISE_Solver | least angle regression and its lasso and stepwise modifications \cite efron2004lars and \cite tibshirani1996lasso |
| LARS | `lars_based.hpp` | TSTAGEWISE_Solver | forward stagewise \cite hastie2007stagewise and \cite tibshirani2015stagewise |
| LARS | `lars_based.hpp` | TENET_Solver, TENETAug_Solver | the elastic net \cite zou2005elasticnet, Gram based and row augmented |
| LARS | `lars_based.hpp` | TIENET_Solver, TIENETAug_Solver | the informed elastic net \cite machkour2023ien, pathwise and row augmented |
| LARS | `lars_based.hpp` | TIKNET_Solver | the Tikhonov net, pathwise |
| OMP | `omp_based.hpp` | TOMP_Solver, TOOLS_Solver | orthogonal matching pursuit and orthogonal least squares \cite pati1993omp, \cite chen1989ols and \cite blumensath2007ompols |
| OMP | `omp_based.hpp` | TGP_Solver, TACGP_Solver | gradient and approximate conjugate gradient pursuit \cite blumensath2008gp |
| OMP | `omp_based.hpp` | TMP_Solver, TNCGMP_Solver | matching pursuit \cite mallat1993mp and non-negative constrained generalized matching pursuit \cite locatello2017mpfw |
| AFS | `afs_based.hpp` | TAFS_Solver | adaptive forward stepwise \cite zhang2026afs |
| CD | `cd_based.hpp` | TCDLASSO_Solver, TCDENET_Solver, TCDIENET_Solver, TCDTIKNET_Solver | pathwise coordinate descent \cite friedman2007pathwise and \cite friedman2010glmnet with exact minimizers at every dummy crossing |

The path solvers (LARS and CD families) walk a regularization path, the greedy solvers (OMP and AFS
families) add one variable per step. The informed solvers take a `groups` vector, the Tikhonov
solvers a sparse penalty matrix.

The headers and names of this page, and two helpers. The first prints index vectors, the second
constructs a solver of a given type as a reader, runs it to three dummies and prints its steps and
active predictors:

\snippet{trimleft} tsolvers.cpp includes

\snippet{trimleft} tsolvers.cpp helpers

## One terminated path

The data are drawn with a fixed seed. `X` and `D` are column-major Eigen matrices, mapped into the
solvers, which share them without a copy.

\snippet{trimleft} tsolvers.cpp data

A solver is constructed on `X`, `D`, `y` and a SolverOptions object, and `executeStep()` runs the
path until `T_stop` dummies are active. The path can be continued by calling `executeStep()` again
with a larger `T_stop`.

\snippet{trimleft} tsolvers.cpp path
\output{tsolvers.path.txt}

Indices are 0-based and the dummies continue after the predictors, so the active set of the enlarged
design runs from 0 to \f$p + L - 1\f$.

\snippet{trimleft} tsolvers.cpp actives
\output{tsolvers.actives.txt}

`getActions()` records the path as signed codes per step in the convention of the R package lars,
`+(j + 1)` for the addition and `-(j + 1)` for the removal of the 0-based variable `j`. `getBeta()`
is the coefficient vector of the last step and `getBetaPathSparse()` the whole path, one sparse step
after the other on the original coefficient scale.

\snippet{trimleft} tsolvers.cpp beta
\output{tsolvers.beta.txt}

## The path as a picture

The coefficient paths of a terminated run show what the solver does. A predictor enters at a step and
its coefficient grows from there, and a dummy that enters marks one of the `T_stop` crossings that end
the path. One solver runs to `T_stop = 1`, writes a snapshot of its path, and continues to `T_stop = 3`,
so the second panel extends the first. The documentation build draws the two snapshots with matplotlib
(`figures/make_figures.py`, run by `build_docs.sh`).

\snippet{trimleft} tsolvers.cpp figure
\output{tsolvers.figure.txt}

\image html tsolvers_path.png "Figure 1: T-LARS paths terminated at T_stop = 1 and continued to T_stop = 3." width=90%

The true predictors enter first and keep growing, further predictors follow with small coefficients,
and the dummies end the path. The second panel repeats the first up to its last step and continues
from there, because the continuation walks the stored path on. The label over each step names the
variable that enters.

## The options

SolverOptions holds every input except the data, the same struct for every solver. A solver reads
the fields it needs and ignores the rest. The fields have defaults, so a designated initializer sets
only what differs.

| Field | Meaning |
| --- | --- |
| `normalize`, `intercept`, `scaling_mode` | center and scale the columns of `X` and `D` in place (`ScalingMode::L2` or `ZSCORE`) and center `y` |
| `precentered` | the caller centered the data: the active set is capped at \f$n - 1\f$ and the intercept counts in the degrees of freedom |
| `tol` | the numerical tolerance of the comparisons (machine epsilon by default) |
| `tie_seed` | the seed of the tie-break engine (0 by default, so exact ties resolve reproducibly) |
| `n_threads` | the thread budget of the solver's own passes (1 serial, at least 1, and a pass splits only above about a million entries of the design) |
| `exch_tie_alpha`, `exch_tie_floor` | the exchangeable-tie band of the greedy solvers |
| `lambda2` | the quadratic penalty of the elastic nets (0 by default) |
| `rho_afs`, `ncgmp_variant`, `tenet_aug_use_lars` | solver-specific switches |
| `kappa_dummy`, `fold_dummies`, `group_penalty`, `cd_*` | the dummy coupling and the coordinate descent budgets |
| `verbose` | print the steps |

\snippet{trimleft} tsolvers.cpp options
\output{tsolvers.options.txt}

## Writers and readers

With `normalize = true` or `intercept = true` (the defaults) a solver is a writer. It centers and
scales the columns of `X` and `D` in place, holds the statistics, and writes the original values
back when `restore()` is called or when the object is destroyed. A second writer on a claimed
matrix throws, because a second in-place transformation would move the data by rounding. `y` is
never written.

\snippet{trimleft} tsolvers.cpp writers
\output{tsolvers.writers.txt}

A pipeline that runs several solvers on one design scales the data once in front with a scaler and
constructs readers, solvers with `normalize = false` and `intercept = false` that never write.

\snippet{trimleft} tsolvers.cpp readers

With the data pre-scaled, every family runs on the same design.

\snippet{trimleft} tsolvers.cpp families
\output{tsolvers.families.txt}

The informed solvers and the Tikhonov solvers take their structure as an extra argument.

\snippet{trimleft} tsolvers.cpp structured
\output{tsolvers.structured.txt}

## Saving and loading

`save()` serializes the state of a solver to a file and the static `load()` constructs a solver on
the same data in the same state. The loaded solver is a reader and continues the path where it was
saved.

\snippet{trimleft} tsolvers.cpp save-load
\output{tsolvers.save-load.txt}

## Memory

`memoryBytes()` reports the heap bytes a solver holds (its copy of `y`, the residual, the
correlations, the Cholesky factor, the coefficient path). `X` and `D` are borrowed and not counted.

\snippet{trimleft} tsolvers.cpp memory
\output{tsolvers.memory.txt}

The last section removes the folder of the saved state.

\snippet{trimleft} tsolvers.cpp cleanup
