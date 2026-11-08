# T-Rex+DA: dependency-aware dummies {#trex_da}

[TOC]

## Why dependent variables need their own dummies

The false discovery rate (FDR) guarantee of the T-Rex selector rests on the exchangeability of the
null variables and the dummies along the solution path. Independent dummies are exchangeable with
independent null variables. When the variables are strongly dependent, a null variable that is
correlated with an active one enters the path earlier than an independent dummy would, so the
dummies no longer mirror the nulls and the false discovery proportion (FDP) estimate becomes
anti-conservative. The dependency-aware T-Rex selector (T-Rex+DA) \cite machkour2025datrex
restores the control by estimating the dependency structure and by deflating the vote within groups
of dependent variables. It has been applied to sparse index tracking with autoregressive and
nearest-neighbour stock dependency models \cite machkour2025indextracking and \cite machkour2024ar1.

## An autoregressive design

The headers and names of this page:

\snippet{trimleft} trex_da.cpp includes

Two small helpers print index vectors and the name of a dependency method:

\snippet{trimleft} trex_da.cpp helpers

The variables follow a first-order autoregressive (AR(1)) process across the columns with
coefficient 0.7, and five of them are active. The data are drawn with a fixed seed, the indices are
0-based, and the design is an Eigen matrix in column-major layout, the layout the core reads without
a copy.

\snippet{trimleft} trex_da.cpp data

The plain selector and the dependency-aware selector on the same data and the same seed. The
dependency-aware control TRexDAControlParameter nests the base control as `trex_ctrl`.

\snippet{trimleft} trex_da.cpp compare
\output{trex_da.compare.txt}

TRexDASelector inherits the getters of TRexSelector (`getSelectedIndices()`, `getFDPHatMat()`,
`getTStop()`, `getDummyMultiplierL()`, ...), and `getDAResult()` returns the `DASelectionResult` with
the record of the dependency estimation.

\snippet{trimleft} trex_da.cpp da-result
\output{trex_da.da-result.txt}

## The dependency methods

`method` of TRexDAControlParameter chooses how the dependency is handled, with the DAMethod enum.

| Method | Dependency model | Parameters |
| --- | --- | --- |
| `BT` (default) | binary tree: hierarchical clustering of the variables at correlation thresholds, the threshold is calibrated on a grid | `hc_linkage`, `hc_grid_length` |
| `NN` | nearest neighbours: a correlation-threshold sweep of the same form | `hc_grid_length` |
| `AR1` | first-order autoregressive dummies | `cor_coef` (estimated when `-2`) |
| `EQUI` | equi-correlated dummies | `cor_coef`, `rho_thr_DA` |

The binary tree method clusters the variables by correlation distance with the linkage `hc_linkage`
(`LinkageMethod::Single` by default), evaluates the FDP estimate at `hc_grid_length` correlation
thresholds (the default `0` lets the core take \f$\min(20, p)\f$ points) and keeps the threshold
`rho_thresh`. The autoregressive and equi-correlated methods draw dummies with the dependency of the
data. `cor_coef = -2` (the constant `AUTO_ESTIMATE_CORRELATION`) estimates the coefficient from the
data and `rho_thr_DA` is the correlation below which the equi-correlated model is skipped.

\snippet{trimleft} trex_da.cpp ar1
\output{trex_da.ar1.txt}

## Prior groups

When the groups of dependent variables are known, pass them as `prior_groups`, a vector of label
vectors of length \f$p\f$ with a group label per variable. The dependency estimation then runs within
the given groups, with hierarchical clustering inside each group, and the calibration grid is built
from the within-group dendrogram heights. The deflation acts on data-driven neighbourhoods inside
the known groups, never on the raw groups themselves. Several levels can be supplied, and their
finest common refinement is the constraint. The result then reports `DAMethod::PRIOR_GROUPS` as the
method applied.

\snippet{trimleft} trex_da.cpp prior-groups
\output{trex_da.prior-groups.txt}

## Greedy solvers under T-Rex+DA

The deflation relies on the spread of occurrences within a group of dependent variables. A greedy
solver (the orthogonal matching pursuit family, T-AFS) picks the top candidate of a tie
deterministically and removes that spread. `exch_tie_alpha` of the solver hyperparameters,
`solver_params` of TRexControlParameter, restores it, since statistically indistinguishable top
candidates within a correlated cluster are then picked uniformly at random per step. The recommended
value for greedy solvers under T-Rex+DA is 0.25. The path solvers ignore the option.

\snippet{trimleft} trex_da.cpp greedy
\output{trex_da.greedy.txt}
