# T-Rex+DA: dependent variables {#trex_da}

[TOC]

## Why dependent variables need a correction

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
| `AR1` | first-order autoregressive: a variable is compared with its neighbours inside a window | `cor_coef` (estimated when `-2`), `rho_thr_DA` |
| `EQUI` | equi-correlated: a variable is compared with all other variables | `cor_coef` (estimated when `-2`), `rho_thr_DA` |

Each method stands for the dependency model of its paper and is made for data that follow that
model. The binary tree method follows the hierarchical graphical model of \cite machkour2025datrex.
It is made for variables that form disjoint groups of highly correlated variables, and it was
studied with one active variable per group. The autoregressive method follows \cite machkour2024ar1.
It is made for variables whose correlation decays with their distance in the column order, as in a
first-order autoregressive process, and it does not assume one active variable per group. The
nearest-neighbour method follows \cite machkour2025indextracking. It is made for overlapping groups
of highly correlated variables, as they are characteristic for stock returns. The equi-correlated
method is made for variables that all share one correlation.

The binary tree method clusters the variables by correlation distance with the linkage `hc_linkage`
(`LinkageMethod::Single` by default), evaluates the FDP estimate at `hc_grid_length` correlation
thresholds (the default `0` lets the core take \f$\min(20, p)\f$ points) and keeps the threshold
`rho_thresh`. The autoregressive and equi-correlated methods draw the same independent dummies as the
other methods and penalize the relative occurrences. The autoregressive method compares a variable
with the neighbours inside a window, which holds the lags whose correlation `|cor_coef|^lag` reaches
`rho_thr_DA`, a threshold in \f$[0, 1)\f$. The equi-correlated method compares a variable with all
other variables and is skipped when `|cor_coef|` does not exceed `rho_thr_DA`. Both methods use the
magnitude of the coefficient. Setting `cor_coef` to `-2` (the constant `AUTO_ESTIMATE_CORRELATION`)
estimates the coefficient from the data, and any other value outside \f$[-1, 1]\f$ gives a warning
and is replaced by the estimate.

\snippet{trimleft} trex_da.cpp ar1
\output{trex_da.ar1.txt}

## The penalty of the last step

All five methods penalize the relative occurrences in the same way. With \f$\Phi_t(j)\f$ the
relative occurrence of variable \f$j\f$ after \f$t\f$ included dummies (step \f$t\f$ of the T-loop)
and \f$T\f$ the stopping time at which the path is evaluated, its last step, the selector works with

\f[
\Phi_t^{\mathrm{DA}}(j) = \Psi_T(j) \, \Phi_t(j), \qquad
\Psi_T(j) = \frac{1}{2 - \min_{j'} |\Phi_T(j) - \Phi_T(j')|}, \qquad t = 1, \dots, T,
\f]

where \f$j'\f$ runs over the variables that \f$j\f$ is compared with, its group, its window or all
other variables. A variable without a group gets \f$\Psi_T(j) = 1/2\f$.

This is a deliberate deviation from the papers. Read literally, they recompute the penalty at every
step, \f$\Psi_t(j)\f$ from \f$\Phi_t\f$, while the selector takes the penalty of the last step for
all steps. The reason is the FDP estimator. It sums over the increments
\f$\Phi_t^{\mathrm{DA}}(j) - \Phi_{t-1}^{\mathrm{DA}}(j)\f$ and counts expected false discoveries
only while no increment is negative. One penalty per variable keeps the increments non-negative. A
penalty per step does not, because the penalty of a variable changes from step to step, and the FDP
estimate then becomes too small.

A Monte Carlo comparison of the two constructions with 1,000 runs per setting on block,
autoregressive and independent designs supports the choice. With the nearest-neighbour method the
per-step penalty exceeded a target FDR of 0.1, with 0.15 on an independent design and 0.13 on a
block design, where the last-step penalty gave 0.00 and 0.03. No setting had the last-step penalty
above the target and the per-step penalty below it. The per-step penalty selected more variables
only where its estimate was invalid, and on autoregressive designs the last-step penalty had the
higher true positive rate. The methods `AR1` and `EQUI` used the per-step penalty up to version
1.0.0, so their selections can differ slightly from earlier versions, while `BT`, `NN` and the prior
groups have always used the penalty of the last step.

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

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of T-Rex+DA under autoregressive, block,
heavy-tailed, grouped, nearest-neighbour and equi-correlated dependency, with their results and a
README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/trex_selector_methods/ad/trex_da),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/trex_selector_methods/ad/trex_da)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/trex_selector_methods/ad/trex_da).