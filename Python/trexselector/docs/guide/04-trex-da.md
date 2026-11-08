---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# T-Rex+DA: dependency-aware dummies

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Why dependent variables need their own dummies

The false discovery rate (FDR) guarantee of the T-Rex selector rests on the exchangeability of the
null variables and the dummies along the solution path. Independent dummies are exchangeable with
independent null variables. When the variables are strongly dependent, a null variable that is
correlated with an active one enters the path earlier than an independent dummy would, so the
dummies no longer mirror the nulls and the false discovery proportion (FDP) estimate becomes
anti-conservative. The dependency-aware T-Rex selector (T-Rex+DA) {cite}`machkour2025datrex`
restores the control by estimating the dependency structure and by deflating the vote within groups
of dependent variables. It has been applied to sparse index tracking with autoregressive and
nearest-neighbour stock dependency models {cite}`machkour2025indextracking,machkour2024ar1`.

## An autoregressive design

The variables follow a first-order autoregressive (AR(1)) process across the columns with
coefficient 0.7, and five of them are active. The data are drawn with a fixed seed, the indices are
0-based, and the design is created in Fortran order, the layout the core reads without a copy.

```{code-cell} python
import numpy as np
import trexselector as tsn
from trexselector.utils import compute_fdp, compute_tpp

rng = np.random.default_rng(329)
n, p, rho = 100, 100, 0.7
X = np.zeros((n, p), order="F")
X[:, 0] = rng.standard_normal(n)
for j in range(1, p):
    X[:, j] = rho * X[:, j - 1] + np.sqrt(1 - rho**2) * rng.standard_normal(n)
true_support = [9, 29, 49, 69, 89]
beta = np.zeros(p)
beta[true_support] = 2.0
y = X @ beta + rng.standard_normal(n)
```

The plain selector and the dependency-aware selector on the same data and the same seed.

```{code-cell} python
ctrl = tsn.TRexControlParameter()
ctrl.K = 10
plain = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=7, verbose=False)
plain.select()

da_ctrl = tsn.TRexDAControlParameter()
da_ctrl.method = tsn.DAMethod.BT
da = tsn.TRexDASelector(X, y, tFDR=0.1, da_control=da_ctrl, trex_control=ctrl, seed=7,
                         verbose=False)
da_result = da.select()

print(f"{'selector':<10} {'selected':>8} {'FDP':>6} {'TPP':>6}")
for name, s in (("T-Rex", plain), ("T-Rex+DA", da)):
    print(f"{name:<10} {len(s.selected_indices):>8} "
          f"{compute_fdp(s.selected_indices, true_support):>6.2f} "
          f"{compute_tpp(s.selected_indices, true_support):>6.2f}")
```

`TRexDASelector` inherits the attributes of `TRexSelector` (`selected_indices`, `fdp_hat_mat`,
`T_stop`, `L`, ...), and its `select()` returns a `DASelectionResult` with the record of the
dependency estimation.

```{code-cell} python
print(da_result.method)
print(da_result.rho_thresh)
print(np.round(da_result.rho_grid.ravel(), 3))
```

## The dependency methods

`method` of `TRexDAControlParameter` chooses how the dependency is handled, with the `DAMethod`
enum.

| Method | Dependency model | Parameters |
| --- | --- | --- |
| `BT` (default) | binary tree: hierarchical clustering of the variables at correlation thresholds, the threshold is calibrated on a grid | `hc_linkage`, `hc_grid_length` |
| `NN` | nearest neighbours: a correlation-threshold sweep of the same form | `hc_grid_length` |
| `AR1` | first-order autoregressive dummies | `cor_coef` (estimated when `-2`) |
| `EQUI` | equi-correlated dummies | `cor_coef`, `rho_thr_DA` |

The binary tree method clusters the variables by correlation distance with the linkage `hc_linkage`
(`LinkageMethod.Single` by default), evaluates the FDP estimate at `hc_grid_length` correlation
thresholds (the default `0` lets the core take $\min(20, p)$ points) and keeps the threshold
`rho_thresh`. The autoregressive and equi-correlated methods draw dummies with the dependency of the
data. `cor_coef = -2` estimates the coefficient from the data and `rho_thr_DA` is the correlation
below which the equi-correlated model is skipped.

```{code-cell} python
ar1_ctrl = tsn.TRexDAControlParameter()
ar1_ctrl.method = tsn.DAMethod.AR1
ar1 = tsn.TRexDASelector(X, y, tFDR=0.1, da_control=ar1_ctrl, trex_control=ctrl, seed=7,
                          verbose=False)
ar1_result = ar1.select()
print(ar1.selected_indices)
print(round(ar1_result.cor_coef, 3))
```

## Prior groups

When the groups of dependent variables are known, pass them as `prior_groups`, a list of integer
sequences of length $p$ with a group label per variable. The dependency estimation then runs within
the given groups, with hierarchical clustering inside each group, and the calibration grid is built
from the within-group dendrogram heights. The deflation acts on data-driven neighbourhoods inside
the known groups, never on the raw groups themselves. Several levels can be supplied, and their
finest common refinement is the constraint. The result then reports `DAMethod.PRIOR_GROUPS` as the
method applied.

```{code-cell} python
blocks = np.repeat(np.arange(1, 11), 10).tolist()
pg_ctrl = tsn.TRexDAControlParameter()
pg_ctrl.prior_groups = [blocks]
pg = tsn.TRexDASelector(X, y, tFDR=0.1, da_control=pg_ctrl, trex_control=ctrl, seed=7,
                         verbose=False)
pg_result = pg.select()
print(pg.selected_indices)
print(pg_result.method)
```

## Greedy solvers under T-Rex+DA

The deflation relies on the spread of occurrences within a group of dependent variables. A greedy
solver (the orthogonal matching pursuit family, T-AFS) picks the top candidate of a tie
deterministically and removes that spread. `exch_tie_alpha` of the solver hyperparameters,
`solver_params` of `TRexControlParameter`, restores it, since statistically indistinguishable top
candidates within a correlated cluster are then picked uniformly at random per step. The recommended
value for greedy solvers under T-Rex+DA is 0.25. The path solvers ignore the option.

```{code-cell} python
greedy_ctrl = tsn.TRexControlParameter()
greedy_ctrl.solver_type = tsn.SolverTypeForTRex.TOMP
greedy_ctrl.K = 10
greedy_ctrl.solver_params.exch_tie_alpha = 0.25
bt_ctrl = tsn.TRexDAControlParameter()
bt_ctrl.method = tsn.DAMethod.BT
greedy = tsn.TRexDASelector(X, y, tFDR=0.1, da_control=bt_ctrl, trex_control=greedy_ctrl, seed=7,
                            verbose=False)
greedy.select()
greedy.selected_indices
```
