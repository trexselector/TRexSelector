---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# The T-solvers

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Terminating forward selection

A T-solver is a forward selection algorithm on the enlarged design $[\boldsymbol X, \boldsymbol D]$
of predictors and dummies that stops as soon as `T_stop` dummies are active. The T-Rex selectors run
one inside every random experiment, and every T-solver is available as a class of its own, for a
single terminated path, for a study of the solvers, or as a building block. The terminated form of
each algorithm follows the T-Rex paper {cite}`machkour2025trex`, the greedy family the terminating
greedy pursuits {cite}`scheidt2026greedy`. The classes live in the four submodules of
`trexselector.tsolvers.ad`.

| Family | Submodule | Classes | Algorithm |
| --- | --- | --- | --- |
| LARS | `lars_based` | `TLARS_Solver`, `TLASSO_Solver`, `TSTEPWISE_Solver` | least angle regression and its lasso and stepwise modifications {cite}`efron2004lars,tibshirani1996lasso` |
| LARS | `lars_based` | `TSTAGEWISE_Solver` | forward stagewise {cite}`hastie2007stagewise,tibshirani2015stagewise` |
| LARS | `lars_based` | `TENET_Solver`, `TENETAug_Solver` | the elastic net {cite}`zou2005elasticnet`, Gram based and row augmented |
| LARS | `lars_based` | `TIENET_Solver`, `TIENETAug_Solver` | the informed elastic net {cite}`machkour2023ien`, pathwise and row augmented |
| LARS | `lars_based` | `TIKNET_Solver` | the Tikhonov net, pathwise |
| OMP | `omp_based` | `TOMP_Solver`, `TOOLS_Solver` | orthogonal matching pursuit and orthogonal least squares {cite}`pati1993omp,chen1989ols,blumensath2007ompols` |
| OMP | `omp_based` | `TGP_Solver`, `TACGP_Solver` | gradient and approximate conjugate gradient pursuit {cite}`blumensath2008gp` |
| OMP | `omp_based` | `TMP_Solver`, `TNCGMP_Solver` | matching pursuit {cite}`mallat1993mp` and non-negative constrained generalized matching pursuit {cite}`locatello2017mpfw` |
| AFS | `afs_based` | `TAFS_Solver` | adaptive forward stepwise {cite}`zhang2026afs` |
| CD | `cd_based` | `TCDLASSO_Solver`, `TCDENET_Solver`, `TCDIENET_Solver`, `TCDTIKNET_Solver` | pathwise coordinate descent {cite}`friedman2007pathwise,friedman2010glmnet` with exact minimizers at every dummy crossing |

The path solvers (LARS and CD families) walk a regularization path, the greedy solvers (OMP and AFS
families) add one variable per step. The informed solvers take a `groups` array, the Tikhonov
solvers a penalty matrix `tikhonov_K`.

## One terminated path

The data are drawn with a fixed seed. `X` and `D` are Fortran-ordered float64 arrays, which the
solvers share without a copy (a C-ordered array is copied first, so an in-place transformation would
not reach the caller's array).

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(31)
n = 100
p = 50
L = 50
true_support = [1, 8, 16, 32]
beta = np.zeros(p)
beta[true_support] = [2.0, -2.0, 1.5, -1.5]
X = np.asfortranarray(rng.standard_normal((n, p)))
D = np.asfortranarray(rng.standard_normal((n, L)))
y = X @ beta + rng.standard_normal(n)
```

A solver is constructed on `X`, `D`, `y` and a `SolverOptions` object, and `executeStep()` runs the
path until `T_stop` dummies are active. The path can be continued by calling `executeStep()` again
with a larger `T_stop`.

```{code-cell} python
from trexselector.tsolvers.ad import SolverOptions
from trexselector.tsolvers.ad.lars_based import TLARS_Solver

solver = TLARS_Solver(X, D, y)
solver.executeStep(T_stop=3)
print("steps:            ", solver.getNumSteps())
print("active predictors:", solver.getActivePredictorIndices())
print("active dummies:   ", solver.getActiveDummyIndices())
```

Indices are 0-based and the dummies continue after the predictors, so the active set of the enlarged
design runs from 0 to $p + L - 1$.

```{code-cell} python
print("actives:          ", solver.getActives())
print("dummy start index:", solver.getDummyStartIndex())
print("first actions:    ", solver.getActions()[:4])
```

`getActions()` records the path as signed codes per step in the convention of the R package lars,
`+(j + 1)` for the addition and `-(j + 1)` for the removal of the 0-based variable `j`. `getBeta()`
is the coefficient vector of the last step and `getBetaPath()` the whole path, one column per step
with the all-zero start.

```{code-cell} python
print("beta path:", solver.getBetaPath().shape)
print("beta on the support:", np.round(solver.getBeta()[true_support], 3))
lam = solver.getLambda()                     # one value per step
k = len(lam)                                 # the other diagnostics start at the all-zero step
path_diag = np.column_stack([lam, solver.getRSS()[-k:], solver.getR2()[-k:],
                             solver.getCp()[-k:], solver.getDoF()[-k:]])
print("lambda, rss, r2, cp, dof per step:")
print(np.round(path_diag[:5], 3))
solver.restore(X, D, y)    # give X and D back before the next solver claims them
```

## The path as a picture

The coefficient paths of a terminated run show what the solver does. A predictor enters at a step and
its coefficient grows from there, and a dummy that enters marks one of the `T_stop` crossings that end
the path. One solver runs to `T_stop = 1`, keeps a snapshot of its path, and continues to `T_stop = 3`,
so the second panel extends the first.

```{code-cell} python
import matplotlib.pyplot as plt


def snapshot(solver):
    """The dense path and the variable that enters at each step (0-based)."""
    path = solver.getBetaPath().copy()                    # (p + L) x (steps + 1)
    entering = [abs(codes[0]) - 1 for codes in solver.getActions()]
    return path, entering


def draw_path(ax, path, entering, title):
    steps = np.arange(path.shape[1])
    for j in np.flatnonzero(np.any(path != 0, axis=1)):
        if j >= p:
            ax.plot(steps, path[j], color="C1", linestyle="--", linewidth=1.2)
        elif j in true_support:
            ax.plot(steps, path[j], color="C0", linewidth=1.8)
        else:
            ax.plot(steps, path[j], color="0.45", linewidth=1.0)
    for step, j in enumerate(entering, start=1):
        if j >= p:
            ax.axvline(step, color="C1", linestyle=":", linewidth=0.9)
    top = ax.secondary_xaxis("top")                       # the variable that enters at each step
    top.set_xticks(range(1, len(entering) + 1))
    top.set_xticklabels([f"D{j - p}" if j >= p else str(j) for j in entering], fontsize=7)
    top.tick_params(length=0)
    top.spines["top"].set_visible(False)
    ax.set_title(title)
    ax.set_xlabel("Step")


s = TLARS_Solver(X, D, y)
s.executeStep(T_stop=1)
first = snapshot(s)
s.executeStep(T_stop=3)
second = snapshot(s)
s.restore(X, D, y)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
draw_path(ax1, *first, "T_stop = 1")
draw_path(ax2, *second, "T_stop = 3, continued")
ax1.set_ylabel("Coefficient")
handles = [plt.Line2D([], [], color="C0", linewidth=1.8, label="true predictor"),
           plt.Line2D([], [], color="0.45", linewidth=1.0, label="other predictor"),
           plt.Line2D([], [], color="C1", linestyle="--", label="dummy"),
           plt.Line2D([], [], color="C1", linestyle=":", linewidth=0.9, label="dummy enters")]
fig.legend(handles=handles, loc="outside center right")
plt.show()
```

The true predictors enter first and keep growing, further predictors follow with small coefficients,
and the dummies end the path. The second panel repeats the first up to its last step and continues
from there, because the continuation walks the stored path on.

## The options

`SolverOptions` holds every input except the data, the same object for every solver. A solver reads
the fields it needs and ignores the rest. It is constructed with keywords, or its attributes are set
one by one.

| Field | Meaning |
| --- | --- |
| `normalize`, `intercept`, `scaling_mode` | center and scale the columns of `X` and `D` in place (`ScalingMode.L2` or `ZSCORE`) and center `y` |
| `precentered` | the caller centered the data: the active set is capped at $n - 1$ and the intercept counts in the degrees of freedom |
| `tol` | the numerical tolerance of the comparisons (machine epsilon by default) |
| `tie_seed` | the seed of the tie-break engine (0 by default, so exact ties resolve reproducibly) |
| `n_threads` | the thread budget of the solver's own passes (1 serial, at least 1, and a pass splits only above about a million entries of the design) |
| `exch_tie_alpha`, `exch_tie_floor` | the exchangeable-tie band of the greedy solvers |
| `lambda2` | the quadratic penalty of the elastic nets (0 by default) |
| `rho_afs`, `ncgmp_variant`, `tenet_aug_use_lars` | solver-specific switches |
| `verbose` | print the steps |

```{code-cell} python
from trexselector.tsolvers.ad.lars_based import TENET_Solver

en = TENET_Solver(X, D, y, SolverOptions(lambda2=0.5))
en.executeStep(T_stop=3)
print("active predictors:", en.getActivePredictorIndices())
en.restore(X, D, y)
```

## Writers and readers

With `normalize = True` or `intercept = True` (the defaults) a solver is a writer. It centers and
scales the columns of `X` and `D` in place, holds the statistics, and writes the original values
back when `restore()` is called or when the object is garbage collected. A second writer on a
claimed array is an error, because a second in-place transformation would move the data by rounding.
`y` is never written.

```{code-cell} python
from trexselector.tsolvers.ad.lars_based import TLASSO_Solver

X_before = X.copy(order="F")
w = TLASSO_Solver(X, D, y)
print("the scaled columns:  ", np.abs(X - X_before).max())
w.restore(X, D, y)
print("the original values: ", np.abs(X - X_before).max())
```

A pipeline that runs several solvers on one design scales the data once in front with a scaler and
constructs readers, solvers with `normalize = False` and `intercept = False` that never write.

```{code-cell} python
from trexselector.ml_methods import LpNormScaler, NormType

sx = LpNormScaler(NormType.L2, center=True)
sx.fit_transform_inplace(X)
sd = LpNormScaler(NormType.L2, center=True)
sd.fit_transform_inplace(D)
y_c = y - y.mean()
reader_opts = SolverOptions(normalize=False, intercept=False, precentered=True)
```

With the data pre-scaled, every family runs on the same design.

```{code-cell} python
from trexselector.tsolvers.ad import lars_based, omp_based, afs_based, cd_based

constructors = {
    "TLARS": lars_based.TLARS_Solver, "TLASSO": lars_based.TLASSO_Solver,
    "TSTEPWISE": lars_based.TSTEPWISE_Solver, "TSTAGEWISE": lars_based.TSTAGEWISE_Solver,
    "TENET": lars_based.TENET_Solver, "TENETAug": lars_based.TENETAug_Solver,
    "TOMP": omp_based.TOMP_Solver, "TOOLS": omp_based.TOOLS_Solver,
    "TGP": omp_based.TGP_Solver, "TACGP": omp_based.TACGP_Solver,
    "TMP": omp_based.TMP_Solver, "TNCGMP": omp_based.TNCGMP_Solver,
    "TAFS": afs_based.TAFS_Solver,
    "TCDLASSO": cd_based.TCDLASSO_Solver, "TCDENET": cd_based.TCDENET_Solver,
}
print(f"{'solver':<11}{'steps':>6}  predictors")
for name, make in constructors.items():
    s = make(X, D, y_c, reader_opts)
    s.executeStep(T_stop=3)
    print(f"{name:<11}{s.getNumSteps():>6}  {s.getActivePredictorIndices()}")
```

The informed solvers and the Tikhonov solvers take their structure as an extra argument.

```{code-cell} python
from trexselector.tsolvers.ad.lars_based import TIENET_Solver
from trexselector.tsolvers.ad.cd_based import TCDTIKNET_Solver

groups = np.repeat(np.arange(5, dtype=np.int32), 10)
structured_opts = SolverOptions(normalize=False, intercept=False, precentered=True, lambda2=0.5)
ien = TIENET_Solver(X, D, y_c, groups, structured_opts)
ien.executeStep(T_stop=3)
print("TIENET:   ", ien.getActivePredictorIndices())
tik = TCDTIKNET_Solver(X, D, y_c, np.eye(p), structured_opts)
tik.executeStep(T_stop=3)
print("TCDTIKNET:", tik.getActivePredictorIndices())
```

## Saving and loading

`save()` serializes the state of a solver to a file and `load()` restores it into a solver
constructed on the same data in the same state. The loaded solver is a reader and continues the path
where it was saved.

```{code-cell} python
import os
import shutil
import tempfile

work_dir = tempfile.mkdtemp()
f = os.path.join(work_dir, "tlars_state.bin")
a = TLARS_Solver(X, D, y_c, reader_opts)
a.executeStep(T_stop=2)
a.save(f)
b = TLARS_Solver(X, D, y_c, reader_opts)
b.load(f)
print("steps after load:", b.getNumSteps())
b.executeStep(T_stop=4)
c = TLARS_Solver(X, D, y_c, reader_opts)
c.executeStep(T_stop=4)
print("same actives:", b.getActives() == c.getActives())
```

## Memory

`getMemoryBytes()` reports the heap bytes a solver holds (its copy of `y`, the residual, the
correlations, the Cholesky factor, the coefficient path). `X` and `D` are borrowed and not counted.

```{code-cell} python
a.getMemoryBytes()
```

```{code-cell} python
:tags: [remove-cell]
sx.inverse_transform_inplace(X)
sd.inverse_transform_inplace(D)
shutil.rmtree(work_dir)
```

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of every T-solver and five performance
studies, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/tsolvers/ad),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/tsolvers/ad) and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/tsolvers/ad).