---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# Working without data copies

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## The rule

The backend works on the caller's data wherever possible. A NumPy array is bound by its memory
address, a memory-mapped file by its mapping, and the preprocessing a method needs is applied in
place and undone afterwards. A design of several gigabytes is therefore held once, by NumPy or by
the file, and not a second time by the backend. This page states the contract that makes it safe.
The bindings hand the memory of a NumPy array {cite}`harris2020numpy` to the C++ core through
pybind11 {cite}`jakob2017pybind11`, as an Eigen map {cite}`guennebaud2010eigen` of the same buffer,
and a memory-mapped file through Boost.Interprocess {cite}`gaztanaga2008interprocess`.

## Pass double matrices

The backend binds double storage in column-major order, the layout Eigen reads without a copy. A
selector converts its input with `np.asfortranarray(X, dtype=np.float64)`, so an integer array, a
single-precision array or a C-ordered array is copied first, once, and a Fortran-ordered array of
doubles passes through untouched. The selector's `X` attribute holds the array the backend works on,
and `np.shares_memory()` tells whether it is the caller's. To keep a large design copy-free, create
it as `np.float64` in Fortran order.

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(1)
Xi = rng.integers(0, 3, size=(60, 20))
print(Xi.dtype)
sel_i = tsn.TRexSelector(Xi, rng.standard_normal(60), tFDR=0.1, seed=2, verbose=False)
print(np.shares_memory(sel_i.X, Xi))      # the integers were converted, a copy
sel_i.restore()
Xf = np.asfortranarray(Xi, dtype=np.float64)   # one explicit conversion, zero-copy from here on
sel_f = tsn.TRexSelector(Xf, rng.standard_normal(60), tFDR=0.1, seed=2, verbose=False)
print(np.shares_memory(sel_f.X, Xf))
sel_f.restore()
```

## Selectors: in place at construction, back at the end

A selector centers and scales `X` in place when it is constructed and writes the original values
back when `select()` returns. In between, `X` holds the working design. `restore()` writes back
early for a selector that will not run, and a selector that is garbage collected without running
writes back in its finalizer.

```{code-cell} python
n, p = 80, 40
X = np.asfortranarray(rng.standard_normal((n, p)))
y = X[:, 0] * 2 + rng.standard_normal(n)
X_before = X.copy()
ctrl = tsn.TRexControlParameter()
ctrl.K = 10
sel = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=2, verbose=False)
print(np.round(np.sqrt((X ** 2).sum(axis=0))[:3], 6))   # the working design: unit columns
sel.select()
print(np.max(np.abs(X - X_before)))                      # the original values are back
```

## T-solvers: writers and readers

A T-solver with `normalize = True` or `intercept = True` is a writer. It transforms `X` and `D` in
place, keeps the statistics, and holds a claim on both matrices until `restore()` or garbage
collection writes them back. A second writer on a claimed matrix is refused, because a second
transformation would move the data by rounding. The `restore()` of a T-solver takes the three arrays
it was constructed on.

```{code-cell} python
from trexselector.tsolvers.ad import SolverOptions
from trexselector.tsolvers.ad.lars_based import TLARS_Solver, TLASSO_Solver

D = np.asfortranarray(rng.standard_normal((n, p)))
w1 = TLARS_Solver(X, D, y)
try:
    TLASSO_Solver(X, D, y)
except RuntimeError as e:
    print("error:", str(e)[:90], "...")
w1.restore(X, D, y)
```

A solver with `normalize = False` and `intercept = False` is a reader. It runs on the data as it
finds them and never writes, so any number of readers can share one design. A pipeline scales the
design once with a scaler and constructs readers.

```{code-cell} python
from trexselector.ml_methods import LpNormScaler, NormType

sx = LpNormScaler(norm_type=NormType.L2, center=True)
sx.fit_transform_inplace(X)
sd = LpNormScaler(norm_type=NormType.L2, center=True)
sd.fit_transform_inplace(D)
y_c = y - y.mean()
opts = SolverOptions(normalize=False, intercept=False, precentered=True)
r1 = TLARS_Solver(X, D, y_c, opts)
r2 = TLASSO_Solver(X, D, y_c, opts)
r1.executeStep(T_stop=2)
r2.executeStep(T_stop=2)
sx.inverse_transform_inplace(X)
sd.inverse_transform_inplace(D)
print(np.max(np.abs(X - X_before)))
```

`y` is never written. The solver centers an internal copy of the response.

## Loading a saved solver

`load()` restores a solver into one constructed on the same data in the same state. Load into a
reader on pre-scaled data, or into a writer constructed with the same options on the same raw data,
which then keeps the duty to write back.

## Mapped files are never read into memory

A memory-mapped matrix is bound through the same interface as a NumPy array. The backend touches the
pages it needs and the file keeps the data, so the working set of a run is the backend's own
buffers, never a copy of the design. The accessor `to_numpy()` of the handle is a view of the
mapping as well, and element access reads from the file, so no accessor copies the design either.

## The low-memory engines

The ridge and elastic-net classes have a `low_memory` engine that avoids every allocation of the
size of `X`.

| Class | `low_memory = False` | `low_memory = True` |
| --- | --- | --- |
| `RidgeCV` | each fold gathers and standardizes its training rows and takes a thin singular value decomposition (SVD), about three to four copies of `X` in flight | the dual (kernel) form, $n \times n$ kernels accumulated from column blocks read in place, memory $O(n^2)$ |
| `ElasticNet` | coordinate descent on a centered and standardized working copy, one $n \times p$ matrix | coordinate descent on `X` itself, no allocation of the size of `X` |
| `ElasticNetCV` | one working copy per fit | every fold on `X` itself with a training mask |

The two engines compute the same estimator on the same grid and agree to rounding. The selectors
that tune $\lambda_2$ by cross-validation (GVS, SPCA, Tikhonov) use these engines, and the elastic-net
tuner of the sparse principal component selector keeps its peak memory under 0.3 times the size of `X`.

```{code-cell} python
from trexselector.ml_methods import RidgeCV

cv_a = RidgeCV()
cv_a.fit(X, y, n_folds=5, n_lambda=30)
cv_b = RidgeCV(low_memory=True)
cv_b.fit(X, y, n_folds=5, n_lambda=30)
print(np.allclose(cv_a.cv_mse(), cv_b.cv_mse()))
```

## Measuring the footprint

`getMemoryBytes()` reports what a T-solver holds on the heap. `X` and `D` are borrowed and not
counted.

```{code-cell} python
r1.getMemoryBytes()
```
