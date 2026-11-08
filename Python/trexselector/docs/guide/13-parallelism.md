---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# Parallelism: who owns the cores

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Two layers of work

A T-Rex run has two layers that could use several cores. The outer layer is the $K$ random
experiments, which are independent of each other and share the design. The inner layer is the work
inside one T-solver: the passes over the design that refresh the correlations, the search for the
next knot of the path, the standardization. The rule of the package is that one layer owns the cores
at a time. Either the experiments run in parallel and every solver runs serially, or the experiments
run in sequence and a solver may split its passes. Nested teams of threads are never formed.

The default is serial. A selector runs its experiments one after the other, every solver runs on one
thread, and a T-solver used on its own runs serially. Parallel work is opt-in through two attributes
of `TRexControlParameter`: `max_outer_threads` for the experiments and `max_inner_threads` for the
solvers' passes.

## The default thread budget of the package

The package is single-threaded by default on every machine. Every parallel region of the backend
runs on an explicit budget, and the methods without an option of their own (the scalers, the
clustering, the ML regressions, `normalize_design()`, the dummy and data generators) take the
package's default budget, which starts at 1. `get_num_threads()` reads that default,
`set_num_threads()` sets it for the process and returns the previous value, and the `n_threads`
argument of those methods overrides it per object. Inside a selector the same kernels take
`max_inner_threads`. `get_max_threads()` reports the threads OpenMP can run on this machine, an
upper bound for a budget. In a build without OpenMP every setting is a no-op and the machine count
is one. The hidden first cell of this page set the default to two, so the first value below is 2.

```{code-cell} python
from trexselector.utils import get_max_threads, get_num_threads, set_num_threads

print(get_num_threads())
print(get_max_threads())
```

The selectors and the solvers do not read the default budget. Their budget comes from the control.

## The solvers' passes

`TRexControlParameter.max_inner_threads = m` gives every solver inside the random experiments a
budget of `m` threads for its passes over the design (the correlation refreshes, the step-length
search, the standardization). The default is 1. A pass is split only when it touches at least a
million entries of the design, because the cost of forming a team of threads exceeds the gain below
that, so small problems run serially whatever the budget says. The seeded generators (the dummies,
the synthetic designs) count their work in the same unit with the measured cost of a random draw
(about 80 entries) and of seeding an engine (up to 32,000), so they split from a few thousand
entries. The budget is live when the experiments run in sequence. Inside parallel experiments the
solvers stay serial. A T-solver used on its own takes the same budget through
`SolverOptions(n_threads=m)`, and `setNumThreads()` changes it on a live solver.

## Parallel random experiments

`TRexControlParameter.parallel_rnd_experiments = True` with `max_outer_threads = m` runs the random
experiments on `m` threads (at most $K$). Inside the parallel experiments every solver runs
serially, so the number of threads in use is `max_outer_threads` and nothing more. Both options are
whole numbers of at least 1. The lean L-loop strategies (`SEEDED`, `PERMUTATION`,
`PERMUTATION_SEEDED`) run their experiments in sequence in every case and give the cores to the
solvers.

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(12)
n, p = 100, 150
beta = np.r_[np.ones(5), np.zeros(p - 5)]
X = np.asfortranarray(rng.standard_normal((n, p)))
y = X @ beta + rng.standard_normal(n)
```

```{code-cell} python
def control(**settings):
    ctrl = tsn.TRexControlParameter()
    ctrl.K = 10
    for name, value in settings.items():
        setattr(ctrl, name, value)
    return ctrl

serial = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=control(), seed=42, verbose=False)
serial.select()
ctrl_outer = control(parallel_rnd_experiments=True, max_outer_threads=2)
parallel = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl_outer, seed=42, verbose=False)
parallel.select()
ctrl_inner = control(max_inner_threads=2)
budget = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl_inner, seed=42, verbose=False)
budget.select()
print(list(serial.selected_indices))
print(list(serial.selected_indices) == list(parallel.selected_indices))
print(list(serial.selected_indices) == list(budget.selected_indices))
```

## Determinism

A seeded selector derives one seed per random experiment from `seed` and the experiment index, so
the dummies of experiment $k$ do not depend on which thread draws them or on how many threads there
are. The selection of a seeded run is therefore the same with experiments in sequence and in
parallel, and the same for any thread count of the inner passes. Every reduction inside the solvers
is summed in a fixed order for the same reason.

## Combining with your own parallel loops

A Monte Carlo study usually parallelizes over trials. Then the trials own the cores, and the package
should run serially inside every worker, which is the default. Two details matter.

- With the `fork` start method of `multiprocessing` the workers are forks of the process and inherit
  its default budget. Keep it at one in the process, or set it to one in the workers with
  `set_num_threads(1)`, so that the forks do not oversubscribe the machine.
- With the `spawn` start method, `concurrent.futures.ProcessPoolExecutor` on macOS and Windows, or
  joblib's loky backend, the workers are fresh interpreters with the default budget of a fresh
  process, one thread. Import the package in the workers and keep the defaults.

```python
from concurrent.futures import ProcessPoolExecutor

from trexselector.utils import compute_fdp, set_num_threads


def trial(mc):
    set_num_threads(1)
    rng = np.random.default_rng(mc)
    X = np.asfortranarray(rng.standard_normal((n, p)))
    y = X @ beta + rng.standard_normal(n)
    sel = tsn.TRexSelector(X, y, tFDR=0.1, verbose=False)   # seed = -1 per trial
    sel.select()
    return compute_fdp(sel.selected_indices, np.flatnonzero(beta).tolist())


with ProcessPoolExecutor(max_workers=4) as pool:
    fdps = list(pool.map(trial, range(40)))
print(np.mean(fdps))
```

## Builds without OpenMP

Every wheel on PyPI bundles the OpenMP runtime of its platform, so an installed wheel threads out of
the box. A build from source takes OpenMP from the toolchain, and the CMake switch
`TREX_WITH_OPENMP=OFF` forces the single-threaded build where OpenMP exists. A single-threaded build
selects the same variables, it only runs the experiments and the passes in sequence.

## The CRAN rule

CRAN checks the R package on at most two cores. The pages of this guide and the tests of the Python
package follow the same rule with the idiom of their hidden first cell, which caps the default
budget for the duration of the page.

```python
old_threads = set_num_threads(min(2, get_max_threads()))
# ... the work ...
set_num_threads(old_threads)
```
