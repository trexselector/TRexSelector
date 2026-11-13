---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# Memory mapping: disk-backed designs

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## A design larger than the memory

A memory-mapped matrix lives in a file on disk, and the operating system pages the parts in use into
memory. The C++ backend binds such a file exactly as it binds a NumPy array, so every selector,
every T-solver, the scalers, the low-memory regression engines and the clustering accept one. That
is how a sparse regression problem with five million variables was solved on a laptop
{cite}`scheidt2023bigtrex`. The file holds the matrix as raw doubles in column-major order, the
layout of Fortran-ordered NumPy arrays and of Eigen, without a header.

The file is mapped through Boost.Interprocess {cite}`gaztanaga2008interprocess`, the library of the
Boost C++ Libraries for memory-mapped files and shared memory. The library binds any buffer, mapped
or not, so the selectors never see where the design lives, as Figure 1 shows.

```{figure} figures/memory_agnostic.png
:width: 90%
:align: center

Figure 1: A design in memory and a design in a file reach the same implementation through one mapped
view.
```

## Creating and opening a mapped matrix

`numpy_to_memmap()` writes an array to a file and returns a handle, a `MemoryMappedMatrix`, and the
constructor of that class opens an existing file with its dimensions and the access mode.
`convert_to_memory_mapped()` writes the file without returning a handle. The handle's `to_numpy()`
is a zero-copy view of the mapping, a Fortran-ordered array of doubles that reads and writes the
file.

```{code-cell} python
import os
import tempfile

import numpy as np
import trexselector as tsn
from trexselector.utils import AccessMode, MemoryMappedMatrix, numpy_to_memmap

rng = np.random.default_rng(5)
n, p = 100, 200
X = np.asfortranarray(rng.standard_normal((n, p)))
workdir = tempfile.mkdtemp()
path = os.path.join(workdir, "design.bin")
mm = numpy_to_memmap(path, X)
print(type(mm).__name__, (mm.rows(), mm.cols()))
print(os.path.getsize(path) == 8 * n * p)
```

The handle supports the element access one needs to inspect and to edit the file, and the view gives
the block access of NumPy. Reading through the view touches the pages of the file, the backend never
copies the design.

```{code-cell} python
view = mm.to_numpy()
print(np.round(view[:2, :3], 3))
print(np.array_equal(view[:2, :3], X[:2, :3]))
mm[0, 0] = 42.0
print(mm[0, 0])
mm[0, 0] = X[0, 0]
print(np.max(np.abs(view - X)))
ro = MemoryMappedMatrix(path, n, p, AccessMode.ReadOnly)
try:
    ro[0, 0] = 0.0
except RuntimeError as e:
    print("error:", e)
```

## Streaming a design into the file

A design that does not fit in memory is written column by column. The file is raw doubles in
column-major order, so writing one column after another with `tofile()` produces it.

```{code-cell} python
path2 = os.path.join(workdir, "streamed.bin")
with open(path2, "wb") as con:
    for j in range(p):
        np.random.default_rng(1000 + j).standard_normal(n).tofile(con)
big = MemoryMappedMatrix(path2, n, p, AccessMode.ReadWrite)
print(np.max(np.abs(big.to_numpy()[:, 2] - np.random.default_rng(1000 + 2).standard_normal(n))))
```

## Selectors on a mapped design

The selectors take the view of the handle where they take `X`. A selector centers and scales the
design in place and restores it when `select()` returns, so the file must be open in `ReadWrite`
mode (the page on [working without data copies](12-no-copies.md) has the contract).

```{code-cell} python
beta = np.r_[np.full(4, 1.5), np.zeros(p - 4)]
y = X @ beta + rng.standard_normal(n)
ctrl = tsn.TRexControlParameter()
ctrl.K = 10
sel_mem = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=3, verbose=False)
sel_mem.select()
sel_map = tsn.TRexSelector(mm.to_numpy(), y, tFDR=0.1, trex_control=ctrl, seed=3, verbose=False)
sel_map.select()
print(list(sel_mem.selected_indices) == list(sel_map.selected_indices))
print(np.max(np.abs(mm.to_numpy() - X)))      # the file is restored up to rounding
```

The same holds for `TRexDASelector`, `TRexGVSSelector`, `TRexScreeningSelector`,
`TRexBiobankScreeningSelector`, `TRexSPCASelector` and `TRexTikhonovSelector`.
`TRexControlParameter.use_memory_mapping = True` additionally keeps the dummy matrices of the random
experiments in mapped files, for the case where $K \cdot L \cdot p$ dummy columns do not fit either.

One file serves one selector at a time. Do not run several selectors on the same file at the same
time under `ScalingMode.L2` or `ScalingMode.ZSCORE`, neither through a second handle in one process
nor in several processes. Each of them rewrites the file while it runs, so it would standardize data
that another one has already changed, and the selections come out wrong. The same holds for
T-solvers with `normalize=True` or `intercept=True`. The next section shows how several processes
use one design at once.

### A read-only file

The default mode `ReadWrite` serves every use and is the one the examples of this guide work with.
The view of a file opened in `ReadOnly` mode is an array that is not writable, and it serves where
nothing is written into it. For the selectors that is `ScalingMode.NONE`, under which a selector
takes the design as supplied. The file then holds the design the selector works on, centered and
scaled beforehand, and `normalize_design()` applies exactly the preprocessing of the selector. Under
`ScalingMode.L2` and `ScalingMode.ZSCORE` the view is refused with an error that names this way out.
The rule reads the write flag of the array, so it holds for every array that is not writable.
Several processes can share one read-only file, because none of them writes into it.

```{code-cell} python
import textwrap

pre = tsn.normalize_design(X, y)
path_pre = os.path.join(workdir, "design_pre.bin")
numpy_to_memmap(path_pre, pre["X"])
pre_view = MemoryMappedMatrix(path_pre, n, p, AccessMode.ReadOnly).to_numpy()
print(pre_view.flags.writeable)
none = tsn.TRexControlParameter()
none.K = 10
none.scaling_mode = tsn.ScalingMode.NONE
sel_ro = tsn.TRexSelector(pre_view, pre["y"], tFDR=0.1, trex_control=none, seed=3, verbose=False)
sel_ro.select()
print(list(sel_ro.selected_indices) == list(sel_mem.selected_indices))
try:
    tsn.TRexSelector(pre_view, pre["y"], trex_control=ctrl, verbose=False)
except ValueError as e:
    print(textwrap.fill(f"error: {e}", 90))
```

`TRexSPCASelector` centers the design in every run, so it needs a writable view under every scaling
mode.

## T-solvers and ML methods

A T-solver binds a memory-mapped `X` and a memory-mapped `D` as it binds an array in memory, the
scalers transform the file in place, the low-memory engines of the regressions read it in place, and
the clustering accepts it as data. The access mode follows the same rule as for the selectors. A
T-solver with `normalize=True` or `intercept=True` and the scalers write into the file and need
`ReadWrite` mode, while a T-solver with both options off and the fits of the regressions only read
it and take the view of a `ReadOnly` file, for `X` and for `D`.

```{code-cell} python
from trexselector.tsolvers.ad.lars_based import TLARS_Solver
from trexselector.ml_methods import ElasticNet, ZScoreScaler

D = np.asfortranarray(rng.standard_normal((n, p)))
s = TLARS_Solver(view, D, y)
s.executeStep(T_stop=2)
print(s.getActivePredictorIndices())
s.restore(view, D, y)
z = ZScoreScaler()
z.fit_transform_inplace(view)
print(round(float(view[:, 0].mean()), 12))
z.inverse_transform_inplace(view)
en = ElasticNet(low_memory=True)
en.fit(view, y, alpha=1.0, n_lambda=20)
print(en.coef().shape)
```

## When it pays

A mapped design pays when the design does not fit into memory, when several solvers share one design
and the working set should be loaded once, and when a design is produced by a streaming process and
read afterwards. For a design that fits comfortably in memory a NumPy array is faster, since no
paging happens. Mapped files are never read into the process heap by the backend, so the memory
footprint of a run is the backend's working data, see `getMemoryBytes()` on the T-solvers.

```{code-cell} python
:tags: [remove-cell]
import shutil
del view, mm, ro, big, s, sel_map, pre_view, sel_ro
shutil.rmtree(workdir)
```

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of memory-mapped designs, with their results
and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/memory_mapping/demo_memory_mapping),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/memory_mapping/demo_memory_mapping)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/memory_mapping/demo_memory_mapping).