---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# T-Rex SPCA: sparse principal components under FDR control

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Sparse loadings with a guarantee

Principal component analysis (PCA) finds the directions of largest variance, and every variable
contributes to every loading vector. Sparse PCA asks for loadings with few non-zero entries, so that
a component can be read as a small set of variables. The T-Rex sparse PCA (T-Rex SPCA)
{cite}`machkour2024spca` decides which variables enter a loading vector with the T-Rex selector. For
each of the $M$ components, ordinary PCA gives the component as a response, the T-Rex selector
selects the variables for it at the target false discovery rate (FDR), and a sparse loading vector
is assembled from the selected set. The selected set of every component is FDR-controlled.

## Data with two sparse factors

Two latent factors, the first loading on variables 0 to 9 and the second on variables 30 to 39
(0-based), plus noise.

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(77)
n, p = 200, 60
Z_true = rng.standard_normal((n, 2))
V_true = np.zeros((p, 2))
V_true[0:10, 0] = 1 / np.sqrt(10)
V_true[30:40, 1] = 1 / np.sqrt(10)
X = (6 * np.outer(Z_true[:, 0], V_true[:, 0]) + 3 * np.outer(Z_true[:, 1], V_true[:, 1])
     + rng.standard_normal((n, p)))
```

## Running the selector

`TRexSPCASelector(X, spca_ctrl, seed)` takes the design and the control, and `select(M, tFDR)`
extracts `M` components at the target FDR. The elastic-net penalty of the inner selector is chosen by
cross-validation by default. `X` is centered in place during `select()` and restored on return. The
base control of the inner selector, with the number of random experiments `K`, is nested in the
group control of the `TRexSPCAControlParameter`.

```{code-cell} python
spca_ctrl = tsn.TRexSPCAControlParameter()
spca_ctrl.gvs_ctrl.trex_ctrl.K = 10

spca = tsn.TRexSPCASelector(X, spca_ctrl=spca_ctrl, seed=9)
spca.select(M=2, tFDR=0.1)
spca.active_sets
```

The results are the scores `Z` ($n \times M$), the sparse loadings `V` ($p \times M$), the selected
sets, and the explained variance.

```{code-cell} python
print("Z:", spca.Z.shape)
print("V:", spca.V.shape)
print("adjusted_ev:  ", np.round(spca.adjusted_ev, 3))
print("cumulative_ev:", np.round(spca.cumulative_ev, 3))
```

```{code-cell} python
:tags: [hide-input]
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
variable = np.arange(spca.V.shape[0])
supports = [range(0, 10), range(30, 40)]
for m, ax in enumerate(axes):
    on = np.isin(variable, list(supports[m]))
    ax.vlines(variable[on], 0, spca.V[on, m], color="C0", linewidth=2, label="true support")
    ax.vlines(variable[~on], 0, spca.V[~on, m], color="C1", linewidth=2, label="other")
    ax.axhline(0, color="0.3", linewidth=0.5)
    ax.set_xlabel("Variable")
    ax.set_ylabel(f"Loading {m + 1}")
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside center right")
fig.suptitle("The two sparse loading vectors. The true supports are variables 0 to 9 and 30 to 39.",
             fontsize=9)
```

Sparse components are not orthogonal in general, so the explained variance is adjusted for the
non-orthogonality by a QR decomposition of the scores. `adjusted_ev` is the marginal share per
component and `cumulative_ev` the cumulative percentage.

## The control

`TRexSPCAControlParameter` sets the inner machinery.

- `mode`: `SPCAMode.ActiveSet` (default) assembles the loading by a ridge regression of the component
  on the selected variables with penalty `lambda2_ridge_loadings`, `SPCAMode.Thresholded` keeps the
  top $|A|$ ordinary loadings and normalizes them.
- `gvs_ctrl.lambda_2` and `gvs_ctrl.lambda2_method`: the elastic-net penalty of the inner group
  selector, fixed or chosen by ridge cross-validation (`CV_1SE_CCD` by default when `lambda_2 < 0`).
  A value of 0 is the degenerate lasso case and gives a warning.
- `en_solver`: `ENSolverType.TENET_AUG` (default) or `TENET`, two forms of the same elastic net.
- `gvs_ctrl.trex_ctrl.K` and `gvs_ctrl.trex_ctrl.max_dummy_multiplier`: the random experiments and
  the dummy cap of the inner selector.

T-Rex SPCA is defined on the elastic net, and `gvs_type = GVSType.IEN` is rejected.

```{code-cell} python
thr_ctrl = tsn.TRexSPCAControlParameter()
thr_ctrl.mode = tsn.SPCAMode.Thresholded
thr_ctrl.gvs_ctrl.trex_ctrl.K = 10

thr = tsn.TRexSPCASelector(X, spca_ctrl=thr_ctrl, seed=9)
thr.select(M=2, tFDR=0.1)
print("active sets:  ", thr.active_sets)
print("cumulative_ev:", np.round(thr.cumulative_ev, 3))
```

## Memory-mapped data

`X` may be the `to_numpy()` view of a `MemoryMappedMatrix` opened in `AccessMode.ReadWrite` (see
[the memory-mapping page](11-memory-mapping.md)). It is centered in the file during `select()` and
restored. `convert_to_memory_mapped()` writes an array to such a file.

```{code-cell} python
import os
import shutil
import tempfile
from trexselector.utils import AccessMode, MemoryMappedMatrix, convert_to_memory_mapped

workdir = tempfile.mkdtemp()
path = os.path.join(workdir, "X.bin")
convert_to_memory_mapped(X, path)

X_map = MemoryMappedMatrix(path, n, p, AccessMode.ReadWrite)
mapped = tsn.TRexSPCASelector(X_map.to_numpy(), spca_ctrl=spca_ctrl, seed=9)
mapped.select(M=2, tFDR=0.1)
print("same active sets:", mapped.active_sets == spca.active_sets)
print("file restored:   ", np.allclose(X_map.to_numpy(), X))

del mapped, X_map
shutil.rmtree(workdir)
```

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of T-Rex SPCA, with their results and a
README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/trex_selector_methods/ad/trex_spca),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/trex_selector_methods/ad/trex_spca)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/trex_selector_methods/ad/trex_spca).