---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# Getting started with trexselector

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## What the package does

The Terminating-Random Experiments (T-Rex) selector selects variables while controlling the false
discovery rate (FDR) at a target level chosen by the user, and it does so with a high true positive rate
(TPR) {cite}`machkour2025trex`. It applies to low- and high-dimensional designs and is tailored to the
high-dimensional problem. It takes a design matrix $\boldsymbol{X}$ with $n$ observations of $p$
variables, a response vector $\boldsymbol{y}$ and the target FDR, and runs $K$ random experiments in
which the real variables compete with computer-generated dummy variables, random columns with finite
mean and variance that act as flagged null variables. In every experiment the dummies are appended to
the design, a forward selection solver walks its path and stops as soon as $T$ dummies have entered, and
the share of experiments in which a variable was selected is its vote. Because a dummy enters as easily
as a null variable, the included dummies give an estimate of the false discoveries, and a calibration
over the number of dummies, the stopping time and the voting level selects as many variables as the FDR
target allows, with a finite-sample guarantee that holds for a broad class of forward selection solvers.
The scheme is not tied to the linear model and scales to millions of variables. The library currently
covers the supervised linear model with a response vector and a design matrix.

The Python package is the Python interface to the TRexSelector C++ library. The same library drives
the R package and the C++ software development kit, so the three surfaces select the same variables
on the same data. Since version 1.0.0 it replaces the pure-Python implementation of the 0.6.x releases
under the same package name [[version 0.6.18](https://pypi.org/project/trexselector/0.6.18/)].

The package currently has five groups of tools.

- **T-Rex selectors.** The core `TRexSelector` and its variants for dependent variables
  (`TRexDASelector`), grouped variables (`TRexGVSSelector`), fast screening (`TRexScreeningSelector`,
  `TRexBiobankScreeningSelector`), sparse principal component analysis (PCA) (`TRexSPCASelector`) and
  a general Tikhonov penalty (`TRexTikhonovSelector`).
- **T-solvers.** The twenty terminating forward selection solvers the selectors run inside their
  random experiments, available as classes of their own under `trexselector.tsolvers.ad`.
- **Machine learning (ML) methods.** Scalers, ridge and elastic-net regression with cross-validation,
  singular value decomposition (SVD), PCA and hierarchical clustering under `trexselector.ml_methods`.
- **Memory-mapped matrices.** Disk-backed designs are accepted everywhere an array is accepted. The
  algorithms are agnostic to the data representation, so the library processes designs larger than
  the random-access memory (RAM) of the machine out of core.
- **Utilities.** The selection metrics (false discovery proportion, true positive proportion,
  precision and recall), the normalizer, and the threading functions under `trexselector.utils`.

## Installation

The released package is on the Python Package Index (PyPI)
[[package page](https://pypi.org/project/trexselector/)]. Install it into the active environment:

```bash
pip install trexselector
```

Binary wheels exist for Linux (x86_64), macOS (arm64) and Windows (x86_64) for Python 3.10 to 3.15,
every wheel with the OpenMP runtime of its platform, so the package runs on several threads out of
the box. The README of the package describes the build from source.

The package builds on NumPy {cite}`harris2020numpy` for its arrays and on pybind11
{cite}`jakob2017pybind11` for the bindings.

## A first run

A small Gaussian data set with $n = 50$ observations, $p = 100$ variables and three active variables,
drawn with a fixed seed:

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(123)
n, p, num_act = 50, 100, 3
beta = np.zeros(p)
beta[:num_act] = 1.0
X = rng.standard_normal((n, p))
y = X @ beta + rng.standard_normal(n)
true_support = np.flatnonzero(beta).tolist()
true_support
```

The selectors are classes. Construct one, run it with `select()` and read the results from its
attributes. A fixed `seed` makes the dummy draws reproducible. The default `seed = -1` draws them from
the hardware entropy source, which is the right choice for Monte Carlo studies, see the T-Rex selector
page.

```{code-cell} python
sel = tsn.TRexSelector(X, y, tFDR=0.1, seed=1234, verbose=False)
sel.select()
sel.selected_indices
```

The two metrics of variable selection compare the selection with the true support. The false
discovery proportion (FDP) is the share of false positives among the selected variables and the true
positive proportion (TPP) is the share of the true support that was found. Their expectations are
the FDR and the TPR.

```{code-cell} python
from trexselector.utils import compute_fdp, compute_tpp

print("FDP:", compute_fdp(sel.selected_indices, true_support))
print("TPP:", compute_tpp(sel.selected_indices, true_support))
```

## The result fields

Every selector exposes the calibration it arrived at.

```{code-cell} python
print("selected:", sel.num_selected)   # number of selected variables
print("T_stop:  ", sel.T_stop)         # stopping time: dummies included when the T-loop stopped
print("L:       ", sel.L)              # calibrated dummy multiplier, the selector used L * p dummies
print("v_thresh:", sel.v_thresh)       # voting level the selection was read at
print("K:       ", sel.K)              # random experiments
```

`fdp_hat_mat` holds the FDP estimate per stopping time (rows) and voting level (columns), `phi_mat` the
relative occurrence of every variable over the $K$ experiments after each stopping time, and
`voting_grid` the voting levels. The T-Rex selector page walks through them.

```{code-cell} python
print(sel.fdp_hat_mat.shape)
print(sel.voting_grid[:6])
print(np.round(sel.fdp_hat_mat[sel.T_stop - 1], 3))
```

## T-Rex variants of the library

| Situation | Selector | Page |
| --- | --- | --- |
| Independent or modestly dependent variables | `TRexSelector` | T-Rex selector |
| Strongly dependent variables (autoregressive, equi-correlated, clustered) | `TRexDASelector` | T-Rex+DA |
| Groups of correlated variables that should enter together | `TRexGVSSelector` | T-Rex GVS |
| Very large $p$, a fast first pass before the full selector | `TRexScreeningSelector` | Screen-T-Rex |
| Many phenotypes against one genotype matrix | `TRexBiobankScreeningSelector` | Screen-T-Rex |
| Sparse principal components | `TRexSPCASelector` | T-Rex SPCA |
| A structured penalty (smoothness, graph Laplacian, groups) | `TRexTikhonovSelector` | T-Rex Tikhonov |

The solver inside the random experiments is set in `TRexControlParameter.solver_type`. T-LARS
(terminating least angle regression) is the default and the solver of the T-Rex paper. T-LASSO adds
the lasso modification, T-ENET the elastic net, the orthogonal matching pursuit (OMP) family and the
adaptive forward stepwise family are greedy and faster per step, and the coordinate descent (CD)
family computes exact penalized minimizers at every dummy crossing. The T-solvers page describes the
families, and every solver can be run on its own.

```{code-cell} python
ctrl = tsn.TRexControlParameter()
ctrl.solver_type = tsn.SolverTypeForTRex.TLASSO
ctrl.K = 10
sel2 = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=1234, verbose=False)
sel2.select()
sel2.selected_indices
```

## Citing the package

The citation metadata of the repository is in its `CITATION.cff`. The methods are published in
{cite}`machkour2025trex` (T-Rex), {cite}`machkour2025datrex` (T-Rex+DA),
{cite}`machkour2022gvs,machkour2023ien` (T-Rex GVS), {cite}`machkour2023screen` (Screen-T-Rex),
{cite}`machkour2024spca` (T-Rex SPCA), {cite}`scheidt2023bigtrex` (memory mapping) and
{cite}`scheidt2026greedy` (the terminating greedy pursuits).
