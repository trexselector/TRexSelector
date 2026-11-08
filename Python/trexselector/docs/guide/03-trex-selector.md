---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# The T-Rex selector

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## FDR and TPR

Let $\widehat{\mathcal{A}} \subseteq \{1, \ldots, p\}$ be the set of selected variables and
$\mathcal{A}$ the set of true active variables. With $1 \vee a = \max\{1, a\}$, the false discovery
rate (FDR) and the true positive rate (TPR) are

$$
\operatorname{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{A}} \setminus \mathcal{A}|}
{1 \vee |\widehat{\mathcal{A}}|}\right],
\qquad
\operatorname{TPR} = \mathbb{E}\left[\frac{|\mathcal{A} \cap \widehat{\mathcal{A}}|}
{1 \vee |\mathcal{A}|}\right].
$$

The quantities inside the expectations are the false discovery proportion (FDP) and the true
positive proportion (TPP) of one selection. `compute_fdp()` and `compute_tpp()` of
`trexselector.utils` compute them from index sets, `compute_fdp_dense()` and `compute_tpp_dense()`
from coefficient vectors. The T-Rex selector controls the FDR at the target `tFDR` and maximizes the
number of selected variables under that constraint {cite}`machkour2025trex`.

## The framework

```{figure} figures/T-Rex_framework.png
:width: 85%
:align: center

Figure 1: The T-Rex framework.
```

The selector runs $K$ random experiments. In experiment $k$ the design $\boldsymbol{X}$ ($n \times
p$) is enlarged by a dummy matrix with $L \cdot p$ columns drawn independently of the data, a
forward selection solver walks the solution path of the enlarged design, and the path is terminated
as soon as $T$ dummies have entered the model. The dummies are a negative control group. A null
variable and a dummy are exchangeable along the path, so the number of included dummies estimates
the number of included null variables.

```{figure} figures/predictor_matrix_with_dummies.png
:width: 65%
:align: center

Figure 2: The enlarged design (predictors and dummies).
```

The candidate sets of the $K$ experiments are fused by a vote. A variable is selected at voting
level $v$ when it appears in more than $v \cdot K$ of the candidate sets. The calibration searches
the voting level and the stopping time $T$ for the largest selection whose FDP estimate stays below
the target. The number of dummies is calibrated too. The L-loop starts at $L = 1$ and raises $L$
until the FDP estimate at the reference voting level falls below the target or the cap
`max_dummy_multiplier` is reached.

## Running the selector

The data are drawn with a fixed seed, $n = 100$ observations, $p = 150$ variables and five active
variables. The design is created in Fortran order, the layout the core reads without a copy (see the
section on the design below). Indices are 0-based.

```{code-cell} python
import numpy as np
import trexselector as tsn
from trexselector.utils import compute_fdp, compute_tpp

rng = np.random.default_rng(123)
n, p = 100, 150
beta = np.zeros(p)
beta[:5] = 1.0
true_support = np.flatnonzero(beta).tolist()
X = np.asfortranarray(rng.standard_normal((n, p)))
y = X @ beta + rng.standard_normal(n)
```

```{code-cell} python
sel = tsn.TRexSelector(X, y, tFDR=0.1, seed=42, verbose=False)
sel.select()
sel.selected_indices
```

### Seeds

`seed` seeds the dummy draws. With a fixed seed a run is reproducible on one platform. The dummies
come from the C++ standard library's random engines, whose Normal distribution differs between
standard library implementations, so a seeded run on macOS, Linux and Windows draws different
dummies and can select differently in the borderline cases. For a Monte Carlo study over many trials
pass `seed=-1` (the default): every trial then draws its own dummies from the hardware entropy
source, while a fixed seed would give every trial the same per-experiment seeds and correlate the
trials.

### The number of random experiments

`K` is the number of random experiments. The T-Rex paper uses $K = 20$, the default. A larger $K$
makes the vote finer, since the voting grid has steps of $1/K$, at a proportional cost.

```{code-cell} python
sel.voting_grid
```

### Early stopping of the T-loop

The T-loop raises the stopping time $T$ one dummy at a time. `use_max_T_stop = True` (default)
bounds it, and `tloop_stagnation_stop` stops the loop when the FDP estimate has not improved for
`tloop_max_stagnant_steps` steps. The stagnation stop defaults to `None`, which resolves to off for
the path solvers (T-LARS, T-LASSO, the elastic nets, the CD family) and to on for the greedy
solvers.

## The result fields

```{code-cell} python
print("T_stop:  ", sel.T_stop)
print("L:       ", sel.L)
print("v_thresh:", sel.v_thresh)
```

`fdp_hat_mat` is the FDP estimate per stopping time (rows, up to `T_stop`) and voting level
(columns), `r_mat` the number of selected variables per stopping time and voting level, and
`phi_mat` the relative occurrence of each variable over the $K$ experiments after each stopping
time. The selection sits where the estimate is below the target and the count is largest. The rows
are 0-based, so the row of the last stopping time is `T_stop - 1`.

```{code-cell} python
fdp = sel.fdp_hat_mat
print(fdp.shape)
print(np.round(fdp[sel.T_stop - 1], 3))
sel.r_mat
```

```{code-cell} python
import matplotlib.pyplot as plt

phi = sel.phi_mat[sel.T_stop - 1]
variable = np.arange(p)
is_true = np.isin(variable, true_support)
fig, ax = plt.subplots(figsize=(6, 4))
ax.vlines(variable[is_true], 0, phi[is_true], colors="C0", linewidth=2, label="true support")
ax.vlines(variable[~is_true], 0, phi[~is_true], colors="0.6", linewidth=2, label="other")
ax.axhline(sel.v_thresh, linestyle="--", color="C1", label="voting threshold")
ax.set_xlabel("Variable")
ax.set_ylabel("Relative occurrence")
ax.legend()
plt.show()
```

Relative occurrences after the last stopping time. The true support stands out, the dashed line is
the voting level the selection was read at.

`phi_prime` is the calibrated occurrence that enters the FDP estimate (at most 1, possibly
negative), `x_means`, `x_l2_norms` and `y_mean` are the centering and scaling constants the selector
applied, and `num_dummies` is the number of dummies of the last run.

## The L-loop strategies

`lloop_strategy` of `TRexControlParameter` says how the dummies of the L-loop are produced, with the
`LLoopStrategy` enum.

| Strategy | Dummies | Storage |
| --- | --- | --- |
| `STANDARD` (default) | fresh matrices at every L step | $K$ matrices |
| `HCONCAT` | the same $K$ matrices grow by $p$ columns per L step | $K$ matrices |
| `SEEDED` | the draws of `STANDARD`, re-derived from the seed at every step | none |
| `PERMUTATION` | permuted columns of the design | $K$ matrices |
| `PERMUTATION_SEEDED` | the permutation re-derived from the seed | none |
| `SKIPL` | no L-loop, the cap `max_dummy_multiplier` is used directly | $K$ matrices |

`SKIPL` skips the calibration of $L$ and runs with the cap at once. Under a large fixed $L$ without
the stagnation stop its FDP estimate is anti-conservative. The expected null inflow per T step
shrinks with $L$, so the T-loop runs deep. There the same few null variables keep entering across
the random experiments, because $X$ and $y$ are the same in all of them, and they collect more
occurrence above the reference voting level than the estimate books as false. The total inflow
itself stays about right. In the Examples demo 04 (signal-to-noise ratio 5, $T$ about 145) the
strategy with twenty times $p$ dummies under T-LARS reached an FDR of 0.133 at a target of 0.1, and
0.026 with the stagnation stop on. The selector warns from `max_dummy_multiplier = 10` on when the
stop is off, and the calibrated strategies are the safe choice.

The seeded strategies store nothing and suit designs where $K \cdot L \cdot p$ columns do not fit in
memory. `SEEDED` selects exactly as `STANDARD` for the same seed, and `PERMUTATION_SEEDED` exactly
as `PERMUTATION`. Fresh and grown dummies give the same FDR. A paired study of 500 trials in nine
cells (signal-to-noise ratio 0.5, 1 and 2 times target FDR 0.05, 0.1 and 0.2, $n = 300$, $p =
1{,}000$, ten active variables, $K = 20$, T-LARS) found the FDR equal within Monte Carlo error in
every cell and the grown form ahead in TPR by 0.6 percentage points on average, because it
calibrates a larger $L$ and runs longer T-loops. The table lists the cells by signal-to-noise ratio
(SNR) and target FDR (tFDR).

| SNR | tFDR | FDR fresh | FDR grown | TPR fresh | TPR grown |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 0.05 | 0.011 | 0.010 | 0.113 | 0.115 |
| 0.5 | 0.10 | 0.041 | 0.041 | 0.253 | 0.266 |
| 0.5 | 0.20 | 0.103 | 0.103 | 0.395 | 0.398 |
| 1.0 | 0.05 | 0.016 | 0.017 | 0.598 | 0.608 |
| 1.0 | 0.10 | 0.038 | 0.038 | 0.738 | 0.738 |
| 1.0 | 0.20 | 0.096 | 0.096 | 0.845 | 0.851 |
| 2.0 | 0.05 | 0.014 | 0.015 | 0.941 | 0.952 |
| 2.0 | 0.10 | 0.036 | 0.035 | 0.983 | 0.982 |
| 2.0 | 0.20 | 0.137 | 0.129 | 0.998 | 0.997 |

```{code-cell} python
for strategy in (tsn.LLoopStrategy.HCONCAT, tsn.LLoopStrategy.STANDARD, tsn.LLoopStrategy.SEEDED):
    ctrl = tsn.TRexControlParameter()
    ctrl.lloop_strategy = strategy
    ctrl.K = 10
    s_sel = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=42, verbose=False)
    s_sel.select()
    print(f"{strategy.name:<10} L = {s_sel.L}  T_stop = {s_sel.T_stop}  "
          f"selected: {' '.join(str(j) for j in s_sel.selected_indices)}")
```

## Dummy distributions

The dummies are standard Normal by default. `DummyDistribution` offers the uniform, Rademacher,
Student t, Laplace, Gumbel, triangle, uniform sphere, Mammen {cite}`mammen1993bootstrap`,
constrained sparse Rademacher and logistic distributions through its named constructors,
`DummyDistribution.normal()`, `DummyDistribution.student_t(df)` and so on. The FDR guarantee rests
on the exchangeability of dummies and null variables along the path, which the standard Normal
provides for the centered and scaled working design.

```{code-cell} python
ctrl = tsn.TRexControlParameter()
ctrl.dummy_distribution = tsn.DummyDistribution.student_t(3.0)
ctrl.K = 10
t_sel = tsn.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl, seed=42, verbose=False)
t_sel.select()
t_sel.selected_indices
```

## Scaling modes

The selector works on a centered and scaled design. `scaling_mode = ScalingMode.L2` (default)
centers `X` and `y` and divides each column of `X` by its L2 norm, `ScalingMode.ZSCORE` divides by
the sample standard deviation instead, and the dummies follow the same convention.
`ScalingMode.NONE` takes `X` and `y` as supplied, for a design that was pre-processed by the caller
or for a dictionary whose columns must not be centered (steering vectors, Fourier atoms).
`normalize_design()` reproduces the selector's own preprocessing, so a pre-normalized design under
`ScalingMode.NONE` selects exactly as the raw design under `ScalingMode.L2`.

```{code-cell} python
d = tsn.normalize_design(X, y, scaling_mode=tsn.ScalingMode.L2)
print(np.round(d["X"].mean(axis=0)[:3], 12))
print(np.linalg.norm(d["X"], axis=0)[:3])
ctrl = tsn.TRexControlParameter()
ctrl.scaling_mode = tsn.ScalingMode.NONE
n_sel = tsn.TRexSelector(d["X"], d["y"], tFDR=0.1, trex_control=ctrl, seed=42, verbose=False)
n_sel.select()
np.array_equal(n_sel.selected_indices, sel.selected_indices)
```

Non-finite entries and constant columns are rejected in every mode with a `ValueError` that names
the column.

## The design is used in place

The selector centers and scales `X` in place at construction and writes the original values back
when `select()` returns, so no copy of the design is made. Between construction and `select()` the
array holds the working design. `restore()` gives it back early when a constructed selector will not
run. The binding hands a Fortran-ordered `float64` array to the core as it is, which is why the
design of this page was created with `np.asfortranarray()`. An array of another layout or type is
converted first, and the selector then works on that converted copy while the caller's array stays
untouched. The page on working without data copies has the full contract.

```{code-cell} python
X_before = X.copy()
r_sel = tsn.TRexSelector(X, y, tFDR=0.1, seed=42, verbose=False)
print(np.abs(X - X_before).max())      # the working design
r_sel.restore()
print(np.abs(X - X_before).max())      # the original values are back
```

## FDR control in a small simulation

A Monte Carlo study over three target levels. Each trial draws its own data and its own dummies.

```{code-cell} python
tFDR_vec = [0.1, 0.2, 0.3]
MC = 20
fdp = np.empty((MC, len(tFDR_vec)))
tpp = np.empty_like(fdp)
for t, target in enumerate(tFDR_vec):
    for mc in range(MC):
        seed = 1000 * (t + 1) + mc + 1
        rng_mc = np.random.default_rng(seed)
        X_mc = np.asfortranarray(rng_mc.standard_normal((n, p)))
        y_mc = X_mc @ beta + rng_mc.standard_normal(n)
        ctrl = tsn.TRexControlParameter()
        ctrl.K = 10
        s = tsn.TRexSelector(X_mc, y_mc, tFDR=target, trex_control=ctrl, seed=seed, verbose=False)
        s.select()
        fdp[mc, t] = compute_fdp(s.selected_indices, true_support)
        tpp[mc, t] = compute_tpp(s.selected_indices, true_support)
for target, f, tp in zip(tFDR_vec, fdp.mean(axis=0), tpp.mean(axis=0)):
    print(f"tFDR {target:.1f}   FDR {f:.3f}   TPR {tp:.3f}")
```

```{code-cell} python
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.5))
ax1.plot(tFDR_vec, fdp.mean(axis=0), "o-", label=f"mean over {MC} trials")
ax1.plot([0, 0.35], [0, 0.35], "--", color="C1", label="FDR equal to the target")
ax1.set_xlim(0.05, 0.35)
ax1.set_ylim(0, 0.5)
ax1.set_xlabel("Target FDR")
ax1.set_ylabel("FDR")
ax2.plot(tFDR_vec, tpp.mean(axis=0), "o-")
ax2.set_xlim(0.05, 0.35)
ax2.set_ylim(0, 1.05)
ax2.set_xlabel("Target FDR")
ax2.set_ylabel("TPR")
fig.legend(*ax1.get_legend_handles_labels(), loc="outside center right")
plt.show()
```

Estimated FDR against the target (left) and the TPR (right), 20 trials per level. The estimated FDR
stays below the target at every level. The simulations of the T-Rex paper {cite}`machkour2025trex`
show the same over a far wider range of settings.
