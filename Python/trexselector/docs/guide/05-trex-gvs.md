---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# T-Rex GVS: grouped variable selection

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Groups of correlated variables

In genomics and finance, variables come in groups of highly correlated members that carry the same
signal. The lasso picks one member of a group and drops the others, while the elastic net
{cite}`zou2005elasticnet` spreads the coefficient over the group and lets the members enter together.
The T-Rex selector for grouped variable selection (GVS) {cite}`machkour2022gvs` runs the terminating
elastic net (T-ENET) inside the random experiments and draws its dummies in groups too, so that the
dummies mirror the group structure of the null variables. The informed elastic net (IEN)
{cite}`machkour2023ien` goes one step further and encodes the known or estimated groups in the
penalty, so that the members of a group share a coefficient.

## Data with groups

Twenty groups of five variables, each group driven by one latent factor. The first two groups are
active. The group labels and the variable indices are 0-based.

```{code-cell} python
import numpy as np
import trexselector as tsn
from trexselector.utils import compute_fdp, compute_tpp

rng = np.random.default_rng(11)
n, p, group_size = 100, 100, 5
groups = np.repeat(np.arange(p // group_size), group_size)
latent = rng.standard_normal((n, groups.max() + 1))
X = 0.9 * latent[:, groups] + np.sqrt(1 - 0.9**2) * rng.standard_normal((n, p))
true_support = np.flatnonzero(np.isin(groups, [0, 1])).tolist()
beta = np.zeros(p)
beta[true_support] = 1.0
y = X @ beta + rng.standard_normal(n)
```

## The elastic net policy

`gvs_type = GVSType.EN` is the default. The selector clusters the variables by correlation (single
linkage, cut at `corr_max`), draws its dummies per cluster and runs the T-ENET solver chosen by
`en_solver` (`ENSolverType.TENET` by default, `TENET_AUG` the row-augmented form and `TCDENET` the
coordinate descent form, which select the same variables). The elastic-net penalty $\lambda_2$ is
fixed here. The group control is a `TRexGVSControlParameter`, the base control a
`TRexControlParameter`, and `select()` returns a `GVSSelectionResult`.

```{code-cell} python
gvs_ctrl = tsn.TRexGVSControlParameter()
gvs_ctrl.gvs_type = tsn.GVSType.EN
gvs_ctrl.lambda_2 = 1.0
ctrl = tsn.TRexControlParameter()
ctrl.K = 10

gvs = tsn.TRexGVSSelector(X, y, tFDR=0.1, gvs_control=gvs_ctrl, trex_control=ctrl,
                          seed=5, verbose=False)
res = gvs.select()
print("selected:", sorted(gvs.selected_indices))
print("FDP:", compute_fdp(gvs.selected_indices, true_support))
print("TPP:", compute_tpp(gvs.selected_indices, true_support))
```

The result records the clustering and the penalty.

```{code-cell} python
print("gvs_type:      ", res.gvs_type.name)
print("max_clusters:  ", res.max_clusters)
print("hc_method_used:", res.hc_method_used)
print("lambda2_used:  ", res.lambda2_used)
labels, counts = np.unique(res.groups_vec, return_counts=True)
print("first clusters:", dict(zip(labels[:5].tolist(), counts[:5].tolist())))
```

## Choosing lambda_2

A negative `lambda_2` (the default `-1`) lets the selector determine the penalty. `lambda2_method`
names the rule. For the elastic net the default is `LambdaSelectionMethod.CV_1SE_CCD`, a ridge
cross-validation by coordinate descent read with the one-standard-error rule, with `CV_MIN_CCD`,
`CV_1SE_SVD` and `CV_MIN_SVD` as alternatives. The one-standard-error rule is preferred because the
minimum rule tends to pin the penalty to the floor of the grid. `cv_n_folds`, `cv_n_lambda` and
`cv_seed` control the cross-validation. `lambda2_cv` of the result is the health record of the
tuner: whether it ran, its seconds, and `n_profile_at_edge`, the number of $\lambda_2$ values whose
profiled $\lambda_1$ optimum sat on a truncation edge of the path (0 is healthy).

```{code-cell} python
cv_ctrl = tsn.TRexGVSControlParameter()
cv_ctrl.gvs_type = tsn.GVSType.EN
cv_ctrl.cv_n_lambda = 100

cv = tsn.TRexGVSSelector(X, y, tFDR=0.1, gvs_control=cv_ctrl, trex_control=ctrl,
                         seed=5, verbose=False)
cv_res = cv.select()
diag = cv_res.lambda2_cv
print("lambda2_used:     ", cv_res.lambda2_used)
print("ran:              ", diag.ran)
print("seconds:          ", round(diag.seconds, 3))
print("n_profile_at_edge:", diag.n_profile_at_edge)
print("grid position:    ", diag.grid_index, "of", diag.grid_size)
print("boundary_pick:    ", diag.boundary_pick)
```

## The informed elastic net

`gvs_type = GVSType.IEN` encodes the groups in the penalty. The groups come from `prior_groups` in
the control (0-based labels of length $p$) or from the correlation clustering. The solver is set in
`TRexControlParameter.solver_type`: `SolverTypeForTRex.TIENET` (default, the native pathwise form),
`TIENET_AUG` (row augmented) or `TCDIENET` (coordinate descent). The IEN rules of `lambda2_method`
(`CV_1SE_IEN_CCD`, the default on this track, `CV_MIN_IEN_CCD`, `CV_1SE_TIK_SVD`, `CV_MIN_TIK_SVD`)
consume the group structure and return $\lambda_2$ on the solver scale. The penalty is left to the
IEN rule here.

```{code-cell} python
ien_ctrl = tsn.TRexGVSControlParameter()
ien_ctrl.gvs_type = tsn.GVSType.IEN
ien_ctrl.prior_groups = groups.tolist()
ien_base = tsn.TRexControlParameter()
ien_base.solver_type = tsn.SolverTypeForTRex.TIENET
ien_base.K = 10

ien = tsn.TRexGVSSelector(X, y, tFDR=0.1, gvs_control=ien_ctrl, trex_control=ien_base,
                          seed=5, verbose=False)
ien_res = ien.select()
print("selected:    ", sorted(ien.selected_indices))
print("lambda2_used:", round(ien_res.lambda2_used, 4))
print("FDP:", compute_fdp(ien.selected_indices, true_support))
print("TPP:", compute_tpp(ien.selected_indices, true_support))
```

With the groups known and the penalty informed, the active groups enter complete.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of T-Rex GVS on block designs of several
kinds, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/trex_selector_methods/ad/trex_gvs),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/trex_selector_methods/ad/trex_gvs)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/trex_selector_methods/ad/trex_gvs).