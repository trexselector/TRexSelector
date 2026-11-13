---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# T-Rex Tikhonov: a general penalty matrix

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## From the elastic net to a penalty matrix

The elastic net adds the ridge penalty $\frac{\lambda_2}{2} \|\boldsymbol\beta\|_2^2$ to the lasso.
Tikhonov regularization in general form {cite}`tikhonov1963regularization,park2018tikhonov` replaces
it by
$\frac{\lambda_2}{2} \|\boldsymbol\Gamma \boldsymbol\beta\|_2^2 =
\frac{\lambda_2}{2} \boldsymbol\beta^\top \boldsymbol K \boldsymbol\beta$
with $\boldsymbol K = \boldsymbol\Gamma^\top \boldsymbol\Gamma$. The matrix encodes what the analyst
knows about the coefficients. $\boldsymbol K = \boldsymbol I$ gives the ridge penalty of the elastic
net, a first-difference operator $\boldsymbol\Gamma$ gives the graph Laplacian of a chain and asks
neighbouring coefficients to be equal (smoothness), a graph Laplacian in general
{cite}`smola2003graphkernels,chung1997spectral` couples the variables along the edges of a graph, a
group-mean matrix reproduces the penalty of the informed elastic net, and a block-diagonal centering
matrix asks the members of a group to share a coefficient. The Tikhonov T-Rex selector runs the
terminating Tikhonov-net solver TIKNET, a pathwise homotopy algorithm, with such a matrix inside the
random experiments.

## A smooth coefficient profile

The active variables are a contiguous block with equal coefficients, the case a chain Laplacian is
made for. The data are drawn with a fixed seed, and the indices are 0-based.

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(2718)
n = 120
p = 60
true_support = list(range(20, 30))
beta = np.zeros(p)
beta[true_support] = 1.0
X = rng.standard_normal((n, p))
y = X @ beta + rng.standard_normal(n)
```

`gamma_to_k()` builds $\boldsymbol K$ from an operator $\boldsymbol\Gamma$ and returns it as a dense
array. The selector converts the matrix to a sparse one internally, so a dense array is accepted
everywhere a penalty matrix is.

```{code-cell} python
gamma = np.diff(np.eye(p), axis=0)    # first differences, (p - 1) x p
K = tsn.gamma_to_k(gamma)
print(K.shape)
print(K[:4, :4])
```

## Running the selector

The Tikhonov control is required, because it carries $\boldsymbol K$. The penalty weight is fixed
here.

```{code-cell} python
from trexselector.utils import compute_fdp, compute_tpp

tik_ctrl = tsn.TRexTikhonovControlParameter()
tik_ctrl.tikhonov_K = K
tik_ctrl.lambda_2 = 1.0
ctrl = tsn.TRexControlParameter()
ctrl.K = 10

tik = tsn.TRexTikhonovSelector(X, y, tFDR=0.1, tik_control=tik_ctrl, trex_control=ctrl,
                               seed=4, verbose=False)
tik.select()
print("selected:", tik.selected_indices)
print("FDP:", compute_fdp(tik.selected_indices, true_support))
print("TPP:", compute_tpp(tik.selected_indices, true_support))
print("lambda_2 used:", tik.lambda2_used)
```

`TRexTikhonovSelector` inherits the attributes of `TRexSelector` and adds `lambda2_used`,
`lambda2_cv_diagnostics` and `lambda2_pilot`.

## Choosing lambda_2

`lambda_2 < 0` (the default) hands the choice to `lambda2_method`, a `TikLambda2Method`. Three rules
exist.

- `PILOT_KNEE` (the default) reads $\lambda_2$ off a few random experiments per candidate of
  `pilot_grid`, by the selector's own principle: the smallest $\lambda_2$ whose pilot selection
  count is within `pilot_knee_tolerance` of the largest, the knee of the profile. `pilot_K`,
  `pilot_T` and `pilot_L_multiplier` size the pilot, `pilot_statistic` names the count
  (`PilotStatistic.MAJORITY`, `CALIBRATED`, `GROUP_MAJORITY`).
- The k-fold cross-validation rules consume $\boldsymbol K$ and return $\lambda_2$ on the solver
  scale: `CV_1SE_IEN_CCD`, `CV_MIN_IEN_CCD`, `CV_1SE_TIK_SVD` and `CV_MIN_TIK_SVD`. `cv_n_folds`,
  `cv_n_lambda1`, `cv_n_lambda2`, `cv_dfmax_ratio` and `cv_seed` control them. They judge a fitted
  estimate, so on structured designs they tend to pick the edges of their grid.
- `tikhonov_lambda2_power_rule()` runs one full selection per grid value with common dummies and
  takes the knee of that profile. It is the reference rule behind the pilot, at the cost of one
  selection per grid value, and it returns every probe's selection.

```{code-cell} python
pilot_ctrl = tsn.TRexTikhonovControlParameter()
pilot_ctrl.tikhonov_K = K
pilot_ctrl.pilot_grid = [0.0, 1.0, 10.0]
pilot_ctrl.pilot_K = 4
pilot_ctrl.pilot_T = 2

pilot = tsn.TRexTikhonovSelector(X, y, tFDR=0.1, tik_control=pilot_ctrl, trex_control=ctrl,
                                 seed=4, verbose=False)
pilot.select()
print("lambda_2 used:", pilot.lambda2_used)
print("grid:  ", pilot.lambda2_pilot.grid)
print("counts:", pilot.lambda2_pilot.counts)
```

```{code-cell} python
rule_ctrl = tsn.TRexTikhonovControlParameter()
rule_ctrl.tikhonov_K = K
profile = tsn.tikhonov_lambda2_power_rule(X, y, tFDR=0.1, tik_control=rule_ctrl, trex_control=ctrl,
                                          seed=4, grid=[0.0, 1.0, 10.0])
print("lambda_2:", profile.lambda2)
print("counts:  ", profile.counts())
```

## Dummies for a structured penalty

The false discovery rate (FDR) guarantee needs dummies that are exchangeable with the null variables
under the penalty. Two mechanisms take care of that. `cluster_dummies = True` (default) draws the
dummies as per-cluster coloured layers on a correlation clustering of the variables (`prior_groups`,
or hierarchical clustering with `corr_max` and `hc_linkage`), and every dummy joins the
$\boldsymbol K$-coupling of the variable it mirrors (the folded coupling). `cluster_dummies = False`
exists for validation runs.

## Group centering

`penalty_structure = TikPenaltyStructure.GROUP_CENTERING` lets the selector build $\boldsymbol K$
itself as the block-diagonal centering penalty $\boldsymbol I - \boldsymbol J / p_m$ over the groups
of `prior_groups`, so the members of a group share a coefficient. The structure needs `prior_groups`
(0-based group ids, one per variable), and the selector refuses it without them. An all-singleton
partition is legitimate, the penalty vanishes and the path is the lasso on per-variable cluster
dummies.

```{code-cell} python
groups = np.repeat(np.arange(6), 10)

gc_ctrl = tsn.TRexTikhonovControlParameter()
gc_ctrl.penalty_structure = tsn.TikPenaltyStructure.GROUP_CENTERING
gc_ctrl.prior_groups = groups
gc_ctrl.lambda_2 = 1.0

gc = tsn.TRexTikhonovSelector(X, y, tFDR=0.1, tik_control=gc_ctrl, trex_control=ctrl,
                              seed=4, verbose=False)
gc.select()
gc.selected_indices
```

## The solver

The selector runs `SolverTypeForTRex.TIKNET`, the pathwise homotopy algorithm of the Tikhonov net. The
neutral default `TLARS` of a fresh control maps to it.
