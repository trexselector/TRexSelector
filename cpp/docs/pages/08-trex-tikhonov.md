# T-Rex+Tikhonov: a general penalty matrix {#trex_tikhonov}

[TOC]

## From the elastic net to a penalty matrix

The elastic net adds the ridge penalty \f$\frac{\lambda_2}{2} \|\boldsymbol\beta\|_2^2\f$ to the
lasso. Tikhonov regularization in general form \cite tikhonov1963regularization and
\cite park2018tikhonov replaces it by
\f$\frac{\lambda_2}{2} \|\boldsymbol\Gamma \boldsymbol\beta\|_2^2 =
\frac{\lambda_2}{2} \boldsymbol\beta^\top \boldsymbol K \boldsymbol\beta\f$
with \f$\boldsymbol K = \boldsymbol\Gamma^\top \boldsymbol\Gamma\f$. The matrix encodes what the
analyst knows about the coefficients. \f$\boldsymbol K = \boldsymbol I\f$ gives the ridge penalty of
the elastic net, a first-difference operator \f$\boldsymbol\Gamma\f$ gives the graph Laplacian of a
chain and asks neighbouring coefficients to be equal (smoothness), a graph Laplacian in general
\cite smola2003graphkernels and \cite chung1997spectral couples the variables along the edges of a
graph, a group-mean matrix reproduces the penalty of the informed elastic net, and a block-diagonal
centering matrix asks the members of a group to share a coefficient. The Tikhonov T-Rex selector
runs the terminating Tikhonov-net solver TSolverTIKNET, a pathwise homotopy algorithm, with such a
matrix inside the random experiments.

The headers and names of this page:

\snippet{trimleft} trex_tikhonov.cpp includes

A small helper prints index vectors:

\snippet{trimleft} trex_tikhonov.cpp helpers

## A smooth coefficient profile

The active variables are a contiguous block with equal coefficients, the case a chain Laplacian is
made for. The data are drawn with a fixed seed, and the indices are 0-based.

\snippet{trimleft} trex_tikhonov.cpp data

`TRexTikhonovSelector::gammaToK()` builds \f$\boldsymbol K\f$ from a sparse operator
\f$\boldsymbol\Gamma\f$ and returns it pruned and compressed. The penalty matrix is an
`Eigen::SparseMatrix<double>` everywhere, in the control and in the solvers.

\snippet{trimleft} trex_tikhonov.cpp penalty
\output{trex_tikhonov.penalty.txt}

## Running the selector

The Tikhonov control is required, because it carries \f$\boldsymbol K\f$. It nests the base control
as `trex_ctrl`, so the number of random experiments and the solver are set there. The penalty weight
is fixed here.

\snippet{trimleft} trex_tikhonov.cpp run
\output{trex_tikhonov.run.txt}

TRexTikhonovSelector inherits the getters of TRexSelector and adds `getLambda2Used()`,
`getLambda2CVDiagnostics()`, `getLambda2PilotResult()` and `getDummyClusters()`. A selector runs once,
as every selector does.

## Choosing lambda_2

`lambda_2 < 0` (the default) hands the choice to `lambda2_method`, a `TikLambda2Method`. Three rules
exist.

- `PILOT_KNEE` (the default) reads \f$\lambda_2\f$ off a few random experiments per candidate of
  `pilot_grid`, by the selector's own principle: the smallest \f$\lambda_2\f$ whose pilot selection
  count is within `pilot_knee_tolerance` of the largest, the knee of the profile. `pilot_K`,
  `pilot_T` and `pilot_L_multiplier` size the pilot, `pilot_statistic` names the count
  (`PilotStatistic::MAJORITY`, `CALIBRATED`, `GROUP_MAJORITY`). `runPilotOnly()` runs the pilot
  without the selection.
- The k-fold cross-validation rules consume \f$\boldsymbol K\f$: `CV_1SE_IEN_CCD`, `CV_MIN_IEN_CCD`,
  `CV_1SE_TIK_SVD` and `CV_MIN_TIK_SVD`. `cv_n_folds`, `cv_n_lambda1`, `cv_n_lambda2`,
  `cv_dfmax_ratio` and `cv_seed` control them. They judge a fitted estimate, so on structured
  designs they tend to pick the edges of their grid.
- `tikhonov_lambda2_power_rule()` runs one full selection per grid value with common dummies and
  takes the knee of that profile. It is the reference rule behind the pilot, at the cost of one
  selection per grid value, and it returns every probe's selection. Its grid and tolerance come in a
  `Lambda2PowerRuleConfig`, whose grid defaults to `default_lambda2_grid()`, the grid of the pilot
  rule. `gvs_lambda2_power_rule()` is the same rule on the T-Rex+GVS selector.

A supplied `lambda_2` is used as it is, on the scale of the standardized columns. A grid, `pilot_grid`
or the grid of a power rule, is stated for columns of unit norm. Under `ScalingMode::ZSCORE` the
columns have the squared norm \f$n - 1\f$, so the rules run at the grid times \f$n - 1\f$, and one
grid means the same shrinkage under every scaling mode. Every rule reports its pick on the scale of a
supplied value, so `getLambda2Used()` can be passed as `lambda_2`.

\snippet{trimleft} trex_tikhonov.cpp pilot
\output{trex_tikhonov.pilot.txt}

\snippet{trimleft} trex_tikhonov.cpp power-rule
\output{trex_tikhonov.power-rule.txt}

## Dummies for a structured penalty

The false discovery rate (FDR) guarantee needs dummies that are exchangeable with the null variables
under the penalty. Two mechanisms take care of that. `cluster_dummies = true` (default) draws the
dummies as per-cluster coloured layers on a correlation clustering of the variables (`prior_groups`,
or hierarchical clustering with `corr_max` and `hc_linkage`), and every dummy joins the
\f$\boldsymbol K\f$-coupling of the variable it mirrors (the folded coupling). `cluster_dummies =
false` exists for validation runs. `getDummyClusters()` reports the clusters of a run.

## Group centering

`penalty_structure = TikPenaltyStructure::GROUP_CENTERING` lets the selector build
\f$\boldsymbol K\f$ itself as the block-diagonal centering penalty
\f$\boldsymbol I - \boldsymbol J / p_m\f$ over the groups of `prior_groups`, so the members of a
group share a coefficient. The structure needs `prior_groups` (0-based group ids, one per variable),
and the selector refuses it without them. An all-singleton partition is legitimate, the penalty
vanishes and the path is the lasso on per-variable cluster dummies.
`TRexTikhonovSelector::buildGroupCenteringK()` returns this matrix for a group assignment.

\snippet{trimleft} trex_tikhonov.cpp group-centering
\output{trex_tikhonov.group-centering.txt}

## The solver

The selector runs `SolverTypeForTRex::TIKNET`, the pathwise homotopy algorithm of the Tikhonov net. It is
the solver of the base control nested in the Tikhonov control, `trex_ctrl`. A fresh
`TRexControl` names `TLARS`, which the constructor refuses as every other solver, so a base
control that is assigned to `trex_ctrl` has to set `solver_type` to `TIKNET`.
