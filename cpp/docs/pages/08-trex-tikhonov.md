# T-Rex Tikhonov: a general penalty matrix {#trex_tikhonov}

[TOC]

## From the elastic net to a penalty matrix

The elastic net adds the ridge penalty \f$\lambda_2 \|\boldsymbol\beta\|_2^2\f$ to the lasso.
Tikhonov regularization in general form \cite tikhonov1963regularization and
\cite park2018tikhonov replaces it by \f$\lambda_2 \|\boldsymbol\Gamma \boldsymbol\beta\|_2^2 =
\lambda_2 \boldsymbol\beta^\top \boldsymbol K \boldsymbol\beta\f$ with \f$\boldsymbol K =
\boldsymbol\Gamma^\top \boldsymbol\Gamma\f$. The matrix encodes what the analyst knows about the
coefficients. \f$\boldsymbol K = \boldsymbol I\f$ is the elastic net, a first-difference operator
\f$\boldsymbol\Gamma\f$ gives the graph Laplacian of a chain and asks neighbouring coefficients to
be equal (smoothness), a graph Laplacian in general \cite smola2003graphkernels and
\cite chung1997spectral couples the variables along the edges of a graph, and a block-diagonal
centering matrix reproduces the penalty of the informed elastic net. The Tikhonov T-Rex selector
runs the terminating Tikhonov-net solvers (TCDTIKNET_Solver by coordinate descent, TIKNET_Solver by
the pathwise algorithm) with such a matrix inside the random experiments.

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
`getLambda2CVDiagnostics()` and `getLambda2PilotResult()`.

## Choosing lambda_2

`lambda_2 < 0` (the default) hands the choice to `lambda2_method`, a `TikLambda2Method`. Three rules
exist.

- `PILOT_KNEE` (the default) reads \f$\lambda_2\f$ off a few random experiments per candidate of
  `pilot_grid`, by the selector's own principle: the smallest \f$\lambda_2\f$ whose pilot selection
  count is within `pilot_knee_tolerance` of the largest, the knee of the profile. `pilot_K`,
  `pilot_T` and `pilot_L_multiplier` size the pilot, `pilot_statistic` names the count
  (`PilotStatistic::MAJORITY`, `CALIBRATED`, `GROUP_MAJORITY`).
- The k-fold cross-validation rules consume \f$\boldsymbol K\f$ and return \f$\lambda_2\f$ on the
  solver scale: `CV_1SE_IEN_CCD`, `CV_MIN_IEN_CCD`, `CV_1SE_TIK_SVD` and `CV_MIN_TIK_SVD`.
  `cv_n_folds`, `cv_n_lambda1`, `cv_n_lambda2`, `cv_dfmax_ratio` and `cv_seed` control them. They
  judge a fitted estimate, so on structured designs they tend to pick the edges of their grid.
- `tikhonov_lambda2_power_rule()` runs one full selection per grid value with common dummies and
  takes the knee of that profile. It is the reference rule behind the pilot, at the cost of one
  selection per grid value, and it returns every probe's selection. Its grid and tolerance come in a
  `Lambda2PowerRuleConfig`.

\snippet{trimleft} trex_tikhonov.cpp pilot
\output{trex_tikhonov.pilot.txt}

\snippet{trimleft} trex_tikhonov.cpp power-rule
\output{trex_tikhonov.power-rule.txt}

## Dummies for a structured penalty

The false discovery rate (FDR) guarantee needs dummies that are exchangeable with the null variables under
the penalty. Two defaults take care of that. `cluster_dummies = true` draws the dummies as per-cluster
coloured layers on a correlation clustering of the variables (`prior_groups`, or hierarchical clustering
with `corr_max` and `hc_linkage`), and `fold_dummy_coupling = true` lets every dummy join the
\f$\boldsymbol K\f$-coupling of the variable it mirrors (the FOLDED coupling, available with
TCDTIKNET_Solver). Both switches exist for validation runs. With \f$\boldsymbol K = \boldsymbol
I\f$, independent dummies and the decoupled dummy block, the selector reduces to the plain elastic
net.

\snippet{trimleft} trex_tikhonov.cpp elastic-net
\output{trex_tikhonov.elastic-net.txt}

## Group centering

`penalty_structure = TikPenaltyStructure::GROUP_CENTERING` lets the selector build
\f$\boldsymbol K\f$ itself as the block-diagonal centering penalty \f$\boldsymbol I - \boldsymbol J
/ p_m\f$ over the dummy clusters, so the members of a group share a coefficient, the penalty of the
informed elastic net. The groups come from `prior_groups` (0-based group ids, one per variable) or
from the clustering. An all-singleton partition is legitimate, the penalty vanishes and the path is
the lasso on per-variable cluster dummies. On this path the TCDTIKNET_Solver runs through the
group-centering hooks of TCDIENET_Solver, so its warnings (sweep budget, certification fallback)
carry the tag `[TCDIENET]`.

\snippet{trimleft} trex_tikhonov.cpp group-centering
\output{trex_tikhonov.group-centering.txt}

## The solver

`solver_type` of the nested base control picks the Tikhonov-family solver:
`SolverTypeForTRex::TCDTIKNET` (the default, the neutral `TLARS` of a fresh control maps to it) or
`TIKNET`, the pathwise algorithm. Both solve the same problem.

\snippet{trimleft} trex_tikhonov.cpp solver
\output{trex_tikhonov.solver.txt}
