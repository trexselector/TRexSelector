# T-Rex GVS: grouped variable selection {#trex_gvs}

[TOC]

## Groups of correlated variables

In genomics and finance, variables come in groups of highly correlated members that carry the same
signal. The lasso picks one member of a group and drops the others, while the elastic net
\cite zou2005elasticnet spreads the coefficient over the group and lets the members enter together.
The T-Rex selector for grouped variable selection (GVS) \cite machkour2022gvs runs the terminating
elastic net (T-ENET) inside the random experiments and draws its dummies in groups too, so that the
dummies mirror the group structure of the null variables. The informed elastic net (IEN)
\cite machkour2023ien goes one step further and encodes the known or estimated groups in the
penalty, so that the members of a group share a coefficient.

## Data with groups

The headers and names of this page, and a helper that prints a selection in sorted order:

\snippet{trimleft} trex_gvs.cpp includes

\snippet{trimleft} trex_gvs.cpp helpers

Twenty groups of five variables, each group driven by one latent factor. The first two groups are
active. The group labels and the variable indices are 0-based.

\snippet{trimleft} trex_gvs.cpp data

## The elastic net policy

`gvs_type = GVSType::EN` is the default. The selector clusters the variables by correlation (single
linkage, cut at `corr_max`), draws its dummies per cluster and runs the T-ENET solver chosen by
`en_solver` (`ENSolverType::TENET` by default, `TENET_AUG` the row-augmented form and `TCDENET` the
coordinate descent form, which select the same variables). The elastic-net penalty \f$\lambda_2\f$
is fixed here. The group control is a TRexGVSControlParameter, which nests the base control as
`trex_ctrl`, and `getGVSResult()` returns the `GVSSelectionResult` after `select()`.

\snippet{trimleft} trex_gvs.cpp en
\output{trex_gvs.en.txt}

The result records the clustering and the penalty.

\snippet{trimleft} trex_gvs.cpp en-result
\output{trex_gvs.en-result.txt}

## Choosing lambda_2

A negative `lambda_2` (the default `-1`) lets the selector determine the penalty. `lambda2_method`
names the rule. For the elastic net the default is `LambdaSelectionMethod::CV_1SE_CCD`, a ridge
cross-validation by coordinate descent read with the one-standard-error rule, with `CV_MIN_CCD`,
`CV_1SE_SVD` and `CV_MIN_SVD` as alternatives. The one-standard-error rule is preferred because the
minimum rule tends to pin the penalty to the floor of the grid. `cv_n_folds`, `cv_n_lambda` and
`cv_seed` control the cross-validation. `getLambda2CVDiagnostics()` (also `lambda2_cv` of the
result) is the health record of the tuner: whether it ran, its seconds, and `n_profile_at_edge`, the
number of \f$\lambda_2\f$ values whose profiled \f$\lambda_1\f$ optimum sat on a truncation edge of
the path (0 is healthy).

\snippet{trimleft} trex_gvs.cpp cv
\output{trex_gvs.cv.txt}

## The informed elastic net

`gvs_type = GVSType::IEN` encodes the groups in the penalty. The groups come from `prior_groups` in
the control (0-based labels of length \f$p\f$) or from the correlation clustering. The solver is set
in `trex_ctrl.solver_type`: `SolverTypeForTRex::TIENET` (default, the native pathwise form),
`TIENET_AUG` (row augmented) or `TCDIENET` (coordinate descent). The IEN rules of `lambda2_method`
(`CV_1SE_IEN_CCD`, the default on this track, `CV_MIN_IEN_CCD`, `CV_1SE_TIK_SVD`, `CV_MIN_TIK_SVD`)
consume the group structure and return \f$\lambda_2\f$ on the solver scale. The penalty is left to
the IEN rule here.

\snippet{trimleft} trex_gvs.cpp ien
\output{trex_gvs.ien.txt}

With the groups known and the penalty informed, the active groups enter complete.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of T-Rex GVS on block designs of several
kinds, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/trex_selector_methods/ad/trex_gvs),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/trex_selector_methods/ad/trex_gvs)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/trex_selector_methods/ad/trex_gvs).