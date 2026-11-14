# TRexSelector 2.0.2

## Changed

* `TRexDASelector` takes the penalty of the methods `"AR1"` and `"EQUI"` from the last step of the
  evaluated path and applies it to every step, as `"BT"`, `"NN"` and the prior groups always did. Up to
  2.0.1 the two methods recomputed the penalty at every step, which is the literal reading of the
  papers. A penalty per step makes the deflated occurrence history non-monotone, and the FDP estimate is
  then too small. In a Monte Carlo comparison of the two constructions for all methods (1,000 runs per
  setting) a penalty per step exceeded a target FDR of 0.1 with the nearest-neighbour method (0.15 and
  0.13 against 0.00 and 0.03 with the last step), and `"AR1"` gained 0.01 to 0.02 in true positive rate
  on AR(1) designs with the last step.
* Selections of `"AR1"` can differ slightly from 2.0.1. Those of `"EQUI"` did not change in any tested
  setting. The help of `TRexDASelector` and the vignette "T-Rex+DA" state the construction, the
  deviation from the papers and the evidence.
* `tsolver_options()` defaults to `lambda2 = 0.1` and `rho_afs = 0.3`, the values a selector hands to its
  solvers. Up to 2.0.1 they were 0 and 1. A stand-alone elastic-net solver built without `lambda2` no
  longer equals its lasso, and a stand-alone `TAFS_Solver` no longer runs forward stepwise. `lambda2 = 0`
  and `rho_afs = 1` give the earlier behaviour. The selectors do not change.

## Fixed

* The help of `trex_da_control()` no longer offers `hc_linkage = "Ward"`, which the selector rejects, and
  the error names the supported linkages. The help of the fields `FDP_hat_array_BT`, `Phi_array_BT` and
  `R_array_BT` of `TRexDASelector` states what the arrays hold. The vignette "T-Rex+DA" is titled
  "T-Rex+DA: dependent variables", because every dependency-aware method draws ordinary independent
  dummies.
* The strategy table of the README and of the vignette "The T-Rex selector" describes `"PERMUTATION"`
  as what it is, one base dummy matrix whose rows are permuted per experiment.
* Four vignettes print the difference between the design before and after a run. Their comments now say
  "restored up to rounding", because the design comes back with differences of the order of 1e-16, and
  the vignette "Working without data copies" says so in its text.
* The help of the fields `max_num_dummies` and `max_T_stop` says what they hold, the cap of the dummy
  multiplier and whether the bound on the stopping time is on. The help of `max_dummy_multiplier` and
  `use_max_T_stop` in `trex_control()` names the bound, `ceiling(n / 2)`.
* The help of `trex_da_control(da_method = )` and the vignette "T-Rex+DA" say which dependency model
  each method is made for.
* The field `da_method_used` of `TRexDASelector` returned `"UNKNOWN"` for a run with prior groups. It
  returns `"PRIOR_GROUPS"`.
* The help of `TRexDASelector` gains the section "The reference threshold". The methods `"BT"`, `"NN"`
  and the prior groups run their loops on one threshold of the grid and search all thresholds for the
  selection, so `fdp_hat_mat`, `phi_mat` and `phi_prime` belong to the reference threshold, while
  `r_mat`, `v_thresh` and the selection belong to `rho_thresh`. The vignette "T-Rex+DA" says the same.
* The help of the fields `rho_thresh`, `rho_grid`, `da_method_used` and `cor_coef` of `TRexDASelector`
  says what they hold, among it that `rho_thresh` is -2 for `"AR1"` and `"EQUI"`. The help of
  `trex_da_control()` says which methods use `hc_linkage` and `hc_grid_length`, how `cor_coef` is
  estimated, and gives the formula of the AR1 window, `ceiling(log(rho_thr_DA) / log(abs(cor_coef)))`
  columns to each side.
* The help of `TRexDASelector` gains the section "FDR limitation of BT, NN and prior groups". At a very
  high correlation within the groups and with no independent predictors the three methods can exceed
  the target FDR with the solver `"TLARS"`. The C++ documentation carried the paragraph before.
* The vignette "T-Rex Tikhonov: a general penalty matrix" and the help of `trex_tikhonov_control()` say
  which penalty is which. The group-mean matrix reproduces the penalty of the informed elastic net, the
  centering matrix of `penalty_structure = "GROUP_CENTERING"` asks the members of a group to share a
  coefficient and needs `prior_groups`, the penalty carries the factor one half of the code, and
  `tikhonov_K` is optional under `"GROUP_CENTERING"`. The vignette no longer runs the deprecated solver
  `"TCDTIKNET"` with the decoupled dummy block.
* Three vignettes say why they copy the design with `X + 0`.
* The help of `tol` in `trex_control()` and in `tsolver_options()` says that the selector and the
  stand-alone solver have different defaults on purpose.
* The vignettes no longer show the deprecated coordinate-descent solver family. Its help pages keep the
  deprecation notice.
* The vignette "T-Rex GVS: grouped variable selection" describes the informed elastic net as its paper
  does. It keeps the grouping effect of the elastic net and is the elastic net when every variable is
  its own group.

## New

* `trex_da_control()` gains `bt_selection_mode` and `trex_gvs_control()` gains `group_labels`. Both
  parameters existed in C++ and Python.

## Deprecated

* The coordinate-descent solver family is deprecated and will be removed in TRexSelector 2.1.0. This
  covers the classes `TCDLASSO_Solver`, `TCDENET_Solver`, `TCDIENET_Solver` and `TCDTIKNET_Solver` with
  their base `TCDSolver_Base`, the solver names `"TCDLASSO"`, `"TCDENET"`, `"TCDIENET"` and
  `"TCDTIKNET"` of `trex_control(solver = )`, and `trex_gvs_control(en_solver = "TCDENET")`. The
  homotopy solvers of the same problems replace them, `TLASSO_Solver`, `TENET_Solver`, `TIENET_Solver`
  and `TIKNET_Solver` with the solver names `"TLASSO"`, `"TENET"`, `"TIENET"` and `"TIKNET"`. The
  family brought no gain over them in speed or in the selections, and its two structured members do
  not converge within a practical sweep budget on every design.
* Constructing one of the four classes raises a deprecation warning, and so does a selector that
  receives one of the solver names. The warning comes once per call and names the replacement.
* The options that only the family reads go with it. These are `cd_lambda_rel_tol`, `cd_tol_certify`,
  `cd_tol_probe`, `cd_gram_cap`, `cd_max_sweeps` and `cd_certified_localization` of `trex_control()`
  and `tsolver_options()`, `kappa_dummy`, `fold_dummies` and `group_penalty` of `tsolver_options()`,
  and `fold_dummy_coupling` of `trex_tikhonov_control()`. Their help says so, and they raise no warning.

# TRexSelector 2.0.1

A patch release. CRAN's daily checks of 2.0.0 reported "Installation failed" on the macOS builders for
R-release (x86_64) and R-oldrel (arm64 and x86_64), which compile with Apple clang 14.0.3. That compiler
accepts C++20 but its standard library has no ranges algorithms, and the backend used `std::ranges::find`
and `std::ranges::sort` at three places. The three calls are `std::find` and `std::sort` now, with identical
results.

On Linux with R 4.6 or R-devel the package compiled without OpenMP, because `R CMD config` no longer serves
`SHLIB_OPENMP_CXXFLAGS` there. The configure script reads the flags from R's Makeconf instead, so those builds
run their parallel regions again. The README is reorganised. It opens with a description of the method, gives
an overview of the selectors, T-solvers, ML methods and utilities, and installs from CRAN first, as the
getting-started vignette does. A new vignette, "Data dimensionality and FDR control", follows the getting-started
vignette with the three regimes of the number of variables against the number of observations in one map,
and an overview vignette, `vignette("TRexSelector")`, introduces the package and lists the pages. The
overview and the thirteen pages share a stylesheet and a navigation header, with a table of contents on
the left that follows the scroll, the list of pages under it, and previous/next links with the copyright
line at the end of every page, a search box over the titles and headings of the pages in the top bar,
and they follow a dark system appearance. Four pages carry figures (the T-Rex framework and a design with
dummies, the regime map, the memory model), five draw their results with ggplot2, and the method pages link
the demos of the Examples repository. The package title reads "T-Rex Selector: High-Dimensional Variable
Selection and FDR Control", and the package entry of `citation("TRexSelector")` names Team T-Rex as its
author. Nothing else changes.

# TRexSelector 2.0.0

The package is rewritten as R bindings to the TRexSelector C++ library. The pure-R implementation of
1.0.0 and its dependencies (tlars, MASS, glmnet, boot, parallel, doParallel, foreach, doRNG) are gone.

## Breaking changes

* The functional interface of 1.0.0 (`trex()`, `FDP()`, `TPP()`, `add_dummies_GVS()` and the other
  helpers) is replaced by R6 classes and control lists:
  `TRexSelector$new(X, y, tFDR, control = trex_control(...))$select()`, with `trex_fdp()` and
  `trex_tpp()` as the metrics.
* R 4.3 or newer and a C++20 compiler are required. OpenMP is optional.

## New

* Selector variants as R6 classes: `TRexDASelector` (dependency-aware), `TRexGVSSelector` (grouped
  variable selection with the elastic net and the informed elastic net), `TRexScreeningSelector` and
  `TRexBiobankScreeningSelector` (screening), `TRexSPCASelector` (sparse PCA) and
  `TRexTikhonovSelector` (a general Tikhonov penalty).
* The twenty terminating solvers as classes of their own (`TLARS_Solver`, `TLASSO_Solver`, ...) with
  `tsolver_options()`, working on the caller's matrices without copying.
* Machine learning utilities: scalers, ridge and elastic-net regression with cross-validation, SVD,
  PCA and hierarchical clustering (the `ML*` classes and `ml_*` functions).
* Memory-mapped matrices (`mem_matrix()`, `mem_matrix_from()`), accepted everywhere a matrix is.
* Twelve vignettes, and `citation("TRexSelector")` lists the package entry with the method papers.
