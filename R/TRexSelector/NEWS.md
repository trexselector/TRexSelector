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
* `hc_linkage` of `trex_da_control()`, `trex_gvs_control()` and `trex_tikhonov_control()` takes the
  linkage in lower case, `"single"`, `"complete"`, `"average"` or `"wpgma"`, as the field
  `hc_method_used` returns it. The capitalized names of 2.0.1 (`"Single"`) are still accepted.
* `TRexScreeningSelector` and `TRexBiobankScreeningSelector` take the `lloop_strategy` of the control.
  Up to 2.0.1 they replaced it by `"STANDARD"` without a word. `"SEEDED"`, `"PERMUTATION"` and
  `"PERMUTATION_SEEDED"` now run, as in C++ and Python, and `"HCONCAT"` and `"SKIPL"` are refused with
  an error, because the screening is one pass at L = p. The default `"STANDARD"` gives the results of
  2.0.1.
* `TRexBiobankScreeningSelector` is quiet by default (`verbose = FALSE`), as in C++ and Python. Up to
  2.0.1 it printed its progress by default.
* `trex_biobank_control()` gains `target_FDR_trex`, the target FDR level of the fallback T-Rex selector
  (default 0.1), as in C++ and Python. `tFDR` of the constructor defaults to `NULL` and still works: a
  supplied value replaces the target of the control.
* `TRexTikhonovSelector` and `trex_tikhonov_lambda2_power_rule()` read their `lambda_2` grid
  (`pilot_grid`, `grid`) as values for columns of unit norm. Under `scaling_mode = "ZSCORE"` they run at
  the grid times `n - 1`, as the cross-validation rules already returned their pick. Up to 2.0.1 the
  grid was used as it was, so under `"ZSCORE"` the default grid meant a shrinkage `n - 1` times weaker
  than under `"L2"`. Selections with `"ZSCORE"` and a grid can differ from 2.0.1, those with the default
  `"L2"` do not. The pick is reported on the working scale, and the pilot record and the result of the
  power rule gain `grid_scale`. A single supplied `lambda_2` is used as it is, as before.
* `TRexSPCASelector` takes `tFDR` in [0, 1], as the other selectors do and as C++ and Python do. Up to
  2.0.1 it refused 0 and 1. An `M` above `min(n, p)` is refused at construction. Up to 2.0.1 that error
  came in `select()`.

## Fixed

* A selector runs once. `select()` gives the design back when it ends, and up to 2.0.1 a second
  `select()` on the same object ran on the design at its original scale, without an error, and returned
  a wrong selection (all variables in a test with 10 true ones among 60). The second call is an error
  now, for `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`, `TRexScreeningSelector` and
  `TRexTikhonovSelector`. Construct a new selector to run again. The fields keep the results of the run.
* A target FDR that is `NaN` is refused. Up to 2.0.1 it was accepted and nothing was selected.
* `TRexTikhonovSelector` checks its control at construction. A penalty matrix that is not finite, not
  symmetric or not positive semi-definite, and a wrong number of folds of a cross-validation rule were
  reported in `select()`, under the name of a solver or a tuner. `lambda_2 = NaN` ended in an internal
  error, and `lambda_2 = Inf` was accepted and selected nothing. All of them are errors of the
  constructor now.
* Under `penalty_structure = "GROUP_CENTERING"` a supplied `tikhonov_K` is compared with the centering
  matrix of `prior_groups`, and another matrix is an error. Up to 2.0.1 it was taken as it was, and a
  cross-validation rule then tuned `lambda_2` on a penalty that the solver did not apply.
* `TRexTikhonovSelector` checks the penalty matrix once. Up to 2.0.1 every solver of a selection
  checked it again with a matrix factorization, at least 110 times in a default run. For a dense matrix
  at p = 2,000 these checks alone took about 73 seconds.
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
* The help of `TRexDASelector` says how the reference threshold is rounded. The package rounds a half
  up, while the 1.x versions rounded it to the even number, so they took one element lower for the grid
  lengths 6, 14, 22, 30 and so on. The default length 20 is not affected.
* `trex_gvs_control()` checks `groups` with the words of R: "groups must number the groups 1 to M
  without gaps". The errors of the core had asked for 0-based labels and had named a missing group
  by its 0-based number. A wrong length of `groups` is reported with the name of the argument.
* The help of `trex_gvs_control()` and of `TRexGVSSelector` is corrected. `corr_max` bounds the
  correlation between clusters with `"Single"` linkage only. A supplied `lambda_2` refers to the scale
  of the standardized columns, not to LARS units, so `scaling_mode = "ZSCORE"` needs the value times
  `n - 1`. `tenet_aug_use_lars` also acts for `gvs_type = "IEN"` with the solver `"TIENET_AUG"`.
  `gvs_type` is the penalty of the grouped selection. The fields `lambda2_used`, `max_clusters`,
  `hc_method_used` and `groups` say what they hold, and the help names the L-loop strategies the
  selector takes.
* The help of `TRexScreeningSelector` and of `trex_screen_control()` is corrected. The fields
  `estimated_correlation`, `beta_mat`, `dummy_betas`, `confidence_level` and `used_bootstrap` say what
  they hold. The variants of `trex_method` are described as dependency-aware filters of the screening.
  A new section says which values the fields inherited from `TRexSelector` have after a screening run,
  among them `fdp_hat_mat`, which stays 0, because the screening estimates no FDP per voting level.
* The help of `TRexBiobankScreeningSelector` states the decision rule per phenotype, names the values
  of `method_used` and says that `use_bootstrap_CI` of the screening control has no effect there,
  because both variants are computed for every phenotype. Its example printed `NULL`, because
  `Gauss_data$y` is a one-column matrix and the selector then returns a list of records. It now passes
  a vector and prints the selection.
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
* The help of `TRexSPCASelector` and of `trex_spca_control()` is corrected. `cumulative_ev` is a
  fraction of the total variance of the centered design, not a percentage, and `adjusted_ev` is a
  variance, given with its definition. In `"Thresholded"` mode the non-zero loadings can be other
  variables than the selected set `active_sets`, for which the FDR control holds. A component without a
  selected variable has a zero column in `V` and `Z`. A supplied `lambda_2` refers to the scale of the
  standardized columns, not to LARS units. The help of the control says that every other parameter of
  the inner T-Rex+GVS runs stays at its default. The short form of the method is T-Rex+SPCA. The
  vignette "T-Rex SPCA" gives `adjusted_ev` and `cumulative_ev` with their right units.

## New

* `trex_da_control()` gains `bt_selection_mode` and `trex_gvs_control()` gains `group_labels`. Both
  parameters existed in C++ and Python.
* `tsolver_options()` gains `tikhonov_validated`. With `TRUE` the solvers `TIKNET_Solver` and
  `TCDTIKNET_Solver` take their penalty matrix as checked by the caller and skip their own check.

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
