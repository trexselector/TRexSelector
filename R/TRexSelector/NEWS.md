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
* `TRexDASelector` runs the greedy solvers with the exchangeable-tie band by default. `exch_tie_alpha`
  of `trex_control()` defaults to -1, which leaves the choice to the selector: `TRexDASelector` runs
  with a band of 0.25, every other selector without the band. A value of 0 switches the band off, and a
  positive value sets its width. Up to 2.0.1 the default was 0. Without the band the selector exceeded
  a target FDR of 0.2 with every greedy solver at a correlation of 0.9 between neighbouring variables,
  with an FDR of 0.21 to 0.27 in a Monte Carlo study on AR(1) designs (1,000 runs per setting), against
  0.13 to 0.17 with the band. The study shows the control up to a correlation of 0.9. At 0.95 the
  target no longer held, for the path solver `"TLARS"` as well. Selections of `TRexDASelector` with a
  greedy solver can differ from 2.0.1, and `exch_tie_alpha = 0` gives the run of 2.0.1.
* The exchangeable-tie band acts in every greedy solver. `"TGP"`, `"TACGP"`, `"TMP"` and `"TNCGMP"`
  ignored `exch_tie_alpha` up to 2.0.1, in a selector and as the stand-alone classes `TGP_Solver`,
  `TACGP_Solver`, `TMP_Solver` and `TNCGMP_Solver`. The help of the option names the seven solvers.
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
* `TRexTikhonovSelector` and `trex_tikhonov_lambda2_power_rule()` run on the homotopy solver `"TIKNET"`
  and on it alone, as C++ and Python do. Up to 2.0.1 the default was the coordinate-descent solver
  `"TCDTIKNET"`, which stops at its sweep budget before it has converged once `lambda_2` times the
  penalty matrix is large. The default pilot rule then picked the largest `lambda_2` of its grid from
  unconverged runs: on the chain Laplacian of the vignette the default run of 2.0.1 selected 33
  variables with an FDP of 0.70 at a target of 0.1. Results change for every run of 2.0.1 that named
  no solver.
* The `control` of `TRexTikhonovSelector` and of `trex_tikhonov_lambda2_power_rule()` defaults to
  `NULL`, which runs the defaults of `trex_control()` with the solver `"TIKNET"`. A control that is
  passed is taken as it is and has to name the solver, `trex_control(solver = "TIKNET", ...)`. Up to
  2.0.1 the default `"TLARS"` of `trex_control()` was replaced by the default solver without a word,
  and every other solver gave a warning and was replaced. Both are errors now, so a call that passes
  `trex_control(K = 10)` needs the solver in it.
* `TRexTikhonovSelector` no longer accepts the coordinate-descent solver `"TCDTIKNET"`, on any of the
  three surfaces. It was not finally developed: once `lambda_2` times K is large it stops at its sweep
  budget before it has converged. A base control that names it is an error. With it goes
  `fold_dummy_coupling = FALSE` of `trex_tikhonov_control()`, the decoupled dummy block that only this
  solver had, which is an error as well. The class `TCDTIKNET_Solver` stays, deprecated, until 2.1.0.
* `TRexGVSSelector` and `trex_gvs_lambda2_power_rule()` take the base control as C++ and Python do.
  `control` defaults to `NULL`, which runs the defaults of `trex_control()` with the solver `"TENET"`. A
  control that is passed is taken as it is and has to name a solver of the track: `"TENET"` or the
  solver of `en_solver` for `gvs_type = "EN"`, and `"TIENET"` or `"TIENET_AUG"` for `gvs_type = "IEN"`.
  Up to 2.0.1 the default `"TLARS"` of `trex_control()` gave the solver of the track without a word, and
  every other solver gave a warning and was replaced. Both are errors now, so a call that passes
  `trex_control(K = 10)` needs the solver in it, and `gvs_type = "IEN"` needs a base control,
  `trex_control(solver = "TIENET", ...)`, where 2.0.1 ran `"TIENET"` on its own. The selections do not
  change.
* `TRexBiobankScreeningSelector` gives every phenotype a seed of its own. With a seed `s` the first
  phenotype runs with `s`, the second with `s + 1` and so on, for its screening and for its fallback
  T-Rex selector. Up to 2.0.1 all phenotypes of a seeded run competed against the same dummies. The
  results of the first phenotype and of a run with one phenotype do not change, those of the later
  phenotypes of a seeded run do.
* `estimated_FDR` of `TRexScreeningSelector` is 1 for an empty selection, as its estimator
  1 / max(1, R) says. Up to 2.0.1 it was 0. `TRexBiobankScreeningSelector` decides as before with its
  default acceptance window of 0.05 to 0.15, and only the two screening estimates that it reports for
  an empty result move from 0 to 1.
* A coordinate-descent solver that stops at its sweep budget `cd_max_sweeps` before it has converged
  ends with an error. Up to 2.0.1 it warned once, and the selection was built on the iterate it had.
  The message names the budget, the violation and the way out, a larger `cd_max_sweeps` or the homotopy
  solver of the same problem.
* A `select()` that ends with an error gives the design back at once. Up to 2.0.1 the design stayed
  standardized until the selector was removed.
* `tloop_max_stagnant_steps = 0` of `trex_control()` is refused, as in C++ and Python. A value of 0 has
  no meaning of its own, and up to 2.0.1 it acted as 1. Use 1 for the same run.
* The errors of the data normalizer count rows and columns from 1, as R does. Up to 2.0.1 a constant
  fifth column was reported as "column 4 (0-based)". The message now says "column 5", an entry that is
  not finite is named as `X[3, 7]`, and the scaling mode is spelled as the argument,
  `scaling_mode = "NONE"`. This holds for the seven selectors, the two power rules and
  `trex_normalize_design()`.
* The package compiles in less time. The Rcpp bindings compile as one unit, so the build has 37 objects
  where 2.0.1 had 48.

## Fixed

* A selector runs once. `select()` gives the design back when it ends, and up to 2.0.1 a second
  `select()` on the same object ran on the design at its original scale, without an error, and returned
  a wrong selection (all variables in a test with 10 true ones among 60). The second call is an error
  now, for `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`, `TRexScreeningSelector` and
  `TRexTikhonovSelector`. Construct a new selector to run again. The fields keep the results of the run.
* A target FDR that is `NaN` is refused. Up to 2.0.1 it was accepted and nothing was selected.
* The AR1 window of `TRexDASelector` follows the magnitude of `cor_coef`. Up to 2.0.1 a supplied
  negative coefficient gave the window of all variables. `rho_thr_DA = 1` is refused for
  `da_method = "AR1"`, where it switched the correction off without a message. Runs with the estimated
  coefficient do not change.
* A supplied `cor_coef` outside [-1, 1] gives a warning and is replaced by the estimate. Up to 2.0.1
  only the value -2 requested the estimate, and a value such as -1.5 was taken as a coefficient.
* `TRexScreeningSelector` checks `R_boot` and `ci_grid_step` of `trex_screen_control()` at construction.
  `R_boot` is at least 2, and `ci_grid_step` lies in (0, 1]. Up to 2.0.1 a step of 0 never ended.
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
* `trex_tikhonov_control()` checks `prior_groups` the same way: "prior_groups must number the groups 1
  to M without gaps". Up to 2.0.1 labels that started at 0 gave "prior_groups must be 0-based
  non-negative integers", and a missing group was named by its 0-based number. A wrong length is
  reported with the name of the argument, by the selector and by the power rule.
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
  a vector and prints the selection. The help also says what `estimated_FDR` of a result is, the
  estimate of the accepted screening variant or, for a fallback result, the target level of the T-Rex
  selector, and that `used_fallback_trex` tells the two apart.
* The help example of `TRexDASelector` runs on an AR(1) design. The example of 2.0.1 ran on independent
  columns and selected nothing.
* The vignette "T-Rex Tikhonov: a general penalty matrix" and the help of `trex_tikhonov_control()` say
  which penalty is which. The group-mean matrix reproduces the penalty of the informed elastic net, the
  centering matrix of `penalty_structure = "GROUP_CENTERING"` asks the members of a group to share a
  coefficient and needs `prior_groups`, the penalty carries the factor one half of the code, and
  `tikhonov_K` is optional under `"GROUP_CENTERING"`. The vignette no longer runs the deprecated solver
  `"TCDTIKNET"` with the decoupled dummy block.
* The help of `TRexTikhonovSelector`, of `trex_tikhonov_control()` and of
  `trex_tikhonov_lambda2_power_rule()` is rewritten after the C++ class. The class gains the sections on
  the penalty, the solver, the dummies, the L-loop strategies and the rules for `lambda_2`. Wrong
  statements are corrected. The class no longer "drives the Tikhonov-net solvers" with TCDTIKNET as an
  equal choice. `lambda_2` and `lambda2_used` were "on the solver scale" and refer to the scale of the
  standardized columns. `lambda2_used` named a supplied value and a cross-validation only, the default
  is the pick of the pilot rule. `prior_groups` was "ignored when `cluster_dummies = FALSE`", although
  the groups also define the penalty under `"GROUP_CENTERING"` and the count `"GROUP_MAJORITY"`.
  `corr_max` bounded the correlation between clusters, which holds for single linkage only. `"ward"`
  linkage is named as rejected. The record `lambda2_pilot` and the result of the power rule say what
  every field holds, among them `curve_soft`, `calibrated`, `group_majority`, `groups_qualifying` and
  `probe_seconds`. The vignette "T-Rex Tikhonov" gains the paragraph on the scale of `lambda_2` and of
  the grids.
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
* The help of `TRexSelector` says what `T_stop` is, the stopping time at which the T-loop stopped, and
  that the selection is taken from an earlier stopping time. It states the rule of the selection: the
  largest number of selected variables among the stopping times before `T_stop` whose FDP estimate does
  not exceed the target. `r_mat` therefore has one row less than `fdp_hat_mat` and `phi_mat`. The help
  of `x_l2_norms` says that the field holds the scale divisors of the scaling mode: the L2 norms under
  `"L2"`, the sample standard deviations under `"ZSCORE"`, ones under `"NONE"`. The vignette "The T-Rex
  selector" explains `T_stop` in the same way.
* The help of `trex_control()` says what `opt_threshold` is (the reference voting level of the L-loop),
  what each L-loop strategy does, that `use_memory_mapping` concerns the dummy matrices and what
  `tloop_max_stagnant_steps` counts. Its note on seeds no longer says that a fixed seed raises the FDR
  of a Monte Carlo study: a seed of 0 or more makes a run reproducible, and a study gives every trial a
  seed of its own or uses -1.
* The help of `trex_dummy_distribution()` names the defaults, among them `s = 0.1` of
  `"ConstrainedSparseRademacher"`, and says that `a` is also the lower bound of `"Triangle"` and `b` its
  mode.
* The help of `TRexScreeningSelector` gains the section "A response without signal". On a response of
  pure noise the screening selected at least one variable in a quarter to a half of the runs in
  simulations, so a small selection is no evidence of a signal, and `estimated_FDR` rests on the
  expected number of false discoveries, about one, and is no bound on the FDP of a run. The help of
  `trex_screen_control()` says that `use_bootstrap_CI = TRUE` mostly selects nothing together with a
  dependency-aware variant. The vignette on the screening says the same and no longer claims that the
  L-loop strategy is forced to `"STANDARD"`.
* The vignette "The T-Rex selector" describes the stagnation stop as it works (the selection no longer
  grows) and says that `scaling_mode = "NONE"` rejects all-zero columns only.
* The references in the help of `TSTAGEWISE_Solver` and `TAFS_Solver` carry the numbers of the two
  articles in their volumes of the Journal of Machine Learning Research. The reference in the help of
  `TNCGMP_Solver` names the Proceedings of Machine Learning Research as its series and carries its pages.

## New

* `trex_screen_control()` takes `trex_method = "TREX_DA_NN"`, a nearest-neighbour variant of the
  screening, as C++ and Python do. A variable is compared with the neighbours inside the AR(1) window
  whose sample correlation reaches `rho_thr_DA`, and a variable without such a neighbour keeps its
  relative occurrence. In simulations (1,000 runs per setting) the variant selected as `"TREX_DA_AR1"`
  at the default threshold and gained power at larger thresholds, at the price of more false
  discoveries. With a threshold at or above the correlation of the data no neighbour is left, and the
  variant selects as `"TREX"`.
* `trex_screen_control()` gains `cor_coef`, `rho_thr_DA` and `n_blocks`, the parameters of the
  dependency-aware variants, as in C++ and Python. They are checked at construction: `rho_thr_DA` lies
  in [0, 1], below 1 for `"TREX_DA_AR1"`, and `n_blocks` between 1 and p for `"TREX_DA_BLOCK_EQUI"`.
* A memory-mapped design opened with `mem_matrix(mode = "readonly")` serves wherever nothing is written
  into it: every selector under `scaling_mode = "NONE"` and every T-solver with `normalize` and
  `intercept` off. Up to 2.0.1 such a design was refused everywhere. `TRexSPCASelector` centers the
  design in every run and refuses it with a message of its own.
* The dummy matrix `D` of a T-solver and the response matrix `Y` of `TRexBiobankScreeningSelector` can
  be memory-mapped files, as the design `X` could.
* `trex_dummy_distribution()` gains `p1` and `p2`, the two probabilities of `"Mammen"`, as in C++ and
  Python. Without them the golden-ratio probabilities apply, as before.
* `trex_da_control()` gains `bt_selection_mode` and `trex_gvs_control()` gains `group_labels`. Both
  parameters existed in C++ and Python.
* `tsolver_options()` gains `tikhonov_validated`. With `TRUE` the solvers `TIKNET_Solver` and
  `TCDTIKNET_Solver` take their penalty matrix as checked by the caller and skip their own check.
* `TRexTikhonovSelector` gains what the C++ class has. The method `run_pilot_only()` runs the `lambda_2`
  pilot without the selection and returns its record. The field `dummy_clusters` reports the clusters
  the dummies were drawn for, with the names `groups`, `max_clusters` and `hc_method_used` of
  `TRexGVSSelector`. The pilot record `lambda2_pilot` gains `dummy_tag`.
* `trex_gvs_lambda2_power_rule()` is the power rule on `TRexGVSSelector`: one selection per value of a
  `lambda_2` grid with common dummies, and the knee of that profile.
* `trex_tikhonov_lambda2_power_rule()` takes a memory-mapped design, as the selectors do. Up to 2.0.1 it
  refused one.
* `trex_default_lambda2_grid()` returns the default `lambda_2` grid of the C++ core. It is the default
  of `pilot_grid` in `trex_tikhonov_control()` and of `grid` in the two power rules, which held the
  numbers themselves up to 2.0.1.
* `trex_build_group_centering_k()` builds the centering matrix of a group assignment, the penalty of
  `penalty_structure = "GROUP_CENTERING"`.

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
  receives one of the solver names. The warning comes once per call and names the replacement. The name
  `"TCDTIKNET"` is the exception, because `TRexTikhonovSelector` refuses it (see "Changed").
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
